#include "lbm/two_phase_solver3d.hpp"

#include "lbm/grid3d.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace lbm {
namespace {

constexpr double min_density = 1.0e-12;
constexpr double pi = 3.141592653589793238462643383279502884;

bool positive_half_link(int cx, int cy, int cz) {
    return cx > 0 || (cx == 0 && cy > 0) || (cx == 0 && cy == 0 && cz > 0);
}

} // namespace

TwoPhaseSolver3D::TwoPhaseSolver3D(
    int nx, int ny, int nz, TwoPhaseConfig3D config)
    : nx_(nx),
      ny_(ny),
      nz_(nz),
      q_(lattice3d(config.lattice_model).q),
      config_(config),
      fa_(detail::checked_grid3d_cells(nx, ny, nz) * static_cast<std::size_t>(q_), 0.0),
      fb_(fa_.size(), 0.0),
      fa_next_(fa_.size(), 0.0),
      fb_next_(fa_.size(), 0.0),
      rho_a_(detail::checked_grid3d_cells(nx, ny, nz), config.rho_low),
      rho_b_(rho_a_.size(), config.rho_high),
      rho_(rho_a_.size(), config.rho_high + config.rho_low),
      porosity_(rho_a_.size(), config.free_flow_porosity),
      solid_(rho_a_.size(), 0) {
    if (nx <= 4 || ny <= 4 || nz <= 4) {
        throw std::invalid_argument(
            "Two-phase 3D grid dimensions must all be greater than 4.");
    }
    detail::validate_collision_parameters3d(
        config_.collision_model, config_.tau_a, config_.mrt);
    detail::validate_collision_parameters3d(
        config_.collision_model, config_.tau_b, config_.mrt);
    if (config_.rho_high <= 0.0 || config_.rho_low <= 0.0) {
        throw std::invalid_argument("Three-dimensional phase densities must be positive.");
    }
    if (config_.droplet_interface_width <= 0.0) {
        throw std::invalid_argument("Three-dimensional droplet interface width must be positive.");
    }
    if (config_.recoloring_strength < 0.0 || config_.recoloring_strength > 1.0) {
        throw std::invalid_argument("Three-dimensional recoloring strength must be in [0, 1].");
    }
    if (config_.bottom_wall_thickness < 1 ||
        config_.bottom_wall_thickness >= ny_ - 2) {
        throw std::invalid_argument("Invalid three-dimensional bottom wall thickness.");
    }
    if (config_.bottom_wall_repulsion_strength < 0.0 ||
        config_.bottom_wall_repulsion_range < 0) {
        throw std::invalid_argument("Three-dimensional wall repulsion must be non-negative.");
    }
    if (config_.free_flow_porosity <= 0.0 || config_.free_flow_porosity > 1.0 ||
        config_.porous_porosity <= 0.0 || config_.porous_porosity > 1.0) {
        throw std::invalid_argument("Three-dimensional porosity must be in (0, 1].");
    }
    if (config_.porous_pore_diameter <= 0.0 || config_.darcy_drag_scale < 0.0 ||
        config_.forchheimer_drag_scale < 0.0) {
        throw std::invalid_argument("Invalid three-dimensional porous-medium parameters.");
    }

    for (auto& component : velocity_) {
        component.assign(rho_.size(), 0.0);
    }
    for (auto& component : force_a_) {
        component.assign(rho_.size(), 0.0);
    }
    for (auto& component : force_b_) {
        component.assign(rho_.size(), 0.0);
    }
}

int TwoPhaseSolver3D::nx() const {
    return nx_;
}

int TwoPhaseSolver3D::ny() const {
    return ny_;
}

int TwoPhaseSolver3D::nz() const {
    return nz_;
}

double TwoPhaseSolver3D::phase_at(int x, int y, int z) const {
    if (is_outside(x, y, z)) {
        throw std::out_of_range("phase_at coordinates are outside the 3D grid.");
    }
    const std::size_t s = scalar_index(x, y, z);
    return (rho_a_[s] - rho_b_[s]) / std::max(rho_[s], min_density);
}

double TwoPhaseSolver3D::porosity_at(int x, int y, int z) const {
    if (is_outside(x, y, z)) {
        throw std::out_of_range("porosity_at coordinates are outside the 3D grid.");
    }
    return porosity_[scalar_index(x, y, z)];
}

bool TwoPhaseSolver3D::solid_at(int x, int y, int z) const {
    if (is_outside(x, y, z)) {
        throw std::out_of_range("solid_at coordinates are outside the 3D grid.");
    }
    return solid_[scalar_index(x, y, z)] != 0;
}

void TwoPhaseSolver3D::initialize_droplet_impact(
    double center_x,
    double center_y,
    double center_z,
    double radius,
    double initial_ux,
    double initial_uy,
    double initial_uz) {
    if (radius <= 0.0) {
        throw std::invalid_argument("Three-dimensional droplet radius must be positive.");
    }
    geometry_initialized_ = false;
    geometry_cells_.clear();

    for (int z = 0; z < nz_; ++z) {
        for (int y = 0; y < ny_; ++y) {
            for (int x = 0; x < nx_; ++x) {
                const std::size_t s = scalar_index(x, y, z);
                const bool wall = x == 0 || x == nx_ - 1 ||
                                  y < config_.bottom_wall_thickness || y == ny_ - 1 ||
                                  z == 0 || z == nz_ - 1;
                solid_[s] = wall ? 1 : 0;
                porosity_[s] = wall ? 0.0 : config_.free_flow_porosity;

                const double dx = static_cast<double>(x) - center_x;
                const double dy = static_cast<double>(y) - center_y;
                const double dz = static_cast<double>(z) - center_z;
                const double signed_distance =
                    radius - std::sqrt(dx * dx + dy * dy + dz * dz);
                const double fraction =
                    wall ? 0.0
                         : 0.5 * (1.0 + std::tanh(
                                            signed_distance /
                                            config_.droplet_interface_width));
                const double rho_a =
                    config_.rho_low +
                    (config_.rho_high - config_.rho_low) * fraction;
                const double rho_b =
                    config_.rho_high -
                    (config_.rho_high - config_.rho_low) * fraction;
                set_equilibrium_cell(
                    x,
                    y,
                    z,
                    rho_a,
                    rho_b,
                    initial_ux * fraction,
                    initial_uy * fraction,
                    initial_uz * fraction);
            }
        }
    }
    compute_macroscopic();
}

void TwoPhaseSolver3D::initialize_geometry_displacement(
    const GeometryMask3D& mask) {
    if (mask.width() != nx_ || mask.height() != ny_ || mask.depth() != nz_) {
        throw std::invalid_argument(
            "3D geometry dimensions must match the two-phase solver.");
    }
    if (!mask.contains(GeometryCell::Inlet) || !mask.contains(GeometryCell::Outlet)) {
        throw std::invalid_argument(
            "Three-dimensional displacement requires an inlet and an outlet.");
    }
    geometry_cells_ = mask.cells();
    geometry_initialized_ = true;

    for (int z = 0; z < nz_; ++z) {
        for (int y = 0; y < ny_; ++y) {
            for (int x = 0; x < nx_; ++x) {
                const std::size_t s = scalar_index(x, y, z);
                const GeometryCell cell = geometry_cells_[s];
                if (cell == GeometryCell::Inlet &&
                    geometry_cells_[scalar_index(1, y, z)] == GeometryCell::Solid) {
                    throw std::invalid_argument(
                        "Each 3D inlet must connect to fluid on its right.");
                }
                if (cell == GeometryCell::Outlet &&
                    geometry_cells_[scalar_index(nx_ - 2, y, z)] ==
                        GeometryCell::Solid) {
                    throw std::invalid_argument(
                        "Each 3D outlet must connect to fluid on its left.");
                }
                solid_[s] = cell == GeometryCell::Solid ? 1 : 0;
                porosity_[s] =
                    solid_[s] ? 0.0
                              : (cell == GeometryCell::Porous
                                     ? config_.porous_porosity
                                     : config_.free_flow_porosity);
                const bool inlet = cell == GeometryCell::Inlet;
                set_equilibrium_cell(
                    x,
                    y,
                    z,
                    inlet ? config_.rho_high : config_.rho_low,
                    solid_[s]
                        ? config_.rho_low
                        : (inlet ? config_.rho_low : config_.rho_high),
                    inlet ? config_.inlet_velocity : 0.0,
                    0.0,
                    0.0);
            }
        }
    }
    apply_geometry_boundaries();
    compute_macroscopic();
}

void TwoPhaseSolver3D::step() {
    if (!geometry_initialized_) {
        throw std::logic_error(
            "initialize_geometry_displacement must be called before open 3D steps.");
    }
    compute_macroscopic();
    compute_forces();
    collide();
    recolor();
    stream();
    apply_geometry_boundaries();
    compute_macroscopic();
}

void TwoPhaseSolver3D::step_closed() {
    compute_macroscopic();
    compute_forces();
    collide();
    recolor();
    stream();
    compute_macroscopic();
}

void TwoPhaseSolver3D::run(int steps) {
    if (steps < 0) {
        throw std::invalid_argument("Steps must be non-negative.");
    }
    for (int i = 0; i < steps; ++i) {
        step();
    }
}

void TwoPhaseSolver3D::run_closed(int steps) {
    if (steps < 0) {
        throw std::invalid_argument("Steps must be non-negative.");
    }
    for (int i = 0; i < steps; ++i) {
        step_closed();
    }
}

TwoPhaseDiagnostics3D TwoPhaseSolver3D::diagnostics() const {
    TwoPhaseDiagnostics3D diagnostics{};
    double weight_sum = 0.0;
    double pore_speed_sum = 0.0;
    std::size_t porous_cells = 0;

    for (int z = 0; z < nz_; ++z) {
        for (int y = 0; y < ny_; ++y) {
            for (int x = 0; x < nx_; ++x) {
                const std::size_t s = scalar_index(x, y, z);
                if (solid_[s]) {
                    continue;
                }
                diagnostics.mass_a += rho_a_[s];
                diagnostics.mass_b += rho_b_[s];
                const double speed =
                    std::sqrt(
                        velocity_[0][s] * velocity_[0][s] +
                        velocity_[1][s] * velocity_[1][s] +
                        velocity_[2][s] * velocity_[2][s]);
                diagnostics.max_speed = std::max(diagnostics.max_speed, speed);
                const double fraction = std::clamp(
                    (rho_a_[s] - config_.rho_low) /
                        std::max(
                            config_.rho_high - config_.rho_low, min_density),
                    0.0,
                    1.0);
                diagnostics.phase_a_centroid_x += x * fraction;
                diagnostics.phase_a_centroid_y += y * fraction;
                diagnostics.phase_a_centroid_z += z * fraction;
                weight_sum += fraction;

                if (porosity_[s] > 0.0 &&
                    porosity_[s] < config_.free_flow_porosity) {
                    const double pore_speed = speed / porosity_[s];
                    pore_speed_sum += pore_speed;
                    diagnostics.porous_max_pore_speed =
                        std::max(diagnostics.porous_max_pore_speed, pore_speed);
                    ++porous_cells;
                }
            }
        }
    }
    if (weight_sum > 0.0) {
        diagnostics.phase_a_centroid_x /= weight_sum;
        diagnostics.phase_a_centroid_y /= weight_sum;
        diagnostics.phase_a_centroid_z /= weight_sum;
    }
    if (porous_cells > 0) {
        diagnostics.porous_mean_pore_speed =
            pore_speed_sum / static_cast<double>(porous_cells);
    }
    diagnostics.interaction_force_balance_x = interaction_force_balance_[0];
    diagnostics.interaction_force_balance_y = interaction_force_balance_[1];
    diagnostics.interaction_force_balance_z = interaction_force_balance_[2];
    return diagnostics;
}

void TwoPhaseSolver3D::write_vtk(const std::string& path) const {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Failed to open two-phase 3D VTK output: " + path);
    }
    out << "# vtk DataFile Version 3.0\n";
    out << "Three-dimensional two-phase lattice Boltzmann flow\n";
    out << "ASCII\nDATASET STRUCTURED_POINTS\n";
    out << "DIMENSIONS " << nx_ << ' ' << ny_ << ' ' << nz_ << '\n';
    out << "ORIGIN 0 0 0\nSPACING 1 1 1\n";
    out << "POINT_DATA " << rho_.size() << '\n';
    out << std::setprecision(17);
    out << "SCALARS phase double 1\nLOOKUP_TABLE default\n";
    for (std::size_t s = 0; s < rho_.size(); ++s) {
        out << (rho_a_[s] - rho_b_[s]) / std::max(rho_[s], min_density)
            << '\n';
    }
    out << "SCALARS porosity double 1\nLOOKUP_TABLE default\n";
    for (double value : porosity_) {
        out << value << '\n';
    }
    out << "SCALARS solid int 1\nLOOKUP_TABLE default\n";
    for (std::uint8_t value : solid_) {
        out << static_cast<int>(value) << '\n';
    }
    out << "VECTORS velocity double\n";
    for (std::size_t s = 0; s < rho_.size(); ++s) {
        out << velocity_[0][s] << ' ' << velocity_[1][s] << ' '
            << velocity_[2][s] << '\n';
    }
}

std::size_t TwoPhaseSolver3D::scalar_index(int x, int y, int z) const {
    return (static_cast<std::size_t>(z) * ny_ + y) * nx_ + x;
}

std::size_t TwoPhaseSolver3D::dist_index(
    int x, int y, int z, int direction) const {
    return scalar_index(x, y, z) * static_cast<std::size_t>(q_) + direction;
}

double TwoPhaseSolver3D::psi(double rho) {
    return 1.0 - std::exp(-std::max(rho, 0.0));
}

void TwoPhaseSolver3D::set_equilibrium_cell(
    int x,
    int y,
    int z,
    double rho_a,
    double rho_b,
    double ux,
    double uy,
    double uz) {
    const std::size_t s = scalar_index(x, y, z);
    rho_a_[s] = rho_a;
    rho_b_[s] = rho_b;
    rho_[s] = rho_a + rho_b;
    velocity_[0][s] = solid_[s] ? 0.0 : ux;
    velocity_[1][s] = solid_[s] ? 0.0 : uy;
    velocity_[2][s] = solid_[s] ? 0.0 : uz;
    for (int q = 0; q < q_; ++q) {
        fa_[dist_index(x, y, z, q)] = detail::equilibrium3d(
            config_.lattice_model,
            q,
            rho_a,
            velocity_[0][s],
            velocity_[1][s],
            velocity_[2][s]);
        fb_[dist_index(x, y, z, q)] = detail::equilibrium3d(
            config_.lattice_model,
            q,
            rho_b,
            velocity_[0][s],
            velocity_[1][s],
            velocity_[2][s]);
    }
}

void TwoPhaseSolver3D::compute_macroscopic() {
    const auto& lattice = lattice3d(config_.lattice_model);
    for (int z = 0; z < nz_; ++z) {
        for (int y = 0; y < ny_; ++y) {
            for (int x = 0; x < nx_; ++x) {
                const std::size_t s = scalar_index(x, y, z);
                if (solid_[s]) {
                    rho_a_[s] = config_.rho_low;
                    rho_b_[s] = config_.rho_low;
                    rho_[s] = rho_a_[s] + rho_b_[s];
                    velocity_[0][s] = 0.0;
                    velocity_[1][s] = 0.0;
                    velocity_[2][s] = 0.0;
                    continue;
                }
                double rho_a = 0.0;
                double rho_b = 0.0;
                std::array<double, 3> momentum{};
                for (int q = 0; q < q_; ++q) {
                    const double fa = fa_[dist_index(x, y, z, q)];
                    const double fb = fb_[dist_index(x, y, z, q)];
                    rho_a += fa;
                    rho_b += fb;
                    const double total = fa + fb;
                    momentum[0] += lattice.cx[q] * total;
                    momentum[1] += lattice.cy[q] * total;
                    momentum[2] += lattice.cz[q] * total;
                }
                if (rho_a <= 0.0 || rho_b <= 0.0 || !std::isfinite(rho_a) ||
                    !std::isfinite(rho_b)) {
                    throw std::runtime_error(
                        "Non-physical density encountered in TwoPhaseSolver3D.");
                }
                rho_a_[s] = rho_a;
                rho_b_[s] = rho_b;
                rho_[s] = rho_a + rho_b;
                for (int axis = 0; axis < 3; ++axis) {
                    velocity_[axis][s] =
                        (momentum[axis] +
                         0.5 * (force_a_[axis][s] + force_b_[axis][s])) /
                        rho_[s];
                }
            }
        }
    }
}

void TwoPhaseSolver3D::compute_forces() {
    const auto& lattice = lattice3d(config_.lattice_model);
    for (int axis = 0; axis < 3; ++axis) {
        std::fill(force_a_[axis].begin(), force_a_[axis].end(), 0.0);
        std::fill(force_b_[axis].begin(), force_b_[axis].end(), 0.0);
        interaction_force_balance_[axis] = 0.0;
    }

    for (int z = 0; z < nz_; ++z) {
        for (int y = 0; y < ny_; ++y) {
            for (int x = 0; x < nx_; ++x) {
                const std::size_t s = scalar_index(x, y, z);
                if (solid_[s]) {
                    continue;
                }
                for (int q = 1; q < q_; ++q) {
                    if (!positive_half_link(
                            lattice.cx[q], lattice.cy[q], lattice.cz[q])) {
                        continue;
                    }
                    const int nx = x + lattice.cx[q];
                    const int ny = y + lattice.cy[q];
                    const int nz = z + lattice.cz[q];
                    if (is_outside(nx, ny, nz)) {
                        continue;
                    }
                    const std::size_t neighbor = scalar_index(nx, ny, nz);
                    if (solid_[neighbor]) {
                        continue;
                    }
                    const double a_to_b =
                        -config_.interaction_strength * lattice.w[q] *
                        psi(rho_a_[s]) * psi(rho_b_[neighbor]);
                    const double b_to_a =
                        -config_.interaction_strength * lattice.w[q] *
                        psi(rho_b_[s]) * psi(rho_a_[neighbor]);
                    const std::array<int, 3> direction{
                        lattice.cx[q], lattice.cy[q], lattice.cz[q]};
                    for (int axis = 0; axis < 3; ++axis) {
                        const double c = static_cast<double>(direction[axis]);
                        force_a_[axis][s] += a_to_b * c;
                        force_b_[axis][neighbor] -= a_to_b * c;
                        force_b_[axis][s] += b_to_a * c;
                        force_a_[axis][neighbor] -= b_to_a * c;
                    }
                }
            }
        }
    }
    for (std::size_t s = 0; s < rho_.size(); ++s) {
        for (int axis = 0; axis < 3; ++axis) {
            interaction_force_balance_[axis] +=
                force_a_[axis][s] + force_b_[axis][s];
        }
    }

    // x 方向开放边界外采用零法向梯度虚拟储液层。
    if (geometry_initialized_) {
        for (int z = 0; z < nz_; ++z) {
            for (int y = 0; y < ny_; ++y) {
                for (int x = 0; x < nx_; ++x) {
                    const std::size_t s = scalar_index(x, y, z);
                    if (solid_[s]) {
                        continue;
                    }
                    for (int q = 1; q < q_; ++q) {
                        const int gx = x + lattice.cx[q];
                        const int gy = y + lattice.cy[q];
                        const int gz = z + lattice.cz[q];
                        if (!is_outside(gx, gy, gz) || gy < 0 || gy >= ny_ ||
                            gz < 0 || gz >= nz_) {
                            continue;
                        }
                        const std::size_t ghost =
                            scalar_index(std::clamp(gx, 0, nx_ - 1), gy, gz);
                        if (solid_[ghost]) {
                            continue;
                        }
                        const double fa =
                            -config_.interaction_strength * lattice.w[q] *
                            psi(rho_a_[s]) * psi(rho_b_[ghost]);
                        const double fb =
                            -config_.interaction_strength * lattice.w[q] *
                            psi(rho_b_[s]) * psi(rho_a_[ghost]);
                        force_a_[0][s] += fa * lattice.cx[q];
                        force_a_[1][s] += fa * lattice.cy[q];
                        force_a_[2][s] += fa * lattice.cz[q];
                        force_b_[0][s] += fb * lattice.cx[q];
                        force_b_[1][s] += fb * lattice.cy[q];
                        force_b_[2][s] += fb * lattice.cz[q];
                    }
                }
            }
        }
    }

    const double wetting_bias =
        config_.wall_adhesion_strength *
        std::cos(config_.contact_angle_degrees * pi / 180.0);
    const std::array<double, 3> body_force{
        config_.body_force_x, config_.body_force_y, config_.body_force_z};

    for (int z = 0; z < nz_; ++z) {
        for (int y = 0; y < ny_; ++y) {
            for (int x = 0; x < nx_; ++x) {
                const std::size_t s = scalar_index(x, y, z);
                if (solid_[s]) {
                    continue;
                }
                std::array<double, 3> wall_normal{};
                for (int q = 1; q < q_; ++q) {
                    const int wx = x + lattice.cx[q];
                    const int wy = y + lattice.cy[q];
                    const int wz = z + lattice.cz[q];
                    if (is_outside(wx, wy, wz) ||
                        !solid_[scalar_index(wx, wy, wz)]) {
                        continue;
                    }
                    wall_normal[0] += lattice.w[q] * lattice.cx[q];
                    wall_normal[1] += lattice.w[q] * lattice.cy[q];
                    wall_normal[2] += lattice.w[q] * lattice.cz[q];
                }
                const double psi_a = psi(rho_a_[s]);
                const double psi_b = psi(rho_b_[s]);
                const double total_rho = std::max(rho_[s], min_density);
                for (int axis = 0; axis < 3; ++axis) {
                    force_a_[axis][s] +=
                        wetting_bias * psi_a * wall_normal[axis];
                    force_b_[axis][s] -=
                        wetting_bias * psi_b * wall_normal[axis];
                    force_a_[axis][s] +=
                        body_force[axis] * rho_a_[s] / total_rho;
                    force_b_[axis][s] +=
                        body_force[axis] * rho_b_[s] / total_rho;
                }

                if (config_.bottom_wall_repulsion_strength > 0.0 &&
                    config_.bottom_wall_repulsion_range > 0) {
                    const int distance = y - config_.bottom_wall_thickness;
                    if (distance >= 0 &&
                        distance < config_.bottom_wall_repulsion_range) {
                        const double fraction =
                            1.0 - static_cast<double>(distance) /
                                      config_.bottom_wall_repulsion_range;
                        const double saturation = rho_a_[s] / total_rho;
                        if (saturation > 0.08) {
                            force_a_[1][s] +=
                                config_.bottom_wall_repulsion_strength *
                                fraction * fraction * saturation * psi_a;
                        }
                    }
                }

                if (porosity_[s] > 0.0 &&
                    porosity_[s] < config_.free_flow_porosity) {
                    const double permeability =
                        permeability_from_porosity(porosity_[s]);
                    const double viscosity =
                        Lattice3DDescriptor::cs2 *
                        (0.5 * (config_.tau_a + config_.tau_b) - 0.5);
                    std::array<double, 3> pore_velocity{};
                    double pore_speed_squared = 0.0;
                    for (int axis = 0; axis < 3; ++axis) {
                        pore_velocity[axis] =
                            velocity_[axis][s] / porosity_[s];
                        pore_speed_squared +=
                            pore_velocity[axis] * pore_velocity[axis];
                    }
                    const double pore_speed = std::sqrt(pore_speed_squared);
                    const double darcy =
                        config_.darcy_drag_scale * viscosity /
                        std::max(permeability, min_density);
                    const double forchheimer =
                        config_.forchheimer_drag_scale /
                        std::sqrt(std::max(permeability, min_density));
                    for (int axis = 0; axis < 3; ++axis) {
                        const double drag =
                            -rho_[s] * (darcy + forchheimer * pore_speed) *
                            pore_velocity[axis];
                        force_a_[axis][s] += drag * rho_a_[s] / total_rho;
                        force_b_[axis][s] += drag * rho_b_[s] / total_rho;
                    }
                }
            }
        }
    }
}

void TwoPhaseSolver3D::collide() {
    for (int z = 0; z < nz_; ++z) {
        for (int y = 0; y < ny_; ++y) {
            for (int x = 0; x < nx_; ++x) {
                const std::size_t s = scalar_index(x, y, z);
                if (solid_[s]) {
                    continue;
                }
                detail::D3Population population_a{};
                detail::D3Population population_b{};
                for (int q = 0; q < q_; ++q) {
                    population_a[q] = fa_[dist_index(x, y, z, q)];
                    population_b[q] = fb_[dist_index(x, y, z, q)];
                }
                const detail::D3Population source_a = detail::guo_source3d(
                    config_.lattice_model,
                    velocity_[0][s],
                    velocity_[1][s],
                    velocity_[2][s],
                    force_a_[0][s],
                    force_a_[1][s],
                    force_a_[2][s]);
                const detail::D3Population source_b = detail::guo_source3d(
                    config_.lattice_model,
                    velocity_[0][s],
                    velocity_[1][s],
                    velocity_[2][s],
                    force_b_[0][s],
                    force_b_[1][s],
                    force_b_[2][s]);
                detail::collide_population3d(
                    config_.lattice_model,
                    population_a,
                    rho_a_[s],
                    velocity_[0][s],
                    velocity_[1][s],
                    velocity_[2][s],
                    config_.tau_a,
                    config_.collision_model,
                    config_.mrt,
                    source_a,
                    1.0 / config_.tau_a);
                detail::collide_population3d(
                    config_.lattice_model,
                    population_b,
                    rho_b_[s],
                    velocity_[0][s],
                    velocity_[1][s],
                    velocity_[2][s],
                    config_.tau_b,
                    config_.collision_model,
                    config_.mrt,
                    source_b,
                    1.0 / config_.tau_b);
                for (int q = 0; q < q_; ++q) {
                    fa_[dist_index(x, y, z, q)] = population_a[q];
                    fb_[dist_index(x, y, z, q)] = population_b[q];
                }
            }
        }
    }
}

void TwoPhaseSolver3D::recolor() {
    if (config_.recoloring_strength <= 0.0) {
        return;
    }
    const auto& lattice = lattice3d(config_.lattice_model);
    std::vector<double> phase(rho_.size(), 0.0);
    for (std::size_t s = 0; s < rho_.size(); ++s) {
        if (!solid_[s]) {
            phase[s] =
                (rho_a_[s] - rho_b_[s]) / std::max(rho_[s], min_density);
        }
    }

    for (int z = 0; z < nz_; ++z) {
        for (int y = 0; y < ny_; ++y) {
            for (int x = 0; x < nx_; ++x) {
                const std::size_t s = scalar_index(x, y, z);
                if (solid_[s]) {
                    continue;
                }
                std::array<double, 3> gradient{};
                for (int q = 1; q < q_; ++q) {
                    const int gx = x + lattice.cx[q];
                    const int gy = y + lattice.cy[q];
                    const int gz = z + lattice.cz[q];
                    if (is_outside(gx, gy, gz)) {
                        continue;
                    }
                    const std::size_t neighbor = scalar_index(gx, gy, gz);
                    if (solid_[neighbor]) {
                        continue;
                    }
                    gradient[0] +=
                        lattice.w[q] * phase[neighbor] * lattice.cx[q];
                    gradient[1] +=
                        lattice.w[q] * phase[neighbor] * lattice.cy[q];
                    gradient[2] +=
                        lattice.w[q] * phase[neighbor] * lattice.cz[q];
                }
                const double gradient_norm =
                    std::sqrt(
                        gradient[0] * gradient[0] +
                        gradient[1] * gradient[1] +
                        gradient[2] * gradient[2]);
                if (gradient_norm <= 1.0e-12) {
                    continue;
                }
                const double total_rho = std::max(rho_[s], min_density);
                const double ratio_a =
                    std::clamp(rho_a_[s] / total_rho, 0.0, 1.0);
                if (ratio_a < 0.04 || ratio_a > 0.96) {
                    continue;
                }
                const double strength =
                    config_.recoloring_strength * ratio_a *
                    (1.0 - ratio_a) * total_rho;
                detail::D3Population total{};
                detail::D3Population recolored{};
                double recolored_mass = 0.0;
                for (int q = 0; q < q_; ++q) {
                    const std::size_t k = dist_index(x, y, z, q);
                    total[q] = fa_[k] + fb_[k];
                    if (total[q] <= min_density) {
                        continue;
                    }
                    const double direction_norm = std::sqrt(
                        static_cast<double>(
                            lattice.cx[q] * lattice.cx[q] +
                            lattice.cy[q] * lattice.cy[q] +
                            lattice.cz[q] * lattice.cz[q]));
                    const double alignment =
                        q == 0
                            ? 0.0
                            : (lattice.cx[q] * gradient[0] +
                               lattice.cy[q] * gradient[1] +
                               lattice.cz[q] * gradient[2]) /
                                  (direction_norm * gradient_norm);
                    recolored[q] = std::clamp(
                        ratio_a * total[q] +
                            strength * lattice.w[q] * alignment,
                        0.0,
                        total[q]);
                    recolored_mass += recolored[q];
                }
                if (recolored_mass > min_density) {
                    const double factor = rho_a_[s] / recolored_mass;
                    for (int q = 0; q < q_; ++q) {
                        recolored[q] =
                            std::clamp(recolored[q] * factor, 0.0, total[q]);
                    }
                }
                double corrected_mass = 0.0;
                for (int q = 0; q < q_; ++q) {
                    corrected_mass += recolored[q];
                }
                const double residual = rho_a_[s] - corrected_mass;
                if (std::abs(residual) > 1.0e-12) {
                    double capacity = 0.0;
                    for (int q = 0; q < q_; ++q) {
                        capacity +=
                            residual > 0.0 ? total[q] - recolored[q]
                                           : recolored[q];
                    }
                    if (capacity > min_density) {
                        for (int q = 0; q < q_; ++q) {
                            const double share =
                                residual > 0.0 ? total[q] - recolored[q]
                                               : recolored[q];
                            recolored[q] += residual * share / capacity;
                        }
                    }
                }
                for (int q = 0; q < q_; ++q) {
                    const std::size_t k = dist_index(x, y, z, q);
                    fa_[k] = std::clamp(recolored[q], 0.0, total[q]);
                    fb_[k] = total[q] - fa_[k];
                }
            }
        }
    }
}

void TwoPhaseSolver3D::stream() {
    const auto& lattice = lattice3d(config_.lattice_model);
    std::fill(fa_next_.begin(), fa_next_.end(), 0.0);
    std::fill(fb_next_.begin(), fb_next_.end(), 0.0);
    for (int z = 0; z < nz_; ++z) {
        for (int y = 0; y < ny_; ++y) {
            for (int x = 0; x < nx_; ++x) {
                if (solid_[scalar_index(x, y, z)]) {
                    continue;
                }
                for (int q = 0; q < q_; ++q) {
                    const int dx = x + lattice.cx[q];
                    const int dy = y + lattice.cy[q];
                    const int dz = z + lattice.cz[q];
                    if (is_outside(dx, dy, dz) ||
                        solid_[scalar_index(dx, dy, dz)]) {
                        const std::size_t target =
                            dist_index(x, y, z, lattice.opposite[q]);
                        fa_next_[target] += fa_[dist_index(x, y, z, q)];
                        fb_next_[target] += fb_[dist_index(x, y, z, q)];
                    } else {
                        const std::size_t target = dist_index(dx, dy, dz, q);
                        fa_next_[target] += fa_[dist_index(x, y, z, q)];
                        fb_next_[target] += fb_[dist_index(x, y, z, q)];
                    }
                }
            }
        }
    }
    fa_.swap(fa_next_);
    fb_.swap(fb_next_);
}

void TwoPhaseSolver3D::apply_geometry_boundaries() {
    if (!geometry_initialized_) {
        return;
    }
    for (int z = 0; z < nz_; ++z) {
        for (int y = 0; y < ny_; ++y) {
            for (int x = 0; x < nx_; ++x) {
                const GeometryCell cell =
                    geometry_cells_[scalar_index(x, y, z)];
                if (cell == GeometryCell::Inlet) {
                    set_equilibrium_cell(
                        x,
                        y,
                        z,
                        config_.rho_high,
                        config_.rho_low,
                        config_.inlet_velocity,
                        0.0,
                        0.0);
                } else if (cell == GeometryCell::Outlet) {
                    for (int q = 0; q < q_; ++q) {
                        fa_[dist_index(x, y, z, q)] =
                            fa_[dist_index(nx_ - 2, y, z, q)];
                        fb_[dist_index(x, y, z, q)] =
                            fb_[dist_index(nx_ - 2, y, z, q)];
                    }
                }
            }
        }
    }
}

double TwoPhaseSolver3D::permeability_from_porosity(double porosity) const {
    const double eps = std::clamp(porosity, 0.05, 1.0);
    const double solid_fraction = std::max(1.0 - eps, 1.0e-6);
    const double diameter = std::max(config_.porous_pore_diameter, 1.0);
    return diameter * diameter * eps * eps * eps /
           (180.0 * solid_fraction * solid_fraction);
}

bool TwoPhaseSolver3D::is_outside(int x, int y, int z) const {
    return x < 0 || x >= nx_ || y < 0 || y >= ny_ || z < 0 || z >= nz_;
}

} // namespace lbm
