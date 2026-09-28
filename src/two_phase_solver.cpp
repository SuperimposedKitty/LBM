#include "lbm/two_phase_solver.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace lbm {
namespace {

constexpr double min_density = 1.0e-12;
constexpr double pi = 3.14159265358979323846;

double local_density(const std::vector<double>& f, int base) {
    double rho = 0.0;
    for (int q = 0; q < D2Q9::q; ++q) {
        rho += f[base + q];
    }
    return rho;
}

} // namespace

TwoPhaseSolver::TwoPhaseSolver(int nx, int ny, TwoPhaseConfig config)
    : nx_(nx),
      ny_(ny),
      config_(config),
      fa_(static_cast<std::size_t>(nx) * ny * D2Q9::q, 0.0),
      fb_(fa_.size(), 0.0),
      fa_next_(fa_.size(), 0.0),
      fb_next_(fa_.size(), 0.0),
      rho_a_(static_cast<std::size_t>(nx) * ny, config.rho_low),
      rho_b_(rho_a_.size(), config.rho_low),
      rho_(rho_a_.size(), config.rho_high + config.rho_low),
      ux_(rho_a_.size(), 0.0),
      uy_(rho_a_.size(), 0.0),
      force_ax_(rho_a_.size(), 0.0),
      force_ay_(rho_a_.size(), 0.0),
      force_bx_(rho_a_.size(), 0.0),
      force_by_(rho_a_.size(), 0.0),
      porosity_(rho_a_.size(), config.free_flow_porosity),
      solid_(rho_a_.size(), 0) {
    if (nx <= 4 || ny <= 4) {
        throw std::invalid_argument("Two-phase grid dimensions must be greater than 4.");
    }
    detail::validate_collision_parameters(
        config_.collision_model, config_.tau_a, config_.mrt);
    detail::validate_collision_parameters(
        config_.collision_model, config_.tau_b, config_.mrt);
    if (config_.rho_high <= 0.0 || config_.rho_low <= 0.0) {
        throw std::invalid_argument("rho_high and rho_low must be positive.");
    }
    if (config_.wall_adhesion_strength < 0.0) {
        throw std::invalid_argument("wall_adhesion_strength must be non-negative.");
    }
    if (config_.bottom_wall_thickness < 1) {
        throw std::invalid_argument("bottom_wall_thickness must be at least 1.");
    }
    if (config_.droplet_interface_width <= 0.0) {
        throw std::invalid_argument("droplet_interface_width must be positive.");
    }
    if (config_.recoloring_strength < 0.0 || config_.recoloring_strength > 1.0) {
        throw std::invalid_argument("recoloring_strength must be in the range [0, 1].");
    }
    if (config_.bottom_wall_repulsion_strength < 0.0) {
        throw std::invalid_argument("bottom_wall_repulsion_strength must be non-negative.");
    }
    if (config_.bottom_wall_repulsion_range < 0) {
        throw std::invalid_argument("bottom_wall_repulsion_range must be non-negative.");
    }
    if (config_.free_flow_porosity <= 0.0 || config_.free_flow_porosity > 1.0 ||
        config_.porous_porosity <= 0.0 || config_.porous_porosity > 1.0) {
        throw std::invalid_argument("porosity values must be in the range (0, 1].");
    }
    if (config_.porous_pore_diameter <= 0.0) {
        throw std::invalid_argument("porous_pore_diameter must be positive.");
    }
    if (config_.darcy_drag_scale < 0.0 || config_.forchheimer_drag_scale < 0.0) {
        throw std::invalid_argument("porous drag scales must be non-negative.");
    }
    if (config_.porous_start_x < 0 && config_.porous_end_x < 0) {
        config_.porous_start_x = 0;
        config_.porous_end_x = 0;
    } else if (config_.porous_start_x < 0 || config_.porous_end_x <= config_.porous_start_x ||
               config_.porous_end_x > nx_) {
        throw std::invalid_argument("porous x range must be empty or inside the domain.");
    }
}

int TwoPhaseSolver::nx() const {
    return nx_;
}

int TwoPhaseSolver::ny() const {
    return ny_;
}

double TwoPhaseSolver::phase_at(int x, int y) const {
    if (is_outside(x, y)) {
        throw std::out_of_range("phase_at coordinates are outside the grid.");
    }
    const int s = scalar_index(x, y);
    // phi 接近 1 表示红色注入相占优，接近 -1 表示蓝色被驱替相占优。
    return (rho_a_[s] - rho_b_[s]) / std::max(rho_[s], min_density);
}

double TwoPhaseSolver::porosity_at(int x, int y) const {
    if (is_outside(x, y)) {
        throw std::out_of_range("porosity_at coordinates are outside the grid.");
    }
    return porosity_[scalar_index(x, y)];
}

bool TwoPhaseSolver::solid_at(int x, int y) const {
    if (is_outside(x, y)) {
        throw std::out_of_range("solid_at coordinates are outside the grid.");
    }
    return solid_[scalar_index(x, y)] != 0;
}

void TwoPhaseSolver::initialize_capillary_displacement() {
    geometry_displacement_initialized_ = false;
    geometry_cells_.clear();
    for (int y = 0; y < ny_; ++y) {
        for (int x = 0; x < nx_; ++x) {
            const int s = scalar_index(x, y);
            // 细管上下边界视为固体壁面，内部初始全部填充蓝色被驱替相。
            solid_[s] = y == 0 || y == ny_ - 1 ? 1 : 0;
            // 中间多孔介质不作为固体堵塞，而是用孔隙度和阻力项表示亚网格孔道。
            porosity_[s] = solid_[s] ? 0.0 :
                (is_porous_column(x) ? config_.porous_porosity : config_.free_flow_porosity);

            const double rho_a = config_.rho_low;
            const double rho_b = solid_[s] ? config_.rho_low : config_.rho_high;
            set_equilibrium_cell(x, y, rho_a, rho_b, 0.0, 0.0);
        }
    }
    apply_inlet_outlet();
    compute_macroscopic();
}

void TwoPhaseSolver::initialize_geometry_displacement(const GeometryMask& mask) {
    if (mask.width() != nx_ || mask.height() != ny_) {
        throw std::invalid_argument("Geometry dimensions must match the two-phase solver grid.");
    }
    if (!mask.contains(GeometryCell::Inlet)) {
        throw std::invalid_argument("Geometry displacement requires at least one inlet I.");
    }
    if (!mask.contains(GeometryCell::Outlet)) {
        throw std::invalid_argument("Geometry displacement requires at least one outlet O.");
    }

    geometry_cells_ = mask.cells();
    for (int y = 0; y < ny_; ++y) {
        for (int x = 0; x < nx_; ++x) {
            const int s = scalar_index(x, y);
            const GeometryCell cell = geometry_cells_[s];
            if (cell == GeometryCell::Inlet &&
                (nx_ < 2 || geometry_cells_[scalar_index(1, y)] == GeometryCell::Solid)) {
                throw std::invalid_argument("Each inlet I must connect to a fluid cell on its right.");
            }
            if (cell == GeometryCell::Outlet &&
                (nx_ < 2 || geometry_cells_[scalar_index(nx_ - 2, y)] == GeometryCell::Solid)) {
                throw std::invalid_argument("Each outlet O must connect to a fluid cell on its left.");
            }
        }
    }

    geometry_displacement_initialized_ = true;
    for (int y = 0; y < ny_; ++y) {
        for (int x = 0; x < nx_; ++x) {
            const int s = scalar_index(x, y);
            const GeometryCell cell = geometry_cells_[s];
            solid_[s] = cell == GeometryCell::Solid ? 1 : 0;
            porosity_[s] = solid_[s]
                               ? 0.0
                               : (cell == GeometryCell::Porous ? config_.porous_porosity
                                                               : config_.free_flow_porosity);

            const bool inlet = cell == GeometryCell::Inlet;
            const double rho_a = inlet ? config_.rho_high : config_.rho_low;
            const double rho_b = solid_[s] ? config_.rho_low
                                           : (inlet ? config_.rho_low : config_.rho_high);
            set_equilibrium_cell(
                x, y, rho_a, rho_b, inlet ? config_.inlet_velocity : 0.0, 0.0);
        }
    }
    apply_geometry_inlet_outlet();
    compute_macroscopic();
}

void TwoPhaseSolver::initialize_droplet_impact(
    double center_x, double center_y, double radius, double initial_ux, double initial_uy) {
    if (radius <= 0.0) {
        throw std::invalid_argument("droplet radius must be positive.");
    }
    if (config_.bottom_wall_thickness >= ny_ - 2) {
        throw std::invalid_argument("bottom_wall_thickness leaves no room for fluid.");
    }

    geometry_displacement_initialized_ = false;
    geometry_cells_.clear();

    const int bottom_wall = config_.bottom_wall_thickness;
    const double interface_width = config_.droplet_interface_width;
    for (int y = 0; y < ny_; ++y) {
        for (int x = 0; x < nx_; ++x) {
            const int s = scalar_index(x, y);
            // 液滴撞击案例采用封闭计算域：底面是被撞击固体表面，其他外边界也反弹。
            solid_[s] = x == 0 || x == nx_ - 1 || y < bottom_wall || y == ny_ - 1 ? 1 : 0;
            porosity_[s] = solid_[s] ? 0.0 : config_.free_flow_porosity;

            const double dx = static_cast<double>(x) - center_x;
            const double dy = static_cast<double>(y) - center_y;
            const double signed_distance = radius - std::sqrt(dx * dx + dy * dy);
            const double droplet_fraction =
                solid_[s] ? 0.0 : 0.5 * (1.0 + std::tanh(signed_distance / interface_width));
            const double rho_a =
                config_.rho_low + (config_.rho_high - config_.rho_low) * droplet_fraction;
            const double rho_b =
                config_.rho_high - (config_.rho_high - config_.rho_low) * droplet_fraction;
            set_equilibrium_cell(
                x, y, rho_a, rho_b, initial_ux * droplet_fraction, initial_uy * droplet_fraction);
        }
    }
    compute_macroscopic();
}

void TwoPhaseSolver::step() {
    // 两相更新：求宏观量 -> 计算相互作用力 -> 带力碰撞 -> 迁移 -> 入口/出口边界。
    compute_macroscopic();
    compute_forces();
    collide();
    recolor();
    stream();
    apply_inlet_outlet();
    compute_macroscopic();
}

void TwoPhaseSolver::step_closed() {
    // 封闭域更新不施加入口/出口，只依赖固壁反弹和体力演化。
    compute_macroscopic();
    compute_forces();
    collide();
    recolor();
    stream();
    compute_macroscopic();
}

void TwoPhaseSolver::run(int steps) {
    if (steps < 0) {
        throw std::invalid_argument("steps must be non-negative.");
    }
    for (int i = 0; i < steps; ++i) {
        step();
    }
    compute_macroscopic();
}

void TwoPhaseSolver::run_closed(int steps) {
    if (steps < 0) {
        throw std::invalid_argument("steps must be non-negative.");
    }
    for (int i = 0; i < steps; ++i) {
        step_closed();
    }
    compute_macroscopic();
}

TwoPhaseDiagnostics TwoPhaseSolver::diagnostics() const {
    TwoPhaseDiagnostics d{};
    double weighted_interface = 0.0;
    double weighted_centroid_y = 0.0;
    double injected_mass = 0.0;
    double porous_pore_speed_sum = 0.0;
    int porous_samples = 0;

    for (int y = 0; y < ny_; ++y) {
        for (int x = 0; x < nx_; ++x) {
            const int s = scalar_index(x, y);
            if (solid_[s]) {
                continue;
            }
            d.mass_a += rho_a_[s];
            d.mass_b += rho_b_[s];
            const double speed = std::sqrt(ux_[s] * ux_[s] + uy_[s] * uy_[s]);
            d.max_speed = std::max(d.max_speed, speed);
            if (porosity_[s] > 0.0 && porosity_[s] < config_.free_flow_porosity) {
                const double pore_speed = speed / porosity_[s];
                porous_pore_speed_sum += pore_speed;
                d.porous_max_pore_speed = std::max(d.porous_max_pore_speed, pore_speed);
                ++porous_samples;
            }

            const double injected_saturation = std::clamp(
                (rho_a_[s] - config_.rho_low) /
                    std::max(config_.rho_high - config_.rho_low, min_density),
                0.0,
                1.0);
            weighted_interface += static_cast<double>(x) * injected_saturation;
            weighted_centroid_y += static_cast<double>(y) * injected_saturation;
            injected_mass += injected_saturation;
        }
    }

    d.interface_x = injected_mass > 0.0 ? weighted_interface / injected_mass : 0.0;
    d.phase_a_centroid_y =
        injected_mass > 0.0 ? weighted_centroid_y / injected_mass : 0.0;
    d.porous_mean_pore_speed =
        porous_samples > 0 ? porous_pore_speed_sum / static_cast<double>(porous_samples) : 0.0;
    d.interaction_force_balance_x = interaction_force_balance_x_;
    d.interaction_force_balance_y = interaction_force_balance_y_;
    return d;
}

void TwoPhaseSolver::write_csv(const std::string& path) const {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Failed to open output file: " + path);
    }

    out << "x,y,rho_a,rho_b,rho,phi,ux,uy,porosity,solid\n";
    out << std::setprecision(17);
    for (int y = 0; y < ny_; ++y) {
        for (int x = 0; x < nx_; ++x) {
            const int s = scalar_index(x, y);
            const double phi = (rho_a_[s] - rho_b_[s]) / std::max(rho_[s], min_density);
            out << x << ',' << y << ',' << rho_a_[s] << ',' << rho_b_[s] << ',' << rho_[s] << ','
                << phi << ',' << ux_[s] << ',' << uy_[s] << ',' << porosity_[s] << ','
                << static_cast<int>(solid_[s]) << '\n';
        }
    }
}

void TwoPhaseSolver::write_vtk(const std::string& path) const {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Failed to open output file: " + path);
    }

    out << "# vtk DataFile Version 3.0\n";
    out << "D2Q9 two-phase capillary displacement\n";
    out << "ASCII\n";
    out << "DATASET STRUCTURED_POINTS\n";
    out << "DIMENSIONS " << nx_ << ' ' << ny_ << " 1\n";
    out << "ORIGIN 0 0 0\n";
    out << "SPACING 1 1 1\n";
    out << "POINT_DATA " << nx_ * ny_ << '\n';
    out << std::setprecision(17);

    out << "SCALARS rho_a double 1\n";
    out << "LOOKUP_TABLE default\n";
    for (double value : rho_a_) {
        out << value << '\n';
    }

    out << "SCALARS rho_b double 1\n";
    out << "LOOKUP_TABLE default\n";
    for (double value : rho_b_) {
        out << value << '\n';
    }

    out << "SCALARS phase double 1\n";
    out << "LOOKUP_TABLE default\n";
    for (std::size_t i = 0; i < rho_.size(); ++i) {
        out << (rho_a_[i] - rho_b_[i]) / std::max(rho_[i], min_density) << '\n';
    }

    out << "SCALARS solid int 1\n";
    out << "LOOKUP_TABLE default\n";
    for (std::uint8_t value : solid_) {
        out << static_cast<int>(value) << '\n';
    }

    out << "SCALARS porosity double 1\n";
    out << "LOOKUP_TABLE default\n";
    for (double value : porosity_) {
        out << value << '\n';
    }

    out << "VECTORS velocity double\n";
    for (std::size_t i = 0; i < rho_.size(); ++i) {
        out << ux_[i] << ' ' << uy_[i] << " 0\n";
    }
}

int TwoPhaseSolver::scalar_index(int x, int y) const {
    return y * nx_ + x;
}

int TwoPhaseSolver::dist_index(int x, int y, int direction) const {
    return (y * nx_ + x) * D2Q9::q + direction;
}

double TwoPhaseSolver::psi(double rho) {
    // Shan-Chen 伪势函数，低密度近似线性，高密度趋于饱和。
    return 1.0 - std::exp(-std::max(rho, 0.0));
}

double TwoPhaseSolver::equilibrium(int direction, double rho, double ux, double uy) {
    // 两个组分共用混合速度，但各自用本组分密度构造平衡分布。
    return detail::equilibrium(direction, rho, ux, uy);
}

void TwoPhaseSolver::compute_macroscopic() {
    for (int y = 0; y < ny_; ++y) {
        for (int x = 0; x < nx_; ++x) {
            const int s = scalar_index(x, y);
            if (solid_[s]) {
                rho_a_[s] = config_.rho_low;
                rho_b_[s] = config_.rho_low;
                rho_[s] = rho_a_[s] + rho_b_[s];
                ux_[s] = 0.0;
                uy_[s] = 0.0;
                continue;
            }

            // 两个组分密度分别由各自分布函数求和；动量按混合物总动量求取。
            double rho_a = 0.0;
            double rho_b = 0.0;
            double momentum_x = 0.0;
            double momentum_y = 0.0;
            for (int q = 0; q < D2Q9::q; ++q) {
                const double fa = fa_[dist_index(x, y, q)];
                const double fb = fb_[dist_index(x, y, q)];
                rho_a += fa;
                rho_b += fb;
                momentum_x += static_cast<double>(D2Q9::cx[q]) * (fa + fb);
                momentum_y += static_cast<double>(D2Q9::cy[q]) * (fa + fb);
            }

            if (rho_a <= 0.0 || rho_b <= 0.0 || !std::isfinite(rho_a) || !std::isfinite(rho_b)) {
                throw std::runtime_error("Non-physical two-phase density encountered.");
            }

            rho_a_[s] = rho_a;
            rho_b_[s] = rho_b;
            rho_[s] = rho_a + rho_b;
            const double total_fx = force_ax_[s] + force_bx_[s];
            const double total_fy = force_ay_[s] + force_by_[s];
            // 半步力修正使速度与带力 LBM 的时间中心一致。
            ux_[s] = (momentum_x + 0.5 * total_fx) / rho_[s];
            uy_[s] = (momentum_y + 0.5 * total_fy) / rho_[s];
        }
    }
}

void TwoPhaseSolver::compute_forces() {
    std::fill(force_ax_.begin(), force_ax_.end(), 0.0);
    std::fill(force_ay_.begin(), force_ay_.end(), 0.0);
    std::fill(force_bx_.begin(), force_bx_.end(), 0.0);
    std::fill(force_by_.begin(), force_by_.end(), 0.0);

    // contact_angle_degrees 通过 cos(theta) 转换为壁面对两相的相反偏好。
    const double contact_angle = config_.contact_angle_degrees * pi / 180.0;
    const double wetting_bias = config_.wall_adhesion_strength * std::cos(contact_angle);

    // 每条无序链路只访问一次，并把两组分之间的作用力等大反向地累加到两个端点。
    constexpr std::array<int, 4> half_directions{1, 2, 5, 6};
    for (int y = 0; y < ny_; ++y) {
        for (int x = 0; x < nx_; ++x) {
            const int s = scalar_index(x, y);
            if (solid_[s]) {
                continue;
            }
            for (int q : half_directions) {
                const int neighbor_x = x + D2Q9::cx[q];
                const int neighbor_y = y + D2Q9::cy[q];
                if (is_outside(neighbor_x, neighbor_y)) {
                    continue;
                }
                const int neighbor = scalar_index(neighbor_x, neighbor_y);
                if (solid_[neighbor]) {
                    continue;
                }

                const double cx = static_cast<double>(D2Q9::cx[q]);
                const double cy = static_cast<double>(D2Q9::cy[q]);
                const double weight = D2Q9::w[q];
                const double a_to_b =
                    -config_.interaction_strength * weight * psi(rho_a_[s]) *
                    psi(rho_b_[neighbor]);
                force_ax_[s] += a_to_b * cx;
                force_ay_[s] += a_to_b * cy;
                force_bx_[neighbor] -= a_to_b * cx;
                force_by_[neighbor] -= a_to_b * cy;

                const double b_to_a =
                    -config_.interaction_strength * weight * psi(rho_b_[s]) *
                    psi(rho_a_[neighbor]);
                force_bx_[s] += b_to_a * cx;
                force_by_[s] += b_to_a * cy;
                force_ax_[neighbor] -= b_to_a * cx;
                force_ay_[neighbor] -= b_to_a * cy;
            }
        }
    }

    interaction_force_balance_x_ = 0.0;
    interaction_force_balance_y_ = 0.0;
    for (std::size_t s = 0; s < rho_.size(); ++s) {
        interaction_force_balance_x_ += force_ax_[s] + force_bx_[s];
        interaction_force_balance_y_ += force_ay_[s] + force_by_[s];
    }

    // 开放边界外侧使用零法向梯度虚拟储液层，补齐缺失的伪势邻居并保持局部各向同性。
    for (int y = 0; y < ny_; ++y) {
        for (int x = 0; x < nx_; ++x) {
            const int s = scalar_index(x, y);
            if (solid_[s]) {
                continue;
            }
            for (int q = 1; q < D2Q9::q; ++q) {
                const int neighbor_x = x + D2Q9::cx[q];
                const int neighbor_y = y + D2Q9::cy[q];
                if (!is_outside(neighbor_x, neighbor_y)) {
                    continue;
                }
                if (neighbor_y < 0 || neighbor_y >= ny_) {
                    continue;
                }
                const int ghost_x = std::clamp(neighbor_x, 0, nx_ - 1);
                const int ghost = scalar_index(ghost_x, neighbor_y);
                if (solid_[ghost]) {
                    continue;
                }
                const double cx = static_cast<double>(D2Q9::cx[q]);
                const double cy = static_cast<double>(D2Q9::cy[q]);
                const double weight = D2Q9::w[q];
                const double force_a =
                    -config_.interaction_strength * weight * psi(rho_a_[s]) *
                    psi(rho_b_[ghost]);
                const double force_b =
                    -config_.interaction_strength * weight * psi(rho_b_[s]) *
                    psi(rho_a_[ghost]);
                force_ax_[s] += force_a * cx;
                force_ay_[s] += force_a * cy;
                force_bx_[s] += force_b * cx;
                force_by_[s] += force_b * cy;
            }
        }
    }

    for (int y = 0; y < ny_; ++y) {
        for (int x = 0; x < nx_; ++x) {
            const int s = scalar_index(x, y);
            if (solid_[s]) {
                continue;
            }

            const double psi_a = psi(rho_a_[s]);
            const double psi_b = psi(rho_b_[s]);

            // 流体-固体黏附力：只统计邻近固体格点，用于体现接触角作用。
            double wall_x = 0.0;
            double wall_y = 0.0;
            for (int q = 1; q < D2Q9::q; ++q) {
                const int neighbor_x = x + D2Q9::cx[q];
                const int neighbor_y = y + D2Q9::cy[q];
                if (is_outside(neighbor_x, neighbor_y)) {
                    continue;
                }
                if (!solid_[scalar_index(neighbor_x, neighbor_y)]) {
                    continue;
                }
                wall_x += D2Q9::w[q] * static_cast<double>(D2Q9::cx[q]);
                wall_y += D2Q9::w[q] * static_cast<double>(D2Q9::cy[q]);
            }

            force_ax_[s] += wetting_bias * psi_a * wall_x;
            force_ay_[s] += wetting_bias * psi_a * wall_y;
            force_bx_[s] -= wetting_bias * psi_b * wall_x;
            force_by_[s] -= wetting_bias * psi_b * wall_y;

            const double total_rho = std::max(rho_[s], min_density);
            if (config_.bottom_wall_repulsion_strength > 0.0 &&
                config_.bottom_wall_repulsion_range > 0) {
                const int distance_from_surface = y - config_.bottom_wall_thickness;
                if (distance_from_surface >= 0 &&
                    distance_from_surface < config_.bottom_wall_repulsion_range) {
                    const double normalized_distance =
                        static_cast<double>(distance_from_surface) /
                        static_cast<double>(config_.bottom_wall_repulsion_range);
                    const double falloff = (1.0 - normalized_distance) *
                                           (1.0 - normalized_distance);
                    const double saturation_a = rho_a_[s] / total_rho;
                    if (saturation_a > 0.08) {
                        const double repulsion =
                            config_.bottom_wall_repulsion_strength * falloff *
                            saturation_a * psi_a;
                        force_ay_[s] += repulsion;
                    }
                }
            }
            force_ax_[s] += config_.body_force_x * rho_a_[s] / total_rho;
            force_bx_[s] += config_.body_force_x * rho_b_[s] / total_rho;
            force_ay_[s] += config_.body_force_y * rho_a_[s] / total_rho;
            force_by_[s] += config_.body_force_y * rho_b_[s] / total_rho;

            if (porosity_[s] > 0.0 && porosity_[s] < config_.free_flow_porosity) {
                const double permeability = permeability_from_porosity(porosity_[s]);
                const double viscosity =
                    D2Q9::cs2 * (0.5 * (config_.tau_a + config_.tau_b) - 0.5);
                const double pore_ux = ux_[s] / porosity_[s];
                const double pore_uy = uy_[s] / porosity_[s];
                const double pore_speed = std::sqrt(pore_ux * pore_ux + pore_uy * pore_uy);
                const double darcy_coeff =
                    config_.darcy_drag_scale * viscosity / std::max(permeability, min_density);
                const double forchheimer_coeff =
                    config_.forchheimer_drag_scale / std::sqrt(std::max(permeability, min_density));
                const double drag_x =
                    -rho_[s] * (darcy_coeff + forchheimer_coeff * pore_speed) * pore_ux;
                const double drag_y =
                    -rho_[s] * (darcy_coeff + forchheimer_coeff * pore_speed) * pore_uy;
                // 多孔介质项是 REV 尺度的 Darcy/Forchheimer 阻力，再按组分占比分配给 A/B。
                force_ax_[s] += drag_x * rho_a_[s] / total_rho;
                force_ay_[s] += drag_y * rho_a_[s] / total_rho;
                force_bx_[s] += drag_x * rho_b_[s] / total_rho;
                force_by_[s] += drag_y * rho_b_[s] / total_rho;
            }
        }
    }
}

void TwoPhaseSolver::collide() {
    for (int y = 0; y < ny_; ++y) {
        for (int x = 0; x < nx_; ++x) {
            const int s = scalar_index(x, y);
            if (solid_[s]) {
                continue;
            }

            detail::D2Q9Population population_a{};
            detail::D2Q9Population population_b{};
            for (int q = 0; q < D2Q9::q; ++q) {
                population_a[q] = fa_[dist_index(x, y, q)];
                population_b[q] = fb_[dist_index(x, y, q)];
            }
            const detail::D2Q9Population source_a = detail::guo_source(
                ux_[s], uy_[s], force_ax_[s], force_ay_[s]);
            const detail::D2Q9Population source_b = detail::guo_source(
                ux_[s], uy_[s], force_bx_[s], force_by_[s]);
            detail::collide_population(
                population_a,
                rho_a_[s],
                ux_[s],
                uy_[s],
                config_.tau_a,
                config_.collision_model,
                config_.mrt,
                source_a,
                1.0 / config_.tau_a);
            detail::collide_population(
                population_b,
                rho_b_[s],
                ux_[s],
                uy_[s],
                config_.tau_b,
                config_.collision_model,
                config_.mrt,
                source_b,
                1.0 / config_.tau_b);
            for (int q = 0; q < D2Q9::q; ++q) {
                fa_[dist_index(x, y, q)] = population_a[q];
                fb_[dist_index(x, y, q)] = population_b[q];
            }
        }
    }
}

void TwoPhaseSolver::recolor() {
    if (config_.recoloring_strength <= 0.0) {
        return;
    }

    std::vector<double> phase(rho_.size(), 0.0);
    for (int y = 0; y < ny_; ++y) {
        for (int x = 0; x < nx_; ++x) {
            const int s = scalar_index(x, y);
            if (solid_[s]) {
                continue;
            }
            phase[s] = (rho_a_[s] - rho_b_[s]) / std::max(rho_[s], min_density);
        }
    }

    for (int y = 0; y < ny_; ++y) {
        for (int x = 0; x < nx_; ++x) {
            const int s = scalar_index(x, y);
            if (solid_[s]) {
                continue;
            }

            double grad_x = 0.0;
            double grad_y = 0.0;
            for (int q = 1; q < D2Q9::q; ++q) {
                const int neighbor_x = x + D2Q9::cx[q];
                const int neighbor_y = y + D2Q9::cy[q];
                if (is_outside(neighbor_x, neighbor_y)) {
                    continue;
                }
                const int neighbor = scalar_index(neighbor_x, neighbor_y);
                if (solid_[neighbor]) {
                    continue;
                }
                grad_x += D2Q9::w[q] * phase[neighbor] * static_cast<double>(D2Q9::cx[q]);
                grad_y += D2Q9::w[q] * phase[neighbor] * static_cast<double>(D2Q9::cy[q]);
            }

            const double grad_norm = std::sqrt(grad_x * grad_x + grad_y * grad_y);
            if (grad_norm <= 1.0e-12) {
                continue;
            }

            const double total_rho = std::max(rho_[s], min_density);
            const double ratio_a = std::clamp(rho_a_[s] / total_rho, 0.0, 1.0);
            const double ratio_b = 1.0 - ratio_a;
            if (ratio_a < 0.04 || ratio_a > 0.96) {
                continue;
            }
            const double strength = config_.recoloring_strength * ratio_a * ratio_b * total_rho;

            std::array<double, D2Q9::q> total_f{};
            std::array<double, D2Q9::q> recolored_a{};
            double recolored_mass_a = 0.0;
            for (int q = 0; q < D2Q9::q; ++q) {
                const int k = dist_index(x, y, q);
                total_f[q] = fa_[k] + fb_[k];
                if (total_f[q] <= min_density) {
                    continue;
                }

                double directional_alignment = 0.0;
                if (q > 0) {
                    const double dir_norm = std::sqrt(
                        static_cast<double>(D2Q9::cx[q] * D2Q9::cx[q] +
                                            D2Q9::cy[q] * D2Q9::cy[q]));
                    directional_alignment =
                        (static_cast<double>(D2Q9::cx[q]) * grad_x +
                         static_cast<double>(D2Q9::cy[q]) * grad_y) /
                        (dir_norm * grad_norm);
                }

                const double correction = strength * D2Q9::w[q] * directional_alignment;
                recolored_a[q] =
                    std::clamp(ratio_a * total_f[q] + correction, 0.0, total_f[q]);
                recolored_mass_a += recolored_a[q];
            }

            if (recolored_mass_a > min_density) {
                const double target_mass_a = rho_a_[s];
                const double factor = target_mass_a / recolored_mass_a;
                for (int q = 0; q < D2Q9::q; ++q) {
                    recolored_a[q] = std::clamp(recolored_a[q] * factor, 0.0, total_f[q]);
                }

                double corrected_mass_a = 0.0;
                for (double value : recolored_a) {
                    corrected_mass_a += value;
                }
                double residual = target_mass_a - corrected_mass_a;
                if (std::abs(residual) > 1.0e-12) {
                    double capacity = 0.0;
                    for (int q = 0; q < D2Q9::q; ++q) {
                        capacity += residual > 0.0 ? total_f[q] - recolored_a[q] : recolored_a[q];
                    }
                    if (capacity > min_density) {
                        for (int q = 0; q < D2Q9::q; ++q) {
                            const double share =
                                residual > 0.0 ? total_f[q] - recolored_a[q] : recolored_a[q];
                            recolored_a[q] += residual * share / capacity;
                        }
                    }
                }
            }

            for (int q = 0; q < D2Q9::q; ++q) {
                const int k = dist_index(x, y, q);
                fa_[k] = std::clamp(recolored_a[q], 0.0, total_f[q]);
                fb_[k] = total_f[q] - fa_[k];
            }
        }
    }
}

void TwoPhaseSolver::stream() {
    std::fill(fa_next_.begin(), fa_next_.end(), 0.0);
    std::fill(fb_next_.begin(), fb_next_.end(), 0.0);

    for (int y = 0; y < ny_; ++y) {
        for (int x = 0; x < nx_; ++x) {
            const int s = scalar_index(x, y);
            if (solid_[s]) {
                continue;
            }

            for (int q = 0; q < D2Q9::q; ++q) {
                const int dst_x = x + D2Q9::cx[q];
                const int dst_y = y + D2Q9::cy[q];
                if (is_outside(dst_x, dst_y) || solid_[scalar_index(dst_x, dst_y)]) {
                    const int opposite = D2Q9::opposite[q];
                    // 撞到固体壁面或计算域外侧时执行反弹，形成无滑移边界。
                    fa_next_[dist_index(x, y, opposite)] += fa_[dist_index(x, y, q)];
                    fb_next_[dist_index(x, y, opposite)] += fb_[dist_index(x, y, q)];
                    continue;
                }

                fa_next_[dist_index(dst_x, dst_y, q)] += fa_[dist_index(x, y, q)];
                fb_next_[dist_index(dst_x, dst_y, q)] += fb_[dist_index(x, y, q)];
            }
        }
    }

    fa_.swap(fa_next_);
    fb_.swap(fb_next_);
}

void TwoPhaseSolver::apply_inlet_outlet() {
    if (geometry_displacement_initialized_) {
        apply_geometry_inlet_outlet();
        return;
    }

    for (int y = 1; y < ny_ - 1; ++y) {
        // 左边界固定为红色注入相入口，持续把 A 相推入细管。
        set_equilibrium_cell(
            0, y, config_.rho_high, config_.rho_low, config_.inlet_velocity, 0.0);

        // 右边界采用零梯度式出口：复制近出口处的密度和受限速度。
        const int source_base = dist_index(nx_ - 2, y, 0);
        const double rho_a = std::max(local_density(fa_, source_base), config_.rho_low);
        const double rho_b = std::max(local_density(fb_, source_base), config_.rho_low);
        double momentum_x = 0.0;
        double momentum_y = 0.0;
        for (int q = 0; q < D2Q9::q; ++q) {
            const double f = fa_[source_base + q] + fb_[source_base + q];
            momentum_x += static_cast<double>(D2Q9::cx[q]) * f;
            momentum_y += static_cast<double>(D2Q9::cy[q]) * f;
        }
        const double rho = rho_a + rho_b;
        const double ux = std::clamp(momentum_x / std::max(rho, min_density), -0.08, 0.08);
        const double uy = std::clamp(momentum_y / std::max(rho, min_density), -0.08, 0.08);
        set_equilibrium_cell(nx_ - 1, y, rho_a, rho_b, ux, uy);
    }
}

void TwoPhaseSolver::apply_geometry_inlet_outlet() {
    for (int y = 0; y < ny_; ++y) {
        for (int x = 0; x < nx_; ++x) {
            const GeometryCell cell = geometry_cells_[scalar_index(x, y)];
            if (cell == GeometryCell::Inlet) {
                set_equilibrium_cell(
                    x, y, config_.rho_high, config_.rho_low, config_.inlet_velocity, 0.0);
                continue;
            }
            if (cell != GeometryCell::Outlet) {
                continue;
            }

            // 出口复制内侧相邻单元的密度和速度，形成近似零梯度开边界。
            const int source_base = dist_index(x - 1, y, 0);
            const double rho_a = std::max(local_density(fa_, source_base), config_.rho_low);
            const double rho_b = std::max(local_density(fb_, source_base), config_.rho_low);
            double momentum_x = 0.0;
            double momentum_y = 0.0;
            for (int q = 0; q < D2Q9::q; ++q) {
                const double f = fa_[source_base + q] + fb_[source_base + q];
                momentum_x += static_cast<double>(D2Q9::cx[q]) * f;
                momentum_y += static_cast<double>(D2Q9::cy[q]) * f;
            }
            const double rho = rho_a + rho_b;
            const double ux =
                std::clamp(momentum_x / std::max(rho, min_density), -0.08, 0.08);
            const double uy =
                std::clamp(momentum_y / std::max(rho, min_density), -0.08, 0.08);
            set_equilibrium_cell(x, y, rho_a, rho_b, ux, uy);
        }
    }
}

void TwoPhaseSolver::set_equilibrium_cell(
    int x, int y, double rho_a, double rho_b, double ux, double uy) {
    const int s = scalar_index(x, y);
    rho_a_[s] = rho_a;
    rho_b_[s] = rho_b;
    rho_[s] = rho_a + rho_b;
    ux_[s] = solid_[s] ? 0.0 : ux;
    uy_[s] = solid_[s] ? 0.0 : uy;

    for (int q = 0; q < D2Q9::q; ++q) {
        fa_[dist_index(x, y, q)] = equilibrium(q, rho_a_[s], ux_[s], uy_[s]);
        fb_[dist_index(x, y, q)] = equilibrium(q, rho_b_[s], ux_[s], uy_[s]);
    }
}

bool TwoPhaseSolver::is_porous_column(int x) const {
    return x >= config_.porous_start_x && x < config_.porous_end_x;
}

double TwoPhaseSolver::permeability_from_porosity(double porosity) const {
    const double eps = std::clamp(porosity, 0.05, 1.0);
    const double solid_fraction = std::max(1.0 - eps, 1.0e-6);
    const double pore_diameter = std::max(config_.porous_pore_diameter, 1.0);
    // Kozeny-Carman 型等效渗透率，把孔隙尺度参数映射到格点尺度阻力。
    return pore_diameter * pore_diameter * eps * eps * eps /
           (180.0 * solid_fraction * solid_fraction);
}

bool TwoPhaseSolver::is_outside(int x, int y) const {
    return x < 0 || x >= nx_ || y < 0 || y >= ny_;
}

FoamSnapshot TwoPhaseSolver::openfoam_snapshot() const {
    FoamSnapshot snapshot;
    snapshot.nx = nx_;
    snapshot.ny = ny_;
    snapshot.nz = 1;
    snapshot.two_dimensional = true;
    snapshot.two_phase = true;
    snapshot.interaction_strength = config_.interaction_strength;
    snapshot.cells.resize(rho_.size());
    for (std::size_t i = 0; i < snapshot.cells.size(); ++i) {
        auto& cell = snapshot.cells[i];
        cell.solid = solid_[i] != 0;
        cell.rho = rho_[i];
        cell.velocity = {ux_[i], uy_[i], 0.0};
        cell.rho_a = rho_a_[i];
        cell.rho_b = rho_b_[i];
        cell.porosity = porosity_[i];
        if (geometry_cells_.size() == snapshot.cells.size()) {
            cell.boundary = geometry_cells_[i] == GeometryCell::Inlet ? 1 :
                            (geometry_cells_[i] == GeometryCell::Outlet ? 2 : 0);
        } else if (!cell.solid) {
            // 细管两端为开放边界；封闭液滴案例在这两列均为固体。
            const auto x = i % static_cast<std::size_t>(nx_);
            cell.boundary = x == 0 ? 1 : (x == static_cast<std::size_t>(nx_ - 1) ? 2 : 0);
        }
    }
    return snapshot;
}

void TwoPhaseSolver::initialize_fields(const FoamSnapshot& input) {
    validate_initial_snapshot(input);
    if (input.nx != nx_ || input.ny != ny_ || input.nz != 1 ||
        input.two_phase != true || input.two_dimensional != true ||
        input.cells.size() != rho_.size())
        throw std::invalid_argument("Imported fields do not match solver dimensions/model.");
    geometry_cells_.resize(input.cells.size());
    geometry_displacement_initialized_ = true;
    for (int z = 0; z < 1; ++z) {
        for (int y = 0; y < ny_; ++y) {
            for (int x = 0; x < nx_; ++x) {
                const auto i = scalar_index(x, y);
                const auto& cell = input.cells[i];
                solid_[i] = cell.solid ? 1 : 0;
                geometry_cells_[i] = cell.solid ? GeometryCell::Solid :
                    (cell.boundary == 1 ? GeometryCell::Inlet :
                     (cell.boundary == 2 ? GeometryCell::Outlet : GeometryCell::Fluid));
                porosity_[i] = cell.solid ? 0.0 : cell.porosity;
                set_equilibrium_cell(x, y, cell.rho_a, cell.rho_b,
                    cell.velocity[0], cell.velocity[1]);
            }
        }
    }
}

} // namespace lbm
