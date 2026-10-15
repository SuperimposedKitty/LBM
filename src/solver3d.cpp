#include "lbm/solver3d.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace lbm {
namespace {

constexpr double pi = 3.141592653589793238462643383279502884;

} // namespace

Solver3D::Solver3D(int nx, int ny, int nz, Solver3DConfig config)
    : grid_(nx, ny, nz, config.lattice_model),
      config_(config) {
    detail::validate_flow_options(config_.flow, false);
    detail::validate_collision_parameters3d(
        config_.collision_model, config_.tau, config_.mrt);
    if (!std::isfinite(config_.initial_rho) || config_.initial_rho <= 0.0) {
        throw std::invalid_argument("Solver3D initial density must be positive and finite.");
    }
    for (int z = 0; z < grid_.nz; ++z) {
        for (int y = 0; y < grid_.ny; ++y) {
            for (int x = 0; x < grid_.nx; ++x) {
                set_equilibrium_cell(x, y, z, config_.initial_rho, 0.0, 0.0, 0.0);
            }
        }
    }
}

Grid3D& Solver3D::grid() {
    return grid_;
}

const Grid3D& Solver3D::grid() const {
    return grid_;
}

void Solver3D::initialize_shear_wave(double amplitude, int mode) {
    if (!std::isfinite(amplitude) || mode <= 0) {
        throw std::invalid_argument("Shear-wave amplitude must be finite and mode positive.");
    }
    masked_flow_initialized_ = false;
    geometry_cells_.clear();
    wall_nodes_.clear();
    origin_={}; spacing_={1,1,1}; units_={}; boundary_diagnostics_={};
    std::fill(grid_.solid.begin(), grid_.solid.end(), std::uint8_t{0});
    for (int z = 0; z < grid_.nz; ++z) {
        for (int y = 0; y < grid_.ny; ++y) {
            const double phase =
                2.0 * pi * mode * static_cast<double>(y) / static_cast<double>(grid_.ny);
            for (int x = 0; x < grid_.nx; ++x) {
                set_equilibrium_cell(
                    x, y, z, config_.initial_rho, amplitude * std::sin(phase), 0.0, 0.0);
            }
        }
    }
}

void Solver3D::initialize_lid_driven_cavity(double lid_velocity) {
    if (!std::isfinite(lid_velocity)) {
        throw std::invalid_argument("Lid velocity must be finite.");
    }
    masked_flow_initialized_ = false;
    geometry_cells_.clear();
    wall_nodes_.clear();
    origin_={}; spacing_={1,1,1}; units_={}; boundary_diagnostics_={};
    for (int z = 0; z < grid_.nz; ++z) {
        for (int y = 0; y < grid_.ny; ++y) {
            for (int x = 0; x < grid_.nx; ++x) {
                const bool wall = x == 0 || x == grid_.nx - 1 || y == 0 ||
                                  y == grid_.ny - 1 || z == 0 || z == grid_.nz - 1;
                const bool lid = y == grid_.ny - 1;
                grid_.solid[grid_.scalar_index(x, y, z)] = wall ? 1 : 0;
                set_equilibrium_cell(
                    x,
                    y,
                    z,
                    config_.initial_rho,
                    lid ? lid_velocity : 0.0,
                    0.0,
                    0.0);
            }
        }
    }
}

void Solver3D::initialize_masked_flow(
    const GeometryMask3D& mask, double inlet_velocity) {
    if (mask.width() != grid_.nx || mask.height() != grid_.ny ||
        mask.depth() != grid_.nz) {
        throw std::invalid_argument("3D geometry dimensions must match the solver grid.");
    }
    if (!std::isfinite(inlet_velocity)) {
        throw std::invalid_argument("Inlet velocity must be finite.");
    }
    if (!mask.contains(GeometryCell::Inlet) || !mask.contains(GeometryCell::Outlet)) {
        throw std::invalid_argument("3D masked flow requires at least one inlet and outlet.");
    }

    geometry_cells_ = mask.cells();
    origin_={}; spacing_={1,1,1}; units_={}; boundary_diagnostics_={};
    wall_nodes_.assign(grid_.rho.size(),{});
    masked_flow_initialized_ = true;
    for (int z = 0; z < grid_.nz; ++z) {
        for (int y = 0; y < grid_.ny; ++y) {
            for (int x = 0; x < grid_.nx; ++x) {
                const std::size_t s = grid_.scalar_index(x, y, z);
                const GeometryCell cell = geometry_cells_[s];
                if (cell == GeometryCell::Inlet &&
                    geometry_cells_[grid_.scalar_index(1, y, z)] == GeometryCell::Solid) {
                    throw std::invalid_argument(
                        "Each 3D inlet must connect to a fluid cell on its right.");
                }
                if (cell == GeometryCell::Outlet &&
                    geometry_cells_[grid_.scalar_index(grid_.nx - 2, y, z)] ==
                        GeometryCell::Solid) {
                    throw std::invalid_argument(
                        "Each 3D outlet must connect to a fluid cell on its left.");
                }
                grid_.solid[s] = cell == GeometryCell::Solid ? 1 : 0;
                wall_nodes_[s].obstacle = grid_.solid[s] && x>0 && x<grid_.nx-1 &&
                    y>0 && y<grid_.ny-1 && z>0 && z<grid_.nz-1;
                set_equilibrium_cell(
                    x,
                    y,
                    z,
                    config_.initial_rho,
                    cell == GeometryCell::Inlet ? inlet_velocity : 0.0,
                    0.0,
                    0.0);
            }
        }
    }
    apply_masked_boundaries(inlet_velocity);
    compute_macroscopic(true);
}

void Solver3D::step() {
    compute_macroscopic(false);
    collide(false);
    stream_periodic();
    compute_macroscopic(false);
}

void Solver3D::step_lid_driven_cavity(double lid_velocity) {
    compute_macroscopic(true, lid_velocity);
    collide(true);
    stream_lid_driven_cavity(lid_velocity);
    compute_macroscopic(true, lid_velocity);
}

void Solver3D::step_masked_flow(double inlet_velocity) {
    if (!masked_flow_initialized_) {
        throw std::logic_error(
            "initialize_masked_flow must be called before step_masked_flow.");
    }
    compute_macroscopic(true);
    collide(true);
    stream_masked_flow();
    apply_masked_boundaries(inlet_velocity);
    compute_macroscopic(true);
}

void Solver3D::run(int steps) {
    if (steps < 0) {
        throw std::invalid_argument("Steps must be non-negative.");
    }
    for (int i = 0; i < steps; ++i) {
        step();
    }
}

void Solver3D::run_lid_driven_cavity(int steps, double lid_velocity) {
    if (steps < 0) {
        throw std::invalid_argument("Steps must be non-negative.");
    }
    for (int i = 0; i < steps; ++i) {
        step_lid_driven_cavity(lid_velocity);
    }
}

Diagnostics3D Solver3D::diagnostics() const {
    Diagnostics3D diagnostics{};
    diagnostics.boundary = boundary_diagnostics_;
    for (std::size_t s = 0; s < grid_.rho.size(); ++s) {
        if (grid_.solid[s]) {
            continue;
        }
        diagnostics.mass += grid_.rho[s];
        diagnostics.max_speed =
            std::max(
                diagnostics.max_speed,
                std::sqrt(
                    grid_.ux[s] * grid_.ux[s] + grid_.uy[s] * grid_.uy[s] +
                    grid_.uz[s] * grid_.uz[s]));
    }
    return diagnostics;
}

void Solver3D::write_vtk(const std::string& path) const {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Failed to open 3D VTK output: " + path);
    }
    out << "# vtk DataFile Version 3.0\n";
    out << "Three-dimensional lattice Boltzmann flow\n";
    out << "ASCII\nDATASET STRUCTURED_POINTS\n";
    out << "DIMENSIONS " << grid_.nx << ' ' << grid_.ny << ' ' << grid_.nz << '\n';
    out << "ORIGIN 0 0 0\nSPACING 1 1 1\n";
    out << "POINT_DATA " << grid_.rho.size() << '\n';
    out << std::setprecision(17);
    out << "SCALARS density double 1\nLOOKUP_TABLE default\n";
    for (double value : grid_.rho) {
        out << value << '\n';
    }
    out << "SCALARS solid int 1\nLOOKUP_TABLE default\n";
    for (std::uint8_t value : grid_.solid) {
        out << static_cast<int>(value) << '\n';
    }
    out << "VECTORS velocity double\n";
    for (std::size_t s = 0; s < grid_.rho.size(); ++s) {
        out << grid_.ux[s] << ' ' << grid_.uy[s] << ' ' << grid_.uz[s] << '\n';
    }
}

void Solver3D::set_equilibrium_cell(
    int x, int y, int z, double rho, double ux, double uy, double uz) {
    const std::size_t s = grid_.scalar_index(x, y, z);
    grid_.rho[s] = rho;
    grid_.ux[s] = ux;
    grid_.uy[s] = uy;
    grid_.uz[s] = uz;
    const auto force = grid_.solid[s] ? std::array<double,3>{} : detail::total_force(config_.flow,rho);
    const auto source = detail::guo_source3d(config_.lattice_model,ux,uy,uz,force[0],force[1],force[2]);
    for (int q = 0; q < grid_.q; ++q) {
        grid_.f[grid_.dist_index(x, y, z, q)] =
            detail::equilibrium3d(config_.lattice_model, q, rho, ux, uy, uz) - 0.5*source[q];
    }
}

void Solver3D::collide(bool fluid_only) {
    for (int z = 0; z < grid_.nz; ++z) {
        for (int y = 0; y < grid_.ny; ++y) {
            for (int x = 0; x < grid_.nx; ++x) {
                const std::size_t s = grid_.scalar_index(x, y, z);
                if (fluid_only && grid_.solid[s]) {
                    continue;
                }
                detail::D3Population population{};
                for (int q = 0; q < grid_.q; ++q) {
                    population[q] = grid_.f[grid_.dist_index(x, y, z, q)];
                }
                const auto force = detail::total_force(config_.flow,grid_.rho[s]);
                const auto source = detail::guo_source3d(config_.lattice_model,
                    grid_.ux[s],grid_.uy[s],grid_.uz[s],force[0],force[1],force[2]);
                detail::collide_population3d(
                    config_.lattice_model,
                    population,
                    grid_.rho[s],
                    grid_.ux[s],
                    grid_.uy[s],
                    grid_.uz[s],
                    config_.tau,
                    config_.collision_model,
                    config_.mrt, source);
                for (int q = 0; q < grid_.q; ++q) {
                    grid_.f[grid_.dist_index(x, y, z, q)] = population[q];
                }
            }
        }
    }
}

void Solver3D::stream_periodic() {
    const auto& lattice = lattice3d(config_.lattice_model);
    std::fill(grid_.f_next.begin(), grid_.f_next.end(), 0.0);
    for (int z = 0; z < grid_.nz; ++z) {
        for (int y = 0; y < grid_.ny; ++y) {
            for (int x = 0; x < grid_.nx; ++x) {
                for (int q = 0; q < grid_.q; ++q) {
                    const int dx = grid_.wrap_x(x + lattice.cx[q]);
                    const int dy = grid_.wrap_y(y + lattice.cy[q]);
                    const int dz = grid_.wrap_z(z + lattice.cz[q]);
                    grid_.f_next[grid_.dist_index(dx, dy, dz, q)] =
                        grid_.f[grid_.dist_index(x, y, z, q)];
                }
            }
        }
    }
    grid_.f.swap(grid_.f_next);
}

void Solver3D::stream_lid_driven_cavity(double lid_velocity) {
    const auto& lattice = lattice3d(config_.lattice_model);
    std::fill(grid_.f_next.begin(), grid_.f_next.end(), 0.0);
    for (int z = 0; z < grid_.nz; ++z) {
        for (int y = 0; y < grid_.ny; ++y) {
            for (int x = 0; x < grid_.nx; ++x) {
                const std::size_t s = grid_.scalar_index(x, y, z);
                if (grid_.solid[s]) {
                    continue;
                }
                for (int q = 0; q < grid_.q; ++q) {
                    const int dx = x + lattice.cx[q];
                    const int dy = y + lattice.cy[q];
                    const int dz = z + lattice.cz[q];
                    if (is_outside(dx, dy, dz) ||
                        grid_.solid[grid_.scalar_index(dx, dy, dz)]) {
                        const bool hits_lid = !is_outside(dx, dy, dz) &&
                                              dy == grid_.ny - 1;
                        const double cu_wall =
                            static_cast<double>(lattice.cx[q]) *
                            (hits_lid ? lid_velocity : 0.0);
                        const double correction =
                            2.0 * lattice.w[q] * grid_.rho[s] * cu_wall /
                            Lattice3DDescriptor::cs2;
                        grid_.f_next[grid_.dist_index(
                            x, y, z, lattice.opposite[q])] +=
                            grid_.f[grid_.dist_index(x, y, z, q)] - correction;
                    } else {
                        grid_.f_next[grid_.dist_index(dx, dy, dz, q)] +=
                            grid_.f[grid_.dist_index(x, y, z, q)];
                    }
                }
            }
        }
    }
    grid_.f.swap(grid_.f_next);
}

void Solver3D::stream_masked_flow() {
    const auto& lattice = lattice3d(config_.lattice_model);
    std::fill(grid_.f_next.begin(), grid_.f_next.end(), 0.0);
    boundary_diagnostics_ = {};
    for (int z = 0; z < grid_.nz; ++z) {
        for (int y = 0; y < grid_.ny; ++y) {
            for (int x = 0; x < grid_.nx; ++x) {
                if (grid_.solid[grid_.scalar_index(x, y, z)]) {
                    continue;
                }
                for (int q = 0; q < grid_.q; ++q) {
                    int dx = x + lattice.cx[q];
                    const int dy = y + lattice.cy[q];
                    int dz = z + lattice.cz[q];
                    if (config_.flow.periodic_x) dx = (dx+grid_.nx)%grid_.nx;
                    if (config_.flow.periodic_z) dz = (dz+grid_.nz)%grid_.nz;
                    if (is_outside(dx, dy, dz) ||
                        grid_.solid[grid_.scalar_index(dx, dy, dz)]) {
                        const auto s = grid_.scalar_index(x,y,z);
                        const int opposite = lattice.opposite[q];
                        const double outgoing = grid_.f[grid_.dist_index(x,y,z,q)];
                        double fraction = 0.5, behind = outgoing, correction = 0;
                        const WallNode* wall = !is_outside(dx,dy,dz) && !wall_nodes_.empty()
                            ? &wall_nodes_[grid_.scalar_index(dx,dy,dz)] : nullptr;
                        if (wall) {
                            correction = 6*lattice.w[q]*grid_.rho[s]*(lattice.cx[q]*wall->velocity[0]+
                                lattice.cy[q]*wall->velocity[1]+lattice.cz[q]*wall->velocity[2]);
                            const auto& here = wall_nodes_[s];
                            if (config_.flow.interpolated_walls && here.distance_known &&
                                wall->distance_known && here.distance>0 && wall->distance<=0)
                                fraction = here.distance/(here.distance-wall->distance);
                        }
                        if (fraction < 0.5) {
                            int bx=x-lattice.cx[q], by=y-lattice.cy[q], bz=z-lattice.cz[q];
                            if (config_.flow.periodic_x) bx=(bx+grid_.nx)%grid_.nx;
                            if (config_.flow.periodic_z) bz=(bz+grid_.nz)%grid_.nz;
                            if (!is_outside(bx,by,bz) && !grid_.solid[grid_.scalar_index(bx,by,bz)])
                                behind=grid_.f[grid_.dist_index(bx,by,bz,q)];
                            else {fraction=0.5; ++boundary_diagnostics_.interpolation_fallbacks;}
                        }
                        const double incoming=detail::reflected_population(outgoing,
                            grid_.f[grid_.dist_index(x,y,z,opposite)],behind,fraction,correction);
                        grid_.f_next[grid_.dist_index(x,y,z,opposite)] += incoming;
                        const bool open_face = is_outside(dx,dy,dz) && !geometry_cells_.empty() &&
                            ((dx<0 && geometry_cells_[s]==GeometryCell::Inlet) ||
                             (dx>=grid_.nx && geometry_cells_[s]==GeometryCell::Outlet));
                        if (!open_face) {
                            boundary_diagnostics_.wall_mass_exchange += incoming-outgoing;
                            const int direction[3]={lattice.cx[q],lattice.cy[q],lattice.cz[q]};
                            for (int d=0; d<3; ++d) {
                                const double f=(outgoing+incoming)*direction[d];
                                boundary_diagnostics_.wall_force[d]+=f;
                                if (wall && wall->obstacle) boundary_diagnostics_.obstacle_force[d]+=f;
                            }
                        }
                    } else {
                        grid_.f_next[grid_.dist_index(dx, dy, dz, q)] +=
                            grid_.f[grid_.dist_index(x, y, z, q)];
                    }
                }
            }
        }
    }
    grid_.f.swap(grid_.f_next);
}

void Solver3D::apply_masked_boundaries(double inlet_velocity) {
    const auto& lattice = lattice3d(config_.lattice_model);
    for (int z = 0; z < grid_.nz; ++z) {
        for (int y = 0; y < grid_.ny; ++y) {
            for (int x = 0; x < grid_.nx; ++x) {
                const std::size_t s = grid_.scalar_index(x, y, z);
                const GeometryCell cell = geometry_cells_[s];
                if (cell!=GeometryCell::Inlet && cell!=GeometryCell::Outlet) continue;
                double before=0;
                for (int q=0; q<grid_.q; ++q) before+=grid_.f[grid_.dist_index(x,y,z,q)];
                if (config_.flow.open_boundary==OpenBoundaryModel::VelocityPressure) {
                    detail::D3Population f{};
                    for (int q=0; q<grid_.q; ++q) f[q]=grid_.f[grid_.dist_index(x,y,z,q)];
                    const auto source=grid_.scalar_index(cell==GeometryCell::Inlet?x+1:x-1,y,z);
                    const double scale=wall_nodes_.empty()?1:wall_nodes_[s].inlet_scale;
                    detail::reconstruct_x_boundary(f,lattice,config_.flow,cell==GeometryCell::Inlet,
                        {inlet_velocity*scale,cell==GeometryCell::Inlet?0:grid_.uy[source],
                                              cell==GeometryCell::Inlet?0:grid_.uz[source]});
                    for (int q=0; q<grid_.q; ++q) grid_.f[grid_.dist_index(x,y,z,q)]=f[q];
                } else if (cell == GeometryCell::Inlet) {
                    const int source_x = 1;
                    double rho = 0.0;
                    double jx = 0.0;
                    double jy = 0.0;
                    double jz = 0.0;
                    for (int q = 0; q < grid_.q; ++q) {
                        const double value =
                            grid_.f[grid_.dist_index(source_x, y, z, q)];
                        rho += value;
                        jx += lattice.cx[q] * value;
                        jy += lattice.cy[q] * value;
                        jz += lattice.cz[q] * value;
                    }
                    if (rho <= 0.0 || !std::isfinite(rho)) {
                        throw std::runtime_error(
                            "Non-physical density next to a 3D inlet.");
                    }
                    for (int q = 0; q < grid_.q; ++q) {
                        const double non_equilibrium =
                            grid_.f[grid_.dist_index(source_x, y, z, q)] -
                            detail::equilibrium3d(
                                config_.lattice_model,
                                q,
                                rho,
                                jx / rho,
                                jy / rho,
                                jz / rho);
                        grid_.f[grid_.dist_index(x, y, z, q)] =
                            detail::equilibrium3d(
                                config_.lattice_model,
                                q,
                                config_.initial_rho,
                                inlet_velocity,
                                0.0,
                                0.0) +
                            non_equilibrium;
                    }
                } else if (cell == GeometryCell::Outlet) {
                    for (int q = 0; q < grid_.q; ++q) {
                        grid_.f[grid_.dist_index(x, y, z, q)] =
                            grid_.f[grid_.dist_index(grid_.nx - 2, y, z, q)];
                    }
                }
                double after=0;
                for (int q=0; q<grid_.q; ++q) after+=grid_.f[grid_.dist_index(x,y,z,q)];
                boundary_diagnostics_.mass_exchange+=after-before;
            }
        }
    }
}

void Solver3D::compute_macroscopic(bool fluid_only, double lid_velocity) {
    const auto& lattice = lattice3d(config_.lattice_model);
    for (int z = 0; z < grid_.nz; ++z) {
        for (int y = 0; y < grid_.ny; ++y) {
            for (int x = 0; x < grid_.nx; ++x) {
                const std::size_t s = grid_.scalar_index(x, y, z);
                if (fluid_only && grid_.solid[s]) {
                    grid_.rho[s] = config_.initial_rho;
                    grid_.ux[s] = y == grid_.ny - 1 ? lid_velocity : 0.0;
                    grid_.uy[s] = 0.0;
                    grid_.uz[s] = 0.0;
                    continue;
                }
                double rho = 0.0;
                double jx = 0.0;
                double jy = 0.0;
                double jz = 0.0;
                for (int q = 0; q < grid_.q; ++q) {
                    const double value = grid_.f[grid_.dist_index(x, y, z, q)];
                    rho += value;
                    jx += lattice.cx[q] * value;
                    jy += lattice.cy[q] * value;
                    jz += lattice.cz[q] * value;
                }
                if (rho <= 0.0 || !std::isfinite(rho)) {
                    throw std::runtime_error(
                        "Non-physical density encountered in Solver3D.");
                }
                grid_.rho[s] = rho;
                const auto force=detail::total_force(config_.flow,rho);
                grid_.ux[s] = (jx+0.5*force[0]) / rho;
                grid_.uy[s] = (jy+0.5*force[1]) / rho;
                grid_.uz[s] = (jz+0.5*force[2]) / rho;
            }
        }
    }
}

bool Solver3D::is_outside(int x, int y, int z) const {
    return x < 0 || x >= grid_.nx || y < 0 || y >= grid_.ny || z < 0 ||
           z >= grid_.nz;
}

FoamSnapshot Solver3D::openfoam_snapshot() const {
    FoamSnapshot snapshot;
    snapshot.nx = grid_.nx;
    snapshot.ny = grid_.ny;
    snapshot.nz = grid_.nz;
    snapshot.two_dimensional = false;
    snapshot.two_phase = false;
    snapshot.origin=origin_; snapshot.spacing=spacing_; snapshot.units=units_;
    snapshot.cells.resize(grid_.rho.size());
    for (std::size_t i = 0; i < snapshot.cells.size(); ++i) {
        auto& cell = snapshot.cells[i];
        cell.solid = grid_.solid[i] != 0;
        cell.rho = grid_.rho[i];
        cell.velocity = {grid_.ux[i], grid_.uy[i], grid_.uz[i]};
        if (wall_nodes_.size()==snapshot.cells.size()) {
            cell.wall_distance=wall_nodes_[i].distance; cell.wall_distance_known=wall_nodes_[i].distance_known;
            cell.obstacle=wall_nodes_[i].obstacle; cell.inlet_scale=wall_nodes_[i].inlet_scale;
            if (cell.solid) cell.velocity=wall_nodes_[i].velocity;
        }
        if (geometry_cells_.size() == snapshot.cells.size()) {
            cell.boundary = geometry_cells_[i] == GeometryCell::Inlet ? 1 :
                            (geometry_cells_[i] == GeometryCell::Outlet ? 2 : 0);
        }
    }
    return snapshot;
}

void Solver3D::initialize_fields(const FoamSnapshot& input) {
    validate_initial_snapshot(input);
    if (input.nx != grid_.nx || input.ny != grid_.ny || input.nz != grid_.nz ||
        input.two_phase != false || input.two_dimensional != false ||
        input.cells.size() != grid_.rho.size())
        throw std::invalid_argument("Imported fields do not match solver dimensions/model.");
    geometry_cells_.resize(input.cells.size());
    origin_=input.origin; spacing_=input.spacing; units_=input.units; boundary_diagnostics_={};
    wall_nodes_.resize(input.cells.size());
    masked_flow_initialized_ = true;
    for (int z = 0; z < grid_.nz; ++z) {
        for (int y = 0; y < grid_.ny; ++y) {
            for (int x = 0; x < grid_.nx; ++x) {
                const auto i = grid_.scalar_index(x, y, z);
                const auto& cell = input.cells[i];
                grid_.solid[i] = cell.solid ? 1 : 0;
                wall_nodes_[i] = {cell.wall_distance,cell.wall_distance_known,cell.obstacle,
                                  cell.solid?cell.velocity:std::array<double,3>{},cell.inlet_scale};
                geometry_cells_[i] = cell.solid ? GeometryCell::Solid :
                    (cell.boundary == 1 ? GeometryCell::Inlet :
                     (cell.boundary == 2 ? GeometryCell::Outlet : GeometryCell::Fluid));
                set_equilibrium_cell(x, y, z, cell.rho, cell.velocity[0], cell.velocity[1], cell.velocity[2]);
            }
        }
    }
}

} // namespace lbm
