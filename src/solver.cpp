#include "lbm/solver.hpp"

#include "lbm/lattice.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace lbm {
namespace {

constexpr double pi = 3.141592653589793238462643383279502884;

} // namespace

Solver::Solver(int nx, int ny, SolverConfig config)
    : grid_(nx, ny),
      config_(config) {
    detail::validate_collision_parameters(
        config_.collision_model, config_.tau, config_.mrt);

    for (int y = 0; y < grid_.ny; ++y) {
        for (int x = 0; x < grid_.nx; ++x) {
            const int s = grid_.scalar_index(x, y);
            grid_.rho[s] = config_.initial_rho;
            grid_.ux[s] = 0.0;
            grid_.uy[s] = 0.0;
            for (int q = 0; q < D2Q9::q; ++q) {
                grid_.f[grid_.dist_index(x, y, q)] =
                    equilibrium(q, grid_.rho[s], grid_.ux[s], grid_.uy[s]);
            }
        }
    }
}

Grid& Solver::grid() {
    return grid_;
}

const Grid& Solver::grid() const {
    return grid_;
}

void Solver::initialize_shear_wave(double amplitude, int mode) {
    if (mode <= 0) {
        throw std::invalid_argument("mode must be positive.");
    }

    masked_flow_initialized_ = false;
    geometry_cells_.clear();

    for (int y = 0; y < grid_.ny; ++y) {
        const double phase = 2.0 * pi * mode * static_cast<double>(y) /
                             static_cast<double>(grid_.ny);
        for (int x = 0; x < grid_.nx; ++x) {
            const int s = grid_.scalar_index(x, y);
            grid_.rho[s] = config_.initial_rho;
            grid_.ux[s] = amplitude * std::sin(phase);
            grid_.uy[s] = 0.0;
            for (int q = 0; q < D2Q9::q; ++q) {
                grid_.f[grid_.dist_index(x, y, q)] =
                    equilibrium(q, grid_.rho[s], grid_.ux[s], grid_.uy[s]);
            }
        }
    }
}

void Solver::initialize_lid_driven_cavity(double lid_velocity) {
    masked_flow_initialized_ = false;
    geometry_cells_.clear();
    for (int y = 0; y < grid_.ny; ++y) {
        for (int x = 0; x < grid_.nx; ++x) {
            const bool is_wall = x == 0 || y == 0 || x == grid_.nx - 1 || y == grid_.ny - 1;
            const bool is_lid = y == grid_.ny - 1;
            const int s = grid_.scalar_index(x, y);

            grid_.solid[s] = is_wall ? 1 : 0;
            grid_.rho[s] = config_.initial_rho;
            grid_.ux[s] = is_lid ? lid_velocity : 0.0;
            grid_.uy[s] = 0.0;

            for (int q = 0; q < D2Q9::q; ++q) {
                grid_.f[grid_.dist_index(x, y, q)] =
                    equilibrium(q, grid_.rho[s], grid_.ux[s], grid_.uy[s]);
            }
        }
    }
}

void Solver::initialize_masked_flow(const GeometryMask& mask, double inlet_velocity) {
    if (mask.width() != grid_.nx || mask.height() != grid_.ny) {
        throw std::invalid_argument("Geometry dimensions must match the solver grid.");
    }
    if (!std::isfinite(inlet_velocity)) {
        throw std::invalid_argument("inlet_velocity must be finite.");
    }
    if (!mask.contains(GeometryCell::Inlet)) {
        throw std::invalid_argument("Masked flow geometry must contain at least one inlet I.");
    }
    if (!mask.contains(GeometryCell::Outlet)) {
        throw std::invalid_argument("Masked flow geometry must contain at least one outlet O.");
    }

    geometry_cells_ = mask.cells();
    for (int y = 0; y < grid_.ny; ++y) {
        for (int x = 0; x < grid_.nx; ++x) {
            const int s = grid_.scalar_index(x, y);
            const GeometryCell cell = geometry_cells_[s];
            if (cell == GeometryCell::Inlet &&
                (grid_.nx < 2 || geometry_cells_[grid_.scalar_index(1, y)] == GeometryCell::Solid)) {
                throw std::invalid_argument("Each inlet I must connect to a fluid cell on its right.");
            }
            if (cell == GeometryCell::Outlet &&
                (grid_.nx < 2 ||
                 geometry_cells_[grid_.scalar_index(grid_.nx - 2, y)] == GeometryCell::Solid)) {
                throw std::invalid_argument("Each outlet O must connect to a fluid cell on its left.");
            }
        }
    }

    masked_flow_initialized_ = true;
    for (int y = 0; y < grid_.ny; ++y) {
        for (int x = 0; x < grid_.nx; ++x) {
            const int s = grid_.scalar_index(x, y);
            const GeometryCell cell = geometry_cells_[s];
            grid_.solid[s] = cell == GeometryCell::Solid ? 1 : 0;
            grid_.rho[s] = config_.initial_rho;
            grid_.ux[s] = cell == GeometryCell::Inlet ? inlet_velocity : 0.0;
            grid_.uy[s] = 0.0;
            for (int q = 0; q < D2Q9::q; ++q) {
                grid_.f[grid_.dist_index(x, y, q)] =
                    equilibrium(q, grid_.rho[s], grid_.ux[s], grid_.uy[s]);
            }
        }
    }

    apply_masked_boundaries(inlet_velocity);
    compute_macroscopic_masked();
}

void Solver::step() {
    // 周期单相更新：宏观量求矩 -> 可配置碰撞 -> 迁移。
    compute_macroscopic();
    collide();
    stream_periodic();
    compute_macroscopic();
}

void Solver::step_lid_driven_cavity(double lid_velocity) {
    // 方腔案例跳过固体壁面内的碰撞，并在迁移阶段执行反弹边界。
    compute_macroscopic_fluid_only(lid_velocity);
    collide_fluid_only();
    stream_lid_driven_cavity(lid_velocity);
    compute_macroscopic_fluid_only(lid_velocity);
}

void Solver::step_masked_flow(double inlet_velocity) {
    if (!masked_flow_initialized_) {
        throw std::logic_error("initialize_masked_flow must be called before step_masked_flow.");
    }
    if (!std::isfinite(inlet_velocity)) {
        throw std::invalid_argument("inlet_velocity must be finite.");
    }

    compute_macroscopic_masked();
    collide_fluid_only();
    stream_masked_flow();
    apply_masked_boundaries(inlet_velocity);
    compute_macroscopic_masked();
}

void Solver::run(int steps) {
    if (steps < 0) {
        throw std::invalid_argument("steps must be non-negative.");
    }
    for (int i = 0; i < steps; ++i) {
        step();
    }
    compute_macroscopic();
}

void Solver::run_lid_driven_cavity(int steps, double lid_velocity) {
    if (steps < 0) {
        throw std::invalid_argument("steps must be non-negative.");
    }
    for (int i = 0; i < steps; ++i) {
        step_lid_driven_cavity(lid_velocity);
    }
    compute_macroscopic_fluid_only(lid_velocity);
}

Diagnostics Solver::diagnostics() const {
    Diagnostics d{};
    for (int y = 0; y < grid_.ny; ++y) {
        for (int x = 0; x < grid_.nx; ++x) {
            const int s = grid_.scalar_index(x, y);
            d.mass += grid_.rho[s];
            const double speed = std::sqrt(grid_.ux[s] * grid_.ux[s] + grid_.uy[s] * grid_.uy[s]);
            d.max_speed = std::max(d.max_speed, speed);
        }
    }
    return d;
}

void Solver::write_csv(const std::string& path) const {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Failed to open output file: " + path);
    }

    out << "x,y,rho,ux,uy\n";
    out << std::setprecision(17);
    for (int y = 0; y < grid_.ny; ++y) {
        for (int x = 0; x < grid_.nx; ++x) {
            const int s = grid_.scalar_index(x, y);
            out << x << ',' << y << ',' << grid_.rho[s] << ',' << grid_.ux[s] << ','
                << grid_.uy[s] << '\n';
        }
    }
}

void Solver::write_vtk(const std::string& path) const {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Failed to open output file: " + path);
    }

    out << "# vtk DataFile Version 3.0\n";
    out << "D2Q9 LBM flow field\n";
    out << "ASCII\n";
    out << "DATASET STRUCTURED_POINTS\n";
    out << "DIMENSIONS " << grid_.nx << ' ' << grid_.ny << " 1\n";
    out << "ORIGIN 0 0 0\n";
    out << "SPACING 1 1 1\n";
    out << "POINT_DATA " << grid_.nx * grid_.ny << '\n';

    out << "SCALARS density double 1\n";
    out << "LOOKUP_TABLE default\n";
    out << std::setprecision(17);
    for (int y = 0; y < grid_.ny; ++y) {
        for (int x = 0; x < grid_.nx; ++x) {
            out << grid_.rho[grid_.scalar_index(x, y)] << '\n';
        }
    }

    out << "SCALARS speed double 1\n";
    out << "LOOKUP_TABLE default\n";
    for (int y = 0; y < grid_.ny; ++y) {
        for (int x = 0; x < grid_.nx; ++x) {
            const int s = grid_.scalar_index(x, y);
            out << std::sqrt(grid_.ux[s] * grid_.ux[s] + grid_.uy[s] * grid_.uy[s]) << '\n';
        }
    }

    out << "VECTORS velocity double\n";
    for (int y = 0; y < grid_.ny; ++y) {
        for (int x = 0; x < grid_.nx; ++x) {
            const int s = grid_.scalar_index(x, y);
            out << grid_.ux[s] << ' ' << grid_.uy[s] << " 0\n";
        }
    }
}

double Solver::equilibrium(int direction, double rho, double ux, double uy) {
    // Maxwell-Boltzmann 平衡分布的低马赫数二阶展开。
    return detail::equilibrium(direction, rho, ux, uy);
}

void Solver::collide() {
    for (int y = 0; y < grid_.ny; ++y) {
        for (int x = 0; x < grid_.nx; ++x) {
            collide_cell(x, y);
        }
    }
}

void Solver::collide_fluid_only() {
    for (int y = 0; y < grid_.ny; ++y) {
        for (int x = 0; x < grid_.nx; ++x) {
            const int s = grid_.scalar_index(x, y);
            if (grid_.solid[s]) {
                continue;
            }
            collide_cell(x, y);
        }
    }
}

void Solver::collide_cell(int x, int y) {
    const int s = grid_.scalar_index(x, y);
    detail::D2Q9Population population{};
    for (int q = 0; q < D2Q9::q; ++q) {
        population[q] = grid_.f[grid_.dist_index(x, y, q)];
    }
    detail::collide_population(
        population,
        grid_.rho[s],
        grid_.ux[s],
        grid_.uy[s],
        config_.tau,
        config_.collision_model,
        config_.mrt);
    for (int q = 0; q < D2Q9::q; ++q) {
        grid_.f[grid_.dist_index(x, y, q)] = population[q];
    }
}

void Solver::stream_periodic() {
    std::fill(grid_.f_next.begin(), grid_.f_next.end(), 0.0);

    // 迁移步骤把碰撞后的分布函数沿各自格子速度送到邻居单元。
    for (int y = 0; y < grid_.ny; ++y) {
        for (int x = 0; x < grid_.nx; ++x) {
            for (int q = 0; q < D2Q9::q; ++q) {
                const int dst_x = grid_.wrap_x(x + D2Q9::cx[q]);
                const int dst_y = grid_.wrap_y(y + D2Q9::cy[q]);
                grid_.f_next[grid_.dist_index(dst_x, dst_y, q)] =
                    grid_.f[grid_.dist_index(x, y, q)];
            }
        }
    }

    grid_.f.swap(grid_.f_next);
}

void Solver::stream_lid_driven_cavity(double lid_velocity) {
    std::fill(grid_.f_next.begin(), grid_.f_next.end(), 0.0);

    for (int y = 0; y < grid_.ny; ++y) {
        for (int x = 0; x < grid_.nx; ++x) {
            const int s = grid_.scalar_index(x, y);
            if (grid_.solid[s]) {
                continue;
            }

            for (int q = 0; q < D2Q9::q; ++q) {
                const int dst_x = x + D2Q9::cx[q];
                const int dst_y = y + D2Q9::cy[q];
                const bool outside = dst_x < 0 || dst_x >= grid_.nx || dst_y < 0 || dst_y >= grid_.ny;

                if (outside || grid_.solid[grid_.scalar_index(dst_x, dst_y)]) {
                    const bool hits_lid = !outside && dst_y == grid_.ny - 1;
                    const double wall_ux = hits_lid ? lid_velocity : 0.0;
                    const double wall_uy = 0.0;
                    const double cu_wall = static_cast<double>(D2Q9::cx[q]) * wall_ux +
                                           static_cast<double>(D2Q9::cy[q]) * wall_uy;
                    const int opposite = D2Q9::opposite[q];
                    // 运动壁面反弹在顶盖方向加入动量修正。
                    const double correction = 2.0 * D2Q9::w[q] * grid_.rho[s] * cu_wall / D2Q9::cs2;
                    grid_.f_next[grid_.dist_index(x, y, opposite)] =
                        grid_.f[grid_.dist_index(x, y, q)] - correction;
                    continue;
                }

                grid_.f_next[grid_.dist_index(dst_x, dst_y, q)] =
                    grid_.f[grid_.dist_index(x, y, q)];
            }
        }
    }

    grid_.f.swap(grid_.f_next);
}

void Solver::stream_masked_flow() {
    std::fill(grid_.f_next.begin(), grid_.f_next.end(), 0.0);

    for (int y = 0; y < grid_.ny; ++y) {
        for (int x = 0; x < grid_.nx; ++x) {
            const int s = grid_.scalar_index(x, y);
            if (grid_.solid[s]) {
                continue;
            }

            for (int q = 0; q < D2Q9::q; ++q) {
                const int dst_x = x + D2Q9::cx[q];
                const int dst_y = y + D2Q9::cy[q];
                const bool outside =
                    dst_x < 0 || dst_x >= grid_.nx || dst_y < 0 || dst_y >= grid_.ny;
                if (outside || grid_.solid[grid_.scalar_index(dst_x, dst_y)]) {
                    grid_.f_next[grid_.dist_index(x, y, D2Q9::opposite[q])] +=
                        grid_.f[grid_.dist_index(x, y, q)];
                    continue;
                }
                grid_.f_next[grid_.dist_index(dst_x, dst_y, q)] +=
                    grid_.f[grid_.dist_index(x, y, q)];
            }
        }
    }

    grid_.f.swap(grid_.f_next);
}

void Solver::apply_masked_boundaries(double inlet_velocity) {
    for (int y = 0; y < grid_.ny; ++y) {
        for (int x = 0; x < grid_.nx; ++x) {
            const int s = grid_.scalar_index(x, y);
            const GeometryCell cell = geometry_cells_[s];
            if (cell == GeometryCell::Inlet) {
                const int source_x = x + 1;
                double source_rho = 0.0;
                double source_momentum_x = 0.0;
                double source_momentum_y = 0.0;
                for (int q = 0; q < D2Q9::q; ++q) {
                    const double value = grid_.f[grid_.dist_index(source_x, y, q)];
                    source_rho += value;
                    source_momentum_x += static_cast<double>(D2Q9::cx[q]) * value;
                    source_momentum_y += static_cast<double>(D2Q9::cy[q]) * value;
                }
                if (source_rho <= 0.0 || !std::isfinite(source_rho)) {
                    throw std::runtime_error("Non-physical density next to masked inlet.");
                }
                const double source_ux = source_momentum_x / source_rho;
                const double source_uy = source_momentum_y / source_rho;
                for (int q = 0; q < D2Q9::q; ++q) {
                    const double non_equilibrium =
                        grid_.f[grid_.dist_index(source_x, y, q)] -
                        equilibrium(q, source_rho, source_ux, source_uy);
                    grid_.f[grid_.dist_index(x, y, q)] =
                        equilibrium(q, config_.initial_rho, inlet_velocity, 0.0) +
                        non_equilibrium;
                }
            } else if (cell == GeometryCell::Outlet) {
                // 复制出口内侧相邻格点的分布函数，形成一阶零梯度出口。
                for (int q = 0; q < D2Q9::q; ++q) {
                    grid_.f[grid_.dist_index(x, y, q)] =
                        grid_.f[grid_.dist_index(x - 1, y, q)];
                }
            }
        }
    }
}

void Solver::compute_macroscopic() {
    for (int y = 0; y < grid_.ny; ++y) {
        for (int x = 0; x < grid_.nx; ++x) {
            const int s = grid_.scalar_index(x, y);
            double rho = 0.0;
            double momentum_x = 0.0;
            double momentum_y = 0.0;

            // 密度和动量分别是分布函数的零阶矩和一阶矩。
            for (int q = 0; q < D2Q9::q; ++q) {
                const double fq = grid_.f[grid_.dist_index(x, y, q)];
                rho += fq;
                momentum_x += static_cast<double>(D2Q9::cx[q]) * fq;
                momentum_y += static_cast<double>(D2Q9::cy[q]) * fq;
            }

            if (rho <= 0.0 || !std::isfinite(rho)) {
                throw std::runtime_error("Non-physical density encountered.");
            }

            grid_.rho[s] = rho;
            grid_.ux[s] = momentum_x / rho;
            grid_.uy[s] = momentum_y / rho;
        }
    }
}

void Solver::compute_macroscopic_fluid_only(double lid_velocity) {
    for (int y = 0; y < grid_.ny; ++y) {
        for (int x = 0; x < grid_.nx; ++x) {
            const int s = grid_.scalar_index(x, y);
            if (grid_.solid[s]) {
                grid_.rho[s] = config_.initial_rho;
                grid_.ux[s] = y == grid_.ny - 1 ? lid_velocity : 0.0;
                grid_.uy[s] = 0.0;
                continue;
            }

            double rho = 0.0;
            double momentum_x = 0.0;
            double momentum_y = 0.0;

            for (int q = 0; q < D2Q9::q; ++q) {
                const double fq = grid_.f[grid_.dist_index(x, y, q)];
                rho += fq;
                momentum_x += static_cast<double>(D2Q9::cx[q]) * fq;
                momentum_y += static_cast<double>(D2Q9::cy[q]) * fq;
            }

            if (rho <= 0.0 || !std::isfinite(rho)) {
                throw std::runtime_error("Non-physical density encountered.");
            }

            grid_.rho[s] = rho;
            grid_.ux[s] = momentum_x / rho;
            grid_.uy[s] = momentum_y / rho;
        }
    }
}

void Solver::compute_macroscopic_masked() {
    for (int y = 0; y < grid_.ny; ++y) {
        for (int x = 0; x < grid_.nx; ++x) {
            const int s = grid_.scalar_index(x, y);
            if (grid_.solid[s]) {
                grid_.rho[s] = config_.initial_rho;
                grid_.ux[s] = 0.0;
                grid_.uy[s] = 0.0;
                continue;
            }

            double rho = 0.0;
            double momentum_x = 0.0;
            double momentum_y = 0.0;
            for (int q = 0; q < D2Q9::q; ++q) {
                const double fq = grid_.f[grid_.dist_index(x, y, q)];
                rho += fq;
                momentum_x += static_cast<double>(D2Q9::cx[q]) * fq;
                momentum_y += static_cast<double>(D2Q9::cy[q]) * fq;
            }
            if (rho <= 0.0 || !std::isfinite(rho)) {
                throw std::runtime_error("Non-physical density encountered in masked flow.");
            }
            grid_.rho[s] = rho;
            grid_.ux[s] = momentum_x / rho;
            grid_.uy[s] = momentum_y / rho;
        }
    }
}

FoamSnapshot Solver::openfoam_snapshot() const {
    FoamSnapshot snapshot;
    snapshot.nx = grid_.nx;
    snapshot.ny = grid_.ny;
    snapshot.nz = 1;
    snapshot.two_dimensional = true;
    snapshot.two_phase = false;
    snapshot.cells.resize(grid_.rho.size());
    for (std::size_t i = 0; i < snapshot.cells.size(); ++i) {
        auto& cell = snapshot.cells[i];
        cell.solid = grid_.solid[i] != 0;
        cell.rho = grid_.rho[i];
        cell.velocity = {grid_.ux[i], grid_.uy[i], 0.0};
        if (geometry_cells_.size() == snapshot.cells.size()) {
            cell.boundary = geometry_cells_[i] == GeometryCell::Inlet ? 1 :
                            (geometry_cells_[i] == GeometryCell::Outlet ? 2 : 0);
        }
    }
    return snapshot;
}

void Solver::initialize_fields(const FoamSnapshot& input) {
    validate_initial_snapshot(input);
    if (input.nx != grid_.nx || input.ny != grid_.ny || input.nz != 1 ||
        input.two_phase != false || input.two_dimensional != true ||
        input.cells.size() != grid_.rho.size())
        throw std::invalid_argument("Imported fields do not match solver dimensions/model.");
    geometry_cells_.resize(input.cells.size());
    masked_flow_initialized_ = true;
    for (int z = 0; z < 1; ++z) {
        for (int y = 0; y < grid_.ny; ++y) {
            for (int x = 0; x < grid_.nx; ++x) {
                const auto i = grid_.scalar_index(x, y);
                const auto& cell = input.cells[i];
                grid_.solid[i] = cell.solid ? 1 : 0;
                geometry_cells_[i] = cell.solid ? GeometryCell::Solid :
                    (cell.boundary == 1 ? GeometryCell::Inlet :
                     (cell.boundary == 2 ? GeometryCell::Outlet : GeometryCell::Fluid));
                grid_.rho[i] = cell.rho;
                grid_.ux[i] = cell.velocity[0];
                grid_.uy[i] = cell.velocity[1];
                for (int q = 0; q < D2Q9::q; ++q)
                    grid_.f[grid_.dist_index(x,y,q)] = equilibrium(q, cell.rho, cell.velocity[0], cell.velocity[1]);
            }
        }
    }
}

} // namespace lbm
