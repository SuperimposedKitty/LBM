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
      config_(config),
      omega_(1.0 / config.tau) {
    if (config_.tau <= 0.5) {
        throw std::invalid_argument("tau must be greater than 0.5 for positive viscosity.");
    }

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

void Solver::step() {
    // 周期单相更新：宏观量求矩 -> BGK 碰撞 -> 迁移。
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
    const double cu = static_cast<double>(D2Q9::cx[direction]) * ux +
                      static_cast<double>(D2Q9::cy[direction]) * uy;
    const double u2 = ux * ux + uy * uy;
    return D2Q9::w[direction] * rho *
           (1.0 + 3.0 * cu + 4.5 * cu * cu - 1.5 * u2);
}

void Solver::collide() {
    for (int y = 0; y < grid_.ny; ++y) {
        for (int x = 0; x < grid_.nx; ++x) {
            const int s = grid_.scalar_index(x, y);
            for (int q = 0; q < D2Q9::q; ++q) {
                const int k = grid_.dist_index(x, y, q);
                const double feq = equilibrium(q, grid_.rho[s], grid_.ux[s], grid_.uy[s]);
                // BGK/SRT 碰撞把每个方向的分布函数松弛到局部平衡态。
                grid_.f[k] -= omega_ * (grid_.f[k] - feq);
            }
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
            for (int q = 0; q < D2Q9::q; ++q) {
                const int k = grid_.dist_index(x, y, q);
                const double feq = equilibrium(q, grid_.rho[s], grid_.ux[s], grid_.uy[s]);
                grid_.f[k] -= omega_ * (grid_.f[k] - feq);
            }
        }
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

} // namespace lbm
