#include "lbm/two_phase_solver.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace lbm {
namespace {

constexpr double min_density = 1.0e-12;

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
      omega_a_(1.0 / config.tau_a),
      omega_b_(1.0 / config.tau_b),
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
      solid_(rho_a_.size(), 0) {
    if (nx <= 4 || ny <= 4) {
        throw std::invalid_argument("Two-phase grid dimensions must be greater than 4.");
    }
    if (config_.tau_a <= 0.5 || config_.tau_b <= 0.5) {
        throw std::invalid_argument("tau_a and tau_b must be greater than 0.5.");
    }
    if (config_.rho_high <= 0.0 || config_.rho_low <= 0.0) {
        throw std::invalid_argument("rho_high and rho_low must be positive.");
    }
    if (config_.initial_interface_x <= 1 || config_.initial_interface_x >= nx_ - 2) {
        config_.initial_interface_x = nx_ / 5;
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
    return (rho_a_[s] - rho_b_[s]) / std::max(rho_[s], min_density);
}

bool TwoPhaseSolver::solid_at(int x, int y) const {
    if (is_outside(x, y)) {
        throw std::out_of_range("solid_at coordinates are outside the grid.");
    }
    return solid_[scalar_index(x, y)] != 0;
}

void TwoPhaseSolver::initialize_capillary_displacement() {
    for (int y = 0; y < ny_; ++y) {
        for (int x = 0; x < nx_; ++x) {
            const int s = scalar_index(x, y);
            solid_[s] = y == 0 || y == ny_ - 1 ? 1 : 0;

            const bool injected_region = x <= config_.initial_interface_x;
            const double rho_a = injected_region ? config_.rho_high : config_.rho_low;
            const double rho_b = injected_region ? config_.rho_low : config_.rho_high;
            set_equilibrium_cell(x, y, rho_a, rho_b, 0.0, 0.0);
        }
    }
    apply_inlet_outlet();
    compute_macroscopic();
}

void TwoPhaseSolver::step() {
    compute_macroscopic();
    compute_forces();
    collide();
    stream();
    apply_inlet_outlet();
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

TwoPhaseDiagnostics TwoPhaseSolver::diagnostics() const {
    TwoPhaseDiagnostics d{};
    double weighted_interface = 0.0;
    double injected_mass = 0.0;

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

            const double volume_fraction_a = rho_a_[s] / std::max(rho_[s], min_density);
            weighted_interface += static_cast<double>(x) * volume_fraction_a;
            injected_mass += volume_fraction_a;
        }
    }

    d.interface_x = injected_mass > 0.0 ? weighted_interface / injected_mass : 0.0;
    return d;
}

void TwoPhaseSolver::write_csv(const std::string& path) const {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Failed to open output file: " + path);
    }

    out << "x,y,rho_a,rho_b,rho,phi,ux,uy,solid\n";
    out << std::setprecision(17);
    for (int y = 0; y < ny_; ++y) {
        for (int x = 0; x < nx_; ++x) {
            const int s = scalar_index(x, y);
            const double phi = (rho_a_[s] - rho_b_[s]) / std::max(rho_[s], min_density);
            out << x << ',' << y << ',' << rho_a_[s] << ',' << rho_b_[s] << ',' << rho_[s] << ','
                << phi << ',' << ux_[s] << ',' << uy_[s] << ',' << static_cast<int>(solid_[s]) << '\n';
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
    return 1.0 - std::exp(-std::max(rho, 0.0));
}

double TwoPhaseSolver::equilibrium(int direction, double rho, double ux, double uy) {
    const double cu = static_cast<double>(D2Q9::cx[direction]) * ux +
                      static_cast<double>(D2Q9::cy[direction]) * uy;
    const double u2 = ux * ux + uy * uy;
    return D2Q9::w[direction] * rho *
           (1.0 + 3.0 * cu + 4.5 * cu * cu - 1.5 * u2);
}

double TwoPhaseSolver::forcing_term(
    int direction, double ux, double uy, double fx, double fy, double omega) {
    const double cx = static_cast<double>(D2Q9::cx[direction]);
    const double cy = static_cast<double>(D2Q9::cy[direction]);
    const double cu = cx * ux + cy * uy;
    const double term_x = (cx - ux) / D2Q9::cs2 + cu * cx / (D2Q9::cs2 * D2Q9::cs2);
    const double term_y = (cy - uy) / D2Q9::cs2 + cu * cy / (D2Q9::cs2 * D2Q9::cs2);
    return (1.0 - 0.5 * omega) * D2Q9::w[direction] * (term_x * fx + term_y * fy);
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

    for (int y = 1; y < ny_ - 1; ++y) {
        for (int x = 0; x < nx_; ++x) {
            const int s = scalar_index(x, y);
            if (solid_[s]) {
                continue;
            }

            double sum_bx = 0.0;
            double sum_by = 0.0;
            double sum_ax = 0.0;
            double sum_ay = 0.0;
            for (int q = 1; q < D2Q9::q; ++q) {
                const int nx = x + D2Q9::cx[q];
                const int ny = y + D2Q9::cy[q];
                if (ny <= 0 || ny >= ny_ - 1) {
                    continue;
                }
                const int wrapped_x = nx < 0 ? 0 : (nx >= nx_ ? nx_ - 1 : nx);
                const int neighbor = scalar_index(wrapped_x, ny);
                if (solid_[neighbor]) {
                    continue;
                }

                sum_bx += D2Q9::w[q] * psi(rho_b_[neighbor]) * static_cast<double>(D2Q9::cx[q]);
                sum_by += D2Q9::w[q] * psi(rho_b_[neighbor]) * static_cast<double>(D2Q9::cy[q]);
                sum_ax += D2Q9::w[q] * psi(rho_a_[neighbor]) * static_cast<double>(D2Q9::cx[q]);
                sum_ay += D2Q9::w[q] * psi(rho_a_[neighbor]) * static_cast<double>(D2Q9::cy[q]);
            }

            const double psi_a = psi(rho_a_[s]);
            const double psi_b = psi(rho_b_[s]);
            force_ax_[s] = -config_.interaction_strength * psi_a * sum_bx;
            force_ay_[s] = -config_.interaction_strength * psi_a * sum_by;
            force_bx_[s] = -config_.interaction_strength * psi_b * sum_ax;
            force_by_[s] = -config_.interaction_strength * psi_b * sum_ay;

            const double total_rho = std::max(rho_[s], min_density);
            force_ax_[s] += config_.body_force_x * rho_a_[s] / total_rho;
            force_bx_[s] += config_.body_force_x * rho_b_[s] / total_rho;
        }
    }
}

void TwoPhaseSolver::collide() {
    for (int y = 1; y < ny_ - 1; ++y) {
        for (int x = 0; x < nx_; ++x) {
            const int s = scalar_index(x, y);
            if (solid_[s]) {
                continue;
            }

            for (int q = 0; q < D2Q9::q; ++q) {
                const int k = dist_index(x, y, q);
                const double feq_a = equilibrium(q, rho_a_[s], ux_[s], uy_[s]);
                const double feq_b = equilibrium(q, rho_b_[s], ux_[s], uy_[s]);
                fa_[k] -= omega_a_ * (fa_[k] - feq_a);
                fb_[k] -= omega_b_ * (fb_[k] - feq_b);
                fa_[k] += forcing_term(q, ux_[s], uy_[s], force_ax_[s], force_ay_[s], omega_a_);
                fb_[k] += forcing_term(q, ux_[s], uy_[s], force_bx_[s], force_by_[s], omega_b_);
            }
        }
    }
}

void TwoPhaseSolver::stream() {
    std::fill(fa_next_.begin(), fa_next_.end(), 0.0);
    std::fill(fb_next_.begin(), fb_next_.end(), 0.0);

    for (int y = 1; y < ny_ - 1; ++y) {
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
    for (int y = 1; y < ny_ - 1; ++y) {
        set_equilibrium_cell(
            0, y, config_.rho_high, config_.rho_low, config_.inlet_velocity, 0.0);

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

bool TwoPhaseSolver::is_outside(int x, int y) const {
    return x < 0 || x >= nx_ || y < 0 || y >= ny_;
}

} // namespace lbm
