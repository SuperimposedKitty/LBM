#pragma once

#include "lbm/lattice.hpp"

#include <cstdint>
#include <string>
#include <vector>

namespace lbm {

struct TwoPhaseConfig {
    // 两个组分各自拥有松弛时间和分布函数。
    double tau_a = 1.0;
    double tau_b = 1.0;
    double rho_high = 1.0;
    double rho_low = 0.02;
    // 正的 Shan-Chen 耦合强度用于促使 A/B 两相分离。
    double interaction_strength = 3.0;
    double inlet_velocity = 0.02;
    double body_force_x = 1.0e-6;
    // 接触角从 A 相一侧度量；小于 90 度时，注入相更容易润湿固壁。
    double contact_angle_degrees = 90.0;
    double wall_adhesion_strength = 0.08;
};

struct TwoPhaseDiagnostics {
    double mass_a = 0.0;
    double mass_b = 0.0;
    double max_speed = 0.0;
    double interface_x = 0.0;
};

class TwoPhaseSolver {
public:
    TwoPhaseSolver(int nx, int ny, TwoPhaseConfig config);

    int nx() const;
    int ny() const;
    double phase_at(int x, int y) const;
    bool solid_at(int x, int y) const;

    void initialize_capillary_displacement();
    void step();
    void run(int steps);

    TwoPhaseDiagnostics diagnostics() const;
    void write_csv(const std::string& path) const;
    void write_vtk(const std::string& path) const;

private:
    int nx_{};
    int ny_{};
    TwoPhaseConfig config_{};
    double omega_a_{};
    double omega_b_{};

    // fa/fb 分别是红/蓝组分的分布函数；rho/ux/uy 是混合物宏观场。
    std::vector<double> fa_;
    std::vector<double> fb_;
    std::vector<double> fa_next_;
    std::vector<double> fb_next_;
    std::vector<double> rho_a_;
    std::vector<double> rho_b_;
    std::vector<double> rho_;
    std::vector<double> ux_;
    std::vector<double> uy_;
    std::vector<double> force_ax_;
    std::vector<double> force_ay_;
    std::vector<double> force_bx_;
    std::vector<double> force_by_;
    std::vector<std::uint8_t> solid_;

    int scalar_index(int x, int y) const;
    int dist_index(int x, int y, int direction) const;

    static double psi(double rho);
    static double equilibrium(int direction, double rho, double ux, double uy);
    static double forcing_term(int direction, double ux, double uy, double fx, double fy, double omega);

    void compute_macroscopic();
    void compute_forces();
    void collide();
    void stream();
    void apply_inlet_outlet();
    void set_equilibrium_cell(int x, int y, double rho_a, double rho_b, double ux, double uy);
    bool is_outside(int x, int y) const;
};

} // namespace lbm
