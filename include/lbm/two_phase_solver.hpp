#pragma once

#include "lbm/openfoam.hpp"

#include "lbm/collision.hpp"
#include "lbm/geometry_mask.hpp"
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
    double body_force_y = 0.0;
    // 接触角从 A 相一侧度量；小于 90 度时，注入相更容易润湿固壁。
    double contact_angle_degrees = 90.0;
    double wall_adhesion_strength = 0.08;
    // 液滴撞击案例可单独加厚底部固壁，并压窄初始相界面以减少过度弥散。
    int bottom_wall_thickness = 1;
    double droplet_interface_width = 2.0;
    // Recoloring 强度用于抑制两相界面数值扩散；0 表示关闭。
    double recoloring_strength = 0.0;
    // 液滴撞击可用短程壁面排斥力表示非润湿表面的弹性反推。
    double bottom_wall_repulsion_strength = 0.0;
    int bottom_wall_repulsion_range = 0;
    // 孔隙度为 1 表示自由流动区；中间多孔区默认孔隙度为 0.3。
    double free_flow_porosity = 1.0;
    double porous_porosity = 0.3;
    int porous_start_x = -1;
    int porous_end_x = -1;
    // 等效孔径和阻力缩放用于把孔隙尺度阻力映射到当前格点尺度。
    double porous_pore_diameter = 40.0;
    double darcy_drag_scale = 0.06;
    double forchheimer_drag_scale = 0.02;
    CollisionModel collision_model = CollisionModel::MRT;
    MrtRelaxationRates mrt{};
};

struct TwoPhaseDiagnostics {
    double mass_a = 0.0;
    double mass_b = 0.0;
    double max_speed = 0.0;
    double interface_x = 0.0;
    double phase_a_centroid_y = 0.0;
    double porous_mean_pore_speed = 0.0;
    double porous_max_pore_speed = 0.0;
    double interaction_force_balance_x = 0.0;
    double interaction_force_balance_y = 0.0;
};

class TwoPhaseSolver {
public:
    TwoPhaseSolver(int nx, int ny, TwoPhaseConfig config);

    int nx() const;
    int ny() const;
    double phase_at(int x, int y) const;
    double porosity_at(int x, int y) const;
    bool solid_at(int x, int y) const;

    void initialize_capillary_displacement();
    void initialize_geometry_displacement(const GeometryMask& mask);
    void initialize_droplet_impact(
        double center_x, double center_y, double radius, double initial_ux, double initial_uy);
    void step();
    void step_closed();
    void run(int steps);
    void run_closed(int steps);

    TwoPhaseDiagnostics diagnostics() const;
    void write_csv(const std::string& path) const;
    void write_vtk(const std::string& path) const;
    FoamSnapshot openfoam_snapshot() const;
    // 用宏观初始场构造平衡分布，不代表非平衡分布的断点续算。
    void initialize_fields(const FoamSnapshot& input);

private:
    int nx_{};
    int ny_{};
    TwoPhaseConfig config_{};

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
    std::vector<double> porosity_;
    std::vector<std::uint8_t> solid_;
    std::vector<GeometryCell> geometry_cells_;
    bool geometry_displacement_initialized_ = false;
    double interaction_force_balance_x_ = 0.0;
    double interaction_force_balance_y_ = 0.0;

    int scalar_index(int x, int y) const;
    int dist_index(int x, int y, int direction) const;

    static double psi(double rho);
    static double equilibrium(int direction, double rho, double ux, double uy);

    void compute_macroscopic();
    void compute_forces();
    void collide();
    void recolor();
    void stream();
    void apply_inlet_outlet();
    void apply_geometry_inlet_outlet();
    void set_equilibrium_cell(int x, int y, double rho_a, double rho_b, double ux, double uy);
    bool is_porous_column(int x) const;
    double permeability_from_porosity(double porosity) const;
    bool is_outside(int x, int y) const;
};

} // namespace lbm
