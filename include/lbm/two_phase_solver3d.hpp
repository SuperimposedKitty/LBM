#pragma once

#include "lbm/collision3d.hpp"
#include "lbm/geometry_mask3d.hpp"

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace lbm {

struct TwoPhaseConfig3D {
    double tau_a = 1.0;
    double tau_b = 1.0;
    double rho_high = 1.0;
    double rho_low = 0.02;
    double interaction_strength = 3.0;
    double inlet_velocity = 0.015;
    double body_force_x = 0.0;
    double body_force_y = -2.0e-6;
    double body_force_z = 0.0;
    double contact_angle_degrees = 100.0;
    double wall_adhesion_strength = 0.08;
    int bottom_wall_thickness = 2;
    double droplet_interface_width = 1.5;
    double recoloring_strength = 0.15;
    double bottom_wall_repulsion_strength = 0.0;
    int bottom_wall_repulsion_range = 0;
    double free_flow_porosity = 1.0;
    double porous_porosity = 0.3;
    double porous_pore_diameter = 24.0;
    double darcy_drag_scale = 0.06;
    double forchheimer_drag_scale = 0.02;
    Lattice3DModel lattice_model = Lattice3DModel::D3Q27;
    CollisionModel collision_model = CollisionModel::MRT;
    MrtRelaxationRates3D mrt{};
};

struct TwoPhaseDiagnostics3D {
    double mass_a = 0.0;
    double mass_b = 0.0;
    double max_speed = 0.0;
    double phase_a_centroid_x = 0.0;
    double phase_a_centroid_y = 0.0;
    double phase_a_centroid_z = 0.0;
    double porous_mean_pore_speed = 0.0;
    double porous_max_pore_speed = 0.0;
    double interaction_force_balance_x = 0.0;
    double interaction_force_balance_y = 0.0;
    double interaction_force_balance_z = 0.0;
};

class TwoPhaseSolver3D {
public:
    TwoPhaseSolver3D(int nx, int ny, int nz, TwoPhaseConfig3D config);

    int nx() const;
    int ny() const;
    int nz() const;
    double phase_at(int x, int y, int z) const;
    double porosity_at(int x, int y, int z) const;
    bool solid_at(int x, int y, int z) const;

    void initialize_droplet_impact(
        double center_x,
        double center_y,
        double center_z,
        double radius,
        double initial_ux,
        double initial_uy,
        double initial_uz);
    void initialize_geometry_displacement(const GeometryMask3D& mask);
    void step();
    void step_closed();
    void run(int steps);
    void run_closed(int steps);
    TwoPhaseDiagnostics3D diagnostics() const;
    void write_vtk(const std::string& path) const;

private:
    int nx_{};
    int ny_{};
    int nz_{};
    int q_{};
    TwoPhaseConfig3D config_{};

    std::vector<double> fa_;
    std::vector<double> fb_;
    std::vector<double> fa_next_;
    std::vector<double> fb_next_;
    std::vector<double> rho_a_;
    std::vector<double> rho_b_;
    std::vector<double> rho_;
    std::array<std::vector<double>, 3> velocity_;
    std::array<std::vector<double>, 3> force_a_;
    std::array<std::vector<double>, 3> force_b_;
    std::vector<double> porosity_;
    std::vector<std::uint8_t> solid_;
    std::vector<GeometryCell> geometry_cells_;
    bool geometry_initialized_ = false;
    std::array<double, 3> interaction_force_balance_{};

    std::size_t scalar_index(int x, int y, int z) const;
    std::size_t dist_index(int x, int y, int z, int direction) const;
    static double psi(double rho);
    void set_equilibrium_cell(
        int x,
        int y,
        int z,
        double rho_a,
        double rho_b,
        double ux,
        double uy,
        double uz);
    void compute_macroscopic();
    void compute_forces();
    void collide();
    void recolor();
    void stream();
    void apply_geometry_boundaries();
    double permeability_from_porosity(double porosity) const;
    bool is_outside(int x, int y, int z) const;
};

} // namespace lbm
