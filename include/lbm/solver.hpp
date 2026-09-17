#pragma once

#include "lbm/openfoam.hpp"

#include "lbm/collision.hpp"
#include "lbm/geometry_mask.hpp"
#include "lbm/grid.hpp"

#include <string>
#include <vector>

namespace lbm {

struct SolverConfig {
    // tau 控制运动黏度：nu = cs2 * (tau - 0.5)。
    double tau = 0.8;
    double initial_rho = 1.0;
    CollisionModel collision_model = CollisionModel::MRT;
    MrtRelaxationRates mrt{};
};

struct Diagnostics {
    double mass = 0.0;
    double max_speed = 0.0;
};

class Solver {
public:
    Solver(int nx, int ny, SolverConfig config);

    Grid& grid();
    const Grid& grid() const;

    void initialize_shear_wave(double amplitude, int mode);
    void initialize_lid_driven_cavity(double lid_velocity);
    void initialize_masked_flow(const GeometryMask& mask, double inlet_velocity);
    void step();
    void step_lid_driven_cavity(double lid_velocity);
    void step_masked_flow(double inlet_velocity);
    void run(int steps);
    void run_lid_driven_cavity(int steps, double lid_velocity);
    Diagnostics diagnostics() const;
    void write_csv(const std::string& path) const;
    void write_vtk(const std::string& path) const;
    FoamSnapshot openfoam_snapshot() const;

private:
    Grid grid_;
    SolverConfig config_;
    std::vector<GeometryCell> geometry_cells_;
    bool masked_flow_initialized_ = false;

    static double equilibrium(int direction, double rho, double ux, double uy);
    void collide();
    void collide_fluid_only();
    void collide_cell(int x, int y);
    void stream_periodic();
    void stream_lid_driven_cavity(double lid_velocity);
    void stream_masked_flow();
    void apply_masked_boundaries(double inlet_velocity);
    void compute_macroscopic();
    void compute_macroscopic_fluid_only(double lid_velocity);
    void compute_macroscopic_masked();
};

} // namespace lbm
