#pragma once

#include "lbm/grid.hpp"

#include <string>

namespace lbm {

struct SolverConfig {
    double tau = 0.8;
    double initial_rho = 1.0;
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
    void step();
    void step_lid_driven_cavity(double lid_velocity);
    void run(int steps);
    void run_lid_driven_cavity(int steps, double lid_velocity);
    Diagnostics diagnostics() const;
    void write_csv(const std::string& path) const;
    void write_vtk(const std::string& path) const;

private:
    Grid grid_;
    SolverConfig config_;
    double omega_;

    static double equilibrium(int direction, double rho, double ux, double uy);
    void collide();
    void collide_fluid_only();
    void stream_periodic();
    void stream_lid_driven_cavity(double lid_velocity);
    void compute_macroscopic();
    void compute_macroscopic_fluid_only(double lid_velocity);
};

} // namespace lbm
