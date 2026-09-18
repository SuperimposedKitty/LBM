#include "lbm/solver.hpp"

#include <iostream>
#include <stdexcept>

int main() {
    const int nx = 128;
    const int ny = 64;
    const int steps = 1000;
    const int output_interval = 50;

    lbm::SolverConfig config;
    config.tau = 0.8;
    config.initial_rho = 1.0;

    lbm::Solver solver(nx, ny, config);
    solver.initialize_shear_wave(0.05, 1);

    const auto before = solver.diagnostics();
    lbm::OpenFoamWriter foam(lbm::openfoam_result_path("periodic_shear"), solver.openfoam_snapshot());
    foam.write(0, solver.openfoam_snapshot());
    std::cout << "OpenFOAM case: " << (foam.directory() / "lbm.foam").string() << '\n';

    for (int step = 1; step <= steps; ++step) {
        solver.step();
        if (step % output_interval == 0 || step == steps) {
            foam.write(step, solver.openfoam_snapshot());
        }
    }


    const auto after = solver.diagnostics();
    std::cout << "D2Q9 periodic shear wave example\n";
    std::cout << "grid: " << nx << " x " << ny << '\n';
    std::cout << "steps: " << steps << '\n';
    std::cout << "output interval: " << output_interval << '\n';
    std::cout << "tau: " << config.tau << '\n';
    std::cout << "initial mass: " << before.mass << '\n';
    std::cout << "final mass: " << after.mass << '\n';
    std::cout << "mass difference: " << after.mass - before.mass << '\n';
    std::cout << "final max speed: " << after.max_speed << '\n';

    return 0;
}
