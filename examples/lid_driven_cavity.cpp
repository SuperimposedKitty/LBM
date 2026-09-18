#include "lbm/solver.hpp"

#include <iostream>
#include <stdexcept>

int main() {
    const int nx = 96;
    const int ny = 96;
    const int steps = 8000;
    const int output_interval = 200;
    const double lid_velocity = 0.06;

    lbm::SolverConfig config;
    config.tau = 0.8;
    config.initial_rho = 1.0;

    lbm::Solver solver(nx, ny, config);
    solver.initialize_lid_driven_cavity(lid_velocity);

    const auto before = solver.diagnostics();
    lbm::OpenFoamWriter foam(lbm::openfoam_result_path("lid_driven_cavity"), solver.openfoam_snapshot());
    foam.write(0, solver.openfoam_snapshot());
    std::cout << "OpenFOAM case: " << (foam.directory() / "lbm.foam").string() << '\n';

    for (int step = 1; step <= steps; ++step) {
        solver.step_lid_driven_cavity(lid_velocity);
        if (step % output_interval == 0 || step == steps) {
            foam.write(step, solver.openfoam_snapshot());
        }
    }


    const auto after = solver.diagnostics();
    const double reynolds = lid_velocity * static_cast<double>(ny - 2) /
                            (lbm::D2Q9::cs2 * (config.tau - 0.5));

    std::cout << "D2Q9 lid-driven cavity example\n";
    std::cout << "grid: " << nx << " x " << ny << '\n';
    std::cout << "steps: " << steps << '\n';
    std::cout << "output interval: " << output_interval << '\n';
    std::cout << "tau: " << config.tau << '\n';
    std::cout << "lid velocity: " << lid_velocity << '\n';
    std::cout << "estimated Reynolds number: " << reynolds << '\n';
    std::cout << "initial mass: " << before.mass << '\n';
    std::cout << "final mass: " << after.mass << '\n';
    std::cout << "mass difference: " << after.mass - before.mass << '\n';
    std::cout << "final max speed: " << after.max_speed << '\n';

    return 0;
}
