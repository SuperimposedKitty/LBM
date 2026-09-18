#include "lbm/two_phase_solver.hpp"

#include <iostream>
#include <stdexcept>

int main() {
    const int nx = 200;
    const int ny = 36;
    const int steps = 18000;
    const int output_interval = 300;

    lbm::TwoPhaseConfig config;
    config.tau_a = 1.0;
    config.tau_b = 1.0;
    config.rho_high = 1.0;
    config.rho_low = 0.02;
    config.interaction_strength = 3.0;
    config.inlet_velocity = 0.065;
    config.body_force_x = 2.5e-5;
    config.contact_angle_degrees = 70.0;
    config.wall_adhesion_strength = 0.08;
    config.porous_start_x = 78;
    config.porous_end_x = 122;
    config.porous_porosity = 0.3;
    config.free_flow_porosity = 1.0;
    config.porous_pore_diameter = 40.0;
    config.darcy_drag_scale = 0.06;
    config.forchheimer_drag_scale = 0.02;

    lbm::TwoPhaseSolver solver(nx, ny, config);
    solver.initialize_capillary_displacement();

    const auto before = solver.diagnostics();
    lbm::OpenFoamWriter foam(lbm::openfoam_result_path("capillary_displacement"), solver.openfoam_snapshot());
    foam.write(0, solver.openfoam_snapshot());
    std::cout << "OpenFOAM case: " << (foam.directory() / "lbm.foam").string() << '\n';

    for (int step = 1; step <= steps; ++step) {
        solver.step();
        if (step % output_interval == 0 || step == steps) {
            foam.write(step, solver.openfoam_snapshot());
        }
    }


    const auto after = solver.diagnostics();
    std::cout << "D2Q9 two-phase capillary displacement example\n";
    std::cout << "grid: " << nx << " x " << ny << '\n';
    std::cout << "steps: " << steps << '\n';
    std::cout << "output interval: " << output_interval << '\n';
    std::cout << "tau_a: " << config.tau_a << '\n';
    std::cout << "tau_b: " << config.tau_b << '\n';
    std::cout << "interaction strength: " << config.interaction_strength << '\n';
    std::cout << "inlet velocity: " << config.inlet_velocity << '\n';
    std::cout << "body force x: " << config.body_force_x << '\n';
    std::cout << "contact angle: " << config.contact_angle_degrees << '\n';
    std::cout << "wall adhesion strength: " << config.wall_adhesion_strength << '\n';
    std::cout << "porous medium x range: [" << config.porous_start_x << ", "
              << config.porous_end_x << ")\n";
    std::cout << "porous porosity: " << config.porous_porosity << '\n';
    std::cout << "porous pore diameter: " << config.porous_pore_diameter << '\n';
    std::cout << "Darcy drag scale: " << config.darcy_drag_scale << '\n';
    std::cout << "Forchheimer drag scale: " << config.forchheimer_drag_scale << '\n';
    std::cout << "initial injected-fluid centroid x: " << before.interface_x << '\n';
    std::cout << "final injected-fluid centroid x: " << after.interface_x << '\n';
    std::cout << "injected-fluid centroid displacement: " << after.interface_x - before.interface_x
              << '\n';
    std::cout << "initial mass A: " << before.mass_a << '\n';
    std::cout << "final mass A: " << after.mass_a << '\n';
    std::cout << "initial mass B: " << before.mass_b << '\n';
    std::cout << "final mass B: " << after.mass_b << '\n';
    std::cout << "final max speed: " << after.max_speed << '\n';
    std::cout << "porous mean pore speed: " << after.porous_mean_pore_speed << '\n';
    std::cout << "porous max pore speed: " << after.porous_max_pore_speed << '\n';

    return 0;
}
