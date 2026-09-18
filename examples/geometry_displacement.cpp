#include "lbm/geometry_mask.hpp"
#include "lbm/two_phase_solver.hpp"

#include <exception>
#include <filesystem>
#include <iostream>
#include <stdexcept>

int main(int argc, char* argv[]) {
    try {
        if (argc > 2) {
            throw std::invalid_argument("Usage: lbm_geometry_displacement.exe [geometry.geom]");
        }
        const std::filesystem::path geometry_file =
            argc == 2 ? std::filesystem::path(argv[1])
                      : lbm::geometry_path("channel_obstacle.geom");
        const lbm::GeometryMask mask = lbm::GeometryMask::load(geometry_file);

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
        config.free_flow_porosity = 1.0;
        config.porous_porosity = 0.3;
        config.porous_pore_diameter = 40.0;
        config.darcy_drag_scale = 0.06;
        config.forchheimer_drag_scale = 0.02;

        lbm::TwoPhaseSolver solver(mask.width(), mask.height(), config);
        solver.initialize_geometry_displacement(mask);
        const auto before = solver.diagnostics();
        lbm::OpenFoamWriter foam(lbm::openfoam_result_path("geometry_displacement"), solver.openfoam_snapshot());
        foam.write(0, solver.openfoam_snapshot());
        std::cout << "OpenFOAM case: " << (foam.directory() / "lbm.foam").string() << '\n';

        for (int step = 1; step <= steps; ++step) {
            solver.step();
            if (step % output_interval == 0 || step == steps) {
                foam.write(step, solver.openfoam_snapshot());
            }
        }


        const auto after = solver.diagnostics();
        std::cout << "D2Q9 imported-geometry two-phase displacement\n";
        std::cout << "geometry: " << geometry_file.string() << '\n';
        std::cout << "grid: " << mask.width() << " x " << mask.height() << '\n';
        std::cout << "steps: " << steps << '\n';
        std::cout << "initial A centroid x: " << before.interface_x << '\n';
        std::cout << "final A centroid x: " << after.interface_x << '\n';
        std::cout << "final max speed: " << after.max_speed << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "geometry displacement error: " << error.what() << '\n';
        return 1;
    }
}
