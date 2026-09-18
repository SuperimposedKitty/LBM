#include "lbm/geometry_mask.hpp"
#include "lbm/solver.hpp"

#include <exception>
#include <filesystem>
#include <iostream>
#include <stdexcept>

int main(int argc, char* argv[]) {
    try {
        if (argc > 2) {
            throw std::invalid_argument("Usage: lbm_geometry_flow.exe [geometry.geom]");
        }
        const std::filesystem::path geometry_file =
            argc == 2 ? std::filesystem::path(argv[1])
                      : lbm::geometry_path("channel_obstacle.geom");
        const lbm::GeometryMask mask = lbm::GeometryMask::load(geometry_file);

        const int steps = 5000;
        const int output_interval = 125;
        const double inlet_velocity = 0.045;
        lbm::SolverConfig config;
        config.tau = 0.8;
        config.initial_rho = 1.0;

        lbm::Solver solver(mask.width(), mask.height(), config);
        solver.initialize_masked_flow(mask, inlet_velocity);
        lbm::OpenFoamWriter foam(lbm::openfoam_result_path("geometry_flow"), solver.openfoam_snapshot());
        foam.write(0, solver.openfoam_snapshot());
        std::cout << "OpenFOAM case: " << (foam.directory() / "lbm.foam").string() << '\n';

        for (int step = 1; step <= steps; ++step) {
            solver.step_masked_flow(inlet_velocity);
            if (step % output_interval == 0 || step == steps) {
                foam.write(step, solver.openfoam_snapshot());
            }
        }


        const auto diagnostics = solver.diagnostics();
        std::cout << "D2Q9 imported-geometry single-phase flow\n";
        std::cout << "geometry: " << geometry_file.string() << '\n';
        std::cout << "grid: " << mask.width() << " x " << mask.height() << '\n';
        std::cout << "steps: " << steps << '\n';
        std::cout << "final max speed: " << diagnostics.max_speed << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "geometry flow error: " << error.what() << '\n';
        return 1;
    }
}
