#include "lbm/geometry_mask.hpp"
#include "lbm/geometry_mask3d.hpp"
#include "lbm/two_phase_solver3d.hpp"

#include <filesystem>
#include <iostream>
#include <stdexcept>

int main(int argc, char** argv) {
    try {
        const std::filesystem::path geometry_file =
            argc > 1 ? std::filesystem::path(argv[1])
                     : lbm::geometry_path("channel_obstacle.geom3d");
        const lbm::GeometryMask3D mask = lbm::GeometryMask3D::load(geometry_file);
        constexpr int steps = 800;
        constexpr int output_interval = 40;

        lbm::TwoPhaseConfig3D config;
        config.lattice_model = lbm::Lattice3DModel::D3Q27;
        config.tau_a = 1.0;
        config.tau_b = 1.0;
        config.interaction_strength = 2.0;
        config.inlet_velocity = 0.008;
        config.body_force_y = 0.0;
        config.contact_angle_degrees = 70.0;
        config.wall_adhesion_strength = 0.02;
        config.recoloring_strength = 0.05;
        config.porous_porosity = 0.3;

        lbm::TwoPhaseSolver3D solver(
            mask.width(), mask.height(), mask.depth(), config);
        solver.initialize_geometry_displacement(mask);
        const auto before = solver.diagnostics();

        lbm::OpenFoamWriter foam(lbm::openfoam_result_path("d3_geometry_displacement"), solver.openfoam_snapshot());
        foam.write(0, solver.openfoam_snapshot());
        std::cout << "OpenFOAM case: " << (foam.directory() / "lbm.foam").string() << '\n';
        for (int step = 1; step <= steps; ++step) {
            solver.step();
            if (step % output_interval == 0 || step == steps) {
                foam.write(step, solver.openfoam_snapshot());
            }
        }


        const auto after = solver.diagnostics();
        std::cout << "D3Q27 imported-geometry two-phase displacement\n";
        std::cout << "geometry: " << geometry_file.string() << '\n';
        std::cout << "initial A centroid x: " << before.phase_a_centroid_x
                  << '\n';
        std::cout << "final A centroid x: " << after.phase_a_centroid_x << '\n';
        std::cout << "final max speed: " << after.max_speed << '\n';
        std::cout << "porous mean pore speed: "
                  << after.porous_mean_pore_speed << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "3D geometry displacement error: " << error.what() << '\n';
        return 1;
    }
}
