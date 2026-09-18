#include "lbm/two_phase_solver3d.hpp"

#include <algorithm>
#include <iostream>
#include <stdexcept>

int main() {
    constexpr int nx = 24;
    constexpr int ny = 25;
    constexpr int nz = 22;
    constexpr int steps = 600;
    constexpr int output_interval = 20;

    lbm::TwoPhaseConfig3D config;
    config.lattice_model = lbm::Lattice3DModel::D3Q27;
    config.collision_model = lbm::CollisionModel::MRT;
    config.tau_a = 0.75;
    config.tau_b = 0.75;
    config.interaction_strength = 3.6;
    config.body_force_y = 0.0;
    config.contact_angle_degrees = 150.0;
    config.wall_adhesion_strength = 0.12;
    config.bottom_wall_thickness = 2;
    config.droplet_interface_width = 1.0;
    config.recoloring_strength = 0.5;
    config.bottom_wall_repulsion_strength = 0.04;
    config.bottom_wall_repulsion_range = 3;

    lbm::TwoPhaseSolver3D solver(nx, ny, nz, config);
    solver.initialize_droplet_impact(
        0.5 * (nx - 1),
        9.0,
        0.5 * (nz - 1),
        5.0,
        0.0,
        -0.1,
        0.0);
    const auto before = solver.diagnostics();

    lbm::OpenFoamWriter foam(lbm::openfoam_result_path("d3_droplet_impact"), solver.openfoam_snapshot());
    foam.write(0, solver.openfoam_snapshot());
    std::cout << "OpenFOAM case: " << (foam.directory() / "lbm.foam").string() << '\n';
    double minimum_center_y = before.phase_a_centroid_y;
    double rebound_center_y = minimum_center_y;
    for (int step = 1; step <= steps; ++step) {
        solver.step_closed();
        const auto current = solver.diagnostics();
        if (current.phase_a_centroid_y < minimum_center_y) {
            minimum_center_y = current.phase_a_centroid_y;
            rebound_center_y = minimum_center_y;
        } else {
            rebound_center_y =
                std::max(rebound_center_y, current.phase_a_centroid_y);
        }
        if (step % output_interval == 0 || step == steps) {
            foam.write(step, solver.openfoam_snapshot());
        }
    }


    const auto after = solver.diagnostics();
    std::cout << "D3Q27 MRT droplet impact\n";
    std::cout << "initial mass A: " << before.mass_a << '\n';
    std::cout << "final mass A: " << after.mass_a << '\n';
    std::cout << "minimum center y: " << minimum_center_y << '\n';
    std::cout << "maximum rebound height: "
              << rebound_center_y - minimum_center_y << '\n';
    std::cout << "final max speed: " << after.max_speed << '\n';
    return 0;
}
