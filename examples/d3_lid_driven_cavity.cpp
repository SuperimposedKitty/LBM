#include "lbm/solver3d.hpp"

#include <iostream>
#include <stdexcept>

int main() {
    constexpr int nx = 30;
    constexpr int ny = 22;
    constexpr int nz = 18;
    constexpr int steps = 240;
    constexpr int output_interval = 20;
    constexpr double lid_velocity = 0.045;

    lbm::Solver3DConfig config;
    config.lattice_model = lbm::Lattice3DModel::D3Q19;
    config.collision_model = lbm::CollisionModel::MRT;
    config.tau = 0.82;
    lbm::Solver3D solver(nx, ny, nz, config);
    solver.initialize_lid_driven_cavity(lid_velocity);

    const auto before = solver.diagnostics();
    lbm::OpenFoamWriter foam(lbm::openfoam_result_path("d3_lid_driven_cavity"), solver.openfoam_snapshot());
    foam.write(0, solver.openfoam_snapshot());
    std::cout << "OpenFOAM case: " << (foam.directory() / "lbm.foam").string() << '\n';
    for (int step = 1; step <= steps; ++step) {
        solver.step_lid_driven_cavity(lid_velocity);
        if (step % output_interval == 0 || step == steps) {
            foam.write(step, solver.openfoam_snapshot());
        }
    }


    const auto after = solver.diagnostics();
    std::cout << "D3Q19 MRT lid-driven cavity\n";
    std::cout << "grid: " << nx << " x " << ny << " x " << nz << '\n';
    std::cout << "initial fluid mass: " << before.mass << '\n';
    std::cout << "final fluid mass: " << after.mass << '\n';
    std::cout << "final max speed: " << after.max_speed << '\n';
    return 0;
}
