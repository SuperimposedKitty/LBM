#include "lbm/geometry_mask.hpp"
#include "lbm/geometry_mask3d.hpp"
#include "lbm/solver3d.hpp"

#include <filesystem>
#include <iostream>
#include <stdexcept>

int main(int argc, char** argv) {
    try {
        const std::filesystem::path geometry_file =
            argc > 1 ? std::filesystem::path(argv[1])
                     : lbm::geometry_path("channel_obstacle.geom3d");
        const lbm::GeometryMask3D mask = lbm::GeometryMask3D::load(geometry_file);
        constexpr int steps = 260;
        constexpr int output_interval = 20;
        constexpr double inlet_velocity = 0.035;

        lbm::Solver3DConfig config;
        config.lattice_model = lbm::Lattice3DModel::D3Q19;
        config.tau = 0.84;
        lbm::Solver3D solver(
            mask.width(), mask.height(), mask.depth(), config);
        solver.initialize_masked_flow(mask, inlet_velocity);

        lbm::OpenFoamWriter foam(lbm::openfoam_result_path("d3_geometry_flow"), solver.openfoam_snapshot());
        foam.write(0, solver.openfoam_snapshot());
        std::cout << "OpenFOAM case: " << (foam.directory() / "lbm.foam").string() << '\n';
        for (int step = 1; step <= steps; ++step) {
            solver.step_masked_flow(inlet_velocity);
            if (step % output_interval == 0 || step == steps) {
                foam.write(step, solver.openfoam_snapshot());
            }
        }

        const lbm::Grid3D& grid = solver.grid();

        const auto diagnostics = solver.diagnostics();
        std::cout << "D3Q19 imported-geometry flow\n";
        std::cout << "geometry: " << geometry_file.string() << '\n';
        std::cout << "grid: " << grid.nx << " x " << grid.ny << " x "
                  << grid.nz << '\n';
        std::cout << "final max speed: " << diagnostics.max_speed << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "3D geometry flow error: " << error.what() << '\n';
        return 1;
    }
}
