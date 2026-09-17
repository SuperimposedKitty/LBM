#include "lbm/geometry_mask.hpp"
#include "lbm/solver.hpp"
#include "lbm/svg_animation.hpp"

#include <cmath>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <iostream>
#include <vector>

namespace {

std::vector<double> capture_speed(const lbm::Solver& solver) {
    const auto& grid = solver.grid();
    std::vector<double> frame;
    frame.reserve(static_cast<std::size_t>(grid.nx) * grid.ny);
    for (int y = 0; y < grid.ny; ++y) {
        for (int x = 0; x < grid.nx; ++x) {
            const int s = grid.scalar_index(x, y);
            frame.push_back(std::sqrt(grid.ux[s] * grid.ux[s] + grid.uy[s] * grid.uy[s]));
        }
    }
    return frame;
}

std::vector<std::uint8_t> capture_solid(const lbm::Solver& solver) {
    const auto& grid = solver.grid();
    return grid.solid;
}

} // namespace

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
        const auto solid = capture_solid(solver);
        std::vector<std::vector<double>> frames;
        frames.reserve(static_cast<std::size_t>(steps / output_interval + 2));
        frames.push_back(capture_speed(solver));
        lbm::OpenFoamWriter foam(lbm::openfoam_result_path("geometry_flow"), solver.openfoam_snapshot());
        foam.write(0, solver.openfoam_snapshot());
        std::cout << "OpenFOAM case: " << (foam.directory() / "lbm.foam").string() << '\n';

        for (int step = 1; step <= steps; ++step) {
            solver.step_masked_flow(inlet_velocity);
            if (step % output_interval == 0 || step == steps) {
                frames.push_back(capture_speed(solver));
                foam.write(step, solver.openfoam_snapshot());
            }
        }

        lbm::ScalarAnimationOptions animation;
        animation.cell = 6;
        animation.fps = 10.0;
        animation.title = "Imported geometry: single-phase obstacle flow";
        animation.footer = "Frames: " + std::to_string(frames.size()) + "; grid: " +
                           std::to_string(mask.width()) + " x " +
                           std::to_string(mask.height()) +
                           "; color: speed; dark gray: solid.";
        animation.color_map = lbm::ColorMap::Sequential;
        const auto animation_path = lbm::result_path("geometry_flow_speed_animation.svg");
        lbm::write_scalar_animation_svg(
            animation_path, frames, solid, mask.width(), mask.height(), animation);

        const auto diagnostics = solver.diagnostics();
        std::cout << "D2Q9 imported-geometry single-phase flow\n";
        std::cout << "geometry: " << geometry_file.string() << '\n';
        std::cout << "grid: " << mask.width() << " x " << mask.height() << '\n';
        std::cout << "steps: " << steps << '\n';
        std::cout << "final max speed: " << diagnostics.max_speed << '\n';
        std::cout << "wrote: " << animation_path.string() << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "geometry flow error: " << error.what() << '\n';
        return 1;
    }
}
