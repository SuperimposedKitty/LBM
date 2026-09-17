#include "lbm/solver.hpp"
#include "lbm/svg_animation.hpp"

#include <cmath>
#include <cstdint>
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
    std::vector<std::uint8_t> solid;
    solid.reserve(static_cast<std::size_t>(grid.nx) * grid.ny);
    for (int y = 0; y < grid.ny; ++y) {
        for (int x = 0; x < grid.nx; ++x) {
            solid.push_back(grid.solid[grid.scalar_index(x, y)]);
        }
    }
    return solid;
}

} // namespace

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
    const auto solid = capture_solid(solver);
    std::vector<std::vector<double>> frames;
    frames.reserve(static_cast<std::size_t>(steps / output_interval + 2));
    frames.push_back(capture_speed(solver));
    lbm::OpenFoamWriter foam(lbm::openfoam_result_path("lid_driven_cavity"), solver.openfoam_snapshot());
    foam.write(0, solver.openfoam_snapshot());
    std::cout << "OpenFOAM case: " << (foam.directory() / "lbm.foam").string() << '\n';

    for (int step = 1; step <= steps; ++step) {
        solver.step_lid_driven_cavity(lid_velocity);
        if (step % output_interval == 0 || step == steps) {
            frames.push_back(capture_speed(solver));
            foam.write(step, solver.openfoam_snapshot());
        }
    }

    lbm::ScalarAnimationOptions animation;
    animation.cell = 5;
    animation.fps = 8.0;
    animation.title = "Lid-driven cavity speed";
    animation.footer = "Frames: " + std::to_string(frames.size()) + "; grid: " +
                       std::to_string(nx) + " x " + std::to_string(ny) +
                       "; color: speed; dark gray: wall.";
    animation.color_map = lbm::ColorMap::Sequential;
    const auto animation_path = lbm::result_path("lid_driven_cavity_speed_animation.svg");
    lbm::write_scalar_animation_svg(animation_path, frames, solid, nx, ny, animation);

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
    std::cout << "wrote: " << animation_path.string() << '\n';

    return 0;
}
