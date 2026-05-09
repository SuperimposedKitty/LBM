#include "lbm/solver.hpp"
#include "lbm/svg_animation.hpp"

#include <cmath>
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

} // namespace

int main() {
    const int nx = 128;
    const int ny = 64;
    const int steps = 1000;
    const int output_interval = 50;

    lbm::SolverConfig config;
    config.tau = 0.8;
    config.initial_rho = 1.0;

    lbm::Solver solver(nx, ny, config);
    solver.initialize_shear_wave(0.05, 1);

    const auto before = solver.diagnostics();
    std::vector<std::vector<double>> frames;
    frames.reserve(static_cast<std::size_t>(steps / output_interval + 2));
    frames.push_back(capture_speed(solver));

    for (int step = 1; step <= steps; ++step) {
        solver.step();
        if (step % output_interval == 0 || step == steps) {
            frames.push_back(capture_speed(solver));
        }
    }

    lbm::ScalarAnimationOptions animation;
    animation.cell = 4;
    animation.fps = 8.0;
    animation.title = "Periodic shear wave speed";
    animation.footer = "Frames: " + std::to_string(frames.size()) + "; grid: " +
                       std::to_string(nx) + " x " + std::to_string(ny) +
                       "; color: speed.";
    animation.color_map = lbm::ColorMap::Sequential;
    const auto animation_path = lbm::result_path("periodic_shear_speed_animation.svg");
    lbm::write_scalar_animation_svg(animation_path, frames, {}, nx, ny, animation);

    const auto after = solver.diagnostics();
    std::cout << "D2Q9 periodic shear wave example\n";
    std::cout << "grid: " << nx << " x " << ny << '\n';
    std::cout << "steps: " << steps << '\n';
    std::cout << "output interval: " << output_interval << '\n';
    std::cout << "tau: " << config.tau << '\n';
    std::cout << "initial mass: " << before.mass << '\n';
    std::cout << "final mass: " << after.mass << '\n';
    std::cout << "mass difference: " << after.mass - before.mass << '\n';
    std::cout << "final max speed: " << after.max_speed << '\n';
    std::cout << "wrote: " << animation_path.string() << '\n';

    return 0;
}
