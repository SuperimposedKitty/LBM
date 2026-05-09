#include "lbm/svg_animation.hpp"
#include "lbm/two_phase_solver.hpp"

#include <cstdint>
#include <iostream>
#include <vector>

namespace {

std::vector<double> capture_phase(const lbm::TwoPhaseSolver& solver) {
    std::vector<double> frame;
    frame.reserve(static_cast<std::size_t>(solver.nx()) * solver.ny());
    for (int y = 0; y < solver.ny(); ++y) {
        for (int x = 0; x < solver.nx(); ++x) {
            frame.push_back(solver.phase_at(x, y));
        }
    }
    return frame;
}

std::vector<std::uint8_t> capture_solid(const lbm::TwoPhaseSolver& solver) {
    std::vector<std::uint8_t> solid;
    solid.reserve(static_cast<std::size_t>(solver.nx()) * solver.ny());
    for (int y = 0; y < solver.ny(); ++y) {
        for (int x = 0; x < solver.nx(); ++x) {
            solid.push_back(solver.solid_at(x, y) ? 1 : 0);
        }
    }
    return solid;
}

} // namespace

int main() {
    const int nx = 180;
    const int ny = 34;
    const int steps = 2500;
    const int output_interval = 125;

    lbm::TwoPhaseConfig config;
    config.tau_a = 1.0;
    config.tau_b = 1.0;
    config.rho_high = 1.0;
    config.rho_low = 0.02;
    config.interaction_strength = 3.0;
    config.inlet_velocity = 0.018;
    config.body_force_x = 8.0e-7;
    config.initial_interface_x = 28;

    lbm::TwoPhaseSolver solver(nx, ny, config);
    solver.initialize_capillary_displacement();

    const auto before = solver.diagnostics();
    const auto solid = capture_solid(solver);
    std::vector<std::vector<double>> frames;
    frames.reserve(static_cast<std::size_t>(steps / output_interval + 2));
    frames.push_back(capture_phase(solver));

    for (int step = 1; step <= steps; ++step) {
        solver.step();
        if (step % output_interval == 0 || step == steps) {
            frames.push_back(capture_phase(solver));
        }
    }

    lbm::ScalarAnimationOptions animation;
    animation.cell = 3;
    animation.fps = 8.0;
    animation.title = "Two-phase capillary displacement";
    animation.footer = "Frames: " + std::to_string(frames.size()) + "; grid: " +
                       std::to_string(nx) + " x " + std::to_string(ny) +
                       "; red: injected phase; blue: displaced phase; dark gray: wall.";
    animation.color_map = lbm::ColorMap::Phase;
    animation.fixed_range = true;
    animation.vmin = -1.0;
    animation.vmax = 1.0;
    const auto animation_path = lbm::result_path("capillary_phase_animation.svg");
    lbm::write_scalar_animation_svg(animation_path, frames, solid, nx, ny, animation);

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
    std::cout << "initial injected-fluid centroid x: " << before.interface_x << '\n';
    std::cout << "final injected-fluid centroid x: " << after.interface_x << '\n';
    std::cout << "injected-fluid centroid displacement: " << after.interface_x - before.interface_x
              << '\n';
    std::cout << "initial mass A: " << before.mass_a << '\n';
    std::cout << "final mass A: " << after.mass_a << '\n';
    std::cout << "initial mass B: " << before.mass_b << '\n';
    std::cout << "final mass B: " << after.mass_b << '\n';
    std::cout << "final max speed: " << after.max_speed << '\n';
    std::cout << "wrote: " << animation_path.string() << '\n';

    return 0;
}
