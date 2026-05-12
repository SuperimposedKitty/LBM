#include "lbm/svg_animation.hpp"
#include "lbm/two_phase_solver.hpp"

#include <algorithm>
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

double red_fraction_in_outlet_band(const lbm::TwoPhaseSolver& solver, int band_width) {
    double fraction_sum = 0.0;
    int samples = 0;
    const int x0 = std::max(0, solver.nx() - 1 - band_width);
    for (int y = 0; y < solver.ny(); ++y) {
        for (int x = x0; x < solver.nx() - 1; ++x) {
            if (solver.solid_at(x, y)) {
                continue;
            }
            fraction_sum += std::clamp(0.5 * (solver.phase_at(x, y) + 1.0), 0.0, 1.0);
            ++samples;
        }
    }
    return samples > 0 ? fraction_sum / static_cast<double>(samples) : 0.0;
}

} // namespace

int main() {
    const int nx = 200;
    const int ny = 36;
    const int max_steps = 18000;
    const int output_interval = 300;
    const int outlet_band_width = 14;
    const double completion_fraction = 0.78;

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

    lbm::TwoPhaseSolver solver(nx, ny, config);
    solver.initialize_capillary_displacement();

    const auto before = solver.diagnostics();
    const auto solid = capture_solid(solver);
    std::vector<std::vector<double>> frames;
    frames.reserve(static_cast<std::size_t>(max_steps / output_interval + 2));
    frames.push_back(capture_phase(solver));

    int completed_steps = 0;
    double outlet_red_fraction = red_fraction_in_outlet_band(solver, outlet_band_width);
    for (int step = 1; step <= max_steps; ++step) {
        solver.step();
        completed_steps = step;
        outlet_red_fraction = red_fraction_in_outlet_band(solver, outlet_band_width);
        const bool completed = outlet_red_fraction >= completion_fraction;
        if (step % output_interval == 0 || completed || step == max_steps) {
            frames.push_back(capture_phase(solver));
        }
        if (completed) {
            break;
        }
    }

    lbm::ScalarAnimationOptions animation;
    animation.cell = 3;
    animation.fps = 10.0;
    animation.title = "Two-phase capillary displacement with wall wetting";
    animation.footer = "Frames: " + std::to_string(frames.size()) + "; grid: " +
                       std::to_string(nx) + " x " + std::to_string(ny) +
                       "; contact angle: " +
                       std::to_string(static_cast<int>(config.contact_angle_degrees)) +
                       " deg; red: injected phase; blue: displaced phase; yellow: interface.";
    animation.color_map = lbm::ColorMap::Phase;
    animation.fixed_range = true;
    animation.vmin = -1.0;
    animation.vmax = 1.0;
    const auto animation_path = lbm::result_path("capillary_phase_animation.svg");
    lbm::write_scalar_animation_svg(animation_path, frames, solid, nx, ny, animation);

    const auto after = solver.diagnostics();
    std::cout << "D2Q9 two-phase capillary displacement example\n";
    std::cout << "grid: " << nx << " x " << ny << '\n';
    std::cout << "max steps: " << max_steps << '\n';
    std::cout << "completed steps: " << completed_steps << '\n';
    std::cout << "output interval: " << output_interval << '\n';
    std::cout << "tau_a: " << config.tau_a << '\n';
    std::cout << "tau_b: " << config.tau_b << '\n';
    std::cout << "interaction strength: " << config.interaction_strength << '\n';
    std::cout << "inlet velocity: " << config.inlet_velocity << '\n';
    std::cout << "body force x: " << config.body_force_x << '\n';
    std::cout << "contact angle: " << config.contact_angle_degrees << '\n';
    std::cout << "wall adhesion strength: " << config.wall_adhesion_strength << '\n';
    std::cout << "outlet red fraction: " << outlet_red_fraction << '\n';
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
