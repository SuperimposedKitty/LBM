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

std::vector<double> capture_porosity(const lbm::TwoPhaseSolver& solver) {
    std::vector<double> porosity;
    porosity.reserve(static_cast<std::size_t>(solver.nx()) * solver.ny());
    for (int y = 0; y < solver.ny(); ++y) {
        for (int x = 0; x < solver.nx(); ++x) {
            porosity.push_back(solver.porosity_at(x, y));
        }
    }
    return porosity;
}

} // namespace

int main() {
    const int nx = 200;
    const int ny = 36;
    const int steps = 18000;
    const int output_interval = 300;

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
    config.porous_start_x = 78;
    config.porous_end_x = 122;
    config.porous_porosity = 0.3;
    config.free_flow_porosity = 1.0;
    config.porous_pore_diameter = 40.0;
    config.darcy_drag_scale = 0.06;
    config.forchheimer_drag_scale = 0.02;

    lbm::TwoPhaseSolver solver(nx, ny, config);
    solver.initialize_capillary_displacement();

    const auto before = solver.diagnostics();
    const auto solid = capture_solid(solver);
    const auto porosity = capture_porosity(solver);
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
    animation.fps = 10.0;
    animation.title = "Two-phase capillary displacement with wall wetting";
    animation.footer = "Frames: " + std::to_string(frames.size()) + "; grid: " +
                       std::to_string(nx) + " x " + std::to_string(ny) +
                       "; contact angle: " +
                       std::to_string(static_cast<int>(config.contact_angle_degrees)) +
                       " deg; porous eps: " + std::to_string(config.porous_porosity) +
                       "; hatched: porous medium.";
    animation.color_map = lbm::ColorMap::Phase;
    animation.fixed_range = true;
    animation.vmin = -1.0;
    animation.vmax = 1.0;
    const auto animation_path = lbm::result_path("capillary_phase_animation.svg");
    lbm::write_scalar_animation_svg(animation_path, frames, solid, nx, ny, animation, porosity);

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
    std::cout << "contact angle: " << config.contact_angle_degrees << '\n';
    std::cout << "wall adhesion strength: " << config.wall_adhesion_strength << '\n';
    std::cout << "porous medium x range: [" << config.porous_start_x << ", "
              << config.porous_end_x << ")\n";
    std::cout << "porous porosity: " << config.porous_porosity << '\n';
    std::cout << "porous pore diameter: " << config.porous_pore_diameter << '\n';
    std::cout << "Darcy drag scale: " << config.darcy_drag_scale << '\n';
    std::cout << "Forchheimer drag scale: " << config.forchheimer_drag_scale << '\n';
    std::cout << "initial injected-fluid centroid x: " << before.interface_x << '\n';
    std::cout << "final injected-fluid centroid x: " << after.interface_x << '\n';
    std::cout << "injected-fluid centroid displacement: " << after.interface_x - before.interface_x
              << '\n';
    std::cout << "initial mass A: " << before.mass_a << '\n';
    std::cout << "final mass A: " << after.mass_a << '\n';
    std::cout << "initial mass B: " << before.mass_b << '\n';
    std::cout << "final mass B: " << after.mass_b << '\n';
    std::cout << "final max speed: " << after.max_speed << '\n';
    std::cout << "porous mean pore speed: " << after.porous_mean_pore_speed << '\n';
    std::cout << "porous max pore speed: " << after.porous_max_pore_speed << '\n';
    std::cout << "wrote: " << animation_path.string() << '\n';

    return 0;
}
