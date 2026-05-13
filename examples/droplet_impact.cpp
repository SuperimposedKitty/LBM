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
    const int nx = 128;
    const int ny = 104;
    const int steps = 6400;
    const int output_interval = 120;
    const double droplet_radius = 15.0;
    const double droplet_center_x = 0.5 * static_cast<double>(nx - 1);
    const double droplet_center_y = 34.0;
    const double initial_droplet_uy = -0.085;

    lbm::TwoPhaseConfig config;
    config.tau_a = 1.0;
    config.tau_b = 1.0;
    config.rho_high = 1.0;
    config.rho_low = 0.02;
    config.interaction_strength = 3.2;
    config.body_force_x = 0.0;
    config.body_force_y = 0.0;
    config.contact_angle_degrees = 165.0;
    config.wall_adhesion_strength = 0.13;

    lbm::TwoPhaseSolver solver(nx, ny, config);
    solver.initialize_droplet_impact(
        droplet_center_x, droplet_center_y, droplet_radius, 0.0, initial_droplet_uy);

    const auto before = solver.diagnostics();
    const auto solid = capture_solid(solver);
    std::vector<std::vector<double>> frames;
    frames.reserve(static_cast<std::size_t>(steps / output_interval + 2));
    frames.push_back(capture_phase(solver));

    double min_centroid_y = before.phase_a_centroid_y;
    int min_centroid_step = 0;
    for (int step = 1; step <= steps; ++step) {
        solver.step_closed();
        const auto current = solver.diagnostics();
        if (current.phase_a_centroid_y < min_centroid_y) {
            min_centroid_y = current.phase_a_centroid_y;
            min_centroid_step = step;
        }
        if (step % output_interval == 0 || step == steps) {
            frames.push_back(capture_phase(solver));
        }
    }

    lbm::ScalarAnimationOptions animation;
    animation.cell = 4;
    animation.fps = 12.0;
    animation.title = "Droplet impact and rebound";
    animation.footer = "Frames: " + std::to_string(frames.size()) + "; grid: " +
                       std::to_string(nx) + " x " + std::to_string(ny) +
                       "; red: droplet; blue: ambient phase; yellow: interface; gray: solid.";
    animation.color_map = lbm::ColorMap::Phase;
    animation.fixed_range = true;
    animation.vmin = -1.0;
    animation.vmax = 1.0;
    const auto animation_path = lbm::result_path("droplet_impact_animation.svg");
    lbm::write_scalar_animation_svg(animation_path, frames, solid, nx, ny, animation);

    const auto after = solver.diagnostics();
    std::cout << "D2Q9 two-phase droplet impact example\n";
    std::cout << "grid: " << nx << " x " << ny << '\n';
    std::cout << "steps: " << steps << '\n';
    std::cout << "output interval: " << output_interval << '\n';
    std::cout << "droplet radius: " << droplet_radius << '\n';
    std::cout << "initial droplet center y: " << before.phase_a_centroid_y << '\n';
    std::cout << "minimum droplet center y: " << min_centroid_y << '\n';
    std::cout << "minimum center step: " << min_centroid_step << '\n';
    std::cout << "final droplet center y: " << after.phase_a_centroid_y << '\n';
    std::cout << "rebound height after minimum: " << after.phase_a_centroid_y - min_centroid_y
              << '\n';
    std::cout << "initial vertical velocity: " << initial_droplet_uy << '\n';
    std::cout << "body force y: " << config.body_force_y << '\n';
    std::cout << "contact angle: " << config.contact_angle_degrees << '\n';
    std::cout << "wall adhesion strength: " << config.wall_adhesion_strength << '\n';
    std::cout << "initial mass A: " << before.mass_a << '\n';
    std::cout << "final mass A: " << after.mass_a << '\n';
    std::cout << "initial mass B: " << before.mass_b << '\n';
    std::cout << "final mass B: " << after.mass_b << '\n';
    std::cout << "final max speed: " << after.max_speed << '\n';
    std::cout << "wrote: " << animation_path.string() << '\n';

    return 0;
}
