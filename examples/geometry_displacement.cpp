#include "lbm/geometry_mask.hpp"
#include "lbm/svg_animation.hpp"
#include "lbm/two_phase_solver.hpp"

#include <cstdint>
#include <exception>
#include <filesystem>
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

int main(int argc, char* argv[]) {
    try {
        if (argc > 2) {
            throw std::invalid_argument("Usage: lbm_geometry_displacement.exe [geometry.geom]");
        }
        const std::filesystem::path geometry_file =
            argc == 2 ? std::filesystem::path(argv[1])
                      : lbm::geometry_path("channel_obstacle.geom");
        const lbm::GeometryMask mask = lbm::GeometryMask::load(geometry_file);

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
        config.free_flow_porosity = 1.0;
        config.porous_porosity = 0.3;
        config.porous_pore_diameter = 40.0;
        config.darcy_drag_scale = 0.06;
        config.forchheimer_drag_scale = 0.02;

        lbm::TwoPhaseSolver solver(mask.width(), mask.height(), config);
        solver.initialize_geometry_displacement(mask);
        const auto before = solver.diagnostics();
        const auto solid = capture_solid(solver);
        const auto porosity = capture_porosity(solver);
        std::vector<std::vector<double>> frames;
        frames.reserve(static_cast<std::size_t>(steps / output_interval + 2));
        frames.push_back(capture_phase(solver));
        lbm::OpenFoamWriter foam(lbm::openfoam_result_path("geometry_displacement"), solver.openfoam_snapshot());
        foam.write(0, solver.openfoam_snapshot());
        std::cout << "OpenFOAM case: " << (foam.directory() / "lbm.foam").string() << '\n';

        for (int step = 1; step <= steps; ++step) {
            solver.step();
            if (step % output_interval == 0 || step == steps) {
                frames.push_back(capture_phase(solver));
                foam.write(step, solver.openfoam_snapshot());
            }
        }

        lbm::ScalarAnimationOptions animation;
        animation.cell = 6;
        animation.fps = 10.0;
        animation.title = "Imported geometry: two-phase displacement";
        animation.footer = "Frames: " + std::to_string(frames.size()) + "; grid: " +
                           std::to_string(mask.width()) + " x " +
                           std::to_string(mask.height()) +
                           "; red: injected A; blue: displaced B; hatched: porous.";
        animation.color_map = lbm::ColorMap::Phase;
        animation.fixed_range = true;
        animation.vmin = -1.0;
        animation.vmax = 1.0;
        const auto animation_path = lbm::result_path("geometry_displacement_animation.svg");
        lbm::write_scalar_animation_svg(
            animation_path,
            frames,
            solid,
            mask.width(),
            mask.height(),
            animation,
            porosity);

        const auto after = solver.diagnostics();
        std::cout << "D2Q9 imported-geometry two-phase displacement\n";
        std::cout << "geometry: " << geometry_file.string() << '\n';
        std::cout << "grid: " << mask.width() << " x " << mask.height() << '\n';
        std::cout << "steps: " << steps << '\n';
        std::cout << "initial A centroid x: " << before.interface_x << '\n';
        std::cout << "final A centroid x: " << after.interface_x << '\n';
        std::cout << "final max speed: " << after.max_speed << '\n';
        std::cout << "wrote: " << animation_path.string() << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "geometry displacement error: " << error.what() << '\n';
        return 1;
    }
}
