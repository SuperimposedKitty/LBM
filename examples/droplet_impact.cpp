#include "lbm/svg_animation.hpp"
#include "lbm/two_phase_solver.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>
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

struct DropletShapeStats {
    double marker_mass = 0.0;
    double marker_centroid_y = 0.0;
    int spread_width = 0;
    int contact_width = 0;
    bool valid = false;
};

struct DropletDiagnosticRecord {
    int step = 0;
    double mass_a = 0.0;
    double centroid_y = 0.0;
    int spread_width = 0;
    int contact_width = 0;
};

DropletShapeStats droplet_shape_stats(
    const lbm::TwoPhaseSolver& solver, int bottom_wall_thickness) {
    constexpr double core_phase_threshold = 0.0;
    constexpr double contact_phase_threshold = -0.35;

    DropletShapeStats stats;
    double weighted_y = 0.0;
    int contact_min_x = solver.nx();
    int contact_max_x = -1;
    const int contact_top_y = std::min(bottom_wall_thickness + 3, solver.ny() - 1);

    for (int y = bottom_wall_thickness; y < solver.ny() - 1; ++y) {
        int row_min_x = solver.nx();
        int row_max_x = -1;
        for (int x = 1; x < solver.nx() - 1; ++x) {
            if (solver.solid_at(x, y)) {
                continue;
            }
            const double phase = solver.phase_at(x, y);
            if (phase > core_phase_threshold) {
                row_min_x = std::min(row_min_x, x);
                row_max_x = std::max(row_max_x, x);
                stats.marker_mass += phase;
                weighted_y += static_cast<double>(y) * phase;
            }
            if (y < contact_top_y && phase > contact_phase_threshold) {
                contact_min_x = std::min(contact_min_x, x);
                contact_max_x = std::max(contact_max_x, x);
            }
        }
        if (row_max_x >= row_min_x) {
            stats.spread_width = std::max(stats.spread_width, row_max_x - row_min_x + 1);
        }
    }

    stats.valid = stats.marker_mass > 0.0;
    stats.marker_centroid_y = stats.valid ? weighted_y / stats.marker_mass : 0.0;
    stats.contact_width = contact_max_x >= contact_min_x ? contact_max_x - contact_min_x + 1 : 0;
    return stats;
}

} // namespace

int main() {
    const int nx = 132;
    const int ny = 112;
    const int steps = 4200;
    const int output_interval = 150;
    const int bottom_wall_thickness = 4;
    const double droplet_radius = 22.0;
    const double droplet_center_x = 0.5 * static_cast<double>(nx - 1);
    const double droplet_center_y = 36.0;
    const double initial_droplet_uy = -0.095;

    lbm::TwoPhaseConfig config;
    config.tau_a = 0.9;
    config.tau_b = 0.9;
    config.rho_high = 1.0;
    config.rho_low = 0.02;
    config.interaction_strength = 3.6;
    config.body_force_x = 0.0;
    config.body_force_y = 0.0;
    config.contact_angle_degrees = 178.0;
    config.wall_adhesion_strength = 0.12;
    config.bottom_wall_thickness = bottom_wall_thickness;
    config.droplet_interface_width = 1.25;
    config.recoloring_strength = 0.35;
    config.bottom_wall_repulsion_strength = 1.2e-3;
    config.bottom_wall_repulsion_range = 12;

    lbm::TwoPhaseSolver solver(nx, ny, config);
    solver.initialize_droplet_impact(
        droplet_center_x, droplet_center_y, droplet_radius, 0.0, initial_droplet_uy);

    const auto before = solver.diagnostics();
    const auto solid = capture_solid(solver);
    std::vector<std::vector<double>> frames;
    frames.reserve(static_cast<std::size_t>(steps / output_interval + 2));
    frames.push_back(capture_phase(solver));
    lbm::OpenFoamWriter foam(lbm::openfoam_result_path("droplet_impact"), solver.openfoam_snapshot());
    foam.write(0, solver.openfoam_snapshot());
    std::cout << "OpenFOAM case: " << (foam.directory() / "lbm.foam").string() << '\n';

    std::vector<DropletDiagnosticRecord> diagnostic_records;
    diagnostic_records.reserve(frames.capacity());
    int max_spread_width = 0;
    int max_spread_step = 0;
    int max_contact_width = 0;
    int max_contact_step = 0;
    bool sampled_shape_warning = false;
    bool sampled_mass_warning = false;
    bool centroid_warning = false;
    const double mass_warning_threshold =
        1.0e-6 * std::max(1.0, std::abs(before.mass_a));
    const auto record_sample = [&](int step, const lbm::TwoPhaseDiagnostics& diagnostics) {
        const auto shape = droplet_shape_stats(solver, bottom_wall_thickness);
        diagnostic_records.push_back(DropletDiagnosticRecord{
            step,
            diagnostics.mass_a,
            diagnostics.phase_a_centroid_y,
            shape.spread_width,
            shape.contact_width,
        });
        if (!shape.valid) {
            sampled_shape_warning = true;
        }
        if (std::abs(diagnostics.mass_a - before.mass_a) > mass_warning_threshold) {
            sampled_mass_warning = true;
        }
        if (!std::isfinite(diagnostics.phase_a_centroid_y)) {
            centroid_warning = true;
        }
        if (shape.spread_width > max_spread_width) {
            max_spread_width = shape.spread_width;
            max_spread_step = step;
        }
        if (shape.contact_width > max_contact_width) {
            max_contact_width = shape.contact_width;
            max_contact_step = step;
        }
    };
    record_sample(0, before);

    double min_centroid_y = before.phase_a_centroid_y;
    int min_centroid_step = 0;
    double max_centroid_after_min = before.phase_a_centroid_y;
    int max_centroid_after_min_step = 0;
    for (int step = 1; step <= steps; ++step) {
        solver.step_closed();
        const auto current = solver.diagnostics();
        if (current.phase_a_centroid_y < min_centroid_y) {
            min_centroid_y = current.phase_a_centroid_y;
            min_centroid_step = step;
            max_centroid_after_min = current.phase_a_centroid_y;
            max_centroid_after_min_step = step;
        } else if (step > min_centroid_step && current.phase_a_centroid_y > max_centroid_after_min) {
            max_centroid_after_min = current.phase_a_centroid_y;
            max_centroid_after_min_step = step;
        }
        if (step % output_interval == 0 || step == steps) {
            frames.push_back(capture_phase(solver));
            foam.write(step, solver.openfoam_snapshot());
            record_sample(step, current);
        }
    }

    lbm::ScalarAnimationOptions animation;
    animation.cell = 4;
    animation.fps = 12.0;
    animation.title = "Droplet impact and rebound";
    animation.footer = "Frames: " + std::to_string(frames.size()) + "; grid: " +
                       std::to_string(nx) + " x " + std::to_string(ny) +
                       "; red: droplet; blue: ambient phase; yellow: interface; gray: solid.";
    animation.color_map = lbm::ColorMap::PhaseHighContrast;
    animation.fixed_range = true;
    animation.vmin = -1.0;
    animation.vmax = 1.0;
    animation.highlight_solid_surface = true;
    animation.solid_surface_y = bottom_wall_thickness;
    animation.solid_surface_thickness = 4;
    const auto animation_path = lbm::result_path("droplet_impact_animation.svg");
    lbm::write_scalar_animation_svg(animation_path, frames, solid, nx, ny, animation);

    const auto after = solver.diagnostics();
    const double final_rebound_height = after.phase_a_centroid_y - min_centroid_y;
    const double max_rebound_height = max_centroid_after_min - min_centroid_y;
    const double final_mass_error = after.mass_a - before.mass_a;
    const double final_mass_relative_error =
        final_mass_error / std::max(1.0, std::abs(before.mass_a));
    std::cout << "D2Q9 two-phase droplet impact example\n";
    std::cout << "grid: " << nx << " x " << ny << '\n';
    std::cout << "steps: " << steps << '\n';
    std::cout << "output interval: " << output_interval << '\n';
    std::cout << "diagnostic samples: " << diagnostic_records.size() << '\n';
    std::cout << "droplet radius: " << droplet_radius << '\n';
    std::cout << "bottom wall thickness: " << bottom_wall_thickness << '\n';
    std::cout << "droplet interface width: " << config.droplet_interface_width << '\n';
    std::cout << "initial droplet center y: " << before.phase_a_centroid_y << '\n';
    std::cout << "minimum droplet center y: " << min_centroid_y << '\n';
    std::cout << "minimum center step: " << min_centroid_step << '\n';
    std::cout << "maximum rebound center y: " << max_centroid_after_min << '\n';
    std::cout << "maximum rebound center step: " << max_centroid_after_min_step << '\n';
    std::cout << "final droplet center y: " << after.phase_a_centroid_y << '\n';
    std::cout << "maximum rebound height after minimum: " << max_rebound_height << '\n';
    std::cout << "final rebound height after minimum: " << final_rebound_height << '\n';
    std::cout << "maximum spread width: " << max_spread_width << '\n';
    std::cout << "maximum spread step: " << max_spread_step << '\n';
    std::cout << "maximum near-wall contact width: " << max_contact_width << '\n';
    std::cout << "maximum near-wall contact step: " << max_contact_step << '\n';
    std::cout << "initial vertical velocity: " << initial_droplet_uy << '\n';
    std::cout << "body force y: " << config.body_force_y << '\n';
    std::cout << "contact angle: " << config.contact_angle_degrees << '\n';
    std::cout << "wall adhesion strength: " << config.wall_adhesion_strength << '\n';
    std::cout << "initial mass A: " << before.mass_a << '\n';
    std::cout << "final mass A: " << after.mass_a << '\n';
    std::cout << "final mass A error: " << final_mass_error << '\n';
    std::cout << "final mass A relative error: " << final_mass_relative_error << '\n';
    std::cout << "initial mass B: " << before.mass_b << '\n';
    std::cout << "final mass B: " << after.mass_b << '\n';
    std::cout << "final max speed: " << after.max_speed << '\n';
    if (sampled_shape_warning) {
        std::cout << "warning: sampled droplet core became invisible in the phase field.\n";
    }
    if (sampled_mass_warning || std::abs(final_mass_error) > mass_warning_threshold) {
        std::cout << "warning: A phase mass changed more than the configured threshold.\n";
    }
    if (centroid_warning || !std::isfinite(after.phase_a_centroid_y)) {
        std::cout << "warning: droplet centroid became non-finite.\n";
    }
    if (max_contact_width == 0) {
        std::cout << "warning: no near-wall contact was captured at output samples.\n";
    }
    if (max_rebound_height <= 0.0) {
        std::cout << "warning: rebound was not detected after the minimum center height.\n";
    }
    std::cout << "wrote: " << animation_path.string() << '\n';

    return 0;
}
