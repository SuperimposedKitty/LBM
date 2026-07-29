#include "lbm/solver3d.hpp"
#include "lbm/svg_slice_animation3d.hpp"

#include <cmath>
#include <cstdint>
#include <iostream>
#include <vector>

namespace {

lbm::SliceFrame3D capture_speed(const lbm::Solver3D& solver) {
    const lbm::Grid3D& grid = solver.grid();
    return lbm::capture_orthogonal_slices(
        grid.nx,
        grid.ny,
        grid.nz,
        grid.nx / 2,
        grid.ny / 2,
        grid.nz / 2,
        [&grid](int x, int y, int z) {
            const std::size_t s = grid.scalar_index(x, y, z);
            return std::sqrt(
                grid.ux[s] * grid.ux[s] + grid.uy[s] * grid.uy[s] +
                grid.uz[s] * grid.uz[s]);
        });
}

} // namespace

int main() {
    constexpr int nx = 30;
    constexpr int ny = 22;
    constexpr int nz = 18;
    constexpr int steps = 240;
    constexpr int output_interval = 20;
    constexpr double lid_velocity = 0.045;

    lbm::Solver3DConfig config;
    config.lattice_model = lbm::Lattice3DModel::D3Q19;
    config.collision_model = lbm::CollisionModel::MRT;
    config.tau = 0.82;
    lbm::Solver3D solver(nx, ny, nz, config);
    solver.initialize_lid_driven_cavity(lid_velocity);

    const auto before = solver.diagnostics();
    std::vector<lbm::SliceFrame3D> frames;
    frames.push_back(capture_speed(solver));
    for (int step = 1; step <= steps; ++step) {
        solver.step_lid_driven_cavity(lid_velocity);
        if (step % output_interval == 0 || step == steps) {
            frames.push_back(capture_speed(solver));
        }
    }

    const lbm::Grid3D& grid = solver.grid();
    const lbm::SliceMasks3D masks = lbm::capture_orthogonal_masks(
        nx,
        ny,
        nz,
        nx / 2,
        ny / 2,
        nz / 2,
        [&grid](int x, int y, int z) {
            return grid.solid[grid.scalar_index(x, y, z)] != 0;
        },
        [](int, int, int) { return 1.0; });
    lbm::ScalarAnimationOptions animation;
    animation.cell = 5;
    animation.fps = 8.0;
    animation.title = "D3Q19 MRT lid-driven cavity";
    animation.footer =
        "Three orthogonal center slices; color: speed; dark gray: solid wall.";
    const auto path = lbm::result_path("d3_lid_driven_cavity_speed_animation.svg");
    lbm::write_orthogonal_slice_animation_svg(
        path, frames, masks, nx, ny, nz, animation);

    const auto after = solver.diagnostics();
    std::cout << "D3Q19 MRT lid-driven cavity\n";
    std::cout << "grid: " << nx << " x " << ny << " x " << nz << '\n';
    std::cout << "initial fluid mass: " << before.mass << '\n';
    std::cout << "final fluid mass: " << after.mass << '\n';
    std::cout << "final max speed: " << after.max_speed << '\n';
    std::cout << "wrote: " << path.string() << '\n';
    return 0;
}
