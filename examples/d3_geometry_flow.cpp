#include "lbm/geometry_mask.hpp"
#include "lbm/geometry_mask3d.hpp"
#include "lbm/solver3d.hpp"
#include "lbm/svg_slice_animation3d.hpp"

#include <cmath>
#include <filesystem>
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

int main(int argc, char** argv) {
    try {
        const std::filesystem::path geometry_file =
            argc > 1 ? std::filesystem::path(argv[1])
                     : lbm::geometry_path("channel_obstacle.geom3d");
        const lbm::GeometryMask3D mask = lbm::GeometryMask3D::load(geometry_file);
        constexpr int steps = 260;
        constexpr int output_interval = 20;
        constexpr double inlet_velocity = 0.035;

        lbm::Solver3DConfig config;
        config.lattice_model = lbm::Lattice3DModel::D3Q19;
        config.tau = 0.84;
        lbm::Solver3D solver(
            mask.width(), mask.height(), mask.depth(), config);
        solver.initialize_masked_flow(mask, inlet_velocity);

        std::vector<lbm::SliceFrame3D> frames;
        frames.push_back(capture_speed(solver));
        for (int step = 1; step <= steps; ++step) {
            solver.step_masked_flow(inlet_velocity);
            if (step % output_interval == 0 || step == steps) {
                frames.push_back(capture_speed(solver));
            }
        }

        const lbm::Grid3D& grid = solver.grid();
        const lbm::SliceMasks3D masks = lbm::capture_orthogonal_masks(
            grid.nx,
            grid.ny,
            grid.nz,
            grid.nx / 2,
            grid.ny / 2,
            grid.nz / 2,
            [&grid](int x, int y, int z) {
                return grid.solid[grid.scalar_index(x, y, z)] != 0;
            },
            [&mask](int x, int y, int z) {
                return mask.cell_at(x, y, z) == lbm::GeometryCell::Porous ? 0.3
                                                                          : 1.0;
            });
        lbm::ScalarAnimationOptions animation;
        animation.cell = 6;
        animation.fps = 8.0;
        animation.title = "D3Q19 imported-voxel obstacle flow";
        animation.footer =
            "Three center slices; dark gray: solid; hatched: porous marker.";
        const auto path =
            lbm::result_path("d3_geometry_flow_speed_animation.svg");
        lbm::write_orthogonal_slice_animation_svg(
            path,
            frames,
            masks,
            grid.nx,
            grid.ny,
            grid.nz,
            animation);

        const auto diagnostics = solver.diagnostics();
        std::cout << "D3Q19 imported-geometry flow\n";
        std::cout << "geometry: " << geometry_file.string() << '\n';
        std::cout << "grid: " << grid.nx << " x " << grid.ny << " x "
                  << grid.nz << '\n';
        std::cout << "final max speed: " << diagnostics.max_speed << '\n';
        std::cout << "wrote: " << path.string() << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "3D geometry flow error: " << error.what() << '\n';
        return 1;
    }
}
