#pragma once

#include "lbm/lattice3d.hpp"

#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

namespace lbm {
namespace detail {

inline std::size_t checked_grid3d_cells(int nx, int ny, int nz) {
    if (nx <= 0 || ny <= 0 || nz <= 0) {
        throw std::invalid_argument("Grid3D dimensions must be positive.");
    }
    const std::size_t sx = static_cast<std::size_t>(nx);
    const std::size_t sy = static_cast<std::size_t>(ny);
    const std::size_t sz = static_cast<std::size_t>(nz);
    if (sx > std::numeric_limits<std::size_t>::max() / sy ||
        sx * sy > std::numeric_limits<std::size_t>::max() / sz) {
        throw std::length_error("Grid3D dimensions overflow addressable memory.");
    }
    return sx * sy * sz;
}

} // namespace detail

struct Grid3D {
    int nx{};
    int ny{};
    int nz{};
    Lattice3DModel lattice_model{Lattice3DModel::D3Q19};
    int q{};
    std::vector<double> f;
    std::vector<double> f_next;
    std::vector<double> rho;
    std::vector<double> ux;
    std::vector<double> uy;
    std::vector<double> uz;
    std::vector<std::uint8_t> solid;

    Grid3D(int width, int height, int depth, Lattice3DModel model)
        : nx(width),
          ny(height),
          nz(depth),
          lattice_model(model),
          q(lattice3d(model).q),
          f(detail::checked_grid3d_cells(width, height, depth) *
                static_cast<std::size_t>(q),
            0.0),
          f_next(f.size(), 0.0),
          rho(detail::checked_grid3d_cells(width, height, depth), 1.0),
          ux(rho.size(), 0.0),
          uy(rho.size(), 0.0),
          uz(rho.size(), 0.0),
          solid(rho.size(), 0) {}

    std::size_t scalar_index(int x, int y, int z) const {
        return (static_cast<std::size_t>(z) * ny + y) * nx + x;
    }

    std::size_t dist_index(int x, int y, int z, int direction) const {
        return scalar_index(x, y, z) * static_cast<std::size_t>(q) + direction;
    }

    int wrap_x(int x) const {
        return x < 0 ? x + nx : (x >= nx ? x - nx : x);
    }

    int wrap_y(int y) const {
        return y < 0 ? y + ny : (y >= ny ? y - ny : y);
    }

    int wrap_z(int z) const {
        return z < 0 ? z + nz : (z >= nz ? z - nz : z);
    }
};

} // namespace lbm
