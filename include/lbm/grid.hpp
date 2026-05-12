#pragma once

#include "lbm/lattice.hpp"

#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace lbm {

struct Grid {
    int nx{};
    int ny{};
    // f 按单元优先顺序存储所有 D2Q9 分布函数：(y * nx + x) * q + direction。
    std::vector<double> f;
    std::vector<double> f_next;
    std::vector<double> rho;
    std::vector<double> ux;
    std::vector<double> uy;
    std::vector<std::uint8_t> solid;

    Grid(int width, int height)
        : nx(width),
          ny(height),
          f(static_cast<std::size_t>(width) * height * D2Q9::q, 0.0),
          f_next(f.size(), 0.0),
          rho(static_cast<std::size_t>(width) * height, 1.0),
          ux(rho.size(), 0.0),
          uy(rho.size(), 0.0),
          solid(rho.size(), 0) {
        if (width <= 0 || height <= 0) {
            throw std::invalid_argument("Grid dimensions must be positive.");
        }
    }

    int scalar_index(int x, int y) const {
        return y * nx + x;
    }

    int dist_index(int x, int y, int direction) const {
        return (y * nx + x) * D2Q9::q + direction;
    }

    // 周期边界案例在边缘回绕，而不是设置固体壁面。
    int wrap_x(int x) const {
        if (x < 0) {
            return x + nx;
        }
        if (x >= nx) {
            return x - nx;
        }
        return x;
    }

    int wrap_y(int y) const {
        if (y < 0) {
            return y + ny;
        }
        if (y >= ny) {
            return y - ny;
        }
        return y;
    }
};

} // namespace lbm
