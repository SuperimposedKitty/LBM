#pragma once

#include <array>
#include <stdexcept>

namespace lbm {

enum class Lattice3DModel {
    D3Q19,
    D3Q27
};

struct Lattice3DDescriptor {
    static constexpr int max_q = 27;
    static constexpr double cs2 = 1.0 / 3.0;

    int q{};
    std::array<int, max_q> cx{};
    std::array<int, max_q> cy{};
    std::array<int, max_q> cz{};
    std::array<double, max_q> w{};
    std::array<int, max_q> opposite{};
};

inline constexpr Lattice3DDescriptor d3q19{
    19,
    {0, 1, -1, 0, 0, 0, 0, 1, -1, 1, -1, 1, -1, 1, -1, 0, 0, 0, 0},
    {0, 0, 0, 1, -1, 0, 0, 1, 1, -1, -1, 0, 0, 0, 0, 1, -1, 1, -1},
    {0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0, 1, 1, -1, -1, 1, 1, -1, -1},
    {1.0 / 3.0,
     1.0 / 18.0,
     1.0 / 18.0,
     1.0 / 18.0,
     1.0 / 18.0,
     1.0 / 18.0,
     1.0 / 18.0,
     1.0 / 36.0,
     1.0 / 36.0,
     1.0 / 36.0,
     1.0 / 36.0,
     1.0 / 36.0,
     1.0 / 36.0,
     1.0 / 36.0,
     1.0 / 36.0,
     1.0 / 36.0,
     1.0 / 36.0,
     1.0 / 36.0,
     1.0 / 36.0},
    {0, 2, 1, 4, 3, 6, 5, 10, 9, 8, 7, 14, 13, 12, 11, 18, 17, 16, 15}};

inline constexpr Lattice3DDescriptor d3q27{
    27,
    {0,
     1,
     -1,
     0,
     0,
     0,
     0,
     1,
     -1,
     1,
     -1,
     1,
     -1,
     1,
     -1,
     0,
     0,
     0,
     0,
     1,
     -1,
     1,
     -1,
     1,
     -1,
     1,
     -1},
    {0,
     0,
     0,
     1,
     -1,
     0,
     0,
     1,
     1,
     -1,
     -1,
     0,
     0,
     0,
     0,
     1,
     -1,
     1,
     -1,
     1,
     1,
     -1,
     -1,
     1,
     1,
     -1,
     -1},
    {0,
     0,
     0,
     0,
     0,
     1,
     -1,
     0,
     0,
     0,
     0,
     1,
     1,
     -1,
     -1,
     1,
     1,
     -1,
     -1,
     1,
     1,
     1,
     1,
     -1,
     -1,
     -1,
     -1},
    {8.0 / 27.0,
     2.0 / 27.0,
     2.0 / 27.0,
     2.0 / 27.0,
     2.0 / 27.0,
     2.0 / 27.0,
     2.0 / 27.0,
     1.0 / 54.0,
     1.0 / 54.0,
     1.0 / 54.0,
     1.0 / 54.0,
     1.0 / 54.0,
     1.0 / 54.0,
     1.0 / 54.0,
     1.0 / 54.0,
     1.0 / 54.0,
     1.0 / 54.0,
     1.0 / 54.0,
     1.0 / 54.0,
     1.0 / 216.0,
     1.0 / 216.0,
     1.0 / 216.0,
     1.0 / 216.0,
     1.0 / 216.0,
     1.0 / 216.0,
     1.0 / 216.0,
     1.0 / 216.0},
    {0,
     2,
     1,
     4,
     3,
     6,
     5,
     10,
     9,
     8,
     7,
     14,
     13,
     12,
     11,
     18,
     17,
     16,
     15,
     26,
     25,
     24,
     23,
     22,
     21,
     20,
     19}};

inline const Lattice3DDescriptor& lattice3d(Lattice3DModel model) {
    switch (model) {
    case Lattice3DModel::D3Q19:
        return d3q19;
    case Lattice3DModel::D3Q27:
        return d3q27;
    }
    throw std::invalid_argument("Unknown three-dimensional lattice model.");
}

} // namespace lbm
