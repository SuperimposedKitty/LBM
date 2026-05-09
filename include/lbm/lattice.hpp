#pragma once

#include <array>

namespace lbm {

struct D2Q9 {
    static constexpr int q = 9;
    static constexpr double cs2 = 1.0 / 3.0;

    static constexpr std::array<int, q> cx{0, 1, 0, -1, 0, 1, -1, -1, 1};
    static constexpr std::array<int, q> cy{0, 0, 1, 0, -1, 1, 1, -1, -1};
    static constexpr std::array<double, q> w{
        4.0 / 9.0,
        1.0 / 9.0,
        1.0 / 9.0,
        1.0 / 9.0,
        1.0 / 9.0,
        1.0 / 36.0,
        1.0 / 36.0,
        1.0 / 36.0,
        1.0 / 36.0};
    static constexpr std::array<int, q> opposite{0, 3, 4, 1, 2, 7, 8, 5, 6};
};

} // namespace lbm

