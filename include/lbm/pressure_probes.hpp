#pragma once

#include "lbm/openfoam.hpp"

namespace lbm {

// 坐标与网格使用相同单位；返回与 OpenFOAM p 相同的压力（SI 为 Pa）。
// 仅单相流体中心包围的双线性/三线性插值，不作固壁或域外外推。
std::vector<double> sample_pressure_probes(
    const FoamSnapshot& snapshot,
    const std::vector<std::array<double, 3>>& positions);

} // namespace lbm
