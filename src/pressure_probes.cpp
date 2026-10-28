#include "lbm/pressure_probes.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace lbm {

std::vector<double> sample_pressure_probes(
    const FoamSnapshot& s, const std::vector<std::array<double, 3>>& positions) {
    if (positions.empty()) return {};
    if (s.two_phase) throw std::invalid_argument("pressureProbes currently require singlePhase");
    s.units.validate();
    const std::array<int, 3> count{s.nx,s.ny,s.nz};
    std::size_t size=1;
    for (int axis=0; axis<3; ++axis) {
        if (count[axis]<=0 || static_cast<std::size_t>(count[axis])>s.cells.size()/size ||
            !std::isfinite(s.origin[axis]) || !std::isfinite(s.spacing[axis]) || s.spacing[axis]<=0)
            throw std::invalid_argument("Invalid pressure probe grid");
        size*=count[axis];
    }
    if (size!=s.cells.size() || (s.two_dimensional && s.nz!=1))
        throw std::invalid_argument("Invalid pressure probe grid shape");
    std::vector<double> result;
    result.reserve(positions.size());
    for (std::size_t probe=0; probe<positions.size(); ++probe) {
        const auto fail=[&](const char* why) {
            throw std::invalid_argument("pressureProbes["+std::to_string(probe)+"] "+why);
        };
        std::array<int, 3> lower{};
        std::array<double, 3> fraction{};
        for (int axis=0; axis<3; ++axis) {
            double coordinate=(positions[probe][axis]-s.origin[axis])/s.spacing[axis];
            // 只吸收坐标换算的舍入误差，不把边界面默默夹到最近流体节点。
            constexpr double tolerance=1e-10;
            if (!std::isfinite(coordinate) || coordinate < -tolerance ||
                coordinate > count[axis]-1+tolerance)
                fail("is outside cell-center bounds (D2Q9 z must equal the mesh midplane)");
            coordinate=std::clamp(coordinate,0.0,static_cast<double>(count[axis]-1));
            if (std::abs(coordinate-std::round(coordinate))<tolerance) coordinate=std::round(coordinate);
            lower[axis]=static_cast<int>(std::floor(coordinate));
            fraction[axis]=coordinate-lower[axis];
        }
        double density_delta=0;
        for (int corner=0; corner<8; ++corner) {
            auto index=lower;
            double weight=1;
            for (int axis=0; axis<3; ++axis) {
                const int high=(corner>>axis)&1;
                index[axis]+=high;
                weight*=high?fraction[axis]:1-fraction[axis];
            }
            if (weight==0) continue;
            const auto id=(static_cast<std::size_t>(index[2])*s.ny+index[1])*s.nx+index[0];
            const auto& cell=s.cells.at(id);
            if (cell.solid) fail("interpolation stencil intersects solid; move the probe into fluid");
            if (!std::isfinite(cell.rho) || cell.rho<=0) fail("has non-finite or non-positive fluid density");
            density_delta+=weight*(cell.rho-1);
        }
        const double pressure=s.units.pressure(1+density_delta);
        if (!std::isfinite(pressure)) fail("pressure is non-finite");
        result.push_back(pressure);
    }
    return result;
}

} // namespace lbm
