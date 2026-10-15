#pragma once

#include "lbm/lattice.hpp"
#include "lbm/lattice3d.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace lbm {

enum class OpenBoundaryModel { Legacy, VelocityPressure };

struct SinglePhaseOptions {
    // 体积力 F 与加速度 a 分开给定，总力密度为 F + rho*a，均使用格子单位。
    std::array<double, 3> force_density{};
    std::array<double, 3> acceleration{};
    OpenBoundaryModel open_boundary = OpenBoundaryModel::Legacy;
    double outlet_density = 1.0;
    bool interpolated_walls = false;
    bool periodic_x = false;
    bool periodic_z = false;
};

struct BoundaryDiagnostics {
    // 本步边界算子引起的质量增量；可与总质量变化直接核对。
    double mass_exchange = 0.0;
    double wall_mass_exchange = 0.0;
    std::array<double, 3> wall_force{};
    std::array<double, 3> obstacle_force{};
    std::size_t interpolation_fallbacks = 0;
};

struct WallNode {
    double distance = 0.0;
    bool distance_known = false;
    bool obstacle = false;
    std::array<double, 3> velocity{};
    double inlet_scale = 1.0;
};

namespace detail {
inline void validate_flow_options(const SinglePhaseOptions& o, bool two_d) {
    if (!std::isfinite(o.outlet_density) || o.outlet_density <= 0)
        throw std::invalid_argument("Outlet density must be positive and finite.");
    for (int d = 0; d < 3; ++d)
        if (!std::isfinite(o.force_density[d]) || !std::isfinite(o.acceleration[d]))
            throw std::invalid_argument("Single-phase force must be finite.");
    if (two_d && (o.force_density[2] != 0 || o.acceleration[2] != 0 || o.periodic_z))
        throw std::invalid_argument("D2Q9 cannot use z forcing/periodicity.");
}

inline std::array<double, 3> total_force(const SinglePhaseOptions& o, double rho) {
    return {o.force_density[0] + rho * o.acceleration[0],
            o.force_density[1] + rho * o.acceleration[1],
            o.force_density[2] + rho * o.acceleration[2]};
}

inline const Lattice3DDescriptor& flow_lattice2d() {
    static const auto descriptor = [] {
        Lattice3DDescriptor l;
        l.q = D2Q9::q;
        for (int q = 0; q < l.q; ++q) {
            l.cx[q] = D2Q9::cx[q]; l.cy[q] = D2Q9::cy[q];
            l.w[q] = D2Q9::w[q]; l.opposite[q] = D2Q9::opposite[q];
        }
        return l;
    }();
    return descriptor;
}

// Zou-He 型未知分布重构；切向修正保证三个动量矩，已知分布保持不变。
template <class Population>
void reconstruct_x_boundary(Population& f, const Lattice3DDescriptor& l,
                            const SinglePhaseOptions& o, bool inlet,
                            std::array<double, 3> u) {
    const int normal = inlet ? 1 : -1;
    double known = 0;
    for (int q = 0; q < l.q; ++q)
        if (l.cx[q] == 0) known += f[q];
        else if (l.cx[q] == -normal) known += 2 * f[q];
    double rho = o.outlet_density;
    if (inlet)
        rho = (known - 0.5 * normal * o.force_density[0]) /
              (1 - normal * u[0] + 0.5 * normal * o.acceleration[0]);
    else {
        const auto force = total_force(o, rho);
        u[0] = (normal * (rho - known) + 0.5 * force[0]) / rho;
    }
    if (!(rho > 0) || !std::isfinite(rho) || !std::isfinite(u[0]) || std::abs(u[0]) >= 0.5)
        throw std::runtime_error("Non-physical velocity/pressure boundary state.");
    const auto force = total_force(o, rho);
    std::array<double, 3> raw{u[0] - force[0] / (2*rho),
                              u[1] - force[1] / (2*rho),
                              u[2] - force[2] / (2*rho)};
    for (int q = 0; q < l.q; ++q)
        if (l.cx[q] == normal)
            f[q] = f[l.opposite[q]] + 6*l.w[q]*rho*
                (l.cx[q]*raw[0] + l.cy[q]*raw[1] + l.cz[q]*raw[2]);
    double jy = 0, jz = 0, sy = 0, sz = 0;
    for (int q = 0; q < l.q; ++q) {
        jy += l.cy[q]*f[q]; jz += l.cz[q]*f[q];
        if (l.cx[q] == normal) {
            sy += l.w[q]*l.cy[q]*l.cy[q];
            sz += l.w[q]*l.cz[q]*l.cz[q];
        }
    }
    for (int q = 0; q < l.q; ++q)
        if (l.cx[q] == normal) {
            if (sy > 0) f[q] += l.w[q]*l.cy[q]*(rho*raw[1]-jy)/sy;
            if (sz > 0) f[q] += l.w[q]*l.cz[q]*(rho*raw[2]-jz)/sz;
        }
}

inline double reflected_population(double outgoing, double opposite, double behind,
                                    double fraction, double wall_correction) {
    if (fraction < 0.5)
        return 2*fraction*outgoing + (1-2*fraction)*behind - wall_correction;
    return (outgoing + (2*fraction-1)*opposite - wall_correction)/(2*fraction);
}
} // namespace detail
} // namespace lbm
