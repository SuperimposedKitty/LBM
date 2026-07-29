#pragma once

#include "lbm/collision.hpp"
#include "lbm/lattice3d.hpp"

#include <array>

namespace lbm {

struct MrtRelaxationRates3D {
    // 体积应力与剩余高阶动理学模态使用独立松弛率，剪切模态仍由 1/tau 控制。
    double bulk = 1.1;
    double kinetic = 1.2;
};

namespace detail {

using D3Population = std::array<double, Lattice3DDescriptor::max_q>;

void validate_collision_parameters3d(
    CollisionModel model, double tau, const MrtRelaxationRates3D& rates);
double equilibrium3d(
    Lattice3DModel model,
    int direction,
    double rho,
    double ux,
    double uy,
    double uz);
D3Population equilibrium_population3d(
    Lattice3DModel model, double rho, double ux, double uy, double uz);
D3Population guo_source3d(
    Lattice3DModel model,
    double ux,
    double uy,
    double uz,
    double fx,
    double fy,
    double fz);
void collide_population3d(
    Lattice3DModel lattice_model,
    D3Population& population,
    double rho,
    double ux,
    double uy,
    double uz,
    double tau,
    CollisionModel collision_model,
    const MrtRelaxationRates3D& rates,
    const D3Population& source = {},
    double momentum_relaxation = 0.0);

} // namespace detail
} // namespace lbm
