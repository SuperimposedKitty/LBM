#pragma once

#include "lbm/lattice.hpp"

#include <array>

namespace lbm {

enum class CollisionModel {
    MRT,
    BGK
};

struct MrtRelaxationRates {
    // 非守恒矩的松弛率；剪切矩仍由 1 / tau 控制。
    double bulk = 1.1;
    double energy = 1.0;
    double energy_flux = 1.2;
};

namespace detail {

using D2Q9Population = std::array<double, D2Q9::q>;

void validate_collision_parameters(
    CollisionModel model, double tau, const MrtRelaxationRates& rates);
double equilibrium(int direction, double rho, double ux, double uy);
D2Q9Population equilibrium_population(double rho, double ux, double uy);
D2Q9Population population_to_moments(const D2Q9Population& population);
D2Q9Population moments_to_population(const D2Q9Population& moments);
D2Q9Population guo_source(double ux, double uy, double fx, double fy);
void collide_population(
    D2Q9Population& population,
    double rho,
    double ux,
    double uy,
    double tau,
    CollisionModel model,
    const MrtRelaxationRates& rates,
    const D2Q9Population& source = {},
    double momentum_relaxation = 0.0);

} // namespace detail
} // namespace lbm
