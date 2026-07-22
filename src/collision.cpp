#include "lbm/collision.hpp"

#include <cmath>
#include <stdexcept>
#include <string>

namespace lbm::detail {
namespace {

void validate_rate(double rate, const char* name) {
    if (!std::isfinite(rate) || rate <= 0.0 || rate >= 2.0) {
        throw std::invalid_argument(std::string(name) + " must be in the range (0, 2).");
    }
}

} // namespace

void validate_collision_parameters(
    CollisionModel model, double tau, const MrtRelaxationRates& rates) {
    if (model != CollisionModel::MRT && model != CollisionModel::BGK) {
        throw std::invalid_argument("Unknown collision model.");
    }
    if (!std::isfinite(tau) || tau <= 0.5) {
        throw std::invalid_argument("tau must be finite and greater than 0.5.");
    }
    if (model == CollisionModel::MRT) {
        validate_rate(rates.bulk, "MRT bulk relaxation rate");
        validate_rate(rates.energy, "MRT energy relaxation rate");
        validate_rate(rates.energy_flux, "MRT energy-flux relaxation rate");
        validate_rate(1.0 / tau, "MRT shear relaxation rate");
    }
}

double equilibrium(int direction, double rho, double ux, double uy) {
    const double cu = static_cast<double>(D2Q9::cx[direction]) * ux +
                      static_cast<double>(D2Q9::cy[direction]) * uy;
    const double u2 = ux * ux + uy * uy;
    return D2Q9::w[direction] * rho *
           (1.0 + 3.0 * cu + 4.5 * cu * cu - 1.5 * u2);
}

D2Q9Population equilibrium_population(double rho, double ux, double uy) {
    D2Q9Population result{};
    for (int q = 0; q < D2Q9::q; ++q) {
        result[q] = equilibrium(q, rho, ux, uy);
    }
    return result;
}

D2Q9Population population_to_moments(const D2Q9Population& f) {
    // Lallemand-Luo 正交矩基，方向顺序与 D2Q9::cx/cy 完全一致。
    return {
        f[0] + f[1] + f[2] + f[3] + f[4] + f[5] + f[6] + f[7] + f[8],
        -4.0 * f[0] - f[1] - f[2] - f[3] - f[4] +
            2.0 * (f[5] + f[6] + f[7] + f[8]),
        4.0 * f[0] - 2.0 * (f[1] + f[2] + f[3] + f[4]) +
            f[5] + f[6] + f[7] + f[8],
        f[1] - f[3] + f[5] - f[6] - f[7] + f[8],
        -2.0 * f[1] + 2.0 * f[3] + f[5] - f[6] - f[7] + f[8],
        f[2] - f[4] + f[5] + f[6] - f[7] - f[8],
        -2.0 * f[2] + 2.0 * f[4] + f[5] + f[6] - f[7] - f[8],
        f[1] - f[2] + f[3] - f[4],
        f[5] - f[6] + f[7] - f[8],
    };
}

D2Q9Population moments_to_population(const D2Q9Population& m) {
    constexpr double one_over_36 = 1.0 / 36.0;
    return {
        (4.0 * m[0] - 4.0 * m[1] + 4.0 * m[2]) * one_over_36,
        (4.0 * m[0] - m[1] - 2.0 * m[2] + 6.0 * m[3] - 6.0 * m[4] +
         9.0 * m[7]) * one_over_36,
        (4.0 * m[0] - m[1] - 2.0 * m[2] + 6.0 * m[5] - 6.0 * m[6] -
         9.0 * m[7]) * one_over_36,
        (4.0 * m[0] - m[1] - 2.0 * m[2] - 6.0 * m[3] + 6.0 * m[4] +
         9.0 * m[7]) * one_over_36,
        (4.0 * m[0] - m[1] - 2.0 * m[2] - 6.0 * m[5] + 6.0 * m[6] -
         9.0 * m[7]) * one_over_36,
        (4.0 * m[0] + 2.0 * m[1] + m[2] + 6.0 * m[3] + 3.0 * m[4] +
         6.0 * m[5] + 3.0 * m[6] + 9.0 * m[8]) * one_over_36,
        (4.0 * m[0] + 2.0 * m[1] + m[2] - 6.0 * m[3] - 3.0 * m[4] +
         6.0 * m[5] + 3.0 * m[6] - 9.0 * m[8]) * one_over_36,
        (4.0 * m[0] + 2.0 * m[1] + m[2] - 6.0 * m[3] - 3.0 * m[4] -
         6.0 * m[5] - 3.0 * m[6] + 9.0 * m[8]) * one_over_36,
        (4.0 * m[0] + 2.0 * m[1] + m[2] + 6.0 * m[3] + 3.0 * m[4] -
         6.0 * m[5] - 3.0 * m[6] - 9.0 * m[8]) * one_over_36,
    };
}

D2Q9Population guo_source(double ux, double uy, double fx, double fy) {
    D2Q9Population source{};
    for (int q = 0; q < D2Q9::q; ++q) {
        const double cx = static_cast<double>(D2Q9::cx[q]);
        const double cy = static_cast<double>(D2Q9::cy[q]);
        const double cu = cx * ux + cy * uy;
        const double term_x =
            (cx - ux) / D2Q9::cs2 + cu * cx / (D2Q9::cs2 * D2Q9::cs2);
        const double term_y =
            (cy - uy) / D2Q9::cs2 + cu * cy / (D2Q9::cs2 * D2Q9::cs2);
        source[q] = D2Q9::w[q] * (term_x * fx + term_y * fy);
    }
    return source;
}

void collide_population(
    D2Q9Population& population,
    double rho,
    double ux,
    double uy,
    double tau,
    CollisionModel model,
    const MrtRelaxationRates& rates,
    const D2Q9Population& source,
    double momentum_relaxation) {
    const D2Q9Population equilibrium_values = equilibrium_population(rho, ux, uy);
    const double shear_rate = 1.0 / tau;

    if (model == CollisionModel::BGK) {
        for (int q = 0; q < D2Q9::q; ++q) {
            population[q] -= shear_rate * (population[q] - equilibrium_values[q]);
            population[q] += (1.0 - 0.5 * shear_rate) * source[q];
        }
        return;
    }

    D2Q9Population moments = population_to_moments(population);
    const D2Q9Population equilibrium_moments = population_to_moments(equilibrium_values);
    const D2Q9Population source_moments = population_to_moments(source);
    const D2Q9Population relaxation{
        0.0,
        rates.bulk,
        rates.energy,
        momentum_relaxation,
        rates.energy_flux,
        momentum_relaxation,
        rates.energy_flux,
        shear_rate,
        shear_rate,
    };

    for (int i = 0; i < D2Q9::q; ++i) {
        moments[i] -= relaxation[i] * (moments[i] - equilibrium_moments[i]);
        moments[i] += (1.0 - 0.5 * relaxation[i]) * source_moments[i];
    }
    population = moments_to_population(moments);
}

} // namespace lbm::detail
