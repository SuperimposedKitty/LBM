#include "lbm/collision3d.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>

namespace lbm::detail {
namespace {

constexpr double cs2 = Lattice3DDescriptor::cs2;
constexpr double cs4 = cs2 * cs2;

struct HermiteModes {
    D3Population density{};
    D3Population momentum{};
    D3Population bulk{};
    D3Population shear{};
    D3Population kinetic{};
};

void validate_rate(double value, const char* name) {
    if (!std::isfinite(value) || value <= 0.0 || value >= 2.0) {
        throw std::invalid_argument(std::string(name) + " must be in the range (0, 2).");
    }
}

HermiteModes decompose_modes(Lattice3DModel model, const D3Population& values) {
    const auto& lattice = lattice3d(model);
    HermiteModes modes{};

    double density = 0.0;
    double jx = 0.0;
    double jy = 0.0;
    double jz = 0.0;
    for (int q = 0; q < lattice.q; ++q) {
        density += values[q];
        jx += static_cast<double>(lattice.cx[q]) * values[q];
        jy += static_cast<double>(lattice.cy[q]) * values[q];
        jz += static_cast<double>(lattice.cz[q]) * values[q];
    }

    D3Population second_order_remainder{};
    for (int q = 0; q < lattice.q; ++q) {
        const double cx = static_cast<double>(lattice.cx[q]);
        const double cy = static_cast<double>(lattice.cy[q]);
        const double cz = static_cast<double>(lattice.cz[q]);
        modes.density[q] = lattice.w[q] * density;
        modes.momentum[q] = lattice.w[q] * (cx * jx + cy * jy + cz * jz) / cs2;
        second_order_remainder[q] = values[q] - modes.density[q] - modes.momentum[q];
    }

    double pxx = 0.0;
    double pyy = 0.0;
    double pzz = 0.0;
    double pxy = 0.0;
    double pxz = 0.0;
    double pyz = 0.0;
    for (int q = 0; q < lattice.q; ++q) {
        const double cx = static_cast<double>(lattice.cx[q]);
        const double cy = static_cast<double>(lattice.cy[q]);
        const double cz = static_cast<double>(lattice.cz[q]);
        pxx += second_order_remainder[q] * cx * cx;
        pyy += second_order_remainder[q] * cy * cy;
        pzz += second_order_remainder[q] * cz * cz;
        pxy += second_order_remainder[q] * cx * cy;
        pxz += second_order_remainder[q] * cx * cz;
        pyz += second_order_remainder[q] * cy * cz;
    }

    const double mean_normal = (pxx + pyy + pzz) / 3.0;
    const double dxx = pxx - mean_normal;
    const double dyy = pyy - mean_normal;
    const double dzz = pzz - mean_normal;
    for (int q = 0; q < lattice.q; ++q) {
        const double cx = static_cast<double>(lattice.cx[q]);
        const double cy = static_cast<double>(lattice.cy[q]);
        const double cz = static_cast<double>(lattice.cz[q]);
        const double qxx = cx * cx - cs2;
        const double qyy = cy * cy - cs2;
        const double qzz = cz * cz - cs2;
        const double scale = lattice.w[q] / (2.0 * cs4);
        modes.bulk[q] = scale * mean_normal * (qxx + qyy + qzz);
        modes.shear[q] =
            scale * (qxx * dxx + qyy * dyy + qzz * dzz +
                     2.0 * (cx * cy * pxy + cx * cz * pxz + cy * cz * pyz));
        modes.kinetic[q] =
            second_order_remainder[q] - modes.bulk[q] - modes.shear[q];
    }
    return modes;
}

double relaxed_mode(
    const HermiteModes& modes,
    int q,
    double momentum_rate,
    double bulk_rate,
    double shear_rate,
    double kinetic_rate,
    bool source) {
    const double factor = source ? 0.5 : 1.0;
    return modes.density[q] +
           (1.0 - factor * momentum_rate) * modes.momentum[q] +
           (1.0 - factor * bulk_rate) * modes.bulk[q] +
           (1.0 - factor * shear_rate) * modes.shear[q] +
           (1.0 - factor * kinetic_rate) * modes.kinetic[q];
}

} // namespace

void validate_collision_parameters3d(
    CollisionModel model, double tau, const MrtRelaxationRates3D& rates) {
    if (model != CollisionModel::MRT && model != CollisionModel::BGK) {
        throw std::invalid_argument("Unknown collision model.");
    }
    if (!std::isfinite(tau) || tau <= 0.5) {
        throw std::invalid_argument("Three-dimensional relaxation time must be greater than 0.5.");
    }
    if (model == CollisionModel::MRT) {
        validate_rate(rates.bulk, "Three-dimensional MRT bulk relaxation rate");
        validate_rate(rates.kinetic, "Three-dimensional MRT kinetic relaxation rate");
    }
}

double equilibrium3d(
    Lattice3DModel model,
    int direction,
    double rho,
    double ux,
    double uy,
    double uz) {
    const auto& lattice = lattice3d(model);
    if (direction < 0 || direction >= lattice.q) {
        throw std::out_of_range("Three-dimensional lattice direction is out of range.");
    }
    const double cu = static_cast<double>(lattice.cx[direction]) * ux +
                      static_cast<double>(lattice.cy[direction]) * uy +
                      static_cast<double>(lattice.cz[direction]) * uz;
    const double u2 = ux * ux + uy * uy + uz * uz;
    return lattice.w[direction] * rho *
           (1.0 + cu / cs2 + 0.5 * cu * cu / cs4 - 0.5 * u2 / cs2);
}

D3Population equilibrium_population3d(
    Lattice3DModel model, double rho, double ux, double uy, double uz) {
    const auto& lattice = lattice3d(model);
    D3Population result{};
    for (int q = 0; q < lattice.q; ++q) {
        result[q] = equilibrium3d(model, q, rho, ux, uy, uz);
    }
    return result;
}

D3Population guo_source3d(
    Lattice3DModel model,
    double ux,
    double uy,
    double uz,
    double fx,
    double fy,
    double fz) {
    const auto& lattice = lattice3d(model);
    D3Population source{};
    for (int q = 0; q < lattice.q; ++q) {
        const double cx = static_cast<double>(lattice.cx[q]);
        const double cy = static_cast<double>(lattice.cy[q]);
        const double cz = static_cast<double>(lattice.cz[q]);
        const double cu = cx * ux + cy * uy + cz * uz;
        const double term_x = (cx - ux) / cs2 + cu * cx / cs4;
        const double term_y = (cy - uy) / cs2 + cu * cy / cs4;
        const double term_z = (cz - uz) / cs2 + cu * cz / cs4;
        source[q] = lattice.w[q] * (term_x * fx + term_y * fy + term_z * fz);
    }
    return source;
}

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
    const D3Population& source,
    double momentum_relaxation) {
    validate_collision_parameters3d(collision_model, tau, rates);
    const auto& lattice = lattice3d(lattice_model);
    const D3Population equilibrium_values =
        equilibrium_population3d(lattice_model, rho, ux, uy, uz);
    const double shear_rate = 1.0 / tau;

    if (collision_model == CollisionModel::BGK) {
        const double force_factor = 1.0 - 0.5 * shear_rate;
        for (int q = 0; q < lattice.q; ++q) {
            population[q] += -shear_rate * (population[q] - equilibrium_values[q]) +
                             force_factor * source[q];
        }
        return;
    }

    D3Population non_equilibrium{};
    for (int q = 0; q < lattice.q; ++q) {
        non_equilibrium[q] = population[q] - equilibrium_values[q];
    }
    const HermiteModes non_equilibrium_modes =
        decompose_modes(lattice_model, non_equilibrium);
    const HermiteModes source_modes = decompose_modes(lattice_model, source);

    for (int q = 0; q < lattice.q; ++q) {
        population[q] =
            equilibrium_values[q] +
            relaxed_mode(
                non_equilibrium_modes,
                q,
                momentum_relaxation,
                rates.bulk,
                shear_rate,
                rates.kinetic,
                false) +
            relaxed_mode(
                source_modes,
                q,
                momentum_relaxation,
                rates.bulk,
                shear_rate,
                rates.kinetic,
                true);
    }
}

} // namespace lbm::detail
