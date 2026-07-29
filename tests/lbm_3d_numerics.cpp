#include "lbm/collision3d.hpp"
#include "lbm/solver3d.hpp"
#include "lbm/two_phase_solver3d.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

constexpr double pi = 3.141592653589793238462643383279502884;

void require(bool condition, const std::string& message) {
    if (!condition) {
        throw std::runtime_error(message);
    }
}

std::array<double, 4> moments(
    lbm::Lattice3DModel model, const lbm::detail::D3Population& population) {
    const auto& lattice = lbm::lattice3d(model);
    std::array<double, 4> result{};
    for (int q = 0; q < lattice.q; ++q) {
        result[0] += population[q];
        result[1] += lattice.cx[q] * population[q];
        result[2] += lattice.cy[q] * population[q];
        result[3] += lattice.cz[q] * population[q];
    }
    return result;
}

void test_lattice_invariants(lbm::Lattice3DModel model) {
    const auto& lattice = lbm::lattice3d(model);
    double weight_sum = 0.0;
    for (int q = 0; q < lattice.q; ++q) {
        weight_sum += lattice.w[q];
        const int opposite = lattice.opposite[q];
        require(opposite >= 0 && opposite < lattice.q, "Invalid opposite direction.");
        require(
            lattice.opposite[opposite] == q, "Three-dimensional opposite map is not involutive.");
        require(
            lattice.cx[opposite] == -lattice.cx[q] &&
                lattice.cy[opposite] == -lattice.cy[q] &&
                lattice.cz[opposite] == -lattice.cz[q],
            "Three-dimensional opposite velocity is inconsistent.");
    }
    require(std::abs(weight_sum - 1.0) < 1.0e-15, "3D lattice weights do not sum to one.");
}

void test_equilibrium_moments(lbm::Lattice3DModel model) {
    constexpr double rho = 1.07;
    constexpr double ux = 0.021;
    constexpr double uy = -0.013;
    constexpr double uz = 0.009;
    const auto equilibrium =
        lbm::detail::equilibrium_population3d(model, rho, ux, uy, uz);
    const auto conserved = moments(model, equilibrium);
    require(std::abs(conserved[0] - rho) < 1.0e-14, "3D equilibrium density mismatch.");
    require(std::abs(conserved[1] - rho * ux) < 1.0e-14, "3D equilibrium x momentum mismatch.");
    require(std::abs(conserved[2] - rho * uy) < 1.0e-14, "3D equilibrium y momentum mismatch.");
    require(std::abs(conserved[3] - rho * uz) < 1.0e-14, "3D equilibrium z momentum mismatch.");
}

void test_collision_conservation(lbm::Lattice3DModel model) {
    const auto& lattice = lbm::lattice3d(model);
    lbm::detail::D3Population population{};
    for (int q = 0; q < lattice.q; ++q) {
        population[q] = lattice.w[q] * (1.0 + 0.015 * std::sin(0.7 * q));
    }
    const auto before = moments(model, population);
    const double rho = before[0];
    const double ux = before[1] / rho;
    const double uy = before[2] / rho;
    const double uz = before[3] / rho;
    lbm::detail::collide_population3d(
        model,
        population,
        rho,
        ux,
        uy,
        uz,
        0.78,
        lbm::CollisionModel::MRT,
        {});
    const auto after = moments(model, population);
    for (int i = 0; i < 4; ++i) {
        require(
            std::abs(after[i] - before[i]) < 2.0e-13,
            "Hermite MRT collision does not conserve a hydrodynamic moment.");
    }
}

double shear_amplitude(const lbm::Solver3D& solver, int mode) {
    const lbm::Grid3D& grid = solver.grid();
    double projection = 0.0;
    for (int z = 0; z < grid.nz; ++z) {
        for (int y = 0; y < grid.ny; ++y) {
            const double basis =
                std::sin(2.0 * pi * mode * static_cast<double>(y) / grid.ny);
            for (int x = 0; x < grid.nx; ++x) {
                projection += grid.ux[grid.scalar_index(x, y, z)] * basis;
            }
        }
    }
    return 2.0 * projection /
           static_cast<double>(grid.nx * grid.ny * grid.nz);
}

void test_shear_decay(lbm::Lattice3DModel model) {
    lbm::Solver3DConfig config;
    config.lattice_model = model;
    config.collision_model = lbm::CollisionModel::MRT;
    config.tau = 0.8;
    lbm::Solver3D solver(10, 24, 8, config);
    constexpr double initial_amplitude = 0.01;
    constexpr int mode = 1;
    constexpr int steps = 80;
    solver.initialize_shear_wave(initial_amplitude, mode);
    solver.run(steps);
    const double measured = std::abs(shear_amplitude(solver, mode));
    const double wave_number = 2.0 * pi * mode / 24.0;
    const double viscosity =
        lbm::Lattice3DDescriptor::cs2 * (config.tau - 0.5);
    const double expected =
        initial_amplitude * std::exp(-viscosity * wave_number * wave_number * steps);
    const double relative_error = std::abs(measured - expected) / expected;
    require(relative_error < 0.04, "Three-dimensional shear-wave viscosity mismatch.");
}

void test_static_droplet() {
    lbm::TwoPhaseConfig3D config;
    config.lattice_model = lbm::Lattice3DModel::D3Q27;
    config.interaction_strength = 2.0;
    config.body_force_x = 0.0;
    config.body_force_y = 0.0;
    config.body_force_z = 0.0;
    config.wall_adhesion_strength = 0.0;
    config.recoloring_strength = 0.05;
    config.bottom_wall_thickness = 1;
    lbm::TwoPhaseSolver3D solver(16, 16, 16, config);
    solver.initialize_droplet_impact(7.5, 8.0, 7.5, 3.5, 0.0, 0.0, 0.0);
    const auto before = solver.diagnostics();
    solver.run_closed(20);
    const auto after = solver.diagnostics();
    const double relative_mass_error =
        std::abs(after.mass_a - before.mass_a) / before.mass_a;
    require(relative_mass_error < 1.0e-10, "D3Q27 droplet component mass drifted.");
    require(std::isfinite(after.max_speed), "D3Q27 droplet speed is not finite.");
    require(
        std::abs(after.interaction_force_balance_x) < 1.0e-12 &&
            std::abs(after.interaction_force_balance_y) < 1.0e-12 &&
            std::abs(after.interaction_force_balance_z) < 1.0e-12,
        "D3Q27 symmetric interaction force does not balance.");
}

} // namespace

int main() {
    try {
        for (lbm::Lattice3DModel model :
             {lbm::Lattice3DModel::D3Q19, lbm::Lattice3DModel::D3Q27}) {
            test_lattice_invariants(model);
            test_equilibrium_moments(model);
            test_collision_conservation(model);
            test_shear_decay(model);
        }
        test_static_droplet();
        std::cout << "All three-dimensional LBM numerical tests passed.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Three-dimensional numerical test failure: " << error.what()
                  << '\n';
        return 1;
    }
}
