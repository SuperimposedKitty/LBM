#include "lbm/collision.hpp"
#include "lbm/solver.hpp"
#include "lbm/two_phase_solver.hpp"

#include <algorithm>
#include <cmath>
#include <exception>
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

double max_difference(
    const lbm::detail::D2Q9Population& a, const lbm::detail::D2Q9Population& b) {
    double result = 0.0;
    for (int q = 0; q < lbm::D2Q9::q; ++q) {
        result = std::max(result, std::abs(a[q] - b[q]));
    }
    return result;
}

void macroscopic(
    const lbm::detail::D2Q9Population& population,
    double& rho,
    double& momentum_x,
    double& momentum_y) {
    rho = 0.0;
    momentum_x = 0.0;
    momentum_y = 0.0;
    for (int q = 0; q < lbm::D2Q9::q; ++q) {
        rho += population[q];
        momentum_x += static_cast<double>(lbm::D2Q9::cx[q]) * population[q];
        momentum_y += static_cast<double>(lbm::D2Q9::cy[q]) * population[q];
    }
}

void test_moment_round_trip() {
    const lbm::detail::D2Q9Population population{
        0.41, 0.08, 0.07, 0.06, 0.09, 0.04, 0.05, 0.03, 0.02};
    const auto moments = lbm::detail::population_to_moments(population);
    const auto restored = lbm::detail::moments_to_population(moments);
    require(max_difference(population, restored) < 1.0e-12,
            "D2Q9 moment transform round trip exceeded tolerance.");
}

void test_mrt_conservation() {
    lbm::detail::D2Q9Population population{
        0.41, 0.08, 0.07, 0.06, 0.09, 0.04, 0.05, 0.03, 0.02};
    double rho_before = 0.0;
    double momentum_x_before = 0.0;
    double momentum_y_before = 0.0;
    macroscopic(population, rho_before, momentum_x_before, momentum_y_before);

    lbm::detail::collide_population(
        population,
        rho_before,
        momentum_x_before / rho_before,
        momentum_y_before / rho_before,
        0.8,
        lbm::CollisionModel::MRT,
        {});

    double rho_after = 0.0;
    double momentum_x_after = 0.0;
    double momentum_y_after = 0.0;
    macroscopic(population, rho_after, momentum_x_after, momentum_y_after);
    require(std::abs(rho_after - rho_before) < 1.0e-12,
            "MRT collision did not conserve density.");
    require(std::abs(momentum_x_after - momentum_x_before) < 1.0e-12 &&
                std::abs(momentum_y_after - momentum_y_before) < 1.0e-12,
            "MRT collision did not conserve momentum.");
}

void test_bgk_reference() {
    const lbm::detail::D2Q9Population initial{
        0.41, 0.08, 0.07, 0.06, 0.09, 0.04, 0.05, 0.03, 0.02};
    auto actual = initial;
    constexpr double rho = 0.85;
    constexpr double ux = 0.025;
    constexpr double uy = -0.015;
    constexpr double tau = 0.82;
    const auto source = lbm::detail::guo_source(ux, uy, 2.0e-5, -3.0e-5);
    const auto equilibrium = lbm::detail::equilibrium_population(rho, ux, uy);
    auto expected = initial;
    const double omega = 1.0 / tau;
    for (int q = 0; q < lbm::D2Q9::q; ++q) {
        expected[q] -= omega * (expected[q] - equilibrium[q]);
        expected[q] += (1.0 - 0.5 * omega) * source[q];
    }

    lbm::detail::collide_population(
        actual, rho, ux, uy, tau, lbm::CollisionModel::BGK, {}, source);
    require(max_difference(actual, expected) < 1.0e-14,
            "BGK compatibility branch differs from the legacy formula.");
}

void test_parameter_validation() {
    lbm::MrtRelaxationRates invalid;
    invalid.bulk = 2.0;
    bool rejected = false;
    try {
        lbm::detail::validate_collision_parameters(
            lbm::CollisionModel::MRT, 0.8, invalid);
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    require(rejected, "Invalid MRT relaxation rate was not rejected.");
}

double shear_amplitude(const lbm::Solver& solver, int mode) {
    const auto& grid = solver.grid();
    double projection = 0.0;
    for (int y = 0; y < grid.ny; ++y) {
        const double basis =
            std::sin(2.0 * pi * mode * static_cast<double>(y) / static_cast<double>(grid.ny));
        for (int x = 0; x < grid.nx; ++x) {
            projection += grid.ux[grid.scalar_index(x, y)] * basis;
        }
    }
    return 2.0 * projection / static_cast<double>(grid.nx * grid.ny);
}

void test_shear_wave_viscosity() {
    constexpr int nx = 64;
    constexpr int ny = 64;
    constexpr int steps = 200;
    constexpr int mode = 1;
    lbm::SolverConfig config;
    config.tau = 0.8;
    config.collision_model = lbm::CollisionModel::MRT;
    lbm::Solver solver(nx, ny, config);
    solver.initialize_shear_wave(0.01, mode);
    const double amplitude_before = shear_amplitude(solver, mode);
    solver.run(steps);
    const double amplitude_after = shear_amplitude(solver, mode);
    const double wave_number = 2.0 * pi * mode / static_cast<double>(ny);
    const double measured_viscosity =
        -std::log(amplitude_after / amplitude_before) /
        (wave_number * wave_number * static_cast<double>(steps));
    const double expected_viscosity = lbm::D2Q9::cs2 * (config.tau - 0.5);
    const double relative_error =
        std::abs(measured_viscosity - expected_viscosity) / expected_viscosity;
    require(relative_error < 0.03, "MRT shear-wave viscosity error exceeded 3 percent.");
}

void test_low_viscosity_finiteness() {
    lbm::SolverConfig config;
    config.tau = 0.505;
    config.collision_model = lbm::CollisionModel::MRT;
    lbm::Solver solver(64, 64, config);
    solver.initialize_shear_wave(0.01, 2);
    solver.run(500);
    for (double rho : solver.grid().rho) {
        require(std::isfinite(rho) && rho > 0.0,
                "Low-viscosity MRT produced a non-physical density.");
    }
    require(std::isfinite(solver.diagnostics().max_speed),
            "Low-viscosity MRT produced a non-finite velocity.");
}

struct DropletResult {
    double mass_error_a = 0.0;
    double mass_error_b = 0.0;
    double max_speed = 0.0;
    double force_balance = 0.0;
};

DropletResult run_static_droplet(lbm::CollisionModel model) {
    lbm::TwoPhaseConfig config;
    config.collision_model = model;
    config.tau_a = 0.9;
    config.tau_b = 0.9;
    config.interaction_strength = 3.0;
    config.body_force_x = 0.0;
    config.body_force_y = 0.0;
    config.wall_adhesion_strength = 0.0;
    config.droplet_interface_width = 1.5;
    lbm::TwoPhaseSolver solver(48, 48, config);
    solver.initialize_droplet_impact(23.5, 24.0, 10.0, 0.0, 0.0);
    const auto before = solver.diagnostics();
    solver.run_closed(300);
    const auto after = solver.diagnostics();
    return {
        std::abs(after.mass_a - before.mass_a) / std::max(1.0, std::abs(before.mass_a)),
        std::abs(after.mass_b - before.mass_b) / std::max(1.0, std::abs(before.mass_b)),
        after.max_speed,
        std::hypot(
            after.interaction_force_balance_x, after.interaction_force_balance_y),
    };
}

void test_two_phase_balance_and_mass() {
    const DropletResult mrt = run_static_droplet(lbm::CollisionModel::MRT);
    const DropletResult bgk = run_static_droplet(lbm::CollisionModel::BGK);
    require(mrt.mass_error_a < 1.0e-10 && mrt.mass_error_b < 1.0e-10,
            "Closed MRT droplet did not conserve component masses.");
    require(mrt.force_balance < 1.0e-12,
            "Pairwise Shan-Chen interaction force is not globally balanced.");
    require(std::isfinite(mrt.max_speed) && std::isfinite(bgk.max_speed),
            "Static droplet produced a non-finite spurious velocity.");
    require(mrt.max_speed <= 1.02 * bgk.max_speed + 1.0e-12,
            "MRT increased static-droplet spurious velocity by more than 2 percent.");
    std::cout << "static droplet max speed: MRT=" << mrt.max_speed
              << ", BGK=" << bgk.max_speed << '\n';
}

} // namespace

int main() {
    try {
        test_moment_round_trip();
        test_mrt_conservation();
        test_bgk_reference();
        test_parameter_validation();
        test_shear_wave_viscosity();
        test_low_viscosity_finiteness();
        test_two_phase_balance_and_mass();
        std::cout << "all LBM numerical tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "LBM numerical test failed: " << error.what() << '\n';
        return 1;
    }
}
