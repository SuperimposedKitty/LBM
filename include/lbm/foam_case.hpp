#pragma once

#include "lbm/openfoam.hpp"
#include "lbm/convergence.hpp"
#include "lbm/solver.hpp"
#include "lbm/solver3d.hpp"
#include "lbm/two_phase_solver.hpp"
#include "lbm/two_phase_solver3d.hpp"

namespace lbm {

// 文件顺序和格点顺序显式映射，不能假定 blockMesh 的单元编号。
struct FoamCase {
    FoamSnapshot initial;
    std::vector<std::size_t> cell_order;
    SolverConfig single;
    Solver3DConfig single3d;
    TwoPhaseConfig two;
    TwoPhaseConfig3D two3d;
    int steps = 100;
    int write_interval = 10;
    double inlet_velocity = 0.0;
    double reference_length = 0.0;
    double max_mach = 0.15;
    double density_warning = 0.01;
    bool parabolic_inlet = false;
    double inlet_ramp_steps = 0.0;
    ConvergenceConfig convergence{};
};

struct FoamRunSummary {
    int steps = 0;
    bool converged = false;
    double velocity_change = 0, density_change = 0, force_range = 0;
    double drag = 0, lift = 0, maximum_mach = 0;
};

FoamCase read_foam_case(const std::filesystem::path &path);
void block_mesh_case(const std::filesystem::path &path);
void set_fields_case(const std::filesystem::path &path);
void solve_foam_case(const std::filesystem::path &path);
FoamRunSummary run_foam_case(const FoamCase& simulation, const char* output_case);
void refine_foam_case(const std::filesystem::path& source,
                      const std::filesystem::path& target, int factor);

} // namespace lbm
