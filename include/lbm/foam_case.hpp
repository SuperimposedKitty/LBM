#pragma once

#include "lbm/openfoam.hpp"
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
};

FoamCase read_foam_case(const std::filesystem::path &path);
void block_mesh_case(const std::filesystem::path &path);
void set_fields_case(const std::filesystem::path &path);
void solve_foam_case(const std::filesystem::path &path);

} // namespace lbm
