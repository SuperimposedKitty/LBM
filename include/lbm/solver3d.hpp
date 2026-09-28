#pragma once

#include "lbm/openfoam.hpp"

#include "lbm/collision3d.hpp"
#include "lbm/geometry_mask3d.hpp"
#include "lbm/grid3d.hpp"

#include <string>
#include <vector>

namespace lbm {

struct Solver3DConfig {
    double tau = 0.8;
    double initial_rho = 1.0;
    Lattice3DModel lattice_model = Lattice3DModel::D3Q19;
    CollisionModel collision_model = CollisionModel::MRT;
    MrtRelaxationRates3D mrt{};
};

struct Diagnostics3D {
    double mass = 0.0;
    double max_speed = 0.0;
};

class Solver3D {
public:
    Solver3D(int nx, int ny, int nz, Solver3DConfig config);

    Grid3D& grid();
    const Grid3D& grid() const;

    void initialize_shear_wave(double amplitude, int mode);
    void initialize_lid_driven_cavity(double lid_velocity);
    void initialize_masked_flow(const GeometryMask3D& mask, double inlet_velocity);
    void step();
    void step_lid_driven_cavity(double lid_velocity);
    void step_masked_flow(double inlet_velocity);
    void run(int steps);
    void run_lid_driven_cavity(int steps, double lid_velocity);
    Diagnostics3D diagnostics() const;
    void write_vtk(const std::string& path) const;
    FoamSnapshot openfoam_snapshot() const;
    // 用宏观初始场构造平衡分布，不代表非平衡分布的断点续算。
    void initialize_fields(const FoamSnapshot& input);

private:
    Grid3D grid_;
    Solver3DConfig config_;
    std::vector<GeometryCell> geometry_cells_;
    bool masked_flow_initialized_ = false;

    void set_equilibrium_cell(int x, int y, int z, double rho, double ux, double uy, double uz);
    void collide(bool fluid_only);
    void stream_periodic();
    void stream_lid_driven_cavity(double lid_velocity);
    void stream_masked_flow();
    void apply_masked_boundaries(double inlet_velocity);
    void compute_macroscopic(bool fluid_only, double lid_velocity = 0.0);
    bool is_outside(int x, int y, int z) const;
};

} // namespace lbm
