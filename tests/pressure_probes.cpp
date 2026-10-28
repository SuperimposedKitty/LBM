#include "lbm/foam_case.hpp"

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {
using Vec=std::array<double,3>;
void check(bool ok,const char* why) {if (!ok) throw std::runtime_error(why);}
template<class F> void rejects(F f,const char* fragment) {
    try {f();} catch (const std::exception& e) {
        check(std::string(e.what()).find(fragment)!=std::string::npos,"Incorrect probe rejection diagnostic");
        return;
    }
    throw std::runtime_error("Invalid pressure probe accepted");
}

void interpolation(bool three_d,bool si) {
    lbm::FoamSnapshot s;
    s.nx=4; s.ny=3; s.nz=three_d?3:1; s.two_dimensional=!three_d;
    s.origin={-0.2,0.7,0.05}; s.spacing={0.1,0.1,three_d?0.1:0.03};
    s.units.si=si; s.units.length_step=0.1; s.units.time_step=0.02;
    s.units.density_scale=1000; s.units.reference_pressure=101325; s.units.kinematic_viscosity=0.01;
    s.cells.resize(s.nx*s.ny*s.nz);
    // 多线性函数应被双线性/三线性插值精确还原，含非零原点与 SI 参考压力。
    auto rho=[](const Vec& p) {return 1+0.01*(p[0]+2*p[1]+3*p[2]+p[0]*p[1]+p[0]*p[1]*p[2]);};
    for (int z=0; z<s.nz; ++z) for (int y=0; y<s.ny; ++y) for (int x=0; x<s.nx; ++x) {
        Vec p{s.origin[0]+x*s.spacing[0],s.origin[1]+y*s.spacing[1],s.origin[2]+z*s.spacing[2]};
        s.cells[(z*s.ny+y)*s.nx+x].rho=rho(p);
    }
    std::vector<Vec> points{s.origin,{{-0.135,0.83,three_d?0.12:0.05}},
                           {{0.1,0.9,three_d?0.25:0.05}}};
    const auto p=lbm::sample_pressure_probes(s,points);
    for (std::size_t i=0; i<p.size(); ++i)
        check(std::abs(p[i]-s.units.pressure(rho(points[i])))<(si?1e-9:1e-14),"Pressure interpolation mismatch");
    auto invalid=points; invalid[0][0]-=0.01;
    rejects([&]{lbm::sample_pressure_probes(s,invalid);},"outside cell-center");
    invalid=points; invalid[0][2]-=0.001;
    rejects([&]{lbm::sample_pressure_probes(s,invalid);},"outside cell-center");
    invalid=points; invalid[0][0]=std::numeric_limits<double>::quiet_NaN();
    rejects([&]{lbm::sample_pressure_probes(s,invalid);},"pressureProbes[0]");
    s.cells[0].solid=true;
    rejects([&]{lbm::sample_pressure_probes(s,{s.origin});},"intersects solid");
    Vec mixed=s.origin; mixed[0]+=0.03;
    rejects([&]{lbm::sample_pressure_probes(s,{mixed});},"intersects solid");
    Vec node=s.origin; node[0]+=s.spacing[0];
    check(lbm::sample_pressure_probes(s,{node}).size()==1,"Zero-weight solid wrongly rejected");
    s.cells[0].solid=false; s.cells[0].rho=std::numeric_limits<double>::infinity();
    rejects([&]{lbm::sample_pressure_probes(s,{s.origin});},"fluid density");
    s.cells[0].rho=0;
    rejects([&]{lbm::sample_pressure_probes(s,{s.origin});},"fluid density");
    s.two_phase=true;
    rejects([&]{lbm::sample_pressure_probes(s,points);},"require singlePhase");
    check(lbm::sample_pressure_probes(s,{}).empty(),"Empty optional probes changed two-phase path");
    s.two_phase=false; s.nx=std::numeric_limits<int>::max();
    rejects([&]{lbm::sample_pressure_probes(s,points);},"grid");
}

void runner(bool three_d) {
    lbm::FoamCase c;
    c.initial.nx=c.initial.ny=6; c.initial.nz=three_d?6:1;
    c.initial.two_dimensional=!three_d; c.initial.cells.resize(36*c.initial.nz);
    c.pressure_probes={{{2.25,2.75,three_d?2.5:0}},{{3,3,three_d?3.0:0}}};
    c.steps=0;
    const char* name=three_d?"pressure_probe_test_3d":"pressure_probe_test_2d";
    const auto initial=lbm::run_foam_case(c,name);
    check(initial.steps==0 && initial.pressures.size()==2 && std::abs(initial.pressures[0]-1.0/3)<1e-14,
          "Initial probes missing");
    c.steps=7; c.write_interval=5; c.convergence.interval=2;
    c.convergence.enabled=true; c.convergence.minimum_steps=2; c.convergence.consecutive_windows=2;
    const auto steady=lbm::run_foam_case(c,name);
    check(steady.converged && steady.steps==6 && steady.pressures.size()==2 &&
          std::abs(steady.pressures[1]-1.0/3)<1e-14,"Final steady probes missing");
    c.convergence.enabled=false; c.steps=3; c.convergence.interval=10;
    const auto final=lbm::run_foam_case(c,name);
    check(final.steps==3 && final.pressures.size()==2,"Partial final probe sample missing");
    for (std::size_t i=0; i<c.initial.cells.size(); ++i)
        c.initial.cells[i].rho=1+0.001*std::cos(0.7*(i%c.initial.nx));
    const auto before=lbm::sample_pressure_probes(c.initial,c.pressure_probes);
    const auto evolving=lbm::run_foam_case(c,name);
    lbm::FoamSnapshot evolved;
    if (three_d) {
        lbm::Solver3D solver(6,6,6,c.single3d); solver.initialize_fields(c.initial);
        for (int i=0; i<c.steps; ++i) solver.step_masked_flow(0);
        evolved=solver.openfoam_snapshot();
    } else {
        lbm::Solver solver(6,6,c.single); solver.initialize_fields(c.initial);
        for (int i=0; i<c.steps; ++i) solver.step_masked_flow(0);
        evolved=solver.openfoam_snapshot();
    }
    const auto expected=lbm::sample_pressure_probes(evolved,c.pressure_probes);
    check(std::abs(evolving.pressures[0]-before[0])>1e-8,"Dynamic probe test did not evolve");
    for (std::size_t i=0; i<expected.size(); ++i)
        check(std::abs(evolving.pressures[i]-expected[i])<1e-14,"Final probes contain stale state");
    c.write_interval=1;
    const auto alternate=lbm::run_foam_case(c,name);
    check(alternate.pressures==evolving.pressures,"Write cadence affected pressure probes");
    c.pressure_probes[0][0]=-100;
    rejects([&]{lbm::run_foam_case(c,name);},"pressureProbes[0]");
    check(std::filesystem::exists(lbm::openfoam_result_path(name)/"3"/"p"),
          "Invalid probe erased previous output");
}
} // namespace

int main() {
    try {
        for (bool three_d:{false,true}) {
            for (bool si:{false,true}) interpolation(three_d,si);
            runner(three_d);
        }
        std::cout << "Pressure probes: interpolation, units, bounds, solids, evolving initial/final samples passed\n";
        return 0;
    } catch (const std::exception& e) {std::cerr << e.what() << '\n'; return 1;}
}
