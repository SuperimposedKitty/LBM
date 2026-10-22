#include "lbm/foam_case.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
lbm::FoamRunSummary simulate(int diameter, int length, double end_time) {
    if (diameter<8 || diameter>32 || diameter%4!=0 || length<6 || length>12 ||
        !std::isfinite(end_time) || end_time<1 || end_time>30)
        throw std::invalid_argument("Use diameter 8..32 (multiple of 4), length 6..12, end-time 1..30 s.");
    // 与 cylinder2d_si 使用相同物理几何/流体/均匀入口；细化时按扩散尺度同步减小 dt。
    lbm::FoamCase c;
    auto& s=c.initial;
    const double dx=0.004/diameter, dt=0.00025*64/(diameter*diameter);
    s.nx=length*diameter; s.ny=3*diameter+2;
    s.origin={0.5*dx,-0.5*dx,0.0005}; s.spacing={dx,dx,0.001};
    s.units.si=true; s.units.length_step=dx; s.units.time_step=dt;
    s.units.density_scale=1000; s.units.kinematic_viscosity=6.666666666666667e-6;
    s.units.reference_length=0.004;
    s.cells.resize(static_cast<std::size_t>(s.nx)*s.ny);
    for (int y=0; y<s.ny; ++y)
        for (int x=0; x<s.nx; ++x) {
            auto& cell=s.cells[y*s.nx+x];
            const double px=(x+0.5)*dx, py=(y-0.5)*dx;
            cell.wall_distance=std::hypot(px-0.008,py-0.006)-0.002;
            cell.wall_distance_known=true;
            cell.obstacle=cell.wall_distance<=0;
            cell.solid=y==0 || y==s.ny-1 || cell.obstacle;
            if (!cell.solid) cell.boundary=x==0?1:(x==s.nx-1?2:0);
        }
    c.single.tau=0.5+3*s.units.kinematic_viscosity*dt/(dx*dx);
    c.single.flow.open_boundary=lbm::OpenBoundaryModel::VelocityPressure;
    c.single.flow.interpolated_walls=true;
    c.inlet_velocity=0.03*dt/dx;
    c.inlet_ramp_steps=0.15/dt;
    c.reference_length=diameter;
    c.max_mach=0.1;
    c.steps=static_cast<int>(std::lround(end_time/dt));
    c.write_interval=static_cast<int>(std::lround(0.4/dt));
    c.convergence.enabled=true;
    c.convergence.interval=static_cast<int>(std::lround(0.1/dt));
    c.convergence.minimum_steps=static_cast<int>(std::lround(2.0/dt));
    c.convergence.velocity_tolerance=1e-4;
    c.convergence.density_tolerance=1e-6;
    c.convergence.force_tolerance=1e-3;
    const auto name="cylinder_validation_d"+std::to_string(diameter)+"_l"+std::to_string(length);
    std::cout << "Validation run " << name << " dx=" << dx << " dt=" << dt << '\n';
    const auto result=lbm::run_foam_case(c,name.c_str());
    std::cout << "VALIDATION D=" << diameter << " L/D=" << length << " Re=18 tau=" << c.single.tau
              << " time=" << result.steps*dt << " converged=" << result.converged
              << " Cd=" << result.drag << " Cl=" << result.lift
              << " dU=" << result.velocity_change << " dRho=" << result.density_change
              << " forceRange=" << result.force_range << " maxMach=" << result.maximum_mach << '\n';
    return result;
}

void study(double end_time) {
    const auto coarse=simulate(8,6,end_time);
    const auto medium=simulate(12,6,end_time);
    const auto fine=simulate(16,6,end_time);
    const auto longer=simulate(12,9,end_time);
    if (!coarse.converged || !medium.converged || !fine.converged || !longer.converged)
        throw std::runtime_error("Study not steady: increase end-time; do not report transient values as converged.");
    const double a=coarse.drag-medium.drag, b=medium.drag-fine.drag;
    std::cout << "STUDY Cd=" << coarse.drag << ',' << medium.drag << ',' << fine.drag
              << " mediumFineRelative=" << std::abs(b/fine.drag)
              << " outletLengthRelative=" << std::abs(longer.drag/medium.drag-1);
    // 非等比加密须使用各自格距比；无正的单调阶数时不伪造 Richardson 外推。
    auto ratio=[](double p) {return (std::pow(2.0,p)-std::pow(4.0/3,p))/(std::pow(4.0/3,p)-1);};
    if (a*b>0 && a/b>ratio(0.1) && a/b<ratio(8)) {
        double lo=0.1, hi=8;
        for (int i=0; i<60; ++i) {
            const double p=(lo+hi)/2;
            if (ratio(p)<a/b) lo=p; else hi=p;
        }
        std::cout << " observedOrder=" << (lo+hi)/2;
    } else std::cout << " observedOrder=unresolved";
    std::cout << " (same velocity/pressure BC; not DFG do-nothing)\n";
}
} // namespace

int main(int argc,char** argv) {
    try {
        std::cout << std::unitbuf << std::setprecision(12);
        int diameter=16, length=6; double end_time=12; bool all=false, geometry_option=false;
        for (int i=1; i<argc; ++i) {
            const std::string arg=argv[i];
            if (arg=="--study") {all=true; continue;}
            if ((arg!="--diameter" && arg!="--length" && arg!="--end-time") || i+1>=argc)
                throw std::invalid_argument("Usage: lbm_cylinder_validation [--study] [--diameter N] [--length N] [--end-time seconds]");
            const std::string value=argv[++i]; std::size_t used=0;
            if (arg=="--end-time") end_time=std::stod(value,&used);
            else {
                const int n=std::stoi(value,&used); (arg=="--diameter"?diameter:length)=n;
                geometry_option=true;
            }
            if (used!=value.size()) throw std::invalid_argument("Invalid numerical argument");
        }
        if (all && geometry_option) throw std::invalid_argument("--study uses fixed diameter/length levels; only --end-time may be combined.");
        if (all) study(end_time); else simulate(diameter,length,end_time);
        return 0;
    } catch (const std::exception& e) {std::cerr << e.what() << '\n'; return 1;}
}
