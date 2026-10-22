#include "lbm/foam_case.hpp"

#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {
void check(bool ok,const char* message) {if (!ok) throw std::runtime_error(message);}

void monitor_tests() {
    lbm::ConvergenceConfig c; c.enabled=true; c.minimum_steps=4; c.consecutive_windows=2;
    lbm::ConvergenceMonitor monitor(c);
    check(!monitor.observe(2,0,0,0,true) && monitor.passes()==0,"Startup window counted");
    check(!monitor.observe(4,0,0,0,true),"Converged without consecutive windows");
    check(!monitor.observe(6,0,0,0.1,true) && monitor.passes()==0,"Force oscillation not detected");
    check(!monitor.observe(8,0,0,0,true),"Counter did not restart");
    check(!monitor.observe(10,0,0.1,0,true),"Density drift ignored");
    check(!monitor.observe(12,std::numeric_limits<double>::quiet_NaN(),0,0,false),"NaN accepted");
    check(!monitor.observe(14,0,0,0,false) && monitor.observe(16,0,0,0,false),"Rest state did not converge");
    c.enabled=false; lbm::ConvergenceMonitor disabled(c);
    for (int i=0; i<10; ++i) check(!disabled.observe(100,0,0,0,false),"Disabled monitor stopped run");
    c.velocity_tolerance=0;
    bool rejected=false; try {c.validate();} catch (const std::invalid_argument&) {rejected=true;}
    check(rejected,"Zero tolerance accepted");
}

void runner_tests(bool three_d) {
    lbm::FoamCase c;
    c.initial.nx=c.initial.ny=6; c.initial.nz=three_d?6:1;
    c.initial.two_dimensional=!three_d;
    c.initial.cells.resize(36*c.initial.nz);
    c.steps=7; c.write_interval=5;
    c.convergence.enabled=true; c.convergence.interval=2;
    c.convergence.minimum_steps=2; c.convergence.consecutive_windows=2;
    const auto name=three_d?"convergence_test_3d":"convergence_test_2d";
    const auto result=lbm::run_foam_case(c,name);
    check(result.converged && result.steps==6,"Rest-state termination step mismatch");
    const auto output=lbm::openfoam_result_path(name);
    check(std::filesystem::exists(output/"0"/"U") && std::filesystem::exists(output/"5"/"U") &&
          std::filesystem::exists(output/"6"/"U"),"Early termination did not save final frame");
    c.write_interval=3;
    const auto alternate=lbm::run_foam_case(c,name);
    check(alternate.converged && alternate.steps==result.steps,
          "Output frequency changed independent convergence monitoring");
    c.write_interval=5;
    c.convergence.enabled=false;
    const auto full=lbm::run_foam_case(c,name);
    check(!full.converged && full.steps==7,"Default endTime changed");
    check(!std::filesystem::exists(output/"6") && std::filesystem::exists(output/"7"/"U"),
          "Rerun mixed old convergence time with full run");
    c.convergence.enabled=true; c.convergence.minimum_steps=0; c.convergence.consecutive_windows=3;
    c.steps=5;
    const auto partial=lbm::run_foam_case(c,name);
    check(!partial.converged && partial.steps==5,"Partial final window counted toward convergence");
    c.convergence.enabled=false;
    c.steps=5; c.write_interval=100; c.max_mach=0.025;
    c.single.flow.acceleration[0]=c.single3d.flow.acceleration[0]=0.01;
    bool rejected=false;
    try {lbm::run_foam_case(c,name);} catch (const std::runtime_error& e) {
        rejected=std::string(e.what()).find("Computed Mach exceeds maxMach")!=std::string::npos;
    }
    check(rejected,"Mach guard missed overspeed between outputs");
}
} // namespace

int main() {
    try {
        monitor_tests(); runner_tests(false); runner_tests(true);
        std::cout << "Convergence windows, default endTime, final output and per-step Mach guard passed\n";
        return 0;
    } catch (const std::exception& e) {std::cerr << e.what() << '\n'; return 1;}
}
