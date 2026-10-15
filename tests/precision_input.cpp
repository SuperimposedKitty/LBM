#include "lbm/foam_case.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>

namespace {
void check(bool ok,const char* why) {if (!ok) throw std::runtime_error(why);}
std::string read(const std::filesystem::path& p) {
    std::ifstream f(p); check(bool(f),"Missing test file");
    return {std::istreambuf_iterator<char>(f),{}};
}
void write(const std::filesystem::path& p,const std::string& text) {
    std::ofstream f(p); f.exceptions(std::ios::badbit|std::ios::failbit); f << text;
}
std::string replace(std::string s,const std::string& from,const std::string& to) {
    const auto i=s.find(from); check(i!=std::string::npos,"Missing replacement");
    s.replace(i,from.size(),to); return s;
}
template<class F> void rejects(F f) {
    bool rejected=false; try {f();} catch (const std::exception&) {rejected=true;}
    check(rejected,"Invalid precision input accepted");
}
void copy_inputs(const std::filesystem::path& source,const std::filesystem::path& target) {
    std::filesystem::create_directories(target);
    for (const auto* dir:{"0","system"})
        std::filesystem::copy(source/dir,target/dir,std::filesystem::copy_options::recursive);
    std::filesystem::create_directories(target/"constant");
    for (const auto* file:{"lbmProperties","lbmGeometry"})
        if (std::filesystem::exists(source/"constant"/file))
            std::filesystem::copy_file(source/"constant"/file,target/"constant"/file);
}
double first_scalar(const std::filesystem::path& path) {
    auto s=read(path); auto start=s.find("internalField nonuniform");
    check(start!=std::string::npos,"Missing nonuniform output");
    std::istringstream input(s.substr(s.find('(',start)+1));
    double v; input>>v; check(bool(input),"Missing scalar value"); return v;
}
} // namespace

int main() {
    try {
        auto root=std::filesystem::current_path()/"precision_input_tests";
        int run=0; while (std::filesystem::exists(root/std::to_string(run))) ++run;
        root/=std::to_string(run);
        const auto source=std::filesystem::path(LBM_CASES_DIR)/"cylinder2d_si";
        const auto coarse=root/"coarse", fine=root/"fine";
        copy_inputs(source,coarse);
        lbm::refine_foam_case(coarse,fine,2);
        lbm::block_mesh_case(coarse); lbm::block_mesh_case(fine);
        auto a=lbm::read_foam_case(coarse), b=lbm::read_foam_case(fine);
        const auto coarse_dict=read(coarse/"system"/"lbmDict");
        write(coarse/"system"/"lbmDict",replace(coarse_dict,"openBoundary velocityPressure;",
              "openBoundary legacy;")+"\ninletProfile parabolic;\n");
        rejects([&]{lbm::read_foam_case(coarse);});
        write(coarse/"system"/"lbmDict",coarse_dict);
        check(a.steps==12000 && b.steps==48000 && b.write_interval==4*a.write_interval,
              "Physical time refinement mismatch");
        check(std::abs(a.single.tau-0.58)<1e-12 && std::abs(a.single.tau-b.single.tau)<1e-12,
              "Viscosity conversion/refinement changed tau");
        check(std::abs(a.inlet_velocity-2*b.inlet_velocity)<1e-12 &&
                  std::abs(a.reference_length*2-b.reference_length)<1e-12,
              "Diffusive scaling changed reference Reynolds number");
        check(std::abs(a.single.flow.outlet_density-1)<1e-12,"SI outlet pressure conversion failed");
        int solids=0;
        for (auto i:a.cell_order) {
            const auto& cell=a.initial.cells[i];
            check(cell.wall_distance_known,"Analytic wall distances absent");
            if (cell.solid) ++solids;
        }
        const double circle_area=3.14159265358979323846*16*16;
        check(std::abs(solids/circle_area-1)<0.02,"Cylinder geometry area mismatch");
        lbm::set_fields_case(coarse);
        check(read(coarse/"0"/"U").find("dimensions [0 1 -1 0 0 0 0]")!=std::string::npos,
              "setFields lost SI dimensions");
        rejects([&]{lbm::refine_foam_case(coarse,root/"bad_refine",2);});
        check(!std::filesystem::exists(root/"bad_refine"),"Invalid refinement wrote partial template");

        const auto initial=root/"initial";
        copy_inputs(source,initial);
        write(initial/"constant"/"lbmProperties",replace(read(initial/"constant"/"lbmProperties"),
              "referencePressure 0;","referencePressure 101325;"));
        auto p=read(initial/"0"/"p");
        p=replace(p,"internalField uniform 0;","internalField uniform 101325.3;");
        p=replace(p,"value uniform 0;","value uniform 101325;");
        write(initial/"0"/"p",p);
        write(initial/"0"/"U",replace(read(initial/"0"/"U"),"internalField uniform (0 0 0);",
                                     "internalField uniform (0.01 0.002 0);"));
        auto dict=read(initial/"system"/"lbmDict");
        dict=replace(dict,"bodyForce (0 0 0);","bodyForce (1 0 0);");
        dict=replace(dict,"acceleration (0 0 0);","acceleration (0.001 0 0);");
        write(initial/"system"/"lbmDict",dict);
        lbm::block_mesh_case(initial);
        auto c=lbm::read_foam_case(initial);
        check(std::abs(c.single.flow.force_density[0]-3.125e-8)<1e-16 &&
                  std::abs(c.single.flow.acceleration[0]-3.125e-8)<1e-16,"SI force conversion failed");
        lbm::Solver solver(c.initial.nx,c.initial.ny,c.single); solver.initialize_fields(c.initial);
        auto snapshot=solver.openfoam_snapshot();
        snapshot.origin=c.initial.origin; snapshot.spacing=c.initial.spacing; snapshot.units=c.initial.units;
        auto output=root/"output";
        {
            lbm::OpenFoamWriter writer(output,snapshot); writer.write(0,snapshot); writer.write(3,snapshot);
        }
        check(std::filesystem::exists(output/"0.0001875"/"p"),"Physical output time missing");
        check(std::abs(first_scalar(output/"0"/"p")-101325.3)<1e-8 &&
                  std::abs(first_scalar(output/"0"/"rho")-1000.225)<1e-8,
              "Physical pressure/density round trip failed");
        check(read(output/"0"/"p").find("dimensions [1 -1 -2 0 0 0 0]")!=std::string::npos,
              "Physical pressure dimensions lost");
        check(std::abs(first_scalar(output/"0"/"speed")-std::sqrt(0.01*0.01+0.002*0.002))<1e-12,
              "Physical velocity round trip failed");
        {
            lbm::OpenFoamWriter writer(output,snapshot); writer.write(0,snapshot);
        }
        check(!std::filesystem::exists(output/"0.0001875"),"Physical old time survived rerun");
        write(initial/"system"/"lbmDict",dict+"\ntau 0.8;\n");
        rejects([&]{lbm::read_foam_case(initial);});
        write(initial/"system"/"lbmDict",dict);
        write(initial/"0"/"p",replace(p,"dimensions [1 -1 -2 0 0 0 0]","dimensions [0 0 0 0 0 0 0]"));
        rejects([&]{lbm::read_foam_case(initial);});
        const auto volume=root/"si3d";
        copy_inputs(std::filesystem::path(LBM_CASES_DIR)/"closed3d",volume);
        write(volume/"constant"/"lbmProperties",
              "FoamFile { version 2.0; format ascii; class dictionary; object lbmProperties; }\n"
              "referenceDensity 1000; kinematicViscosity 16; referencePressure 0;\n");
        auto rho3d=read(volume/"0"/"rho");
        rho3d=replace(rho3d,"dimensions [0 0 0 0 0 0 0]","dimensions [1 -3 0 0 0 0 0]");
        rho3d=replace(rho3d,"internalField uniform 1;","internalField uniform 1000;");
        write(volume/"0"/"rho",rho3d);
        auto u3d=read(volume/"0"/"U");
        u3d=replace(u3d,"dimensions [0 0 0 0 0 0 0]","dimensions [0 1 -1 0 0 0 0]");
        u3d=replace(u3d,"internalField uniform (0 0 0);","internalField uniform (0.01 -0.002 0.003);");
        write(volume/"0"/"U",u3d);
        write(volume/"system"/"lbmDict",replace(read(volume/"system"/"lbmDict"),"tau 0.8;",""));
        auto control3d=read(volume/"system"/"controlDict");
        control3d=replace(control3d,"deltaT 1;","deltaT 0.001;");
        control3d=replace(control3d,"endTime 160;","endTime 0.004;");
        write(volume/"system"/"controlDict",control3d);
        lbm::block_mesh_case(volume);
        auto three=lbm::read_foam_case(volume);
        check(std::abs(three.single3d.tau-0.8)<1e-12,"SI 3D viscosity conversion failed");
        lbm::Solver3D flow3d(three.initial.nx,three.initial.ny,three.initial.nz,three.single3d);
        flow3d.initialize_fields(three.initial);
        const auto output3d=flow3d.openfoam_snapshot();
        check(output3d.units==three.initial.units && output3d.spacing==three.initial.spacing,
              "SI snapshot metadata lost");
        const auto fluid_id=three.cell_order.front();
        check(std::abs(output3d.cells[fluid_id].velocity[0]-2.5e-5)<1e-12,
              "SI 3D initial velocity conversion failed");
        flow3d.step_masked_flow(0);
        check(std::isfinite(flow3d.diagnostics().max_speed),"SI 3D step failed");
        std::cout << "SI fields, forcing, curved geometry, refinement and output times passed\n";
        return 0;
    } catch (const std::exception& e) {std::cerr << e.what() << '\n'; return 1;}
}
