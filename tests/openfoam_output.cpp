#include "lbm/openfoam.hpp"
#include "lbm/solver.hpp"
#include "lbm/solver3d.hpp"
#include "lbm/two_phase_solver.hpp"
#include "lbm/two_phase_solver3d.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <regex>
#include <stdexcept>

namespace {
void require(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}
std::string read(const std::filesystem::path& p) {
    std::ifstream in(p);
    require(bool(in), "Missing output file");
    return {std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
}
std::vector<double> numbers(std::string text) {
    static const std::regex number(R"([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)");
    std::vector<double> out;
    for (std::sregex_iterator i(text.begin(), text.end(), number), end; i!=end; ++i)
        out.push_back(std::stod(i->str()));
    return out;
}
std::vector<double> mesh_list(const std::filesystem::path& p) {
    auto text = read(p);
    return numbers(text.substr(text.find('}')+1));
}
std::vector<double> field(const std::filesystem::path& p) {
    auto text = read(p);
    auto begin = text.find("internalField");
    auto values = numbers(text.substr(begin, text.find("boundaryField")-begin));
    require(!values.empty(), "Empty field");
    return values;
}
using Vec = std::array<double, 3>;
Vec sub(Vec a, Vec b) { return {a[0]-b[0],a[1]-b[1],a[2]-b[2]}; }
Vec cross(Vec a, Vec b) { return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]}; }
double dot(Vec a, Vec b) { return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }

void verify(const std::filesystem::path& path, const lbm::FoamSnapshot& s) {
    const auto mesh = path / "constant" / "polyMesh";
    auto points = mesh_list(mesh / "points"), faces = mesh_list(mesh / "faces");
    auto owner = mesh_list(mesh / "owner"), neighbour = mesh_list(mesh / "neighbour");
    int count = 0;
    for (const auto& c : s.cells) if (!c.solid) ++count;
    require(points.size()==1+3*points[0], "Point count mismatch");
    require(faces.size()==1+5*faces[0], "Face count mismatch");
    require(owner.size()==1+faces[0], "Owner count mismatch");
    require(neighbour.size()==1+neighbour[0], "Neighbour count mismatch");
    std::vector<double> volumes(count);
    std::vector<int> cell_faces(count);
    for (int i=0; i<faces[0]; ++i) {
        require(faces[1+i*5]==4, "Not a quad");
        Vec p[4], center{};
        for (int k=0; k<4; ++k) {
            int id=static_cast<int>(faces[2+i*5+k]);
            require(id>=0 && id<points[0], "Bad point label");
            p[k]={points[1+3*id],points[2+3*id],points[3+3*id]};
            for (int d=0; d<3; ++d) center[d]+=p[k][d]/4;
        }
        auto normal = cross(sub(p[1],p[0]),sub(p[3],p[0]));
        require(std::abs(dot(normal,normal)-1)<1e-12, "Bad face area");
        int own=static_cast<int>(owner[1+i]);
        require(own>=0 && own<count, "Bad owner");
        double volume=dot(center,normal)/3;
        volumes[own]+=volume;
        ++cell_faces[own];
        if (i<neighbour[0]) {
            int nei=static_cast<int>(neighbour[1+i]);
            require(nei>own && nei<count, "Internal face order");
            volumes[nei]-=volume;
            ++cell_faces[nei];
        }
    }
    for (int i=0; i<count; ++i) {
        require(cell_faces[i]==6, "Cell is not closed");
        require(std::abs(volumes[i]-1)<1e-12, "Wrong cell volume/orientation");
    }
    const auto bounds=read(mesh / "boundary");
    const std::regex patch(R"((\w+)\s*\{\s*type\s+(\w+);\s*nFaces\s+(\d+);\s*startFace\s+(\d+);)");
    int next=static_cast<int>(neighbour[0]);
    for (std::sregex_iterator i(bounds.begin(),bounds.end(),patch), end; i!=end; ++i) {
        require(std::stoi((*i)[4])==next, "Non-contiguous boundary");
        next+=std::stoi((*i)[3]);
        const auto u=read(path / "0" / "U");
        require(u.find((*i)[1].str()+"\n{")!=std::string::npos, "Missing field patch");
    }
    require(next==faces[0], "Boundary face coverage");
    require((bounds.find("type empty")!=std::string::npos)==s.two_dimensional, "Wrong empty patch");
    auto rho=field(path/"0"/"rho"), pressure=field(path/"0"/"p"), velocity=field(path/"0"/"U");
    require(rho[0]==count && rho.size()==count+1, "Scalar field count");
    require(velocity[0]==count && velocity.size()==3*count+1, "Vector field count");
    int j=0;
    for (const auto& c : s.cells) {
        if (c.solid) continue;
        require(std::abs(rho[1+j]-c.rho)<1e-14, "Density mapping");
        require(std::abs(pressure[1+j]-c.rho/3)<1e-14, "Pressure mapping");
        for (int d=0; d<3; ++d) require(std::abs(velocity[1+3*j+d]-c.velocity[d])<1e-14, "Velocity mapping");
        ++j;
    }
    if (s.two_phase) {
        auto a=field(path/"0"/"alpha.A"), b=field(path/"0"/"alpha.B");
        auto porosity=field(path/"0"/"porosity");
        auto bulk=field(path/"0"/"pBulk");
        int k=0;
        for (const auto& c : s.cells) {
            if (c.solid) continue;
            ++k;
            require(std::abs(a[k]+b[k]-1)<1e-14, "Saturation sum");
            require(std::abs(a[k]-c.rho_a/(c.rho_a+c.rho_b))<1e-14, "Saturation mapping");
            require(std::abs(porosity[k]-c.porosity)<1e-14, "Porosity mapping");
            const double expected=(c.rho+s.interaction_strength*(1-std::exp(-c.rho_a))*(1-std::exp(-c.rho_b)))/3;
            require(std::abs(bulk[k]-expected)<1e-14,"Bulk pressure mapping");
        }
    }
}
} // namespace

int main() {
    try {
        const auto root=std::filesystem::current_path()/"openfoam_test_output";
        // 每次测试使用独立目录，不删除既有文件。
        int run=0;
        while (std::filesystem::exists(root/std::to_string(run))) ++run;
        const auto base=root/std::to_string(run);
        auto check=[&](const char* name, const lbm::FoamSnapshot& s) {
            lbm::OpenFoamWriter writer(base/name,s);
            writer.write(0,s); writer.write(10,s);
            verify(base/name,s);
            require(std::filesystem::exists(base/name/"10"/"U"), "Missing time series");
            auto bad=s;
            bad.cells[0].solid=!bad.cells[0].solid;
            bool rejected=false;
            try { writer.write(20,bad); } catch (const std::invalid_argument&) { rejected=true; }
            require(rejected,"Changed topology accepted");
            require(!std::filesystem::exists(base/name/"20"),"Invalid time written");
            bad=s;
            for (auto& c : bad.cells) if (!c.solid) c.rho=std::numeric_limits<double>::quiet_NaN();
            rejected=false;
            try { writer.write(30,bad); } catch (const std::invalid_argument&) { rejected=true; }
            require(rejected,"Non-finite field accepted");
            require(!std::filesystem::exists(base/name/"30"),"Non-finite field written");
        };
        lbm::Solver flow(9,8,{});
        flow.initialize_lid_driven_cavity(0.03); flow.step_lid_driven_cavity(0.03);
        check("flow2d",flow.openfoam_snapshot());
        lbm::Solver3D flow3d(8,7,6,{});
        flow3d.initialize_lid_driven_cavity(0.03); flow3d.step_lid_driven_cavity(0.03);
        check("flow3d",flow3d.openfoam_snapshot());
        lbm::TwoPhaseSolver two(12,10,{});
        two.initialize_capillary_displacement(); two.step();
        check("two2d",two.openfoam_snapshot());
        lbm::TwoPhaseSolver3D two3d(12,12,10,{});
        two3d.initialize_droplet_impact(6,7,5,2,0,-0.01,0); two3d.step_closed();
        check("two3d",two3d.openfoam_snapshot());
        lbm::FoamSnapshot obstacle;
        obstacle.nx=4; obstacle.ny=3; obstacle.cells.resize(12);
        obstacle.two_phase=true; obstacle.interaction_strength=3;
        for (int i=0;i<12;++i) {
            auto& c=obstacle.cells[i]; c.rho_a=0.7; c.rho_b=0.3;
            c.boundary=i%4==0?1:(i%4==3?2:0); c.porosity=0.3;
        }
        obstacle.cells[5].solid=true;
        check("obstacle",obstacle);
        // 重跑使用同一目录，较短的新计算不能残留旧时间步或旧两相场。
        const auto repeat = base / "overwrite";
        std::filesystem::create_directories(repeat / "old_run_0");
        { std::ofstream note(repeat / "notes.txt"); note << "keep"; }
        { std::ofstream legacy(repeat / "old_run_0" / "archive.txt"); }
        {
            lbm::OpenFoamWriter writer(repeat,obstacle);
            writer.write(0,obstacle); writer.write(100,obstacle);
        }
        auto replacement=flow.openfoam_snapshot();
        {
            lbm::OpenFoamWriter writer(repeat,replacement);
            require(!std::filesystem::exists(repeat / "100"),"Old time survived overwrite");
            writer.write(0,replacement); writer.write(5,replacement);
        }
        require(!std::filesystem::exists(repeat / "0" / "alpha.A"),"Old phase field survived overwrite");
        require(read(repeat / "notes.txt")=="keep","User file deleted");
        require(std::filesystem::exists(repeat / "old_run_0" / "archive.txt"),"Legacy run deleted");
        verify(repeat,replacement);
        require(lbm::openfoam_result_path("test_case").filename()=="test_case","Run suffix remains");
        { std::ofstream note(repeat / "5" / "user.txt"); note << "keep"; }
        bool refused=false;
        try { lbm::OpenFoamWriter writer(repeat,replacement); }
        catch (const std::runtime_error&) { refused=true; }
        require(refused,"Unknown time file removed");
        require(std::filesystem::exists(repeat / "0" / "U"),"Cleanup occurred before preflight");
        std::cout << "OpenFOAM topology, volumes, fields, snapshots and validation passed: " << base << '\n';
        return 0;
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
