#include "lbm/openfoam.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <map>
#include <set>
#include <stdexcept>
#include <string>

namespace lbm {
namespace {
std::ofstream file(const std::filesystem::path& path) {
    std::ofstream out;
    out.exceptions(std::ios::failbit | std::ios::badbit);
    out.open(path);
    out.imbue(std::locale::classic());
    out << std::setprecision(17);
    return out;
}

void header(std::ostream& out, const char* cls, const std::string& location,
            const std::string& object) {
    out << "FoamFile\n{\n version 2.0;\n format ascii;\n class " << cls
        << ";\n location \"" << location << "\";\n object " << object << ";\n}\n\n";
}

std::size_t cell_count(const FoamSnapshot& s) {
    if (s.nx <= 0 || s.ny <= 0 || s.nz <= 0 || (s.two_dimensional && s.nz != 1))
        throw std::invalid_argument("Invalid OpenFOAM grid dimensions.");
    std::size_t count = 1;
    for (int n : {s.nx, s.ny, s.nz}) {
        if (count > static_cast<std::size_t>(std::numeric_limits<int>::max()) / n)
            throw std::length_error("OpenFOAM grid exceeds label range.");
        count *= n;
    }
    if (count != s.cells.size()) throw std::invalid_argument("OpenFOAM cell count mismatch.");
    return count;
}

// 只清理本导出器清单中的时间步，逐文件删除；保留用户文件和旧版运行编号目录。
void prepare_output(const std::filesystem::path& directory) {
    const auto manifest = directory / ".lbm-times";
    const std::set<std::string> fields{
        "U", "rho", "p", "porosity", "rhoA", "rhoB", "alpha.A", "alpha.B", "phase", "pBulk"};
    auto check_path = [&](const std::filesystem::path& path) {
        auto current = directory;
        for (const auto& part : path.lexically_relative(directory)) {
            current /= part;
            if (std::filesystem::is_symlink(std::filesystem::symlink_status(current)))
                throw std::runtime_error("Refusing linked OpenFOAM output path: " + current.string());
        }
    };
    check_path(manifest);
    std::set<std::string> times;
    const bool managed = std::filesystem::exists(manifest);
    if (managed) {
        std::ifstream in(manifest);
        std::string version, time;
        std::getline(in, version);
        if (version != "LBM_OPENFOAM_TIMES_V1") throw std::runtime_error("Invalid LBM output manifest.");
        while (std::getline(in, time)) {
            if (time.empty() || time.find_first_not_of("0123456789") != std::string::npos)
                throw std::runtime_error("Invalid time in LBM output manifest.");
            times.insert(time);
        }
        if (in.bad()) throw std::runtime_error("Cannot read LBM output manifest.");
    }
    for (const auto& relative : {"constant", "system", "lbm.foam"}) {
        const auto path = directory / relative;
        check_path(path);
        if (!managed && std::filesystem::exists(path))
            throw std::runtime_error("Output is not managed by LBM: " + path.string());
    }
    check_path(directory / "constant" / "polyMesh");
    for (const auto* name : {"points", "faces", "owner", "neighbour", "boundary"})
        check_path(directory / "constant" / "polyMesh" / name);
    check_path(directory / "system" / "controlDict");
    if (std::filesystem::exists(directory)) {
        for (const auto& entry : std::filesystem::directory_iterator(directory)) {
            const auto name = entry.path().filename().string();
            if (!name.empty() && name.find_first_not_of("0123456789") == std::string::npos && times.count(name) == 0)
                throw std::runtime_error("Unmanaged time directory: " + entry.path().string());
        }
    }
    // 先检查所有目标，再删除，避免未知文件导致清理到一半才报错。
    for (const auto& time : times) {
        const auto path = directory / time;
        check_path(path);
        if (!std::filesystem::exists(path)) continue;
        for (const auto& entry : std::filesystem::directory_iterator(path)) {
            check_path(entry.path());
            if (!entry.is_regular_file() || fields.count(entry.path().filename().string()) == 0)
                throw std::runtime_error("Unmanaged file in previous time directory: " + entry.path().string());
        }
    }
    for (const auto& time : times) {
        const auto path = directory / time;
        if (!std::filesystem::exists(path)) continue;
        for (const auto& name : fields) std::filesystem::remove(path / name);
        std::filesystem::remove(path); // 仅删除已清空的时间目录，不递归。
    }
    std::filesystem::create_directories(directory);
    auto out = file(manifest);
    out << "LBM_OPENFOAM_TIMES_V1\n";
    out.close();
}
} // namespace

std::filesystem::path openfoam_result_path(const char* case_name) {
    const std::string name = case_name ? case_name : "";
    if (name.empty() || name == "." || name == ".." ||
        name.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-") != std::string::npos)
        throw std::invalid_argument("Invalid OpenFOAM case name.");
    return std::filesystem::path(LBM_RESULT_DIR) / "openfoam" / name;
}

void validate_initial_snapshot(const FoamSnapshot& s) {
    cell_count(s);
    bool any=false;
    for(std::size_t i=0;i<s.cells.size();++i) {
        const auto& c=s.cells[i];
        for(double v:{c.rho,c.rho_a,c.rho_b,c.velocity[0],c.velocity[1],c.velocity[2],c.porosity})
            if(!std::isfinite(v))throw std::invalid_argument("Non-finite initial field");
        if(c.rho<=0 || c.rho_a<0 || c.rho_b<0 || c.boundary<0 || c.boundary>2)
            throw std::invalid_argument("Invalid initial density/boundary");
        if(c.solid)continue;
        any=true;
        if(c.porosity<=0 || c.porosity>1 || (s.two_dimensional&&c.velocity[2]!=0))
            throw std::invalid_argument("Invalid initial porosity/velocity");
        if(s.two_phase && std::abs(c.rho-c.rho_a-c.rho_b)>1e-10*std::max(1.0,c.rho))
            throw std::invalid_argument("Initial component densities do not sum to rho");
        const auto x=i%static_cast<std::size_t>(s.nx);
        if(c.boundary==1 && (x!=0 || s.nx<2 || s.cells[i+1].solid))
            throw std::invalid_argument("Initial inlet must be at x=0 with a fluid neighbour");
        if(c.boundary==2 && (x!=static_cast<std::size_t>(s.nx-1) || s.nx<2 || s.cells[i-1].solid))
            throw std::invalid_argument("Initial outlet must be at right x face with a fluid neighbour");
    }
    if(!any)throw std::invalid_argument("Initial mesh contains no fluid cells");
}

void OpenFoamWriter::validate(const FoamSnapshot& s) const {
    cell_count(s);
    if (s.nx != topology_.nx || s.ny != topology_.ny || s.nz != topology_.nz ||
        s.two_dimensional != topology_.two_dimensional || s.two_phase != topology_.two_phase ||
        s.interaction_strength != topology_.interaction_strength ||
        s.origin != topology_.origin || s.spacing != topology_.spacing)
        throw std::invalid_argument("OpenFOAM snapshot topology/model changed.");
    if (!std::isfinite(s.interaction_strength)) throw std::invalid_argument("Invalid interaction strength.");
    for (int d = 0; d < 3; ++d)
        if (!std::isfinite(s.origin[d]) || !std::isfinite(s.spacing[d]) || s.spacing[d] <= 0)
            throw std::invalid_argument("Invalid mesh coordinates.");
    for (std::size_t i = 0; i < s.cells.size(); ++i) {
        const auto& c = s.cells[i];
        if (c.solid != topology_.cells[i].solid || c.boundary != topology_.cells[i].boundary)
            throw std::invalid_argument("OpenFOAM solid/boundary mask changed.");
        if (c.solid) continue;
        for (double value : {c.velocity[0], c.velocity[1], c.velocity[2], c.rho,
                             c.rho_a, c.rho_b, c.porosity})
            if (!std::isfinite(value)) throw std::invalid_argument("Non-finite OpenFOAM field.");
        if (c.rho <= 0 || c.porosity <= 0 || c.porosity > 1 || c.rho_a < 0 || c.rho_b < 0 ||
            (s.two_phase && c.rho_a + c.rho_b <= 0))
            throw std::invalid_argument("Invalid OpenFOAM density/porosity.");
    }
}

OpenFoamWriter::OpenFoamWriter(const std::filesystem::path& directory, const FoamSnapshot& initial)
    : directory_(directory), topology_(initial),
      patches_{{"walls", "wall", {}}, {"inlet", "patch", {}}, {"outlet", "patch", {}},
               {"outer", "patch", {}}, {"frontAndBack", "empty", {}}} {
    validate(initial);
    const auto& s = initial;
    std::vector<int> labels(s.cells.size(), -1);
    for (std::size_t i = 0; i < s.cells.size(); ++i) {
        if (!s.cells[i].solid) {
            labels[i] = static_cast<int>(fluid_.size());
            fluid_.push_back(i);
        }
    }
    if (fluid_.empty()) throw std::invalid_argument("OpenFOAM mesh has no fluid cells.");
    struct Face { std::array<int, 4> vertices; int owner; int neighbour; };
    std::vector<Face> internal;
    std::array<std::vector<Face>, 5> boundary;
    std::map<std::array<int, 3>, int> point_ids;
    std::vector<std::array<int, 3>> points;
    auto point = [&](int x, int y, int z) {
        std::array<int, 3> key{x, y, z};
        auto entry = point_ids.emplace(key, static_cast<int>(points.size()));
        if (entry.second) points.push_back(key);
        return entry.first->second;
    };
    // 各面顶点从单元外侧看逆时针排列，内部面法向从 owner 指向 neighbour。
    const int delta[6][3] = {{-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}};
    const int corners[6][4][3] = {
        {{0,0,0},{0,0,1},{0,1,1},{0,1,0}}, {{1,0,0},{1,1,0},{1,1,1},{1,0,1}},
        {{0,0,0},{1,0,0},{1,0,1},{0,0,1}}, {{0,1,0},{0,1,1},{1,1,1},{1,1,0}},
        {{0,0,0},{0,1,0},{1,1,0},{1,0,0}}, {{0,0,1},{1,0,1},{1,1,1},{0,1,1}}};
    for (std::size_t id : fluid_) {
        int x = static_cast<int>(id % s.nx);
        int y = static_cast<int>((id / s.nx) % s.ny);
        int z = static_cast<int>(id / s.nx / s.ny);
        for (int d = 0; d < 6; ++d) {
            int xx=x+delta[d][0], yy=y+delta[d][1], zz=z+delta[d][2];
            bool outside = xx<0 || xx>=s.nx || yy<0 || yy>=s.ny || zz<0 || zz>=s.nz;
            int other = outside ? -1 : labels[(static_cast<std::size_t>(zz)*s.ny+yy)*s.nx+xx];
            if (other >= 0 && other < labels[id]) continue;
            Face face{{}, labels[id], other};
            for (int k=0; k<4; ++k)
                face.vertices[k] = point(x+corners[d][k][0], y+corners[d][k][1], z+corners[d][k][2]);
            if (other >= 0) internal.push_back(face);
            else {
                int patch = !outside ? 0 : 3;
                if (outside && s.two_dimensional && d >= 4) patch = 4;
                else if (outside && d == 0 && s.cells[id].boundary == 1) patch = 1;
                else if (outside && d == 1 && s.cells[id].boundary == 2) patch = 2;
                boundary[patch].push_back(face);
                patches_[patch].owners.push_back(face.owner);
            }
        }
    }
    prepare_output(directory_);
    const auto mesh = directory_ / "constant" / "polyMesh";
    std::filesystem::create_directories(mesh);
    auto pts = file(mesh / "points"); header(pts, "vectorField", "constant/polyMesh", "points");
    pts << points.size() << "\n(\n";
    for (auto p : points) pts << '(' << s.origin[0]+s.spacing[0]*(p[0]-0.5) << ' '
        << s.origin[1]+s.spacing[1]*(p[1]-0.5) << ' '
        << s.origin[2]+s.spacing[2]*(p[2]-0.5) << ")\n";
    pts << ")\n"; pts.close();
    std::size_t count = internal.size();
    for (const auto& b : boundary) count += b.size();
    auto faces = file(mesh / "faces"); header(faces, "faceList", "constant/polyMesh", "faces");
    auto owner = file(mesh / "owner"); header(owner, "labelList", "constant/polyMesh", "owner");
    auto neighbour = file(mesh / "neighbour"); header(neighbour, "labelList", "constant/polyMesh", "neighbour");
    faces << count << "\n(\n"; owner << count << "\n(\n";
    neighbour << internal.size() << "\n(\n";
    auto emit = [&](const Face& f) {
        faces << "4(" << f.vertices[0] << ' ' << f.vertices[1] << ' ' << f.vertices[2] << ' ' << f.vertices[3] << ")\n";
        owner << f.owner << '\n';
    };
    for (const auto& f : internal) { emit(f); neighbour << f.neighbour << '\n'; }
    for (const auto& b : boundary) for (const auto& f : b) emit(f);
    faces << ")\n"; owner << ")\n"; neighbour << ")\n";
    faces.close(); owner.close(); neighbour.close();
    auto bounds = file(mesh / "boundary"); header(bounds, "polyBoundaryMesh", "constant/polyMesh", "boundary");
    int nonempty = 0;
    for (const auto& p : patches_) if (!p.owners.empty()) ++nonempty;
    bounds << nonempty << "\n(\n";
    std::size_t start = internal.size();
    for (const auto& p : patches_) {
        if (p.owners.empty()) continue;
        bounds << p.name << "\n{\n type " << p.type << ";\n nFaces " << p.owners.size()
               << ";\n startFace " << start << ";\n}\n";
        start += p.owners.size();
    }
    bounds << ")\n"; bounds.close();
    std::filesystem::create_directories(directory_ / "system");
    auto control = file(directory_ / "system" / "controlDict");
    header(control, "dictionary", "system", "controlDict");
    control << "application lbmPostProcess;\nstartFrom startTime;\nstartTime 0;\nstopAt endTime;\n"
               "endTime 0;\ndeltaT 1;\nwriteControl timeStep;\nwriteInterval 1;\n"
               "writeFormat ascii;\nwritePrecision 17;\nrunTimeModifiable false;\n";
    control.close();
    auto marker = file(directory_ / "lbm.foam"); marker.close();
}

void OpenFoamWriter::write(int step, const FoamSnapshot& s) const {
    if (step < 0) throw std::invalid_argument("OpenFOAM time step must be nonnegative.");
    validate(s);
    const auto time = std::to_string(step);
    const auto path = directory_ / time;
    // 写场之前登记时间步，中途中断产生的部分文件也能在下次运行清理。
    std::ofstream manifest;
    manifest.exceptions(std::ios::failbit | std::ios::badbit);
    manifest.open(directory_ / ".lbm-times", std::ios::app);
    manifest << time << '\n';
    manifest.close();
    std::filesystem::create_directories(path);
    auto field = [&](const char* name, bool vector, auto value) {
        auto out = file(path / name);
        header(out, vector ? "volVectorField" : "volScalarField", time, name);
        // 未配置 SI 映射，所有输出量明确标为无量纲格子单位。
        out << "dimensions [0 0 0 0 0 0 0];\ninternalField nonuniform List<"
            << (vector ? "vector" : "scalar") << ">\n" << fluid_.size() << "\n(\n";
        for (auto i : fluid_) { value(out, s.cells[i]); out << '\n'; }
        out << ");\nboundaryField\n{\n";
        for (const auto& p : patches_) {
            if (p.owners.empty()) continue;
            out << p.name << "\n{\n type " << (std::string(p.type)=="empty" ? "empty" : "calculated") << ";\n";
            if (std::string(p.type)!="empty") {
                out << " value nonuniform List<" << (vector ? "vector" : "scalar") << ">\n"
                    << p.owners.size() << "\n(\n";
                for (int i : p.owners) { value(out, s.cells[fluid_[i]]); out << '\n'; }
                out << ");\n";
            }
            out << "}\n";
        }
        out << "}\n"; out.close();
    };
    field("U", true, [](auto& out, const FoamCell& c) {
        out << '(' << c.velocity[0] << ' ' << c.velocity[1] << ' ' << c.velocity[2] << ')';
    });
    field("rho", false, [](auto& out, const FoamCell& c) { out << c.rho; });
    field("p", false, [](auto& out, const FoamCell& c) { out << c.rho/3.0; });
    field("porosity", false, [](auto& out, const FoamCell& c) { out << c.porosity; });
    if (s.two_phase) {
        field("rhoA", false, [](auto& out, const FoamCell& c) { out << c.rho_a; });
        field("rhoB", false, [](auto& out, const FoamCell& c) { out << c.rho_b; });
        field("alpha.A", false, [](auto& out, const FoamCell& c) { out << c.rho_a/(c.rho_a+c.rho_b); });
        field("alpha.B", false, [](auto& out, const FoamCell& c) { out << c.rho_b/(c.rho_a+c.rho_b); });
        field("phase", false, [](auto& out, const FoamCell& c) { out << (c.rho_a-c.rho_b)/(c.rho_a+c.rho_b); });
        field("pBulk", false, [&](auto& out, const FoamCell& c) {
            out << (c.rho + s.interaction_strength*(1-std::exp(-c.rho_a))*(1-std::exp(-c.rho_b)))/3.0;
        });
    }
}
} // namespace lbm
