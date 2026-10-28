#include "lbm/foam_case.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <chrono>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <type_traits>

namespace lbm {
namespace {
using Tokens = std::vector<std::string>;
using Dict = std::map<std::string, Tokens>;
using Vec = std::array<double, 3>;
void need(bool ok, const std::string &why) {
    if (!ok)
        throw std::runtime_error(why);
}
bool near(double a, double b) {
    return std::abs(a - b) <= 1e-8 * std::max({1.0, std::abs(a), std::abs(b)});
}
bool same_position(double a, double b, double spacing) {
    return std::abs(a - b) <= 1e-8 * std::abs(spacing);
}

// 词法解析保留列表与字典层次，支持两种注释和带引号名称；不执行宏或代码。
Tokens tokenize(const std::filesystem::path &path) {
    std::ifstream in(path, std::ios::binary);
    need(bool(in), "Cannot open " + path.string());
    std::string text{std::istreambuf_iterator<char>(in), {}};
    Tokens out;
    for (std::size_t i = 0; i < text.size();) {
        if (std::isspace(static_cast<unsigned char>(text[i]))) {
            ++i;
            continue;
        }
        if (text.compare(i, 2, "//") == 0) {
            auto e = text.find('\n', i);
            i = e == std::string::npos ? text.size() : e;
            continue;
        }
        if (text.compare(i, 2, "/*") == 0) {
            auto e = text.find("*/", i + 2);
            need(e != std::string::npos, "Unclosed comment: " + path.string());
            i = e + 2;
            continue;
        }
        need(text[i] != '#' && text[i] != '$' && text[i] != '\0',
             "Unsupported directive, macro or binary data: " + path.string());
        if (text[i] == '"') {
            auto e = text.find('"', i + 1);
            need(e != std::string::npos, "Unclosed quoted word");
            out.push_back(text.substr(i + 1, e - i - 1));
            i = e + 1;
            continue;
        }
        if (std::string("{}()[];").find(text[i]) != std::string::npos) {
            out.push_back(text.substr(i++, 1));
            continue;
        }
        auto begin = i;
        while (i < text.size() && text.compare(i, 2, "//") != 0 && text.compare(i, 2, "/*") != 0 &&
               !std::isspace(static_cast<unsigned char>(text[i])) &&
               std::string("{}()[];\"#$").find(text[i]) == std::string::npos)
            ++i;
        need(i > begin, "Unsupported token");
        out.push_back(text.substr(begin, i - begin));
    }
    return out;
}
struct Cursor {
    Tokens t;
    std::size_t i = 0;
    explicit Cursor(Tokens tokens) : t(std::move(tokens)) {}
    std::string get() {
        need(i < t.size(), "Unexpected end of OpenFOAM file");
        return t[i++];
    }
    bool take(const std::string &s) {
        if (i < t.size() && t[i] == s) {
            ++i;
            return true;
        }
        return false;
    }
    void expect(const std::string &s) { need(get() == s, "Expected token " + s); }
    double scalar() {
        auto s = get();
        std::size_t n = 0;
        double v = std::stod(s, &n);
        need(n == s.size() && std::isfinite(v), "Invalid scalar " + s);
        return v;
    }
    int label() {
        double v = scalar();
        need(v >= 0 && v <= 100000000 && std::floor(v) == v, "Invalid/nonintegral label");
        return static_cast<int>(v);
    }
    Vec vec() {
        expect("(");
        Vec v{scalar(), scalar(), scalar()};
        expect(")");
        return v;
    }
    Tokens group(const std::string &open, const std::string &close) {
        expect(open);
        Tokens out;
        int depth = 1;
        while (depth) {
            auto v = get();
            if (v == open)
                ++depth;
            if (v == close)
                --depth;
            if (depth)
                out.push_back(v);
        }
        return out;
    }
    Dict dict() {
        Dict d;
        while (i < t.size()) {
            auto key = get();
            if (key == ";")
                continue;
            Tokens value;
            if (i < t.size() && t[i] == "{") {
                value = group("{", "}");
                take(";");
            } else {
                int paren = 0, bracket = 0, brace = 0;
                while (true) {
                    auto v = get();
                    if (v == ";" && !paren && !bracket && !brace)
                        break;
                    if (v == "(")
                        ++paren;
                    if (v == ")")
                        --paren;
                    if (v == "[")
                        ++bracket;
                    if (v == "]")
                        --bracket;
                    if (v == "{")
                        ++brace;
                    if (v == "}")
                        --brace;
                    need(paren >= 0 && bracket >= 0 && brace >= 0, "Unbalanced dictionary");
                    value.push_back(v);
                }
            }
            need(d.emplace(key, std::move(value)).second, "Duplicate entry: " + key);
        }
        return d;
    }
};
void skip_header(Cursor &c) {
    c.expect("FoamFile");
    auto h = Cursor(c.group("{", "}")).dict();
    need(h.count("format") && h.at("format") == Tokens{"ascii"},
         "Only ASCII OpenFOAM files are supported");
}
Dict dictionary(const std::filesystem::path &p) {
    Cursor c(tokenize(p));
    skip_header(c);
    return c.dict();
}
const Tokens &required(const Dict &d, const std::string &key) {
    auto it = d.find(key);
    need(it != d.end(), "Missing entry: " + key);
    return it->second;
}
std::string word(const Dict &d, const std::string &key, const std::string &fallback = "") {
    auto it = d.find(key);
    if (it == d.end())
        return fallback;
    need(it->second.size() == 1, "Expected word: " + key);
    return it->second[0];
}
double number(const Dict &d, const std::string &key, double fallback) {
    auto it = d.find(key);
    if (it == d.end())
        return fallback;
    Cursor c(it->second);
    auto v = c.scalar();
    need(c.i == c.t.size(), "Expected scalar: " + key);
    return v;
}
void known(const Dict &d, const std::set<std::string> &keys) {
    for (const auto &e : d)
        need(keys.count(e.first) != 0, "Unsupported entry: " + e.first);
}
int integer(double value, const char *key, bool positive = false) {
    need(value >= (positive ? 1 : 0) && value <= 100000000 && std::floor(value) == value,
         std::string("Invalid integer: ") + key);
    return static_cast<int>(value);
}
std::vector<Vec> values(const Tokens &t, int components, std::size_t count) {
    Cursor c(t);
    std::vector<Vec> out;
    auto one = [&]() { return components == 3 ? c.vec() : Vec{c.scalar(), 0, 0}; };
    if (c.take("uniform"))
        out.assign(count, one());
    else {
        c.expect("nonuniform");
        c.expect(components == 3 ? "List<vector>" : "List<scalar>");
        need(static_cast<std::size_t>(c.label()) == count, "Field list size does not match mesh");
        c.expect("(");
        for (std::size_t i = 0; i < count; ++i)
            out.push_back(one());
        c.expect(")");
    }
    need(c.i == c.t.size(), "Extra field tokens");
    return out;
}
struct Field {
    Dict dict;
    std::vector<Vec> data;
    Dict patches;
};
Field read_field(const std::filesystem::path &p, int components, std::size_t count,
                 bool physical = false) {
    Cursor c(tokenize(p));
    c.expect("FoamFile");
    auto h = Cursor(c.group("{", "}")).dict();
    need(word(h, "format") == "ascii", "Only ASCII fields supported");
    need(word(h, "class") == (components == 3 ? "volVectorField" : "volScalarField"),
         "Wrong field class: " + p.string());
    auto d = c.dict();
    Tokens dimensions{"[", "0", "0", "0", "0", "0", "0", "0", "]"};
    const auto name = word(h,"object");
    if (physical && name=="U") dimensions={"[","0","1","-1","0","0","0","0","]"};
    if (physical && name=="rho") dimensions={"[","1","-3","0","0","0","0","0","]"};
    if (physical && name=="p") dimensions={"[","1","-1","-2","0","0","0","0","]"};
    need(required(d, "dimensions") == dimensions, "Unexpected field dimensions: " + p.string());
    return {d, values(required(d, "internalField"), components, count),
            Cursor(required(d, "boundaryField")).dict()};
}

template<class Transform> void transform_field(Field& field, int components, Transform transform) {
    for (auto& v : field.data) v=transform(v);
    for (auto& entry : field.patches) {
        auto b=Cursor(entry.second).dict();
        if (word(b,"type")!="fixedValue") continue;
        const auto v=transform(values(required(b,"value"),components,1)[0]);
        Tokens result{"type","fixedValue",";","value","uniform"};
        if (components==3) result.push_back("(");
        for (int d=0; d<components; ++d) {
            std::ostringstream out;
            out << std::setprecision(17) << v[d];
            result.push_back(out.str());
        }
        if (components==3) result.push_back(")");
        result.push_back(";");
        entry.second=std::move(result);
    }
}
std::vector<int> labels(const std::filesystem::path &p) {
    Cursor c(tokenize(p));
    skip_header(c);
    auto n = c.label();
    c.expect("(");
    std::vector<int> v;
    for (int i = 0; i < n; ++i)
        v.push_back(c.label());
    c.expect(")");
    need(c.i == c.t.size(), "Extra labels");
    return v;
}
struct Patch {
    std::string name, type;
    int start = 0, count = 0;
};
struct Mesh {
    std::vector<Vec> points;
    std::vector<std::array<int, 4>> faces;
    std::vector<int> owner, neighbour;
    std::vector<Patch> patches;
    std::vector<Vec> centers;
    Vec spacing{};
    bool two_d = false;
};
Vec face_center(const Mesh &m, int f) {
    Vec c{};
    for (auto p : m.faces[f])
        for (int d = 0; d < 3; ++d)
            c[d] += m.points[p][d] / 4;
    return c;
}
Mesh mesh_read(const std::filesystem::path &path) {
    Mesh m;
    const auto base = path / "constant" / "polyMesh";
    Cursor p(tokenize(base / "points"));
    skip_header(p);
    int np = p.label();
    p.expect("(");
    for (int i = 0; i < np; ++i)
        m.points.push_back(p.vec());
    p.expect(")");
    Cursor f(tokenize(base / "faces"));
    skip_header(f);
    int nf = f.label();
    f.expect("(");
    for (int i = 0; i < nf; ++i) {
        need(f.label() == 4, "Only quadrilateral Cartesian faces supported");
        f.expect("(");
        std::array<int, 4> face{};
        std::set<int> unique;
        for (auto &v : face) {
            v = f.label();
            need(v < np, "Bad point index");
            unique.insert(v);
        }
        need(unique.size() == 4, "Repeated face vertex");
        f.expect(")");
        m.faces.push_back(face);
    }
    f.expect(")");
    m.owner = labels(base / "owner");
    m.neighbour = labels(base / "neighbour");
    need(m.owner.size() == m.faces.size() && !m.owner.empty() &&
             m.neighbour.size() <= m.owner.size(),
         "Invalid face addressing");
    const int nc = *std::max_element(m.owner.begin(), m.owner.end()) + 1;
    need(nc <= nf, "Invalid cell labels");
    std::vector<std::set<int>> vertices(nc);
    std::vector<int> counts(nc);
    for (int i = 0; i < nf; ++i) {
        ++counts[m.owner[i]];
        vertices[m.owner[i]].insert(m.faces[i].begin(), m.faces[i].end());
        if (i < static_cast<int>(m.neighbour.size())) {
            int n = m.neighbour[i];
            need(n < nc && n > m.owner[i], "Invalid neighbour ordering");
            ++counts[n];
            vertices[n].insert(m.faces[i].begin(), m.faces[i].end());
        }
    }
    m.centers.resize(nc);
    for (int i = 0; i < nc; ++i) {
        need(counts[i] == 6 && vertices[i].size() == 8,
             "Only conforming hexahedral cells supported");
        Vec lo{1e100, 1e100, 1e100}, hi{-1e100, -1e100, -1e100};
        for (int v : vertices[i])
            for (int d = 0; d < 3; ++d) {
                lo[d] = std::min(lo[d], m.points[v][d]);
                hi[d] = std::max(hi[d], m.points[v][d]);
            }
        for (int d = 0; d < 3; ++d) {
            auto h = hi[d] - lo[d];
            need(h > 1e-10, "Degenerate cell");
            if (i == 0)
                m.spacing[d] = h;
            else
                need(near(h / m.spacing[d], 1),
                     "Graded/nonuniform meshes are not supported by this LBM");
            m.centers[i][d] = (lo[d] + hi[d]) / 2;
            for (int v : vertices[i])
                need(same_position(m.points[v][d], lo[d], h) ||
                         same_position(m.points[v][d], hi[d], h),
                     "Skewed/curved cells not supported");
        }
    }
    Cursor b(tokenize(base / "boundary"));
    skip_header(b);
    int nb = b.label();
    b.expect("(");
    int next = static_cast<int>(m.neighbour.size());
    std::set<std::string> names;
    for (int i = 0; i < nb; ++i) {
        Patch patch;
        patch.name = b.get();
        need(names.insert(patch.name).second, "Duplicate patch");
        auto d = Cursor(b.group("{", "}")).dict();
        patch.type = word(d, "type");
        need(patch.type == "wall" || patch.type == "patch" || patch.type == "empty",
             "Unsupported mesh patch type: " + patch.type);
        patch.start = integer(number(d, "startFace", -1), "startFace");
        patch.count = integer(number(d, "nFaces", -1), "nFaces");
        need(patch.start == next && patch.count <= nf - next, "Invalid patch range");
        next += patch.count;
        m.two_d = m.two_d || (patch.type == "empty" && patch.count > 0);
        m.patches.push_back(patch);
    }
    b.expect(")");
    need(next == nf, "Boundary coverage incomplete");
    need(near(m.spacing[0] / m.spacing[1], 1) && (m.two_d || near(m.spacing[0] / m.spacing[2], 1)),
         "LBM requires equal spacing in active dimensions");
    // 检查面几何、朝向与内面邻接，拒绝折叠单元或拓扑相邻但几何不相邻的网格。
    std::vector<std::set<int>> cell_sides(nc);
    for (int i = 0; i < nf; ++i) {
        auto c = face_center(m, i);
        int axis = -1;
        for (int d = 0; d < 3; ++d)
            if (near(std::abs(c[d] - m.centers[m.owner[i]][d]) / m.spacing[d], 0.5)) {
                need(axis < 0, "Nonplanar face");
                axis = d;
            }
        need(axis >= 0, "Face not on cell surface");
        const int side = 2 * axis + (c[axis] > m.centers[m.owner[i]][axis] ? 1 : 0);
        need(cell_sides[m.owner[i]].insert(side).second, "Duplicate cell face");
        if (i < static_cast<int>(m.neighbour.size()))
            need(cell_sides[m.neighbour[i]].insert(side ^ 1).second, "Duplicate neighbour face");
        for (int v : m.faces[i])
            need(same_position(m.points[v][axis], c[axis], m.spacing[axis]), "Nonplanar face");
        auto a = m.points[m.faces[i][0]], v = m.points[m.faces[i][1]], w = m.points[m.faces[i][3]];
        for (int d = 0; d < 3; ++d) {
            v[d] -= a[d];
            w[d] -= a[d];
        }
        Vec normal{v[1] * w[2] - v[2] * w[1], v[2] * w[0] - v[0] * w[2], v[0] * w[1] - v[1] * w[0]};
        need(normal[axis] * (c[axis] - m.centers[m.owner[i]][axis]) > 0,
             "Inverted face orientation");
        if (i < static_cast<int>(m.neighbour.size()))
            for (int d = 0; d < 3; ++d)
                need(same_position(m.centers[m.neighbour[i]][d],
                                   d == axis ? 2 * c[d] - m.centers[m.owner[i]][d]
                                             : m.centers[m.owner[i]][d],
                                   m.spacing[d]),
                     "Non-Cartesian internal neighbour");
    }
    for (const auto &sides : cell_sides)
        need(sides.size() == 6, "Cell surface is not closed");
    return m;
}

struct AnalyticSolid {
    std::string type;
    Vec centre{}, axis{0,0,1}, lo{}, hi{};
    double radius=0;
    double distance(Vec p) const {
        if (type=="box") {
            double outside=0, inside=-1e100;
            for (int d=0; d<3; ++d) {
                const double v=std::max(lo[d]-p[d],p[d]-hi[d]);
                outside+=std::max(v,0.0)*std::max(v,0.0);
                inside=std::max(inside,v);
            }
            return std::sqrt(outside)+std::min(inside,0.0);
        }
        double squared=0, parallel=0;
        for (int d=0; d<3; ++d) {p[d]-=centre[d]; squared+=p[d]*p[d]; parallel+=p[d]*axis[d];}
        if (type=="cylinder") squared-=parallel*parallel;
        return std::sqrt(std::max(squared,0.0))-radius;
    }
};

std::vector<AnalyticSolid> read_geometry(const std::filesystem::path& path) {
    const auto file=path/"constant"/"lbmGeometry";
    if (!std::filesystem::exists(file)) return {};
    auto dict=dictionary(file);
    known(dict,{"solids"});
    Cursor cursor(required(dict,"solids"));
    cursor.expect("(");
    std::vector<AnalyticSolid> solids;
    while (!cursor.take(")")) {
        AnalyticSolid shape;
        shape.type=cursor.get();
        auto d=Cursor(cursor.group("{","}")).dict();
        if (shape.type=="box") {
            known(d,{"min","max"});
            Cursor a(required(d,"min")),b(required(d,"max"));
            shape.lo=a.vec(); shape.hi=b.vec();
            need(a.i==a.t.size() && b.i==b.t.size(),"Invalid box coordinates");
            for (int k=0; k<3; ++k) need(shape.hi[k]>shape.lo[k],"Invalid solid box");
        } else {
            need(shape.type=="sphere" || shape.type=="cylinder","Unknown analytic solid");
            known(d,{"centre","radius","axis"});
            Cursor p(required(d,"centre")); shape.centre=p.vec();
            need(p.i==p.t.size(),"Invalid solid centre");
            shape.radius=number(d,"radius",-1);
            need(shape.radius>0,"Solid radius must be positive");
            if (d.count("axis")) {
                need(shape.type=="cylinder","Only cylinders use an axis");
                Cursor a(d.at("axis")); shape.axis=a.vec();
                need(a.i==a.t.size(),"Invalid cylinder axis");
                double norm=std::sqrt(shape.axis[0]*shape.axis[0]+shape.axis[1]*shape.axis[1]+shape.axis[2]*shape.axis[2]);
                need(norm>0,"Cylinder axis cannot be zero");
                for (auto& v:shape.axis) v/=norm;
            }
        }
        solids.push_back(shape);
    }
    need(cursor.i==cursor.t.size(),"Extra geometry tokens");
    return solids;
}

void configure(FoamCase &c, const std::filesystem::path &path) {
    auto d = dictionary(path / "system" / "lbmDict");
    known(d, {"model", "lattice", "collision", "tau", "tauA", "tauB", "interactionStrength",
              "contactAngle", "wallAdhesion", "recoloring", "bodyForce", "poreDiameter",
              "darcyDrag", "forchheimerDrag", "acceleration", "openBoundary", "wallBoundary",
              "inletProfile", "maxMach", "densityWarning", "referenceLength", "inletRampTime", "convergence",
              "pressureProbes"});
    const auto model = word(d, "model", "singlePhase"), lattice = word(d, "lattice", "D2Q9"),
               collision = word(d, "collision", "MRT");
    need(model == "singlePhase" || model == "twoPhase", "Unknown LBM model");
    need(lattice == "D2Q9" || lattice == "D3Q19" || lattice == "D3Q27", "Unknown lattice");
    need(collision == "MRT" || collision == "BGK", "Unknown collision model");
    c.initial.two_phase = model == "twoPhase";
    c.initial.two_dimensional = lattice == "D2Q9";
    if (d.count("pressureProbes")) {
        Cursor probes(d.at("pressureProbes"));
        probes.expect("(");
        while (!probes.take(")")) c.pressure_probes.push_back(probes.vec());
        need(probes.i==probes.t.size(),"Extra pressureProbes tokens");
        need(!c.initial.two_phase || c.pressure_probes.empty(),"pressureProbes currently require singlePhase");
    }
    c.single.collision_model = c.single3d.collision_model = c.two.collision_model =
        c.two3d.collision_model = collision == "MRT" ? CollisionModel::MRT : CollisionModel::BGK;
    c.single3d.lattice_model = c.two3d.lattice_model =
        lattice == "D3Q27" ? Lattice3DModel::D3Q27 : Lattice3DModel::D3Q19;
    c.single.tau = c.single3d.tau = number(d, "tau", 0.8);
    const auto open=word(d,"openBoundary","legacy"), wall=word(d,"wallBoundary","halfway"),
               profile=word(d,"inletProfile","uniform");
    need(open=="legacy" || open=="velocityPressure","Unknown openBoundary");
    need(wall=="halfway" || wall=="interpolated","Unknown wallBoundary");
    need(profile=="uniform" || profile=="parabolic","Unknown inletProfile");
    need(!c.initial.two_phase || (open=="legacy" && wall=="halfway" && profile=="uniform"),
         "Precision boundary options currently require singlePhase");
    need(profile=="uniform" || c.initial.two_dimensional,"Parabolic inlet currently requires D2Q9");
    need(profile=="uniform" || open=="velocityPressure","Parabolic inlet requires velocityPressure");
    c.single.flow.open_boundary=c.single3d.flow.open_boundary=
        open=="legacy"?OpenBoundaryModel::Legacy:OpenBoundaryModel::VelocityPressure;
    c.single.flow.interpolated_walls=c.single3d.flow.interpolated_walls=wall=="interpolated";
    c.parabolic_inlet=profile=="parabolic";
    c.max_mach=number(d,"maxMach",0.15);
    c.density_warning=number(d,"densityWarning",0.01);
    c.reference_length=number(d,"referenceLength",0);
    need(c.max_mach>0 && c.max_mach<0.5 && c.density_warning>0 && c.reference_length>=0,
         "Invalid precision diagnostic controls");
    c.two.tau_a = c.two3d.tau_a = number(d, "tauA", 1);
    c.two.tau_b = c.two3d.tau_b = number(d, "tauB", 1);
    c.initial.interaction_strength = c.two.interaction_strength = c.two3d.interaction_strength =
        number(d, "interactionStrength", 3);
    c.two.contact_angle_degrees = c.two3d.contact_angle_degrees = number(d, "contactAngle", 90);
    c.two.wall_adhesion_strength = c.two3d.wall_adhesion_strength = number(d, "wallAdhesion", 0.08);
    c.two.recoloring_strength = c.two3d.recoloring_strength = number(d, "recoloring", 0);
    c.two.porous_pore_diameter = c.two3d.porous_pore_diameter = number(d, "poreDiameter", 40);
    c.two.darcy_drag_scale = c.two3d.darcy_drag_scale = number(d, "darcyDrag", 0.06);
    c.two.forchheimer_drag_scale = c.two3d.forchheimer_drag_scale =
        number(d, "forchheimerDrag", 0.02);
    Vec force{};
    if (d.count("bodyForce")) {
        Cursor v(d.at("bodyForce"));
        force = v.vec();
        need(v.i == v.t.size(), "Invalid bodyForce");
    }
    c.single.flow.force_density=c.single3d.flow.force_density=force;
    Vec acceleration{};
    if (d.count("acceleration")) {
        Cursor a(d.at("acceleration")); acceleration=a.vec();
        need(a.i==a.t.size(),"Invalid acceleration");
    }
    need(!c.initial.two_phase || acceleration==Vec{},"Two-phase acceleration option is not implemented");
    c.single.flow.acceleration=c.single3d.flow.acceleration=acceleration;
    need(!c.initial.two_dimensional || force[2] == 0, "D2Q9 cannot use z force");
    c.two.body_force_x = c.two3d.body_force_x = force[0];
    c.two.body_force_y = c.two3d.body_force_y = force[1];
    c.two3d.body_force_z = force[2];
    auto& units=c.initial.units;
    const auto physical=path/"constant"/"lbmProperties";
    if (std::filesystem::exists(physical)) {
        auto p=dictionary(physical);
        known(p,{"referenceDensity","kinematicViscosity","referencePressure","referenceLength"});
        need(!c.initial.two_phase,"SI conversion currently requires singlePhase");
        need(!d.count("tau"),"SI tau is derived from viscosity and deltaT; remove tau");
        units.si=true;
        units.density_scale=number(p,"referenceDensity",-1);
        units.kinematic_viscosity=number(p,"kinematicViscosity",-1);
        units.reference_pressure=number(p,"referencePressure",0);
        units.reference_length=number(p,"referenceLength",c.reference_length);
    }
    auto control = dictionary(path / "system" / "controlDict");
    units.time_step=number(control,"deltaT",1);
    need(number(control, "startTime", 0) == 0 && (units.si || units.time_step == 1) &&
             word(control, "startFrom", "startTime") == "startTime",
         "Use startTime 0; lattice mode requires deltaT 1");
    units.validate();
    c.inlet_ramp_steps=number(d,"inletRampTime",0)/units.time_step;
    need(c.inlet_ramp_steps>=0,"inletRampTime cannot be negative");
    need(!c.initial.two_phase || c.inlet_ramp_steps==0,"Inlet ramp currently requires singlePhase");
    need(word(control, "writeControl", "timeStep") == "timeStep", "Use writeControl timeStep");
    need(word(control, "stopAt", "endTime") == "endTime", "Use stopAt endTime");
    need(word(control, "writeFormat", "ascii") == "ascii", "Only ASCII output is supported");
    const double requested_steps=number(control,"endTime",100)/units.time_step;
    need(near(requested_steps,std::round(requested_steps)),"endTime must be a multiple of deltaT");
    c.steps = integer(std::round(requested_steps), "endTime/deltaT");
    c.write_interval = integer(number(control, "writeInterval", 10), "writeInterval", true);
    if (d.count("convergence")) {
        auto settings=Cursor(d.at("convergence")).dict();
        known(settings,{"enabled","interval","minTime","consecutiveWindows","velocityTolerance",
                        "densityTolerance","forceTolerance"});
        const auto enabled=word(settings,"enabled","false");
        need(enabled=="true" || enabled=="false","convergence.enabled must be true or false");
        auto& s=c.convergence;
        s.enabled=enabled=="true";
        s.interval=integer(number(settings,"interval",0),"convergence.interval");
        const double minimum=number(settings,"minTime",0)/units.time_step;
        need(minimum>=0 && minimum<=c.steps,"convergence.minTime must be within the simulation");
        s.minimum_steps=integer(std::ceil(minimum),"convergence.minTime/deltaT");
        s.consecutive_windows=integer(number(settings,"consecutiveWindows",5),"convergence.consecutiveWindows",true);
        s.velocity_tolerance=number(settings,"velocityTolerance",s.velocity_tolerance);
        s.density_tolerance=number(settings,"densityTolerance",s.density_tolerance);
        s.force_tolerance=number(settings,"forceTolerance",s.force_tolerance);
        s.validate();
        need(!s.enabled || !c.initial.two_phase,"Steady-state termination currently requires singlePhase");
    }
}
} // namespace

FoamCase read_foam_case(const std::filesystem::path &path) {
    FoamCase c;
    configure(c, path);
    const auto mesh = mesh_read(path);
    auto &s = c.initial;
    need(s.two_dimensional == mesh.two_d, "Mesh empty patches do not match selected lattice");
    auto& units=s.units;
    if (units.si) {
        units.length_step=mesh.spacing[0]; units.validate();
        c.single.tau=c.single3d.tau=0.5+3*units.kinematic_viscosity*units.time_step/
                                              (units.length_step*units.length_step);
        const double force_scale=units.time_step*units.time_step/(units.density_scale*units.length_step);
        const double acceleration_scale=units.time_step*units.time_step/units.length_step;
        for (int d=0; d<3; ++d) {
            c.single.flow.force_density[d]*=force_scale;
            c.single.flow.acceleration[d]*=acceleration_scale;
        }
        c.single3d.flow=c.single.flow;
    }
    auto u = read_field(path / "0" / "U", 3, mesh.centers.size(),units.si);
    transform_field(u,3,[&](Vec v) {for(auto& a:v) a/=units.velocity_scale(); return v;});
    const bool has_rho = std::filesystem::exists(path / "0" / "rho");
    auto rho = read_field(path / "0" / (has_rho ? "rho" : "p"), 1, mesh.centers.size(),units.si);
    transform_field(rho,1,[&](Vec v) {
        v[0]=has_rho?v[0]/(units.si?units.density_scale:1):units.rho_from_pressure(v[0]); return v;
    });
    if (has_rho && std::filesystem::exists(path / "0" / "p")) {
        auto pressure = read_field(path / "0" / "p", 1, mesh.centers.size(),units.si);
        for (std::size_t i = 0; i < rho.data.size(); ++i)
            need(near(rho.data[i][0], units.rho_from_pressure(pressure.data[i][0])),
                 "Initial p and rho disagree with the reference equation of state");
        for (const auto& entry:rho.patches) {
            const auto rb=Cursor(entry.second).dict();
            const auto pb=Cursor(required(pressure.patches,entry.first)).dict();
            need(word(rb,"type")==word(pb,"type"),"p/rho boundary types disagree");
            if (word(rb,"type")=="fixedValue")
                need(near(values(required(rb,"value"),1,1)[0][0],
                          units.rho_from_pressure(values(required(pb,"value"),1,1)[0][0])),
                     "p/rho boundary values disagree");
        }
    }
    Field alpha, porosity, solid;
    if (s.two_phase)
        alpha = read_field(path / "0" / "alpha.A", 1, mesh.centers.size());
    bool porous = std::filesystem::exists(path / "0" / "porosity");
    if (porous)
        porosity = read_field(path / "0" / "porosity", 1, mesh.centers.size());
    const bool has_solid = std::filesystem::exists(path / "0" / "solid");
    if (has_solid)
        solid = read_field(path / "0" / "solid", 1, mesh.centers.size());
    Vec lo{1e100, 1e100, 1e100}, hi{-1e100, -1e100, -1e100};
    for (auto center : mesh.centers)
        for (int d = 0; d < 3; ++d) {
            lo[d] = std::min(lo[d], center[d]);
            hi[d] = std::max(hi[d], center[d]);
        }
    if (s.two_dimensional)
        need(same_position(lo[2], hi[2], mesh.spacing[2]), "2D mesh must have exactly one z layer");
    bool open = false;
    for (const auto &p : mesh.patches)
        if (p.count && p.type == "patch")
            open = true;
    std::array<int, 3> pad{open ? 0 : 1, 1, s.two_dimensional ? 0 : 1}, n{};
    std::size_t total = 1;
    for (int d = 0; d < 3; ++d) {
        n[d] =
            integer(std::round((hi[d] - lo[d]) / mesh.spacing[d]), "grid extent") + 1 + 2 * pad[d];
        need(n[d] <= 2000000 && total <= 2000000 / static_cast<std::size_t>(n[d]),
             "Imported dense lattice exceeds 2 million cells");
        total *= n[d];
        s.origin[d] = lo[d] - pad[d] * mesh.spacing[d];
        s.spacing[d] = mesh.spacing[d];
    }
    s.nx = n[0];
    s.ny = n[1];
    s.nz = n[2];
    need(s.nx > 4 && s.ny > 4 && (s.two_dimensional || s.nz > 4),
         "Mesh is too small (need at least 5 lattice nodes per active axis)");
    s.cells.resize(total);
    for (auto &cell : s.cells) {
        cell.solid = true;
        cell.porosity = 0;
        cell.rho_a = cell.rho_b = 0.01;
    }
    auto index = [&](const Vec &center) {
        std::array<int, 3> k{};
        for (int d = 0; d < 3; ++d) {
            auto v = (center[d] - s.origin[d]) / s.spacing[d];
            need(near(v, std::round(v)), "Cell centres are not on a uniform lattice");
            k[d] = static_cast<int>(std::round(v));
            need(k[d] >= 0 && k[d] < n[d], "Invalid mapped cell");
        }
        return (static_cast<std::size_t>(k[2]) * s.ny + k[1]) * s.nx + k[0];
    };
    for (std::size_t i = 0; i < mesh.centers.size(); ++i) {
        auto id = index(mesh.centers[i]);
        auto &cell = s.cells[id];
        need(cell.solid, "Duplicate cell centre");
        c.cell_order.push_back(id);
        cell.solid = false;
        cell.velocity = u.data[i];
        cell.rho = rho.data[i][0];
        cell.porosity = porous ? porosity.data[i][0] : 1;
        need(cell.rho > 0 && cell.porosity > 0 && cell.porosity <= 1, "Invalid density/porosity");
        need(s.two_phase || near(cell.porosity, 1), "Single-phase porous drag is not implemented");
        need(!s.two_dimensional || cell.velocity[2] == 0, "D2Q9 initial U.z must be zero");
        need(std::sqrt(cell.velocity[0] * cell.velocity[0] + cell.velocity[1] * cell.velocity[1] +
                       cell.velocity[2] * cell.velocity[2]) < 0.2,
             "Initial velocity must be low Mach (<0.2 lattice units)");
        if (s.two_phase) {
            double a = alpha.data[i][0];
            need(a >= 0 && a <= 1, "alpha.A outside [0,1]");
            cell.rho_a = cell.rho * a;
            cell.rho_b = cell.rho * (1 - a);
        }
    }
    if (has_solid) {
        for (std::size_t i = 0; i < mesh.centers.size(); ++i) {
            const double value = solid.data[i][0];
            need(near(value, 0) || near(value, 1), "solid must contain only 0 or 1");
            s.cells[c.cell_order[i]].solid = near(value, 1);
            s.cells[c.cell_order[i]].obstacle = near(value, 1);
        }
        need(std::any_of(s.cells.begin(), s.cells.end(),
                         [](const FoamCell &cell) { return !cell.solid; }),
             "solid field leaves no fluid cells");
    }
    const auto geometry=read_geometry(path);
    need(!c.single.flow.interpolated_walls || !geometry.empty(),
         "Interpolated walls require constant/lbmGeometry");
    for (std::size_t id=0; id<s.cells.size() && !geometry.empty(); ++id) {
        Vec point{s.origin[0]+(id%s.nx)*s.spacing[0],
                  s.origin[1]+((id/s.nx)%s.ny)*s.spacing[1],
                  s.origin[2]+(id/s.nx/s.ny)*s.spacing[2]};
        double distance=1e100;
        for (const auto& shape:geometry) distance=std::min(distance,shape.distance(point));
        s.cells[id].wall_distance=distance;
        s.cells[id].wall_distance_known=true;
        if (distance<=0) {s.cells[id].solid=true; s.cells[id].obstacle=true;}
    }
    c.reference_length=units.si?units.reference_length:c.reference_length;
    if (c.reference_length==0)
        c.reference_length=geometry.empty()?(hi[1]-lo[1]+mesh.spacing[1]):
                          (geometry.front().radius>0?2*geometry.front().radius:geometry.front().hi[1]-geometry.front().lo[1]);
    c.reference_length/=mesh.spacing[0];
    bool inlet_seen = false, outlet_seen = false;
    double inlet_rho = 1, inlet_alpha = 1;
    double outlet_rho=1;
    const bool precision=c.single.flow.open_boundary==OpenBoundaryModel::VelocityPressure && !s.two_phase;
    for (const auto &patch : mesh.patches) {
        auto bc = Cursor(required(u.patches, patch.name)).dict();
        auto rb = Cursor(required(rho.patches, patch.name)).dict();
        Dict ab;
        if (s.two_phase)
            ab = Cursor(required(alpha.patches, patch.name)).dict();
        auto type = word(bc, "type"), rt = word(rb, "type");
        if (has_solid) {
            auto sb = Cursor(required(solid.patches, patch.name)).dict();
            need(word(sb, "type") == (patch.type == "empty" ? "empty" : "zeroGradient"),
                 "solid boundary must be empty or zeroGradient");
        }
        for (int f = patch.start; f < patch.start + patch.count; ++f) {
            auto center = mesh.centers[mesh.owner[f]], fc = face_center(mesh, f);
            auto &cell = s.cells[c.cell_order[mesh.owner[f]]];
            if (patch.type == "empty") {
                need(type == "empty" && rt == "empty" &&
                         (!s.two_phase || word(ab, "type") == "empty"),
                     "empty mesh/field mismatch");
                need(near(std::abs(fc[2] - center[2]) / mesh.spacing[2], 0.5),
                     "Only z-normal empty faces supported");
                continue;
            }
            need(!s.two_dimensional || same_position(fc[2], center[2], mesh.spacing[2]),
                 "2D front/back must be empty");
            // 禁止内部挡板：这种零厚度边界无法只用固体体素表示。
            Vec across = center;
            for (int d = 0; d < 3; ++d)
                across[d] = 2 * fc[d] - center[d];
            bool within = true;
            for (int d = 0; d < 3; ++d)
                within = within && (across[d] >= lo[d] - 1e-8 * mesh.spacing[d] &&
                                    across[d] <= hi[d] + 1e-8 * mesh.spacing[d]);
            if (within)
                need(s.cells[index(across)].solid,
                     "Internal baffles are not supported; use solid voxels");
            if (patch.type == "wall") {
                need(type == "noSlip" ||
                         (type == "fixedValue" && values(required(bc, "value"), 3, 1)[0] == Vec{}),
                     "Only stationary no-slip walls supported");
                need(rt == "zeroGradient" && (!s.two_phase || word(ab, "type") == "zeroGradient"),
                     "Wall scalar fields require zeroGradient; wetting uses lbmDict");
            } else if (type == "fixedValue") {
                need(!cell.solid, "Inlet cannot be marked solid");
                need(same_position(center[0], lo[0], mesh.spacing[0]) && fc[0] < center[0] &&
                         same_position(fc[1], center[1], mesh.spacing[1]) &&
                         same_position(fc[2], center[2], mesh.spacing[2]),
                     "Velocity inlet must be on left x face");
                auto v = values(required(bc, "value"), 3, 1)[0];
                need(v[1] == 0 && v[2] == 0 && v[0] >= 0 && v[0] < 0.2,
                     "Inlet must have uniform +x velocity <0.2");
                need(rt == (precision?"zeroGradient":"fixedValue"),
                     "Inlet rho: velocityPressure uses zeroGradient; legacy uses fixedValue");
                double r = precision?cell.rho:values(required(rb, "value"), 1, 1)[0][0];
                need(r > 0, "Invalid inlet density");
                double a = 1;
                if (s.two_phase) {
                    need(word(ab, "type") == "fixedValue", "Inlet alpha.A requires fixedValue");
                    a = values(required(ab, "value"), 1, 1)[0][0];
                    need(a > 0 && a < 1,
                         "Two-phase inlet alpha.A must lie strictly between 0 and 1");
                }
                if (inlet_seen)
                    need(near(c.inlet_velocity, v[0]) && near(inlet_rho, r) && near(inlet_alpha, a),
                         "All inlet patches must share uniform values");
                inlet_seen = true;
                c.inlet_velocity = v[0];
                inlet_rho = r;
                inlet_alpha = a;
                cell.boundary = 1;
                if (c.parabolic_inlet) {
                    const double eta=(center[1]-lo[1]+0.5*mesh.spacing[1])/
                                      (hi[1]-lo[1]+mesh.spacing[1]);
                    cell.inlet_scale=6*eta*(1-eta);
                }
            } else {
                need(!cell.solid, "Outlet cannot be marked solid");
                need(type == "zeroGradient" && rt == (precision?"fixedValue":"zeroGradient") &&
                         (!s.two_phase || word(ab, "type") == "zeroGradient"),
                     "Unsupported outlet field condition");
                need(same_position(center[0], hi[0], mesh.spacing[0]) && fc[0] > center[0] &&
                         same_position(fc[1], center[1], mesh.spacing[1]) &&
                         same_position(fc[2], center[2], mesh.spacing[2]),
                     "Outlet must be on right x face");
                cell.boundary = 2;
                if (precision) {
                    double r=values(required(rb,"value"),1,1)[0][0];
                    need(r>0 && (!outlet_seen || near(r,outlet_rho)),"Outlet densities must be positive and uniform");
                    outlet_rho=r;
                }
                outlet_seen = true;
            }
        }
    }
    need(!open || (inlet_seen && outlet_seen),
         "Open flow requires both left inlet and right outlet");
    for (std::size_t id : c.cell_order)
        if (s.cells[id].boundary) {
            auto adjacent = s.cells[id].boundary == 1 ? id + 1 : id - 1;
            need(adjacent < s.cells.size() && !s.cells[adjacent].solid,
                 "Inlet/outlet needs adjacent fluid cell");
        }
    c.single.initial_rho = c.single3d.initial_rho = precision?outlet_rho:inlet_rho;
    c.single.flow.outlet_density=c.single3d.flow.outlet_density=outlet_rho;
    c.two.inlet_velocity = c.two3d.inlet_velocity = c.inlet_velocity;
    c.two.rho_high = c.two3d.rho_high = inlet_rho * inlet_alpha;
    c.two.rho_low = c.two3d.rho_low = open ? inlet_rho * (1 - inlet_alpha) : 0.02;
    sample_pressure_probes(c.initial,c.pressure_probes);
    return c;
}

namespace {
struct FieldReport {
    double mass=0, rho_min=1e100, rho_max=0, speed=0, energy=0, change=0, inlet=0, outlet=0;
    double divergence_rms=0;
    double density_change=0;
};

FieldReport field_report(const FoamSnapshot& s, const FoamSnapshot* previous=nullptr) {
    FieldReport r;
    double delta=0, norm=0, divergence2=0, density_delta=0, density_norm=0;
    std::size_t fluid=0;
    for (std::size_t i=0; i<s.cells.size(); ++i) {
        const auto& cell=s.cells[i];
        if (cell.solid) continue;
        ++fluid;
        double divergence=0;
        for (int d=0; d<3; ++d) divergence+=detail::velocity_derivative(s,i,d,d);
        divergence2+=divergence*divergence;
        r.mass+=cell.rho;
        density_norm+=cell.rho*cell.rho;
        if (previous) {
            const double diff=cell.rho-previous->cells[i].rho;
            density_delta+=diff*diff;
        }
        r.rho_min=std::min(r.rho_min,cell.rho); r.rho_max=std::max(r.rho_max,cell.rho);
        double speed2=0;
        for (int d=0; d<3; ++d) {
            speed2+=cell.velocity[d]*cell.velocity[d];
            if (previous) {
                const double diff=cell.velocity[d]-previous->cells[i].velocity[d];
                delta+=diff*diff;
            }
        }
        norm+=speed2; r.speed=std::max(r.speed,std::sqrt(speed2));
        r.energy+=0.5*cell.rho*speed2;
        if (cell.boundary==1) r.inlet+=cell.rho*cell.velocity[0];
        if (cell.boundary==2) r.outlet+=cell.rho*cell.velocity[0];
    }
    r.change=previous?std::sqrt(delta/std::max(norm,1e-30)):0;
    r.density_change=previous?std::sqrt(density_delta/std::max(density_norm,1e-30)):0;
    r.divergence_rms=s.units.velocity_scale()*std::sqrt(divergence2/std::max<std::size_t>(fluid,1));
    return r;
}

template <class SolverType>
FoamRunSummary run_case(SolverType &solver, const FoamCase &c, const std::filesystem::path &output) {
    c.convergence.validate();
    need(c.steps>=0 && c.write_interval>0,"Invalid run steps/output interval");
    need(std::isfinite(c.inlet_ramp_steps) && c.inlet_ramp_steps>=0 &&
         std::isfinite(c.max_mach) && c.max_mach>0 && c.max_mach<0.5 &&
         std::isfinite(c.reference_length) && c.reference_length>=0 &&
         std::isfinite(c.inlet_velocity) && c.inlet_velocity>=0,"Invalid run controls");
    solver.initialize_fields(c.initial);
    auto snapshot = [&]() {
        auto s = solver.openfoam_snapshot();
        s.origin = c.initial.origin;
        s.spacing = c.initial.spacing;
        s.units = c.initial.units;
        return s;
    };
    constexpr bool single=std::is_same_v<SolverType,Solver> || std::is_same_v<SolverType,Solver3D>;
    need(single || !c.convergence.enabled,"Steady-state termination currently requires singlePhase");
    auto previous=snapshot();
    const auto initial_pressures=sample_pressure_probes(previous,c.pressure_probes);
    auto last=field_report(previous);
    const double reference_speed=c.inlet_velocity>0?c.inlet_velocity:last.speed;
    const double reference_rho=single?(previous.two_dimensional?c.single.initial_rho:c.single3d.initial_rho):1;
    const double area=previous.two_dimensional?c.reference_length:c.reference_length*c.reference_length;
    const std::size_t fluid=std::count_if(previous.cells.begin(),previous.cells.end(),
                                        [](const FoamCell& cell){return !cell.solid;});
    std::cout << "units=" << (previous.units.si?"SI":"lattice") << " cells=" << fluid
              << " steps=" << c.steps << " deltaT=" << previous.units.time_step;
    if constexpr (single) {
        const double tau=previous.two_dimensional?c.single.tau:c.single3d.tau;
        std::cout << " tau=" << tau << " reference Re=" <<
            reference_speed*c.reference_length/((tau-0.5)/3);
    }
    const int monitor_interval=c.convergence.interval>0?c.convergence.interval:c.write_interval;
    const bool check_force=reference_speed>0 && area>0 && std::any_of(previous.cells.begin(),previous.cells.end(),
        [](const FoamCell& cell){return cell.solid && cell.obstacle;});
    const double force_scale=check_force?0.5*reference_rho*reference_speed*reference_speed*area:1;
    auto settings=c.convergence;
    settings.minimum_steps=std::max(settings.minimum_steps,
        static_cast<int>(std::min(static_cast<double>(c.steps),std::ceil(c.inlet_ramp_steps))));
    ConvergenceMonitor convergence(settings);
    FoamRunSummary summary;
    summary.pressures=initial_pressures;
    auto print_probes=[&](int step, const FoamSnapshot& state) {
        summary.pressures=sample_pressure_probes(state,c.pressure_probes);
        if (summary.pressures.empty()) return;
        std::cout << "PROBES step=" << step << " time=" << step*state.units.time_step
                  << " pressureUnit=" << (state.units.si?"Pa":"lattice") << std::setprecision(12);
        for (std::size_t i=0; i<summary.pressures.size(); ++i)
            std::cout << " p[" << i << "]=" << summary.pressures[i];
        if (summary.pressures.size()==2)
            std::cout << " deltaP01=" << summary.pressures[0]-summary.pressures[1];
        std::cout << '\n';
    };
    summary.maximum_mach=last.speed*std::sqrt(3.0);
    int previous_step=0;
    double drag_min=0, drag_max=0, lift_min=0, lift_max=0;
    std::cout << "\nDiagnostics: lattice mass/rho/Mach; dU/dRho use monitorInterval=" << monitor_interval
              << "; forceRange is the full-window Cd/Cl range\n";
    if constexpr (single)
        need(last.speed*std::sqrt(3.0)<=c.max_mach,"Initial Mach exceeds maxMach");
    print_probes(0,previous);
    auto io_start=std::chrono::steady_clock::now();
    OpenFoamWriter writer(output, previous);
    writer.write(0, previous);
    double io_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-io_start).count();
    double compute_seconds=0, boundary_exchange=0, wall_exchange=0;
    for (int step = 1; step <= c.steps; ++step) {
        auto tick=std::chrono::steady_clock::now();
        if constexpr (std::is_same_v<SolverType, Solver> || std::is_same_v<SolverType, Solver3D>)
            solver.step_masked_flow(c.inlet_velocity*(c.inlet_ramp_steps>0?
                0.5*(1-std::cos(3.14159265358979323846*std::min(1.0,step/c.inlet_ramp_steps))):1));
        else
            solver.step();
        if constexpr (single) {
            // 每步扫描宏观场，避免写盘间隔之间的超限被漏检；不复制全场。
            const auto& grid=solver.grid();
            double speed2=0;
            for (std::size_t i=0; i<grid.rho.size(); ++i) if (!grid.solid[i]) {
                double value=grid.ux[i]*grid.ux[i]+grid.uy[i]*grid.uy[i];
                if constexpr (std::is_same_v<SolverType,Solver3D>) value+=grid.uz[i]*grid.uz[i];
                if (!std::isfinite(value)) throw std::runtime_error("Non-finite single-phase velocity");
                speed2=std::max(speed2,value);
            }
            const double mach=std::sqrt(3*speed2);
            summary.maximum_mach=std::max(summary.maximum_mach,mach);
            need(mach<=c.max_mach,"Computed Mach exceeds maxMach; reduce deltaT/velocity");
            const auto& d=solver.boundary_diagnostics();
            boundary_exchange+=d.mass_exchange;
            wall_exchange+=d.wall_mass_exchange;
            summary.drag=d.obstacle_force[0]/force_scale;
            summary.lift=d.obstacle_force[1]/force_scale;
            drag_min=std::min(drag_min,summary.drag); drag_max=std::max(drag_max,summary.drag);
            lift_min=std::min(lift_min,summary.lift); lift_max=std::max(lift_max,summary.lift);
        }
        compute_seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-tick).count();
        summary.steps=step;
        if (step % monitor_interval == 0 || step == c.steps) {
            auto current=snapshot();
            const auto report=field_report(current,&previous);
            const double mach=report.speed*std::sqrt(3.0);
            summary.velocity_change=report.change;
            summary.density_change=report.density_change;
            summary.force_range=std::max(drag_max-drag_min,lift_max-lift_min);
            summary.converged=step-previous_step==monitor_interval &&
                convergence.observe(previous_step,report.change,report.density_change,summary.force_range,check_force);
            std::cout << "step=" << step << " time=" << step*current.units.time_step
                      << " mass=" << std::setprecision(12) << report.mass
                      << " rho=[" << report.rho_min << ',' << report.rho_max << ']'
                      << " Ma=" << mach << " dU=" << report.change << " dRho=" << report.density_change
                      << " forceRange=" << summary.force_range << " steadyWindows=" << convergence.passes()
                      << " divRMS=" << report.divergence_rms << " energyLU=" << report.energy
                      << " flux(in,out)=" << report.inlet << ',' << report.outlet;
            if constexpr (single) {
                need(mach<=c.max_mach,"Computed Mach exceeds maxMach; reduce deltaT/velocity");
                const auto d=solver.diagnostics();
                const double balance=report.mass-last.mass-boundary_exchange-wall_exchange;
                std::cout << " massBalance=" << balance << " wallMass=" << wall_exchange
                          << " obstacleForce=(" << d.boundary.obstacle_force[0] << ','
                          << d.boundary.obstacle_force[1] << ',' << d.boundary.obstacle_force[2] << ')';
                if (reference_speed>0 && area>0) {
                    const double scale=0.5*reference_rho*reference_speed*reference_speed*area;
                    std::cout << " Cd_ref=" << d.boundary.obstacle_force[0]/scale
                              << " Cl_ref=" << d.boundary.obstacle_force[1]/scale;
                }
                std::cout << " interpolationFallbacks=" << d.boundary.interpolation_fallbacks;
                if (std::max(std::abs(report.rho_max-reference_rho),
                             std::abs(report.rho_min-reference_rho))/reference_rho>c.density_warning)
                    std::cout << " WARNING: density variation exceeds densityWarning";
            }
            std::cout << '\n';
            print_probes(step,current);
            previous=std::move(current); last=report;
            previous_step=step;
            drag_min=drag_max=summary.drag; lift_min=lift_max=summary.lift;
            boundary_exchange=wall_exchange=0;
        }
        if (step%c.write_interval==0 || step==c.steps || summary.converged) {
            tick=std::chrono::steady_clock::now();
            writer.write(step,previous_step==step?previous:snapshot());
            io_seconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-tick).count();
        }
        if (summary.converged) break;
    }
    std::cout << "computeSeconds=" << compute_seconds << " outputSeconds=" << io_seconds
              << " fluidMLUPS=" << (compute_seconds>0?fluid*static_cast<double>(summary.steps)/compute_seconds/1e6:0) << '\n';
    std::cout << "termination=" << (summary.converged?"converged":"endTime") << " maxMach=" << summary.maximum_mach << '\n';
    std::cout << "Completed " << summary.steps << " steps. OpenFOAM: " << (output / "lbm.foam") << '\n';
    return summary;
}
} // namespace

void solve_foam_case(const std::filesystem::path &path) {
    const auto c = read_foam_case(path);
    auto input = std::filesystem::canonical(path);
    while (input.has_relative_path() && input.filename().empty())
        input = input.parent_path();
    const auto output = openfoam_result_path(input.filename().string().c_str());
    need(std::filesystem::weakly_canonical(path) != std::filesystem::weakly_canonical(output),
         "Input case must not be its output directory");
    run_foam_case(c,input.filename().string().c_str());
}

FoamRunSummary run_foam_case(const FoamCase& c, const char* output_case) {
    const auto& s=c.initial;
    const auto output=openfoam_result_path(output_case);
    if (s.two_phase) {
        if (s.two_dimensional) {
            TwoPhaseSolver solver(s.nx, s.ny, c.two);
            return run_case(solver, c, output);
        } else {
            TwoPhaseSolver3D solver(s.nx, s.ny, s.nz, c.two3d);
            return run_case(solver, c, output);
        }
    } else {
        if (s.two_dimensional) {
            Solver solver(s.nx, s.ny, c.single);
            return run_case(solver, c, output);
        } else {
            Solver3D solver(s.nx, s.ny, s.nz, c.single3d);
            return run_case(solver, c, output);
        }
    }
}

void block_mesh_case(const std::filesystem::path &path) {
    auto d = dictionary(path / "system" / "blockMeshDict");
    known(d, {"convertToMeters", "scale", "vertices", "blocks", "edges", "boundary",
              "mergePatchPairs"});
    need(!(d.count("scale") && d.count("convertToMeters")),
         "Specify only scale or convertToMeters");
    double scale = number(d, "scale", number(d, "convertToMeters", 1));
    need(scale > 0, "Invalid mesh scale");
    Cursor v(required(d, "vertices"));
    v.expect("(");
    std::vector<Vec> vertices;
    while (!v.take(")"))
        vertices.push_back(v.vec());
    need(vertices.size() == 8 && v.i == v.t.size(),
         "Native blockMesh supports exactly eight vertices");
    for (auto &p : vertices)
        for (auto &x : p)
            x *= scale;
    Cursor b(required(d, "blocks"));
    b.expect("(");
    b.expect("hex");
    b.expect("(");
    std::array<int, 8> ids{};
    for (auto &id : ids) {
        id = b.label();
        need(id < 8, "Bad block vertex");
    }
    b.expect(")");
    b.expect("(");
    std::array<int, 3> n{b.label(), b.label(), b.label()};
    b.expect(")");
    b.expect("simpleGrading");
    auto grading = b.vec();
    need(grading == Vec{1, 1, 1}, "Native blockMesh requires simpleGrading (1 1 1)");
    b.expect(")");
    need(b.i == b.t.size(), "Native blockMesh supports one block; external blockMesh may generate "
                            "multiple uniform blocks");
    Vec lo = vertices[ids[0]], hi = vertices[ids[6]];
    const int xyz[8][3] = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
                           {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
    for (int i = 0; i < 8; ++i)
        for (int a = 0; a < 3; ++a)
            need(same_position(vertices[ids[i]][a], xyz[i][a] ? hi[a] : lo[a], hi[a] - lo[a]),
                 "Block must be axis-aligned with standard hex ordering");
    for (const auto *key : {"edges", "mergePatchPairs"})
        if (d.count(key))
            need(d.at(key) == Tokens{"(", ")"},
                 std::string("Nonempty ") + key + " is not supported");
    const int face_ids[6][4] = {{0, 4, 7, 3}, {1, 2, 6, 5}, {0, 1, 5, 4},
                                {3, 7, 6, 2}, {0, 3, 2, 1}, {4, 5, 6, 7}};
    std::array<std::string, 6> sides{};
    Cursor patches(required(d, "boundary"));
    patches.expect("(");
    while (!patches.take(")")) {
        auto name = patches.get();
        auto pd = Cursor(patches.group("{", "}")).dict();
        known(pd, {"type", "faces"});
        const auto type = word(pd, "type");
        need((name == "walls" && type == "wall") ||
                 ((name == "inlet" || name == "outlet") && type == "patch") ||
                 (name == "frontAndBack" && type == "empty"),
             "Native mesh patches: walls, inlet, outlet, frontAndBack only");
        Cursor faces(required(pd, "faces"));
        faces.expect("(");
        while (!faces.take(")")) {
            faces.expect("(");
            std::set<int> face;
            for (int k = 0; k < 4; ++k)
                face.insert(faces.label());
            faces.expect(")");
            int side = -1;
            for (int q = 0; q < 6; ++q) {
                std::set<int> expected;
                for (int k : face_ids[q])
                    expected.insert(ids[k]);
                if (face == expected)
                    side = q;
            }
            need(side >= 0 && sides[side].empty(), "Invalid/repeated block boundary face");
            sides[side] = name;
        }
        need(faces.i == faces.t.size(), "Extra patch faces tokens");
    }
    for (const auto &side : sides)
        need(!side.empty(), "All six block faces must be assigned");
    const bool two_d = sides[4] == "frontAndBack" && sides[5] == "frontAndBack";
    const bool open = sides[0] == "inlet" && sides[1] == "outlet";
    need(open || (sides[0] == "walls" && sides[1] == "walls"),
         "Use left inlet/right outlet or closed walls");
    need(sides[2] == "walls" && sides[3] == "walls" &&
             (two_d || (sides[4] == "walls" && sides[5] == "walls")),
         "Invalid block boundary orientation");
    need(!two_d || n[2] == 1, "2D requires one z cell");
    FoamSnapshot s;
    s.two_dimensional = two_d;
    std::array<int, 3> pad{open ? 0 : 1, 1, two_d ? 0 : 1};
    std::size_t count = 1;
    for (int a = 0; a < 3; ++a) {
        need(n[a] > 0 && hi[a] > lo[a], "Invalid block dimensions");
        s.spacing[a] = (hi[a] - lo[a]) / n[a];
        s.origin[a] = lo[a] + (0.5 - pad[a]) * s.spacing[a];
        n[a] += 2 * pad[a];
        need(count <= 2000000 / static_cast<std::size_t>(n[a]),
             "Native mesh exceeds 2 million cells");
        count *= n[a];
    }
    need(near(s.spacing[0] / s.spacing[1], 1) && (two_d || near(s.spacing[0] / s.spacing[2], 1)),
         "LBM requires equal active spacing");
    s.nx = n[0];
    s.ny = n[1];
    s.nz = n[2];
    s.cells.resize(count);
    for (int z = 0; z < s.nz; ++z)
        for (int y = 0; y < s.ny; ++y)
            for (int x = 0; x < s.nx; ++x) {
                auto &c = s.cells[(static_cast<std::size_t>(z) * s.ny + y) * s.nx + x];
                c.solid = x < pad[0] || x >= s.nx - pad[0] || y < pad[1] || y >= s.ny - pad[1] ||
                          z < pad[2] || z >= s.nz - pad[2];
                if (!c.solid && open)
                    c.boundary = x == 0 ? 1 : (x == s.nx - 1 ? 2 : 0);
            }
    // 临时网格独立于用户的 0/，网格重建不会覆盖手工修改的初始场。
    const auto work = path / ".lbm-mesh-work";
    OpenFoamWriter writer(work, s);
    const auto target = path / "constant" / "polyMesh";
    std::filesystem::create_directories(target);
    for (const auto *name : {"points", "faces", "owner", "neighbour", "boundary"})
        std::filesystem::copy_file(work / "constant" / "polyMesh" / name, target / name,
                                   std::filesystem::copy_options::overwrite_existing);
    std::cout << "Wrote mesh: " << target << '\n';
}

void set_fields_case(const std::filesystem::path &path) {
    const auto mesh = mesh_read(path);
    const bool physical=std::filesystem::exists(path/"constant"/"lbmProperties");
    auto d = dictionary(path / "system" / "setFieldsDict");
    known(d, {"defaultFieldValues", "regions"});
    std::map<std::string, Field> fields;
    std::map<std::string, int> components;
    auto assign = [&](const Tokens &tokens, const std::vector<std::size_t> &selected) {
        Cursor c(tokens);
        c.expect("(");
        while (!c.take(")")) {
            auto type = c.get(), name = c.get();
            need(type == "volScalarFieldValue" || type == "volVectorFieldValue",
                 "Unsupported setFields value type");
            need(name == "U" || name == "rho" || name == "p" || name == "alpha.A" ||
                     name == "porosity" || name == "solid",
                 "Unsupported setFields field: " + name);
            int comp = type == "volVectorFieldValue" ? 3 : 1;
            need((name == "U") == (comp == 3), "Wrong setFields field type");
            Vec value = comp == 3 ? c.vec() : Vec{c.scalar(), 0, 0};
            if (!fields.count(name)) {
                fields.emplace(name, read_field(path / "0" / name, comp, mesh.centers.size(),physical));
                components[name] = comp;
            }
            for (auto i : selected)
                fields.at(name).data[i] = value;
            c.take(";");
        }
        need(c.i == c.t.size(), "Extra setFields values");
    };
    std::vector<std::size_t> all;
    for (std::size_t i = 0; i < mesh.centers.size(); ++i)
        all.push_back(i);
    assign(required(d, "defaultFieldValues"), all);
    Cursor regions(required(d, "regions"));
    regions.expect("(");
    while (!regions.take(")")) {
        auto type = regions.get();
        auto region = Cursor(regions.group("{", "}")).dict();
        std::vector<std::size_t> selected;
        if (type == "boxToCell") {
            known(region, {"box", "fieldValues"});
            Cursor box(required(region, "box"));
            auto lo = box.vec(), hi = box.vec();
            need(box.i == box.t.size(), "Invalid region box");
            for (int a = 0; a < 3; ++a)
                need(lo[a] <= hi[a], "Reversed region box");
            for (std::size_t i : all) {
                auto p = mesh.centers[i];
                if (p[0] >= lo[0] && p[0] <= hi[0] && p[1] >= lo[1] && p[1] <= hi[1] &&
                    p[2] >= lo[2] && p[2] <= hi[2])
                    selected.push_back(i);
            }
        } else if (type == "sphereToCell") {
            known(region, {"centre", "radius", "fieldValues"});
            Cursor center(required(region, "centre"));
            auto p = center.vec();
            need(center.i == center.t.size(), "Invalid sphere centre");
            double r = number(region, "radius", -1);
            need(r > 0, "Invalid sphere radius");
            for (std::size_t i : all) {
                double dist = 0;
                for (int a = 0; a < 3; ++a)
                    dist += std::pow(mesh.centers[i][a] - p[a], 2);
                if (dist <= r * r)
                    selected.push_back(i);
            }
        } else
            throw std::runtime_error("Only boxToCell and sphereToCell regions supported");
        assign(required(region, "fieldValues"), selected);
    }
    need(regions.i == regions.t.size(), "Extra regions tokens");
    // 维持输入网格的 cell 编号和原边界字典，不采用求解器的重排编号。
    for (const auto &entry : fields) {
        const auto &name = entry.first;
        const auto &field = entry.second;
        const int comp = components[name];
        std::ofstream out;
        out.exceptions(std::ios::badbit | std::ios::failbit);
        out.open(path / "0" / name);
        out << std::setprecision(17);
        out << "FoamFile\n{ version 2.0; format ascii; class "
            << (comp == 3 ? "volVectorField" : "volScalarField") << "; location \"0\"; object "
            << name << "; }\n";
        out << "dimensions [";
        const auto& dims=required(field.dict,"dimensions");
        for (std::size_t k=1; k+1<dims.size(); ++k) out << (k==1?"":" ") << dims[k];
        out << "];\ninternalField nonuniform List<"
            << (comp == 3 ? "vector" : "scalar") << ">\n"
            << field.data.size() << "\n(\n";
        for (auto value : field.data) {
            if (comp == 3)
                out << '(' << value[0] << ' ' << value[1] << ' ' << value[2] << ')';
            else
                out << value[0];
            out << '\n';
        }
        out << ");\nboundaryField\n{\n";
        for (const auto &patch : field.patches) {
            out << '"' << patch.first << "\"\n{\n";
            for (const auto &token : patch.second)
                out << token << ' ';
            out << "\n}\n";
        }
        out << "}\n";
        out.close();
    }
    std::cout << "Updated initial fields in " << (path / "0") << '\n';
}
void refine_foam_case(const std::filesystem::path& source,
                      const std::filesystem::path& target, int factor) {
    need(factor>=2 && factor<=8,"Refinement factor must be 2..8");
    need(!std::filesystem::exists(target),"Refinement target already exists");
    need(std::filesystem::exists(source/"constant"/"lbmProperties"),"Refinement requires an SI case");
    auto mesh=dictionary(source/"system"/"blockMeshDict");
    Cursor block(required(mesh,"blocks"));
    block.expect("("); block.expect("hex"); block.group("(",")");
    block.expect("(");
    const auto count_begin=block.i;
    int nx=block.label(), ny=block.label(), nz=block.label();
    block.expect(")"); block.expect("simpleGrading"); const auto grading=block.vec(); block.expect(")");
    need(block.i==block.t.size() && grading==Vec{1,1,1},"Refinement requires one uniform hex block");
    auto config=dictionary(source/"system"/"lbmDict");
    need(word(config,"model","singlePhase")=="singlePhase","SI refinement currently requires singlePhase");
    const bool two_d=word(config,"lattice","D2Q9")=="D2Q9";
    need(nx>0 && ny>0 && nz>0 && (!two_d || nz==1),"Invalid refinement cell counts");
    need(nx<=2000000/factor && ny<=2000000/factor && nz<=2000000/factor,"Refined extent too large");
    const auto refined_cells=(static_cast<long long>(nx)*factor+2)*
                             (static_cast<long long>(ny)*factor+2)*
                             (two_d?1:static_cast<long long>(nz)*factor+2);
    need(refined_cells<=2000000,"Refined bounding grid exceeds 2000000 cells including wall padding");
    block.t[count_begin]=std::to_string(nx*factor);
    block.t[count_begin+1]=std::to_string(ny*factor);
    block.t[count_begin+2]=std::to_string(nz*(two_d?1:factor));
    mesh["blocks"]=block.t;
    auto control=dictionary(source/"system"/"controlDict");
    const double dt=number(control,"deltaT",1)/(factor*factor);
    const int interval=integer(number(control,"writeInterval",10)*factor*factor,"refined writeInterval",true);
    std::ostringstream dt_text; dt_text << std::setprecision(17) << dt;
    control["deltaT"]={dt_text.str()}; control["writeInterval"]={std::to_string(interval)};
    if (config.count("convergence")) {
        auto settings=Cursor(config.at("convergence")).dict();
        if (settings.count("interval")) settings["interval"]={std::to_string(integer(
            number(settings,"interval",0)*factor*factor,"refined convergence interval"))};
        Tokens tokens;
        for (const auto& entry:settings) {
            tokens.push_back(entry.first);
            tokens.insert(tokens.end(),entry.second.begin(),entry.second.end());
            tokens.push_back(";");
        }
        config["convergence"]=std::move(tokens);
    }
    // 初始均匀场可保持物理值不变；非均匀列表必须显式重采样，不能直接复制。
    std::vector<std::filesystem::path> fields;
    for (const auto& entry:std::filesystem::directory_iterator(source/"0")) {
        need(entry.is_regular_file(),"Unexpected initial-field directory");
        const auto name=entry.path().filename().string();
        need(name=="U" || name=="rho" || name=="p" || name=="porosity" || name=="solid",
             "Unsupported refinement field: "+name);
        auto f=read_field(entry.path(),name=="U"?3:1,1,true);
        need(required(f.dict,"internalField").front()=="uniform","Refinement needs uniform 0/ templates; regenerate setFields afterward");
        fields.push_back(entry.path());
    }
    auto write_dictionary=[&](const char* name,const Dict& dict) {
        auto file=target/"system"/name;
        std::ofstream out(file); out.exceptions(std::ios::badbit|std::ios::failbit);
        out << "FoamFile\n{version 2.0; format ascii; class dictionary; object " << name << ";}\n";
        for (const auto& entry:dict) {
            out << entry.first << ' ';
            if (entry.first=="convergence") out << "{ ";
            for (const auto& token:entry.second) out << token << ' ';
            out << (entry.first=="convergence"?"}\n":";\n");
        }
    };
    std::filesystem::create_directories(target/"0");
    std::filesystem::create_directories(target/"system");
    std::filesystem::create_directories(target/"constant");
    for (const auto& file:fields) std::filesystem::copy_file(file,target/"0"/file.filename());
    for (const auto* file:{"setFieldsDict"})
        if (std::filesystem::exists(source/"system"/file))
            std::filesystem::copy_file(source/"system"/file,target/"system"/file);
    for (const auto* file:{"lbmProperties","lbmGeometry"})
        if (std::filesystem::exists(source/"constant"/file))
            std::filesystem::copy_file(source/"constant"/file,target/"constant"/file);
    write_dictionary("blockMeshDict",mesh); write_dictionary("controlDict",control);
    write_dictionary("lbmDict",config);
    std::cout << "Refined SI template: " << target << "; dx/" << factor << ", deltaT/"
              << factor*factor << "; physical fields/endTime preserved\n";
}
} // namespace lbm
