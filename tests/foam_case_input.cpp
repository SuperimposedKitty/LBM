#include "lbm/foam_case.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>

namespace {
void check(bool ok, const char *why) {
    if (!ok)
        throw std::runtime_error(why);
}
std::string read(const std::filesystem::path &p) {
    std::ifstream f(p);
    check(bool(f), "Cannot read test file");
    return {std::istreambuf_iterator<char>(f), {}};
}
void write(const std::filesystem::path &p, const std::string &s) {
    std::ofstream f(p);
    f << s;
    check(bool(f), "Cannot write test file");
}
std::string replaced(std::string s, const std::string &from, const std::string &to) {
    auto pos = s.find(from);
    check(pos != std::string::npos, "Missing test replacement");
    s.replace(pos, from.size(), to);
    return s;
}
template <class F> void rejects(F run) {
    bool failed = false;
    try {
        run();
    } catch (const std::exception &) {
        failed = true;
    }
    check(failed, "Invalid case was accepted");
}
template <class S> void verify_initial(S &solver, const lbm::FoamCase &c) {
    solver.initialize_fields(c.initial);
    const auto s = solver.openfoam_snapshot();
    for (auto i : c.cell_order) {
        check(std::abs(s.cells[i].rho - c.initial.cells[i].rho) < 1e-12, "Initial density ignored");
        for (int d = 0; d < 3; ++d)
            check(std::abs(s.cells[i].velocity[d] - c.initial.cells[i].velocity[d]) < 1e-12,
                  "Initial velocity ignored");
        if (s.two_phase)
            check(std::abs(s.cells[i].rho_a - c.initial.cells[i].rho_a) < 1e-12,
                  "Initial phase ignored");
    }
}

// 反转所有 cell 标签，重定向内部面，使输入顺序与格点顺序彻底不同。
void reverse_cells(const std::filesystem::path &path, int cells) {
    const auto base = path / "constant" / "polyMesh";
    auto list = [&](const char *name) {
        auto text = read(base / name);
        std::istringstream in(text.substr(text.find('}') + 1));
        int n;
        char ch;
        in >> n >> ch;
        std::vector<int> v(n);
        for (auto &x : v)
            in >> x;
        return v;
    };
    auto owner = list("owner"), neighbour = list("neighbour");
    auto text = read(base / "faces");
    auto header = text.substr(0, text.find('}') + 1);
    std::istringstream in(text.substr(text.find('}') + 1));
    int n;
    char ch;
    in >> n >> ch;
    std::vector<std::array<int, 4>> faces(n);
    for (auto &f : faces) {
        int size;
        in >> size >> ch;
        check(size == 4, "Expected quad");
        for (auto &v : f)
            in >> v;
        in >> ch;
    }
    for (std::size_t i = 0; i < owner.size(); ++i) {
        owner[i] = cells - 1 - owner[i];
        if (i < neighbour.size()) {
            neighbour[i] = cells - 1 - neighbour[i];
            std::swap(owner[i], neighbour[i]);
            std::reverse(faces[i].begin(), faces[i].end());
        }
    }
    std::ostringstream out;
    out << header << '\n' << n << "\n(\n";
    for (auto f : faces)
        out << "4(" << f[0] << ' ' << f[1] << ' ' << f[2] << ' ' << f[3] << ")\n";
    out << ")\n";
    write(base / "faces", out.str());
    auto save = [&](const char *name, const std::vector<int> &v) {
        auto raw = read(base / name);
        std::ostringstream s;
        s << raw.substr(0, raw.find('}') + 1) << '\n' << v.size() << "\n(\n";
        for (auto x : v)
            s << x << '\n';
        s << ")\n";
        write(base / name, s.str());
    };
    save("owner", owner);
    save("neighbour", neighbour);
}
} // namespace

int main() {
    try {
        auto root = std::filesystem::current_path() / "foam_case_input_tests";
        int run = 0;
        while (std::filesystem::exists(root / std::to_string(run)))
            ++run;
        root /= std::to_string(run);
        std::filesystem::create_directories(root);
        for (const auto *name : {"channel2d", "droplet2d", "displacement3d", "closed3d"}) {
            auto source = std::filesystem::path(LBM_CASES_DIR) / name;
            auto path = root / (std::string("input_test_") + name);
            std::filesystem::create_directories(path);
            // 只复制版本控制中的输入，忽略用户可能生成的网格。
            for (const auto *dir : {"0", "system"})
                std::filesystem::copy(source / dir, path / dir,
                                      std::filesystem::copy_options::recursive);
            const auto original_u = read(path / "0" / "U");
            if (std::string(name) == "closed3d") {
                auto mesh_dict = read(path / "system" / "blockMeshDict");
                write(path / "system" / "blockMeshDict",
                      replaced(mesh_dict, "convertToMeters 1;", "convertToMeters 0.000001;"));
            }
            lbm::block_mesh_case(path);
            check(read(path / "0" / "U") == original_u, "blockMesh overwrote initial field");
            // 手工修改 0/U 后无需重编译，初始化必须保留非零的两个速度分量。
            if (std::string(name) == "channel2d") {
                write(path / "0" / "U", replaced(original_u, "internalField uniform (0 0 0)",
                                                 "internalField uniform (0.012 0.004 0)"));
                auto edited = lbm::read_foam_case(path);
                check(edited.initial.cells[edited.cell_order[0]].velocity[0] == 0.012,
                      "Edited U ignored");
                lbm::Solver solver(edited.initial.nx, edited.initial.ny, edited.single);
                verify_initial(solver, edited);
                write(path / "0" / "U", original_u);
            }
            if (std::string(name) == "channel2d") {
                reverse_cells(path, 24 * 12);
                auto assignment = read(path / "system" / "setFieldsDict");
                write(path / "system" / "setFieldsDict",
                      replaced(assignment, "regions\n(\n);",
                               "regions\n( boxToCell { box (0 0 0) (4 12 1); "
                               "fieldValues (volVectorFieldValue U (0.01 0 0)); } );"));
            }
            lbm::set_fields_case(path);
            auto c = lbm::read_foam_case(path);
            const auto &s = c.initial;
            if (std::string(name) == "closed3d")
                check(std::abs(s.spacing[0] - 1e-6) < 1e-15, "Scaled coordinates lost");
            check(!c.cell_order.empty(), "No imported cells");
            if (std::string(name) == "channel2d") {
                for (auto id : c.cell_order) {
                    const auto x = id % static_cast<std::size_t>(s.nx);
                    check(std::abs(s.cells[id].velocity[0] - (x < 4 ? 0.01 : 0.0)) < 1e-12,
                          "Nonuniform values do not follow permuted cell labels");
                }
            }
            // p 可单独作为密度输入，但与 rho 同时存在时必须一致。
            if (std::string(name) == "channel2d") {
                const auto density = read(path / "0" / "rho");
                auto pressure = read(source / "0" / "rho");
                pressure = replaced(pressure, "object rho;", "object p;");
                pressure = replaced(pressure, "internalField uniform 1;",
                                    "internalField uniform 0.333333333333333333;");
                pressure =
                    replaced(pressure, "value uniform 1;", "value uniform 0.333333333333333333;");
                write(path / "0" / "p", pressure);
                std::filesystem::rename(path / "0" / "rho", path / "0" / "rho.saved");
                auto from_p = lbm::read_foam_case(path);
                check(std::abs(from_p.initial.cells[from_p.cell_order[0]].rho - 1) < 1e-12,
                      "rho=3p conversion failed");
                std::filesystem::rename(path / "0" / "rho.saved", path / "0" / "rho");
                write(path / "0" / "p",
                      replaced(pressure, "internalField uniform 0.333333333333333333;",
                               "internalField uniform 0.2;"));
                rejects([&] { lbm::read_foam_case(path); });
                std::filesystem::remove(path / "0" / "p");
            }
            if (std::string(name) == "channel2d")
                check(c.cell_order.front() > c.cell_order.back(), "Cell permutation not preserved");
            if (std::string(name) == "channel2d") {
                auto obstacle_path = root / "input_test_obstacle";
                std::filesystem::create_directories(obstacle_path);
                for (const auto *dir : {"0", "system"})
                    std::filesystem::copy(source / dir, obstacle_path / dir,
                                          std::filesystem::copy_options::recursive);
                auto solid = s;
                const auto blocked = static_cast<std::size_t>(s.nx * 5 + 10);
                solid.cells[blocked].solid = true;
                lbm::OpenFoamWriter exporter(root / "obstacle_mesh", solid);
                auto mesh_target = obstacle_path / "constant" / "polyMesh";
                std::filesystem::create_directories(mesh_target);
                for (const auto *file : {"points", "faces", "owner", "neighbour", "boundary"})
                    std::filesystem::copy_file(root / "obstacle_mesh" / "constant" / "polyMesh" /
                                                   file,
                                               mesh_target / file);
                auto imported = lbm::read_foam_case(obstacle_path);
                check(imported.cell_order.size() + 1 == c.cell_order.size() &&
                          imported.initial.cells[blocked].solid,
                      "Imported solid obstacle lost");
                lbm::Solver solver(imported.initial.nx, imported.initial.ny, imported.single);
                solver.initialize_fields(imported.initial);
                solver.step_masked_flow(imported.inlet_velocity);
                check(std::isfinite(solver.diagnostics().max_speed), "Obstacle step failed");
            }
            if (s.two_phase) {
                double min_a = 1, max_a = 0;
                for (auto i : c.cell_order) {
                    auto a = s.cells[i].rho_a / s.cells[i].rho;
                    min_a = std::min(min_a, a);
                    max_a = std::max(max_a, a);
                }
                if (std::string(name) == "droplet2d")
                    check(max_a - min_a > 0.9, "Sphere setFields failed");
            }
            if (std::string(name) == "displacement3d") {
                int porous = 0;
                for (auto i : c.cell_order)
                    if (s.cells[i].porosity < 0.5)
                        ++porous;
                check(porous == 4 * 12 * 8, "Porous box mapping failed");
            }
            if (s.two_phase) {
                if (s.two_dimensional) {
                    lbm::TwoPhaseSolver solver(s.nx, s.ny, c.two);
                    verify_initial(solver, c);
                } else {
                    lbm::TwoPhaseSolver3D solver(s.nx, s.ny, s.nz, c.two3d);
                    verify_initial(solver, c);
                }
            } else {
                if (s.two_dimensional) {
                    lbm::Solver solver(s.nx, s.ny, c.single);
                    verify_initial(solver, c);
                } else {
                    lbm::Solver3D solver(s.nx, s.ny, s.nz, c.single3d);
                    verify_initial(solver, c);
                }
            }
            const auto u = read(path / "0" / "U"), rho = read(path / "0" / "rho"),
                       points = read(path / "constant" / "polyMesh" / "points");
            write(path / "system" / "controlDict",
                  replaced(read(path / "system" / "controlDict"), "endTime 200;", "endTime 4;"));
            lbm::solve_foam_case(path);
            check(read(path / "0" / "U") == u && read(path / "0" / "rho") == rho &&
                      read(path / "constant" / "polyMesh" / "points") == points,
                  "Solver modified input case");
            const auto output = lbm::openfoam_result_path(path.filename().string().c_str());
            check(std::filesystem::exists(output / "0" / "U") &&
                      std::filesystem::exists(output / "4" / "U"),
                  "Initial/final output missing");
            write(path / "0" / "rho",
                  replaced(rho, "dimensions [0 0 0 0 0 0 0]", "dimensions [1 -3 0 0 0 0 0]"));
            rejects([&] { lbm::read_foam_case(path); });
            write(path / "0" / "rho", rho);
            auto mesh_dict = read(path / "system" / "blockMeshDict");
            write(path / "system" / "blockMeshDict",
                  replaced(mesh_dict, "simpleGrading (1 1 1)", "simpleGrading (2 1 1)"));
            rejects([&] { lbm::block_mesh_case(path); });
            auto u_bad = replaced(u, "fixedValue", "slip");
            write(path / "0" / "U", u_bad);
            rejects([&] { lbm::read_foam_case(path); });
            write(path / "0" / "U", u);
            std::cout << "PASS case " << name << '\n';
        }
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
