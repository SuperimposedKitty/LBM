#pragma once

#include <array>
#include <filesystem>
#include <vector>

namespace lbm {

// 全部数值使用格子单位；单元中心仍位于整数格点坐标。
struct FoamCell {
    bool solid = false;
    int boundary = 0; // 0: 普通流体，1: 入口，2: 出口
    std::array<double, 3> velocity{};
    double rho = 1.0;
    double rho_a = 0.0;
    double rho_b = 0.0;
    double porosity = 1.0;
};

struct FoamSnapshot {
    int nx = 0, ny = 0, nz = 1;
    bool two_dimensional = true;
    bool two_phase = false;
    double interaction_strength = 0.0;
    std::vector<FoamCell> cells;
};

// 原生 OpenFOAM 后处理案例，不包含可用于续算的有限体积求解器设置。
class OpenFoamWriter {
public:
    OpenFoamWriter(const std::filesystem::path& directory, const FoamSnapshot& initial);
    void write(int step, const FoamSnapshot& snapshot) const;
    const std::filesystem::path& directory() const { return directory_; }

private:
    struct Patch {
        const char* name;
        const char* type;
        std::vector<int> owners;
    };
    std::filesystem::path directory_;
    FoamSnapshot topology_;
    std::vector<std::size_t> fluid_;
    std::vector<Patch> patches_;
    void validate(const FoamSnapshot& snapshot) const;
};

// 每次运行创建独立目录，避免旧时间步混入新计算，也不删除既有结果。
std::filesystem::path openfoam_result_path(const char* case_name);

} // namespace lbm
