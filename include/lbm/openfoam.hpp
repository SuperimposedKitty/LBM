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
    // 格点 (0,0,0) 的坐标和网格步长；场值仍使用格子单位。
    std::array<double, 3> origin{0.0, 0.0, 0.0};
    std::array<double, 3> spacing{1.0, 1.0, 1.0};
};

// 文件导入与公开初始化接口共用校验，失败时不修改求解器状态。
void validate_initial_snapshot(const FoamSnapshot& snapshot);

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

// 返回固定案例目录；构造 writer 时清理上次由本导出器登记的时间步。
std::filesystem::path openfoam_result_path(const char* case_name);

} // namespace lbm
