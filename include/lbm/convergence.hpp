#pragma once

#include <cmath>
#include <initializer_list>
#include <stdexcept>

namespace lbm {
struct ConvergenceConfig {
    bool enabled = false;
    int interval = 0; // 0 表示沿用输出间隔，显式设置后不受写盘频率影响。
    int minimum_steps = 0;
    int consecutive_windows = 5;
    double velocity_tolerance = 1e-5;
    double density_tolerance = 1e-6;
    double force_tolerance = 1e-3;

    void validate() const {
        if (interval<0 || minimum_steps<0 || consecutive_windows<2)
            throw std::invalid_argument("Invalid convergence intervals/windows.");
        for (double v:{velocity_tolerance,density_tolerance,force_tolerance})
            if (!std::isfinite(v) || v<=0)
                throw std::invalid_argument("Convergence tolerances must be finite and positive.");
    }
};

class ConvergenceMonitor {
public:
    explicit ConvergenceMonitor(ConvergenceConfig config):config_(config) { config_.validate(); }
    bool observe(int window_start, double velocity_change, double density_change,
                 double force_range, bool check_force) {
        // 窗口必须整体越过最早检查时间，启动阶段和失败样本都会清空连续计数。
        const bool pass=config_.enabled && window_start>=config_.minimum_steps &&
            std::isfinite(velocity_change) && velocity_change>=0 && velocity_change<=config_.velocity_tolerance &&
            std::isfinite(density_change) && density_change>=0 && density_change<=config_.density_tolerance &&
            (!check_force || (std::isfinite(force_range) && force_range>=0 && force_range<=config_.force_tolerance));
        passes_=pass?passes_+1:0;
        return passes_>=config_.consecutive_windows;
    }
    int passes() const { return passes_; }
private:
    ConvergenceConfig config_;
    int passes_=0;
};
} // namespace lbm
