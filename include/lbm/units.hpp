#pragma once

#include <cmath>
#include <initializer_list>
#include <stdexcept>

namespace lbm {
// 内部始终使用格子变量；仅在文件输入和输出时转换单相 SI 量纲。
struct UnitScale {
    bool si = false;
    double length_step = 1.0;
    double time_step = 1.0;
    double density_scale = 1.0;
    double reference_pressure = 0.0;
    double kinematic_viscosity = 0.0;
    double reference_length = 0.0;

    void validate() const {
        for (double v : {length_step, time_step, density_scale})
            if (!std::isfinite(v) || v <= 0)
                throw std::invalid_argument("Unit scales must be finite and positive.");
        if (!std::isfinite(reference_pressure) || !std::isfinite(kinematic_viscosity) ||
            !std::isfinite(reference_length) || reference_length < 0 ||
            (si && kinematic_viscosity <= 0))
            throw std::invalid_argument("Invalid physical properties.");
    }
    double velocity_scale() const { return si ? length_step/time_step : 1.0; }
    double pressure(double rho) const {
        return si ? reference_pressure + density_scale*velocity_scale()*velocity_scale()*(rho-1)/3
                  : rho/3;
    }
    double rho_from_pressure(double p) const {
        return si ? 1 + 3*(p-reference_pressure)/(density_scale*velocity_scale()*velocity_scale())
                  : 3*p;
    }
    bool operator==(const UnitScale& b) const {
        return si==b.si && length_step==b.length_step && time_step==b.time_step &&
            density_scale==b.density_scale && reference_pressure==b.reference_pressure &&
            kinematic_viscosity==b.kinematic_viscosity && reference_length==b.reference_length;
    }
    bool operator!=(const UnitScale& b) const { return !(*this==b); }
};
} // namespace lbm
