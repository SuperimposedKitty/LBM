#pragma once

#include "lbm/geometry_mask.hpp"

#include <filesystem>
#include <vector>

namespace lbm {

class GeometryMask3D {
public:
    static GeometryMask3D load(const std::filesystem::path& path);

    int width() const;
    int height() const;
    int depth() const;
    GeometryCell cell_at(int x, int y, int z) const;
    bool contains(GeometryCell cell) const;
    const std::vector<GeometryCell>& cells() const;

private:
    int width_ = 0;
    int height_ = 0;
    int depth_ = 0;
    std::vector<GeometryCell> cells_;
};

} // namespace lbm
