#pragma once

#include <filesystem>
#include <string>
#include <vector>

namespace lbm {

enum class GeometryCell {
    Fluid,
    Solid,
    Inlet,
    Outlet,
    Porous
};

class GeometryMask {
public:
    static GeometryMask load(const std::filesystem::path& path);

    int width() const;
    int height() const;
    GeometryCell cell_at(int x, int y) const;
    bool contains(GeometryCell cell) const;
    const std::vector<GeometryCell>& cells() const;

private:
    int width_ = 0;
    int height_ = 0;
    std::vector<GeometryCell> cells_;
};

std::filesystem::path geometry_path(const std::string& filename);

} // namespace lbm
