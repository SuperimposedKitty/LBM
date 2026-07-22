#include "lbm/geometry_mask.hpp"

#include <algorithm>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace lbm {
namespace {

GeometryCell parse_cell(char symbol, std::size_t line, std::size_t column) {
    switch (symbol) {
    case '.':
        return GeometryCell::Fluid;
    case '#':
        return GeometryCell::Solid;
    case 'I':
        return GeometryCell::Inlet;
    case 'O':
        return GeometryCell::Outlet;
    case 'P':
        return GeometryCell::Porous;
    default:
        std::ostringstream message;
        message << "Invalid geometry character '" << symbol << "' at line " << line
                << ", column " << column << ". Allowed characters are . # I O P.";
        throw std::runtime_error(message.str());
    }
}

} // namespace

GeometryMask GeometryMask::load(const std::filesystem::path& path) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("Failed to open geometry file: " + path.string());
    }

    std::vector<std::string> rows;
    std::vector<std::size_t> source_lines;
    std::string line;
    std::size_t line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        if (line.empty()) {
            continue;
        }
        rows.push_back(line);
        source_lines.push_back(line_number);
    }

    if (rows.empty()) {
        throw std::runtime_error("Geometry file contains no non-empty rows: " + path.string());
    }

    const std::size_t width = rows.front().size();
    if (width == 0) {
        throw std::runtime_error("Geometry rows must not be empty: " + path.string());
    }

    GeometryMask mask;
    mask.width_ = static_cast<int>(width);
    mask.height_ = static_cast<int>(rows.size());
    mask.cells_.resize(width * rows.size());

    bool has_fluid = false;
    for (std::size_t file_y = 0; file_y < rows.size(); ++file_y) {
        if (rows[file_y].size() != width) {
            std::ostringstream message;
            message << "Geometry row width mismatch at line " << source_lines[file_y]
                    << ": expected " << width << " characters, got "
                    << rows[file_y].size() << '.';
            throw std::runtime_error(message.str());
        }

        const int solver_y = mask.height_ - 1 - static_cast<int>(file_y);
        for (std::size_t x = 0; x < width; ++x) {
            const GeometryCell cell = parse_cell(rows[file_y][x], source_lines[file_y], x + 1);
            if (cell == GeometryCell::Inlet && x != 0) {
                std::ostringstream message;
                message << "Inlet I must be on the left boundary at line "
                        << source_lines[file_y] << ", column " << x + 1 << '.';
                throw std::runtime_error(message.str());
            }
            if (cell == GeometryCell::Outlet && x + 1 != width) {
                std::ostringstream message;
                message << "Outlet O must be on the right boundary at line "
                        << source_lines[file_y] << ", column " << x + 1 << '.';
                throw std::runtime_error(message.str());
            }

            mask.cells_[static_cast<std::size_t>(solver_y) * width + x] = cell;
            has_fluid = has_fluid || cell != GeometryCell::Solid;
        }
    }

    if (!has_fluid) {
        throw std::runtime_error("Geometry must contain at least one fluid cell: " + path.string());
    }
    return mask;
}

int GeometryMask::width() const {
    return width_;
}

int GeometryMask::height() const {
    return height_;
}

GeometryCell GeometryMask::cell_at(int x, int y) const {
    if (x < 0 || x >= width_ || y < 0 || y >= height_) {
        throw std::out_of_range("Geometry coordinates are outside the mask.");
    }
    return cells_[static_cast<std::size_t>(y) * width_ + x];
}

bool GeometryMask::contains(GeometryCell cell) const {
    return std::find(cells_.begin(), cells_.end(), cell) != cells_.end();
}

const std::vector<GeometryCell>& GeometryMask::cells() const {
    return cells_;
}

std::filesystem::path geometry_path(const std::string& filename) {
#ifdef LBM_GEOMETRY_DIR
    return std::filesystem::path(LBM_GEOMETRY_DIR) / filename;
#else
    return std::filesystem::path("geometry") / filename;
#endif
}

} // namespace lbm
