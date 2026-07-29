#include "lbm/geometry_mask3d.hpp"

#include <algorithm>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

namespace lbm {
namespace {

GeometryCell parse_cell3d(char symbol, int z, int row, int column) {
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
        message << "Invalid 3D geometry character '" << symbol << "' in slice " << z
                << ", row " << row << ", column " << column
                << ". Allowed characters are . # I O P.";
        throw std::runtime_error(message.str());
    }
}

std::string read_nonempty_line(std::ifstream& input, std::size_t& line_number) {
    std::string line;
    while (std::getline(input, line)) {
        ++line_number;
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        if (!line.empty()) {
            return line;
        }
    }
    return {};
}

} // namespace

GeometryMask3D GeometryMask3D::load(const std::filesystem::path& path) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("Failed to open 3D geometry file: " + path.string());
    }

    std::size_t line_number = 0;
    const std::string header = read_nonempty_line(input, line_number);
    std::istringstream header_stream(header);
    std::string magic;
    GeometryMask3D mask;
    if (!(header_stream >> magic >> mask.width_ >> mask.height_ >> mask.depth_) ||
        magic != "LBM_GEOMETRY_3D") {
        throw std::runtime_error(
            "3D geometry header must be: LBM_GEOMETRY_3D <nx> <ny> <nz>.");
    }
    std::string trailing;
    if (header_stream >> trailing) {
        throw std::runtime_error("Unexpected data after 3D geometry dimensions.");
    }
    if (mask.width_ <= 1 || mask.height_ <= 1 || mask.depth_ <= 1) {
        throw std::runtime_error("3D geometry dimensions must all be greater than 1.");
    }

    const std::size_t width = static_cast<std::size_t>(mask.width_);
    const std::size_t height = static_cast<std::size_t>(mask.height_);
    const std::size_t depth = static_cast<std::size_t>(mask.depth_);
    if (width > std::numeric_limits<std::size_t>::max() / height ||
        width * height > std::numeric_limits<std::size_t>::max() / depth) {
        throw std::runtime_error("3D geometry dimensions overflow addressable memory.");
    }
    const std::size_t cell_count = width * height * depth;
    mask.cells_.resize(cell_count);
    bool has_fluid = false;

    for (int file_z = 0; file_z < mask.depth_; ++file_z) {
        const std::string slice_header = read_nonempty_line(input, line_number);
        std::istringstream slice_stream(slice_header);
        std::string slice_magic;
        int slice_index = -1;
        if (!(slice_stream >> slice_magic >> slice_index) || slice_magic != "SLICE" ||
            slice_index != file_z) {
            std::ostringstream message;
            message << "Expected 'SLICE " << file_z << "' near source line " << line_number << '.';
            throw std::runtime_error(message.str());
        }
        if (slice_stream >> trailing) {
            throw std::runtime_error("Unexpected data after a 3D SLICE index.");
        }

        for (int file_y = 0; file_y < mask.height_; ++file_y) {
            const std::string row = read_nonempty_line(input, line_number);
            if (static_cast<int>(row.size()) != mask.width_) {
                std::ostringstream message;
                message << "3D geometry row width mismatch at line " << line_number
                        << ": expected " << mask.width_ << " characters, got " << row.size()
                        << '.';
                throw std::runtime_error(message.str());
            }

            const int solver_y = mask.height_ - 1 - file_y;
            for (int x = 0; x < mask.width_; ++x) {
                const GeometryCell cell =
                    parse_cell3d(row[static_cast<std::size_t>(x)], file_z, file_y + 1, x + 1);
                if (cell == GeometryCell::Inlet && x != 0) {
                    throw std::runtime_error("3D inlet I must be on the x-min boundary.");
                }
                if (cell == GeometryCell::Outlet && x != mask.width_ - 1) {
                    throw std::runtime_error("3D outlet O must be on the x-max boundary.");
                }
                const std::size_t index =
                    (static_cast<std::size_t>(file_z) * mask.height_ + solver_y) * mask.width_ + x;
                mask.cells_[index] = cell;
                has_fluid = has_fluid || cell != GeometryCell::Solid;
            }
        }
    }

    if (!read_nonempty_line(input, line_number).empty()) {
        throw std::runtime_error("Unexpected data after the final 3D geometry slice.");
    }
    if (!has_fluid) {
        throw std::runtime_error("3D geometry must contain at least one fluid cell.");
    }
    return mask;
}

int GeometryMask3D::width() const {
    return width_;
}

int GeometryMask3D::height() const {
    return height_;
}

int GeometryMask3D::depth() const {
    return depth_;
}

GeometryCell GeometryMask3D::cell_at(int x, int y, int z) const {
    if (x < 0 || x >= width_ || y < 0 || y >= height_ || z < 0 || z >= depth_) {
        throw std::out_of_range("Coordinates are outside the 3D geometry mask.");
    }
    return cells_[(static_cast<std::size_t>(z) * height_ + y) * width_ + x];
}

bool GeometryMask3D::contains(GeometryCell cell) const {
    return std::find(cells_.begin(), cells_.end(), cell) != cells_.end();
}

const std::vector<GeometryCell>& GeometryMask3D::cells() const {
    return cells_;
}

} // namespace lbm
