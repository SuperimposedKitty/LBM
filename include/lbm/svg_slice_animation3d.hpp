#pragma once

#include "lbm/svg_animation.hpp"

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace lbm {

struct SliceFrame3D {
    std::vector<double> xy;
    std::vector<double> xz;
    std::vector<double> yz;
};

struct SliceMasks3D {
    std::vector<std::uint8_t> xy_solid;
    std::vector<std::uint8_t> xz_solid;
    std::vector<std::uint8_t> yz_solid;
    std::vector<double> xy_porosity;
    std::vector<double> xz_porosity;
    std::vector<double> yz_porosity;
};

template <typename ValueAt>
SliceFrame3D capture_orthogonal_slices(
    int nx,
    int ny,
    int nz,
    int slice_x,
    int slice_y,
    int slice_z,
    ValueAt value_at) {
    SliceFrame3D frame;
    frame.xy.resize(static_cast<std::size_t>(nx) * ny);
    frame.xz.resize(static_cast<std::size_t>(nx) * nz);
    frame.yz.resize(static_cast<std::size_t>(ny) * nz);
    for (int y = 0; y < ny; ++y) {
        for (int x = 0; x < nx; ++x) {
            frame.xy[static_cast<std::size_t>(y) * nx + x] =
                value_at(x, y, slice_z);
        }
    }
    for (int z = 0; z < nz; ++z) {
        for (int x = 0; x < nx; ++x) {
            frame.xz[static_cast<std::size_t>(z) * nx + x] =
                value_at(x, slice_y, z);
        }
    }
    for (int y = 0; y < ny; ++y) {
        for (int z = 0; z < nz; ++z) {
            frame.yz[static_cast<std::size_t>(y) * nz + z] =
                value_at(slice_x, y, z);
        }
    }
    return frame;
}

template <typename SolidAt, typename PorosityAt>
SliceMasks3D capture_orthogonal_masks(
    int nx,
    int ny,
    int nz,
    int slice_x,
    int slice_y,
    int slice_z,
    SolidAt solid_at,
    PorosityAt porosity_at) {
    SliceMasks3D masks;
    masks.xy_solid.resize(static_cast<std::size_t>(nx) * ny);
    masks.xz_solid.resize(static_cast<std::size_t>(nx) * nz);
    masks.yz_solid.resize(static_cast<std::size_t>(ny) * nz);
    masks.xy_porosity.resize(masks.xy_solid.size());
    masks.xz_porosity.resize(masks.xz_solid.size());
    masks.yz_porosity.resize(masks.yz_solid.size());
    for (int y = 0; y < ny; ++y) {
        for (int x = 0; x < nx; ++x) {
            const std::size_t s = static_cast<std::size_t>(y) * nx + x;
            masks.xy_solid[s] = solid_at(x, y, slice_z) ? 1 : 0;
            masks.xy_porosity[s] = porosity_at(x, y, slice_z);
        }
    }
    for (int z = 0; z < nz; ++z) {
        for (int x = 0; x < nx; ++x) {
            const std::size_t s = static_cast<std::size_t>(z) * nx + x;
            masks.xz_solid[s] = solid_at(x, slice_y, z) ? 1 : 0;
            masks.xz_porosity[s] = porosity_at(x, slice_y, z);
        }
    }
    for (int y = 0; y < ny; ++y) {
        for (int z = 0; z < nz; ++z) {
            const std::size_t s = static_cast<std::size_t>(y) * nz + z;
            masks.yz_solid[s] = solid_at(slice_x, y, z) ? 1 : 0;
            masks.yz_porosity[s] = porosity_at(slice_x, y, z);
        }
    }
    return masks;
}

inline void write_slice_cells3d(
    std::ofstream& out,
    const std::vector<double>& values,
    const std::vector<std::uint8_t>& solid,
    const std::vector<double>& porosity,
    int width_cells,
    int height_cells,
    int origin_x,
    int origin_y,
    int cell,
    double vmin,
    double span,
    ColorMap color_map) {
    for (int y = 0; y < height_cells; ++y) {
        for (int x = 0; x < width_cells; ++x) {
            const std::size_t s = static_cast<std::size_t>(y) * width_cells + x;
            const int px = origin_x + x * cell;
            const int py = origin_y + (height_cells - 1 - y) * cell;
            const std::string fill =
                !solid.empty() && solid[s]
                    ? "#2f343b"
                    : scalar_color(color_map, (values[s] - vmin) / span);
            out << "    <rect x=\"" << px << "\" y=\"" << py
                << "\" width=\"" << cell << "\" height=\"" << cell
                << "\" fill=\"" << fill << "\"/>\n";
            if (!porosity.empty() && porosity[s] > 0.0 && porosity[s] < 0.999 &&
                (solid.empty() || !solid[s])) {
                out << "    <rect x=\"" << px << "\" y=\"" << py
                    << "\" width=\"" << cell << "\" height=\"" << cell
                    << "\" fill=\"url(#porousHatch3d)\" fill-opacity=\"0.55\"/>\n";
            }
        }
    }
}

inline void write_orthogonal_slice_animation_svg(
    const std::filesystem::path& path,
    const std::vector<SliceFrame3D>& frames,
    const SliceMasks3D& masks,
    int nx,
    int ny,
    int nz,
    const ScalarAnimationOptions& options) {
    if (nx <= 0 || ny <= 0 || nz <= 0 || frames.empty()) {
        throw std::invalid_argument(
            "Three-dimensional slice animation requires positive dimensions and frames.");
    }
    const std::size_t xy_size = static_cast<std::size_t>(nx) * ny;
    const std::size_t xz_size = static_cast<std::size_t>(nx) * nz;
    const std::size_t yz_size = static_cast<std::size_t>(ny) * nz;
    for (const SliceFrame3D& frame : frames) {
        if (frame.xy.size() != xy_size || frame.xz.size() != xz_size ||
            frame.yz.size() != yz_size) {
            throw std::invalid_argument("A 3D animation slice has an invalid size.");
        }
    }
    if (masks.xy_solid.size() != xy_size || masks.xz_solid.size() != xz_size ||
        masks.yz_solid.size() != yz_size) {
        throw std::invalid_argument("Three-dimensional solid slice masks have invalid sizes.");
    }

    double vmin = options.vmin;
    double vmax = options.vmax;
    if (!options.fixed_range) {
        vmin = std::numeric_limits<double>::max();
        vmax = std::numeric_limits<double>::lowest();
        for (const SliceFrame3D& frame : frames) {
            for (double value : frame.xy) {
                vmin = std::min(vmin, value);
                vmax = std::max(vmax, value);
            }
            for (double value : frame.xz) {
                vmin = std::min(vmin, value);
                vmax = std::max(vmax, value);
            }
            for (double value : frame.yz) {
                vmin = std::min(vmin, value);
                vmax = std::max(vmax, value);
            }
        }
    }
    const double span = vmax > vmin ? vmax - vmin : 1.0;
    const int cell = std::max(1, options.cell);
    const int margin = 24;
    const int gap = 28;
    const int title_height = 46;
    const int label_height = 22;
    const int xy_width = nx * cell;
    const int xz_width = nx * cell;
    const int yz_width = nz * cell;
    const int plot_height = std::max(ny, nz) * cell;
    const int width = margin * 2 + xy_width + xz_width + yz_width + gap * 2;
    const int height = margin + title_height + label_height + plot_height + 40;
    const double fps = options.fps > 0.0 ? options.fps : 8.0;
    const double duration = static_cast<double>(frames.size()) / fps;
    const int xy_x = margin;
    const int xz_x = xy_x + xy_width + gap;
    const int yz_x = xz_x + xz_width + gap;
    const int plot_y = margin + title_height + label_height;

    std::filesystem::create_directories(path.parent_path());
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Failed to open 3D animation output: " + path.string());
    }
    out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    out << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << width
        << "\" height=\"" << height << "\" viewBox=\"0 0 " << width << ' '
        << height << "\" role=\"img\" aria-label=\"" << options.title << "\">\n";
    out << "  <rect width=\"100%\" height=\"100%\" fill=\"#f7f8fa\"/>\n";
    out << "  <defs><pattern id=\"porousHatch3d\" width=\"6\" height=\"6\""
           " patternUnits=\"userSpaceOnUse\"><path d=\"M0 6 L6 0\""
           " stroke=\"#1f2328\" stroke-opacity=\"0.42\" stroke-width=\"1\"/>"
           "</pattern></defs>\n";
    out << "  <text x=\"" << margin << "\" y=\"" << margin + 18
        << "\" font-family=\"Arial, sans-serif\" font-size=\"18\""
           " font-weight=\"700\" fill=\"#1f2328\">"
        << options.title << "</text>\n";
    out << "  <g font-family=\"Arial, sans-serif\" font-size=\"13\""
           " font-weight=\"600\" fill=\"#343a40\">\n";
    out << "    <text x=\"" << xy_x << "\" y=\"" << plot_y - 8
        << "\">XY center slice</text>\n";
    out << "    <text x=\"" << xz_x << "\" y=\"" << plot_y - 8
        << "\">XZ center slice</text>\n";
    out << "    <text x=\"" << yz_x << "\" y=\"" << plot_y - 8
        << "\">YZ center slice</text>\n";
    out << "  </g>\n";

    for (std::size_t i = 0; i < frames.size(); ++i) {
        out << "  <g shape-rendering=\"crispEdges\" display=\""
            << (i == 0 ? "inline" : "none") << "\">\n";
        write_slice_cells3d(
            out,
            frames[i].xy,
            masks.xy_solid,
            masks.xy_porosity,
            nx,
            ny,
            xy_x,
            plot_y,
            cell,
            vmin,
            span,
            options.color_map);
        write_slice_cells3d(
            out,
            frames[i].xz,
            masks.xz_solid,
            masks.xz_porosity,
            nx,
            nz,
            xz_x,
            plot_y,
            cell,
            vmin,
            span,
            options.color_map);
        write_slice_cells3d(
            out,
            frames[i].yz,
            masks.yz_solid,
            masks.yz_porosity,
            nz,
            ny,
            yz_x,
            plot_y,
            cell,
            vmin,
            span,
            options.color_map);
        out << "    <animate attributeName=\"display\" values=\"";
        write_visibility_values(out, i, frames.size());
        out << "\" dur=\"" << std::fixed << std::setprecision(6) << duration
            << "s\" repeatCount=\"indefinite\"/>\n";
        out << "  </g>\n";
    }

    out << "  <g fill=\"none\" stroke=\"#24292f\" stroke-width=\"1\">";
    out << "<rect x=\"" << xy_x << "\" y=\"" << plot_y << "\" width=\""
        << xy_width << "\" height=\"" << ny * cell << "\"/>";
    out << "<rect x=\"" << xz_x << "\" y=\"" << plot_y << "\" width=\""
        << xz_width << "\" height=\"" << nz * cell << "\"/>";
    out << "<rect x=\"" << yz_x << "\" y=\"" << plot_y << "\" width=\""
        << yz_width << "\" height=\"" << ny * cell << "\"/></g>\n";
    out << "  <text x=\"" << margin << "\" y=\"" << height - 14
        << "\" font-family=\"Arial, sans-serif\" font-size=\"12\""
           " fill=\"#57606a\">";
    if (!options.footer.empty()) {
        out << options.footer;
    } else {
        out << "Frames: " << frames.size() << "; grid: " << nx << " x " << ny
            << " x " << nz;
    }
    out << "</text>\n</svg>\n";
}

} // namespace lbm
