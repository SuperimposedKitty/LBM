#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace lbm {

enum class ColorMap {
    Sequential,
    Phase,
    PhaseHighContrast
};

struct ScalarAnimationOptions {
    int cell = 4;
    double fps = 8.0;
    std::string title = "LBM animation";
    std::string footer;
    ColorMap color_map = ColorMap::Sequential;
    bool fixed_range = false;
    double vmin = 0.0;
    double vmax = 1.0;
    bool highlight_solid_surface = false;
    int solid_surface_y = 1;
    int solid_surface_thickness = 3;
    std::string solid_surface_color = "#111827";
};

inline std::filesystem::path result_path(const std::string& filename) {
#ifdef LBM_RESULT_DIR
    return std::filesystem::path(LBM_RESULT_DIR) / filename;
#else
    return std::filesystem::path("result") / filename;
#endif
}

inline double svg_lerp(double a, double b, double t) {
    return a + (b - a) * t;
}

inline std::string svg_rgb_hex(int r, int g, int b) {
    std::ostringstream out;
    out << '#' << std::hex << std::setfill('0') << std::setw(2) << std::clamp(r, 0, 255)
        << std::setw(2) << std::clamp(g, 0, 255) << std::setw(2) << std::clamp(b, 0, 255);
    return out.str();
}

inline std::string sequential_color(double t) {
    constexpr std::array<std::array<int, 3>, 5> stops{{
        {{35, 55, 115}},
        {{35, 105, 165}},
        {{35, 150, 135}},
        {{170, 190, 80}},
        {{245, 210, 80}},
    }};

    t = std::clamp(t, 0.0, 1.0);
    const double scaled = t * static_cast<double>(stops.size() - 1);
    const int i = std::min(static_cast<int>(scaled), static_cast<int>(stops.size()) - 2);
    const double local = scaled - static_cast<double>(i);
    const int r = static_cast<int>(std::round(svg_lerp(stops[i][0], stops[i + 1][0], local)));
    const int g = static_cast<int>(std::round(svg_lerp(stops[i][1], stops[i + 1][1], local)));
    const int b = static_cast<int>(std::round(svg_lerp(stops[i][2], stops[i + 1][2], local)));
    return svg_rgb_hex(r, g, b);
}

inline std::string phase_color(double t) {
    // 相场色带采用蓝-黄-红；phi 接近 0 的界面用黄色突出显示。
    constexpr std::array<std::array<int, 3>, 5> stops{{
        {{31, 78, 121}},
        {{104, 166, 190}},
        {{247, 201, 72}},
        {{224, 122, 95}},
        {{150, 45, 56}},
    }};

    t = std::clamp(t, 0.0, 1.0);
    const double scaled = t * static_cast<double>(stops.size() - 1);
    const int i = std::min(static_cast<int>(scaled), static_cast<int>(stops.size()) - 2);
    const double local = scaled - static_cast<double>(i);
    const int r = static_cast<int>(std::round(svg_lerp(stops[i][0], stops[i + 1][0], local)));
    const int g = static_cast<int>(std::round(svg_lerp(stops[i][1], stops[i + 1][1], local)));
    const int b = static_cast<int>(std::round(svg_lerp(stops[i][2], stops[i + 1][2], local)));
    return svg_rgb_hex(r, g, b);
}

inline std::string phase_high_contrast_color(double t) {
    // 液滴撞击用更亮的红色和更深的蓝色，避免红色液滴在动画中显得发暗或消失。
    constexpr std::array<std::array<int, 3>, 5> stops{{
        {{18, 52, 104}},
        {{62, 127, 190}},
        {{255, 214, 74}},
        {{244, 112, 82}},
        {{218, 38, 52}},
    }};

    t = std::clamp(t, 0.0, 1.0);
    const double scaled = t * static_cast<double>(stops.size() - 1);
    const int i = std::min(static_cast<int>(scaled), static_cast<int>(stops.size()) - 2);
    const double local = scaled - static_cast<double>(i);
    const int r = static_cast<int>(std::round(svg_lerp(stops[i][0], stops[i + 1][0], local)));
    const int g = static_cast<int>(std::round(svg_lerp(stops[i][1], stops[i + 1][1], local)));
    const int b = static_cast<int>(std::round(svg_lerp(stops[i][2], stops[i + 1][2], local)));
    return svg_rgb_hex(r, g, b);
}

inline std::string scalar_color(ColorMap color_map, double t) {
    if (color_map == ColorMap::Phase) {
        return phase_color(t);
    }
    if (color_map == ColorMap::PhaseHighContrast) {
        return phase_high_contrast_color(t);
    }
    return sequential_color(t);
}

inline void write_visibility_values(std::ofstream& out, std::size_t active, std::size_t count) {
    for (std::size_t i = 0; i < count; ++i) {
        out << (i == active ? "inline" : "none") << ';';
    }
}

inline void write_scalar_animation_svg(
    const std::filesystem::path& path,
    const std::vector<std::vector<double>>& frames,
    const std::vector<std::uint8_t>& solid,
    int nx,
    int ny,
    const ScalarAnimationOptions& options,
    const std::vector<double>& porosity = {}) {
    if (nx <= 0 || ny <= 0) {
        throw std::invalid_argument("Animation grid dimensions must be positive.");
    }
    if (frames.empty()) {
        throw std::runtime_error("No frames available for animation output.");
    }
    const std::size_t cell_count = static_cast<std::size_t>(nx) * ny;
    if (!solid.empty() && solid.size() != cell_count) {
        throw std::invalid_argument("Solid mask size does not match the animation grid.");
    }
    if (!porosity.empty() && porosity.size() != cell_count) {
        throw std::invalid_argument("Porosity mask size does not match the animation grid.");
    }
    for (const auto& frame : frames) {
        if (frame.size() != cell_count) {
            throw std::invalid_argument("Frame size does not match the animation grid.");
        }
    }

    double vmin = options.vmin;
    double vmax = options.vmax;
    if (!options.fixed_range) {
        vmin = std::numeric_limits<double>::max();
        vmax = std::numeric_limits<double>::lowest();
        for (const auto& frame : frames) {
            for (std::size_t s = 0; s < frame.size(); ++s) {
                if (!solid.empty() && solid[s]) {
                    continue;
                }
                vmin = std::min(vmin, frame[s]);
                vmax = std::max(vmax, frame[s]);
            }
        }
        if (vmin > vmax) {
            vmin = 0.0;
            vmax = 1.0;
        }
    }
    const double span = vmax > vmin ? vmax - vmin : 1.0;

    const int cell = std::max(1, options.cell);
    const double fps = options.fps > 0.0 ? options.fps : 8.0;
    const int margin = 24;
    const int plot_w = nx * cell;
    const int plot_h = ny * cell;
    const int width = plot_w + 48;
    const int height = plot_h + 78;
    const double duration = static_cast<double>(frames.size()) / fps;

    std::filesystem::create_directories(path.parent_path());
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("Failed to open animation output file: " + path.string());
    }

    out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    out << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << width << "\" height=\""
        << height << "\" viewBox=\"0 0 " << width << ' ' << height
        << "\" role=\"img\" aria-label=\"" << options.title << "\">\n";
    out << "  <rect width=\"100%\" height=\"100%\" fill=\"#f7f8fa\"/>\n";
    out << "  <text x=\"" << margin << "\" y=\"" << margin
        << "\" font-family=\"Arial, sans-serif\" font-size=\"18\" font-weight=\"700\""
        << " fill=\"#1f2328\">" << options.title << "</text>\n";
    if (!porosity.empty()) {
        out << "  <defs>\n";
        out << "    <pattern id=\"porousHatch\" width=\"6\" height=\"6\""
            << " patternUnits=\"userSpaceOnUse\">\n";
        out << "      <path d=\"M0 6 L6 0\" stroke=\"#1f2328\" stroke-opacity=\"0.42\""
            << " stroke-width=\"1\"/>\n";
        out << "    </pattern>\n";
        out << "  </defs>\n";
    }

    for (std::size_t frame_index = 0; frame_index < frames.size(); ++frame_index) {
        out << "  <g shape-rendering=\"crispEdges\" display=\""
            << (frame_index == 0 ? "inline" : "none") << "\">\n";
        const auto& frame = frames[frame_index];
        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {
                const std::size_t s = static_cast<std::size_t>(y) * nx + x;
                const int px = x * cell + margin;
                const int py = (ny - 1 - y) * cell + margin + 30;
                const double t = (frame[s] - vmin) / span;
                const std::string fill =
                    (!solid.empty() && solid[s]) ? "#2f343b" : scalar_color(options.color_map, t);
                out << "    <rect x=\"" << px << "\" y=\"" << py << "\" width=\"" << cell
                    << "\" height=\"" << cell << "\" fill=\"" << fill << "\"/>\n";
                if (!porosity.empty() && porosity[s] > 0.0 && porosity[s] < 0.999 &&
                    (solid.empty() || !solid[s])) {
                    // 多孔介质叠加斜线纹理，底色仍保留两相相场。
                    out << "    <rect x=\"" << px << "\" y=\"" << py << "\" width=\"" << cell
                        << "\" height=\"" << cell
                        << "\" fill=\"url(#porousHatch)\" fill-opacity=\"0.55\"/>\n";
                }
            }
        }
        out << "    <animate attributeName=\"display\" values=\"";
        write_visibility_values(out, frame_index, frames.size());
        out << "\" dur=\"" << std::fixed << std::setprecision(6) << duration
            << "s\" repeatCount=\"indefinite\"/>\n";
        out << "  </g>\n";
    }

    if (options.highlight_solid_surface) {
        const int surface_y = std::clamp(options.solid_surface_y, 0, ny);
        const int line_y = margin + 30 + (ny - surface_y) * cell;
        out << "  <line x1=\"" << margin << "\" y1=\"" << line_y << "\" x2=\""
            << margin + plot_w << "\" y2=\"" << line_y << "\" stroke=\""
            << options.solid_surface_color << "\" stroke-width=\""
            << std::max(1, options.solid_surface_thickness)
            << "\" stroke-linecap=\"square\"/>\n";
    }

    out << "  <rect x=\"" << margin << "\" y=\"" << margin + 30 << "\" width=\"" << plot_w
        << "\" height=\"" << plot_h
        << "\" fill=\"none\" stroke=\"#24292f\" stroke-width=\"1\"/>\n";
    out << "  <text x=\"" << margin << "\" y=\"" << height - 14
        << "\" font-family=\"Arial, sans-serif\" font-size=\"12\" fill=\"#57606a\">";
    if (!options.footer.empty()) {
        out << options.footer;
    } else {
        out << "Frames: " << frames.size() << "; grid: " << nx << " x " << ny
            << "; range: " << std::setprecision(6) << vmin << " to " << vmax;
    }
    out << "</text>\n";
    out << "</svg>\n";
}

} // namespace lbm
