#!/usr/bin/env python3
"""Create a dependency-free animated SVG from LBM CSV frame files."""

import argparse
import csv
import glob
import math
from pathlib import Path
from typing import Dict, List


def lerp(a: float, b: float, t: float) -> float:
    return a + (b - a) * t


def sequential_ramp(t: float) -> str:
    stops = [
        (35, 55, 115),
        (35, 105, 165),
        (35, 150, 135),
        (170, 190, 80),
        (245, 210, 80),
    ]
    t = max(0.0, min(1.0, t))
    scaled = t * (len(stops) - 1)
    i = min(int(scaled), len(stops) - 2)
    local = scaled - i
    r = round(lerp(stops[i][0], stops[i + 1][0], local))
    g = round(lerp(stops[i][1], stops[i + 1][1], local))
    b = round(lerp(stops[i][2], stops[i + 1][2], local))
    return "#%02x%02x%02x" % (r, g, b)


def phase_ramp(t: float) -> str:
    stops = [
        (31, 78, 121),
        (104, 166, 190),
        (244, 246, 247),
        (224, 122, 95),
        (150, 45, 56),
    ]
    t = max(0.0, min(1.0, t))
    scaled = t * (len(stops) - 1)
    i = min(int(scaled), len(stops) - 2)
    local = scaled - i
    r = round(lerp(stops[i][0], stops[i + 1][0], local))
    g = round(lerp(stops[i][1], stops[i + 1][1], local))
    b = round(lerp(stops[i][2], stops[i + 1][2], local))
    return "#%02x%02x%02x" % (r, g, b)


def color_for(field: str, t: float) -> str:
    if field in {"phi", "phase"}:
        return phase_ramp(t)
    return sequential_ramp(t)


def read_rows(path: Path) -> List[Dict[str, float]]:
    rows = []
    with path.open(newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            parsed: Dict[str, float] = {}
            for key, value in row.items():
                try:
                    parsed[key] = float(value)
                except (TypeError, ValueError):
                    continue
            if "x" not in parsed or "y" not in parsed:
                raise ValueError("CSV rows in %s must include numeric x and y columns" % path)
            if "ux" in parsed and "uy" in parsed:
                parsed["speed"] = math.hypot(parsed["ux"], parsed["uy"])
            rows.append(parsed)
    if not rows:
        raise ValueError("No rows found in %s" % path)
    return rows


def frame_rects(rows: List[Dict[str, float]], field: str, cell: int, nx: int, ny: int,
                vmin: float, vmax: float, margin: int) -> str:
    span = vmax - vmin if vmax > vmin else 1.0
    parts = []
    for row in rows:
        t = (row[field] - vmin) / span
        x = int(row["x"]) * cell + margin
        y = (ny - 1 - int(row["y"])) * cell + margin + 30
        fill = "#2f343b" if row.get("solid", 0.0) >= 0.5 else color_for(field, t)
        parts.append(
            '<rect x="%d" y="%d" width="%d" height="%d" fill="%s"/>'
            % (x, y, cell, cell, fill)
        )
    return "".join(parts)


def write_animation(
    input_glob: str,
    output: Path,
    field: str,
    cell: int,
    fps: float,
    vmin_arg: float = None,
    vmax_arg: float = None,
) -> None:
    paths = [Path(p) for p in sorted(glob.glob(input_glob))]
    if not paths:
        raise ValueError("No frame files matched: %s" % input_glob)

    frames = [read_rows(path) for path in paths]
    if field not in frames[0][0]:
        available = ", ".join(sorted(frames[0][0].keys()))
        raise ValueError("Field %r not found. Available numeric fields: %s" % (field, available))

    nx = max(int(row["x"]) for row in frames[0]) + 1
    ny = max(int(row["y"]) for row in frames[0]) + 1
    values = [row[field] for frame in frames for row in frame if row.get("solid", 0.0) < 0.5]
    if not values:
        values = [row[field] for frame in frames for row in frame]
    if field in {"phi", "phase"} and vmin_arg is None and vmax_arg is None:
        vmin = -1.0
        vmax = 1.0
    else:
        vmin = min(values) if vmin_arg is None else vmin_arg
        vmax = max(values) if vmax_arg is None else vmax_arg

    margin = 24
    plot_w = nx * cell
    plot_h = ny * cell
    width = plot_w + 48
    height = plot_h + 78
    frame_duration = 1.0 / fps
    total_duration = frame_duration * len(frames)

    groups = []
    for i, rows in enumerate(frames):
        begin = i * frame_duration
        visible_values = []
        for j in range(len(frames)):
            visible_values.append("inline" if j == i else "none")
        values_attr = ";".join(visible_values) + ";"
        groups.append(
            '<g shape-rendering="crispEdges" display="%s">%s'
            '<animate attributeName="display" values="%s" dur="%.6fs" repeatCount="indefinite"/>'
            '</g>'
            % (
                "inline" if i == 0 else "none",
                frame_rects(rows, field, cell, nx, ny, vmin, vmax, margin),
                values_attr,
                total_duration,
            )
        )

    svg = '''<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-label="Animated D2Q9 LBM {field} field">
  <rect width="100%" height="100%" fill="#f7f8fa"/>
  <text x="{margin}" y="{margin}" font-family="Arial, sans-serif" font-size="18" font-weight="700" fill="#1f2328">D2Q9 LBM {field} animation</text>
  {groups}
  <rect x="{margin}" y="{top}" width="{plot_w}" height="{plot_h}" fill="none" stroke="#24292f" stroke-width="1"/>
  <text x="{margin}" y="{footer}" font-family="Arial, sans-serif" font-size="12" fill="#57606a">Frames: {count}; grid: {nx} x {ny}; global range: {vmin:.6g} to {vmax:.6g}</text>
</svg>
'''.format(
        width=width,
        height=height,
        field=field,
        margin=margin,
        groups="\n  ".join(groups),
        top=margin + 30,
        plot_w=plot_w,
        plot_h=plot_h,
        footer=height - 14,
        count=len(frames),
        nx=nx,
        ny=ny,
        vmin=vmin,
        vmax=vmax,
    )

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(svg, encoding="utf-8")
    print("wrote: %s" % output)


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate an animated SVG from CSV frames.")
    parser.add_argument("--input-glob", default="output/frames/periodic_shear_*.csv")
    parser.add_argument("--output", default="output/periodic_shear_speed_animation.svg")
    parser.add_argument("--field", default="speed", help="Numeric CSV field to visualize")
    parser.add_argument("--cell", type=int, default=4)
    parser.add_argument("--fps", type=float, default=8.0)
    parser.add_argument("--vmin", type=float, default=None, help="Optional lower color-map bound")
    parser.add_argument("--vmax", type=float, default=None, help="Optional upper color-map bound")
    args = parser.parse_args()

    write_animation(
        args.input_glob,
        Path(args.output),
        args.field,
        max(1, args.cell),
        args.fps,
        args.vmin,
        args.vmax,
    )


if __name__ == "__main__":
    main()
