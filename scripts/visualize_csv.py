#!/usr/bin/env python3
"""Create a dependency-free SVG heatmap from an LBM CSV field."""

import argparse
import csv
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
    return f"#{r:02x}{g:02x}{b:02x}"


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
    return f"#{r:02x}{g:02x}{b:02x}"


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
                raise ValueError(f"CSV rows in {path} must include numeric x and y columns")
            if "ux" in parsed and "uy" in parsed:
                parsed["speed"] = math.hypot(parsed["ux"], parsed["uy"])
            rows.append(parsed)
    if not rows:
        raise ValueError(f"No rows found in {path}")
    return rows


def write_svg(
    rows: List[Dict[str, float]],
    field: str,
    out_path: Path,
    cell: int,
    vmin_arg: float = None,
    vmax_arg: float = None,
) -> None:
    if field not in rows[0]:
        available = ", ".join(sorted(rows[0].keys()))
        raise ValueError(f"Field {field!r} not found. Available numeric fields: {available}")

    nx = max(int(r["x"]) for r in rows) + 1
    ny = max(int(r["y"]) for r in rows) + 1
    values = [r[field] for r in rows if r.get("solid", 0.0) < 0.5]
    if not values:
        values = [r[field] for r in rows]
    if field in {"phi", "phase"} and vmin_arg is None and vmax_arg is None:
        vmin = -1.0
        vmax = 1.0
    else:
        vmin = min(values) if vmin_arg is None else vmin_arg
        vmax = max(values) if vmax_arg is None else vmax_arg
    span = vmax - vmin if vmax > vmin else 1.0

    plot_w = nx * cell
    plot_h = ny * cell
    legend_w = 72
    margin = 24
    width = plot_w + legend_w + margin * 2
    height = plot_h + margin * 2 + 44

    rects = []
    for r in rows:
        value = r[field]
        t = (value - vmin) / span
        x = int(r["x"]) * cell + margin
        y = (ny - 1 - int(r["y"])) * cell + margin + 28
        fill = "#2f343b" if r.get("solid", 0.0) >= 0.5 else color_for(field, t)
        rects.append(f'<rect x="{x}" y="{y}" width="{cell}" height="{cell}" fill="{fill}"/>')

    legend_x = margin + plot_w + 24
    legend_y = margin + 28
    legend = []
    for i in range(plot_h):
        t = 1.0 - i / max(1, plot_h - 1)
        legend.append(
            f'<rect x="{legend_x}" y="{legend_y + i}" width="18" height="1" fill="{color_for(field, t)}"/>'
        )

    svg = f'''<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-label="LBM {field} heatmap">
  <rect width="100%" height="100%" fill="#f7f8fa"/>
  <text x="{margin}" y="{margin}" font-family="Arial, sans-serif" font-size="18" font-weight="700" fill="#1f2328">D2Q9 LBM {field} heatmap</text>
  <g shape-rendering="crispEdges">
    {''.join(rects)}
  </g>
  <rect x="{margin}" y="{margin + 28}" width="{plot_w}" height="{plot_h}" fill="none" stroke="#24292f" stroke-width="1"/>
  <g shape-rendering="crispEdges">
    {''.join(legend)}
  </g>
  <text x="{legend_x + 26}" y="{legend_y + 10}" font-family="Consolas, monospace" font-size="12" fill="#24292f">{vmax:.6g}</text>
  <text x="{legend_x + 26}" y="{legend_y + plot_h}" font-family="Consolas, monospace" font-size="12" fill="#24292f">{vmin:.6g}</text>
  <text x="{margin}" y="{height - 14}" font-family="Arial, sans-serif" font-size="12" fill="#57606a">Input grid: {nx} x {ny}; field: {field}; generated from CSV.</text>
</svg>
'''
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(svg, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate an SVG heatmap from LBM CSV output.")
    parser.add_argument("--input", default="output/periodic_shear.csv", help="Input CSV path")
    parser.add_argument("--output", default="output/periodic_shear_speed.svg", help="Output SVG path")
    parser.add_argument("--field", default="speed", help="Numeric CSV field to visualize")
    parser.add_argument("--cell", type=int, default=6, help="SVG cell size in pixels")
    parser.add_argument("--vmin", type=float, default=None, help="Optional lower color-map bound")
    parser.add_argument("--vmax", type=float, default=None, help="Optional upper color-map bound")
    args = parser.parse_args()

    rows = read_rows(Path(args.input))
    write_svg(rows, args.field, Path(args.output), max(1, args.cell), args.vmin, args.vmax)
    print(f"wrote: {args.output}")


if __name__ == "__main__":
    main()
