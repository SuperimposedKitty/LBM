#!/usr/bin/env python3
"""Create an animated SVG with tracer particles advected by the LBM velocity field."""

import argparse
import csv
import glob
import math
from pathlib import Path
from typing import Dict, List, Tuple


def lerp(a: float, b: float, t: float) -> float:
    return a + (b - a) * t


def ramp(t: float) -> str:
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


def read_rows(path: Path) -> List[Dict[str, float]]:
    rows = []
    with path.open(newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            ux = float(row["ux"])
            uy = float(row["uy"])
            rows.append(
                {
                    "x": int(row["x"]),
                    "y": int(row["y"]),
                    "ux": ux,
                    "uy": uy,
                    "speed": math.hypot(ux, uy),
                }
            )
    if not rows:
        raise ValueError("No rows found in %s" % path)
    return rows


def velocity_at(rows: List[Dict[str, float]], nx: int, ny: int, x: float, y: float) -> Tuple[float, float]:
    ix = int(round(x)) % nx
    iy = int(round(y)) % ny
    row = rows[iy * nx + ix]
    return row["ux"], row["uy"]


def integrate_tracers(frames: List[List[Dict[str, float]]], nx: int, ny: int,
                      cols: int, rows_count: int, gain: float) -> List[List[Tuple[float, float]]]:
    tracers = []
    for j in range(rows_count):
        y = (j + 0.5) * ny / rows_count
        for i in range(cols):
            x = (i + 0.5) * nx / cols
            tracers.append((x, y))

    paths = []
    for x0, y0 in tracers:
        x = x0
        y = y0
        path = []
        for frame in frames:
            path.append((x, y))
            ux, uy = velocity_at(frame, nx, ny, x, y)
            x = (x + ux * gain) % nx
            y = (y + uy * gain) % ny
        paths.append(path)
    return paths


def background_rects(rows: List[Dict[str, float]], nx: int, ny: int, cell: int,
                     margin: int, vmin: float, vmax: float) -> str:
    span = vmax - vmin if vmax > vmin else 1.0
    parts = []
    for row in rows:
        t = (row["speed"] - vmin) / span
        x = int(row["x"]) * cell + margin
        y = (ny - 1 - int(row["y"])) * cell + margin + 34
        parts.append(
            '<rect x="%d" y="%d" width="%d" height="%d" fill="%s"/>'
            % (x, y, cell, cell, ramp(t))
        )
    return "".join(parts)


def polyline(points: List[Tuple[float, float]], cell: int, margin: int, ny: int) -> str:
    coords = []
    for x, y in points:
        px = margin + x * cell
        py = margin + 34 + (ny - y) * cell
        coords.append("%.2f,%.2f" % (px, py))
    return " ".join(coords)


def write_svg(input_glob: str, output: Path, cell: int, fps: float,
              tracer_cols: int, tracer_rows: int, gain: float) -> None:
    paths = [Path(p) for p in sorted(glob.glob(input_glob))]
    if not paths:
        raise ValueError("No frame files matched: %s" % input_glob)

    frames = [read_rows(path) for path in paths]
    nx = max(int(row["x"]) for row in frames[0]) + 1
    ny = max(int(row["y"]) for row in frames[0]) + 1
    speeds = [row["speed"] for frame in frames for row in frame]
    vmin = min(speeds)
    vmax = max(speeds)

    margin = 24
    plot_w = nx * cell
    plot_h = ny * cell
    width = plot_w + 48
    height = plot_h + 84
    frame_duration = 1.0 / fps
    total_duration = frame_duration * len(frames)

    bg_groups = []
    for i, frame in enumerate(frames):
        visible = ["inline" if j == i else "none" for j in range(len(frames))]
        bg_groups.append(
            '<g shape-rendering="crispEdges" display="%s">%s'
            '<animate attributeName="display" values="%s;" dur="%.6fs" repeatCount="indefinite"/>'
            '</g>'
            % (
                "inline" if i == 0 else "none",
                background_rects(frame, nx, ny, cell, margin, vmin, vmax),
                ";".join(visible),
                total_duration,
            )
        )

    tracer_paths = integrate_tracers(frames, nx, ny, tracer_cols, tracer_rows, gain)
    tracers = []
    for path in tracer_paths:
        points = polyline(path, cell, margin, ny)
        x0 = margin + path[0][0] * cell
        y0 = margin + 34 + (ny - path[0][1]) * cell
        tracers.append(
            '<circle cx="%.2f" cy="%.2f" r="2.2" fill="#ffffff" stroke="#111827" stroke-width="0.7">'
            '<animateMotion dur="%.6fs" repeatCount="indefinite" path="M %s"/>'
            '</circle>'
            % (x0, y0, total_duration, points.replace(" ", " L "))
        )

    svg = '''<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}" role="img" aria-label="D2Q9 LBM tracer particle animation">
  <rect width="100%" height="100%" fill="#f7f8fa"/>
  <text x="{margin}" y="{margin}" font-family="Arial, sans-serif" font-size="18" font-weight="700" fill="#1f2328">D2Q9 LBM tracer animation</text>
  <text x="{margin}" y="{sub}" font-family="Arial, sans-serif" font-size="12" fill="#57606a">White particles are advected by the computed velocity field; background color is speed.</text>
  {background}
  <rect x="{margin}" y="{top}" width="{plot_w}" height="{plot_h}" fill="none" stroke="#24292f" stroke-width="1"/>
  <g>
    {tracers}
  </g>
  <text x="{margin}" y="{footer}" font-family="Arial, sans-serif" font-size="12" fill="#57606a">Frames: {count}; grid: {nx} x {ny}; particle advection gain: {gain}</text>
</svg>
'''.format(
        width=width,
        height=height,
        margin=margin,
        sub=margin + 18,
        background="\n  ".join(bg_groups),
        top=margin + 34,
        plot_w=plot_w,
        plot_h=plot_h,
        tracers="\n    ".join(tracers),
        footer=height - 14,
        count=len(frames),
        nx=nx,
        ny=ny,
        gain=gain,
    )

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(svg, encoding="utf-8")
    print("wrote: %s" % output)


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate tracer-particle animation from CSV frames.")
    parser.add_argument("--input-glob", default="output/frames/periodic_shear_*.csv")
    parser.add_argument("--output", default="output/periodic_shear_tracers.svg")
    parser.add_argument("--cell", type=int, default=5)
    parser.add_argument("--fps", type=float, default=8.0)
    parser.add_argument("--tracer-cols", type=int, default=18)
    parser.add_argument("--tracer-rows", type=int, default=9)
    parser.add_argument("--gain", type=float, default=18.0)
    args = parser.parse_args()

    write_svg(
        args.input_glob,
        Path(args.output),
        max(1, args.cell),
        args.fps,
        max(1, args.tracer_cols),
        max(1, args.tracer_rows),
        args.gain,
    )


if __name__ == "__main__":
    main()

