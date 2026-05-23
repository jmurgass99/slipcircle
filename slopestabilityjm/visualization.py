"""
visualization.py — Matplotlib rendering for slope geometry and slip surfaces.

Functions
---------
plot_slope(slope, ax=None)
    Draw ground profile and shade soil layers.
plot_circle(surface, slope, ax=None)
    Overlay a circular slip surface.
plot_result(slope, surface, result, slices=None, reinforcements=None, ax=None)
    Full plot with FS label, water table, and optional slice/reinforcement
    annotation.

Matplotlib is imported lazily so the rest of the library remains usable
without it installed.
"""
from __future__ import annotations
from math import sqrt
import numpy as np

from .geometry import Slope, CircularSurface


_LAYER_COLORS = ["#d9c39d", "#a87f50", "#c1b39b", "#8c6f4a", "#bfae8c"]


def _get_ax(ax):
    if ax is None:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(9, 5))
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.grid(True, alpha=0.3)
    return ax


def plot_slope(slope: Slope, ax=None):
    ax = _get_ax(ax)
    xs = [p[0] for p in slope.profile]
    ys = [p[1] for p in slope.profile]

    # Determine plotting depth: 1.5*height below the toe.
    H = max(slope.height, 1.0)
    y_min = min(ys) - 1.5 * H
    x_min, x_max = min(xs), max(xs)

    # Shade soil layers (top -> bottom) using cumulative layer_depths.
    y_top_layer = max(ys)
    for i, (layer, depth) in enumerate(zip(slope.layers, slope.layer_depths)):
        y_bot_layer = max(ys) - depth
        color = _LAYER_COLORS[i % len(_LAYER_COLORS)]
        # Polygon: ground surface on top, horizontal line at y_bot_layer beneath.
        poly_x = list(xs) + [x_max, x_min]
        poly_y = list(ys) + [y_bot_layer, y_bot_layer]
        if i > 0:
            # Clip from above by previous layer bottom
            prev_y = max(ys) - slope.layer_depths[i - 1]
            poly_x = [x_min, x_max, x_max, x_min]
            poly_y = [prev_y, prev_y, y_bot_layer, y_bot_layer]
        ax.fill(poly_x, poly_y, color=color, alpha=0.55,
                label=f"{layer.name}  c'={layer.cohesion} φ'={layer.friction_angle}")

    # Ground profile line
    ax.plot(xs, ys, "k-", lw=2)

    # Water table
    if slope.water_table_y is not None:
        ax.axhline(slope.water_table_y, color="#1f77b4", ls="--", lw=1.2,
                   label=f"Water table  y={slope.water_table_y} m")

    ax.set_ylim(y_min, max(ys) + 0.2 * H)
    ax.set_xlim(x_min - 0.05 * (x_max - x_min), x_max + 0.05 * (x_max - x_min))
    return ax


def plot_circle(surface: CircularSurface, slope: Slope, ax=None,
                color="red", lw=1.8, label=None):
    ax = _get_ax(ax)
    cx, cy, R = surface.cx, surface.cy, surface.R

    # Only draw the arc inside the slope.
    xs = [p[0] for p in slope.profile]
    x_lo = max(cx - R, min(xs))
    x_hi = min(cx + R, max(xs))
    if x_hi <= x_lo:
        return ax
    x = np.linspace(x_lo, x_hi, 200)
    y = cy - np.sqrt(np.maximum(R ** 2 - (x - cx) ** 2, 0.0))
    ax.plot(x, y, color=color, lw=lw, label=label)
    ax.plot([cx], [cy], color=color, marker="+", ms=10)
    return ax


def plot_result(slope: Slope, surface: CircularSurface, result: dict,
                slices=None, reinforcements=None, ax=None, title=None):
    ax = plot_slope(slope, ax=ax)
    FS = result.get("FS")
    method = result.get("method", "LEM")
    lab = f"{method}: FS = {FS}"
    plot_circle(surface, slope, ax=ax, label=lab)

    if slices:
        for sl in slices:
            ax.plot([sl["x_left"], sl["x_left"]],
                    [sl["y_bot"], sl["y_top"]],
                    color="gray", lw=0.4, alpha=0.7)
        # last right boundary
        sl = slices[-1]
        ax.plot([sl["x_right"], sl["x_right"]],
                [sl["y_bot"], sl["y_top"]],
                color="gray", lw=0.4, alpha=0.7)

    if reinforcements:
        from .reinforcement import GeosyntheticLayer, SoilNail
        for r in reinforcements:
            if isinstance(r, GeosyntheticLayer):
                ax.plot([r.x_left, r.x_right], [r.y, r.y],
                        color="green", lw=2.5)
            elif isinstance(r, SoilNail):
                tx, ty = r.tip
                ax.plot([r.x_head, tx], [r.y_head, ty], color="purple", lw=2)

    ax.legend(loc="lower right", fontsize=8)
    ax.set_title(title or f"Slope stability — FS = {FS}")
    return ax
