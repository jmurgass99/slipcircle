"""
search.py — Grid search for the critical circular failure surface.
"""
from __future__ import annotations
import numpy as np
from dataclasses import dataclass
from .geometry import Slope, CircularSurface
from .bishop import BishopAnalysis
from .janbu import JanbuAnalysis
from .janbu_generalized import JanbuGeneralizedAnalysis
from .spencer import SpencerAnalysis
from .morgenstern_price import MorgensternPriceAnalysis


@dataclass
class CriticalSurfaceSearch:
    """
    Grid search for the minimum FS circular failure surface.

    Searches over a grid of circle centers and radii to find the
    critical surface with the lowest factor of safety.

    Parameters
    ----------
    slope : Slope
    method : str  'bishop' (default)
    n_slices : int  slices for stability calculation
    n_cx, n_cy : int  grid points in x and y for circle centers
    n_R : int  number of radii to try per center
    """
    slope: Slope
    method: str = "bishop"
    n_slices: int = 20
    n_cx: int = 10
    n_cy: int = 10
    n_R: int = 8
    max_depth_below_toe: float = 0.5  # multiple of slope height; limits base failures
    require_toe_exit: bool = True     # circle must exit at or near the toe

    def search(self, verbose: bool = False) -> dict:
        """
        Perform grid search over circle centers and radii.

        The search domain is:
          - Centers: x from toe to crest, y from crest to 1.5×height above crest
          - Radii: from H/2 to 2H

        Returns
        -------
        dict with critical FS, critical surface parameters, and all results.
        """
        toe_x, toe_y = self.slope.toe
        crest_x, crest_y = self.slope.crest
        H = self.slope.height

        # Grid of centers: extend laterally beyond toe/crest by 0.5H, and
        # vertically from just above crest to 2H above. Critical circles
        # for typical slopes have centres above the crest.
        margin = 0.5 * max(H, 1.0)
        cx_range = np.linspace(toe_x - margin, crest_x + margin, self.n_cx)
        cy_range = np.linspace(crest_y + 0.2 * H, crest_y + 2.0 * max(H, 1.0), self.n_cy)
        R_range  = np.linspace(max(H * 0.6, 1.0), H * 2.5, self.n_R)

        method_key = self.method.strip().lower()
        if method_key == "bishop":
            analyzer = BishopAnalysis(self.slope, n_slices=self.n_slices)
            method_label = "Bishop Simplified (grid search)"
        elif method_key == "spencer":
            analyzer = SpencerAnalysis(self.slope, n_slices=self.n_slices)
            method_label = "Spencer (grid search)"
        elif method_key in ("morgenstern_price", "morgenstern-price", "mp"):
            analyzer = MorgensternPriceAnalysis(self.slope, n_slices=self.n_slices)
            method_label = "Morgenstern-Price (grid search)"
        elif method_key == "janbu":
            analyzer = JanbuAnalysis(self.slope, n_slices=self.n_slices)
            method_label = "Janbu Simplified (grid search)"
        elif method_key in ("janbu_generalized", "janbu-generalized"):
            analyzer = JanbuGeneralizedAnalysis(self.slope, n_slices=self.n_slices)
            method_label = "Janbu Generalized (grid search)"
        else:
            raise ValueError(
                f"Unsupported search method {self.method!r}. "
                "Use 'bishop', 'spencer', 'morgenstern_price', 'janbu' or 'janbu_generalized'."
            )

        best_FS  = float("inf")
        best_surf = None
        results   = []
        n_valid   = 0

        from .slices import slices_from_circle
        for cx in cx_range:
            for cy in cy_range:
                for R in R_range:
                    surface = CircularSurface(cx, cy, R)
                    # Discard circles whose centre is below ground or that do
                    # not span enough of the slope (trivial shallow slips).
                    if cy <= self.slope.ground_y(cx) + 0.1:
                        continue
                    # Limit base failures: bottom of circle not deeper than
                    # max_depth_below_toe * H below toe elevation.
                    if (cy - R) < toe_y - self.max_depth_below_toe * H:
                        continue
                    sl = slices_from_circle(self.slope, surface, self.n_slices)
                    if len(sl) < 0.5 * self.n_slices:
                        continue
                    span = sl[-1]["x_right"] - sl[0]["x_left"]
                    if span < 0.5 * H:
                        continue
                    # Optional: require the circle to pass near the toe point
                    # (within ~5% H tolerance). This restricts the search to
                    # toe failures, where most LEM critical circles for
                    # homogeneous c'-φ' slopes are found.
                    if self.require_toe_exit:
                        d_toe = ((cx - toe_x) ** 2 + (cy - toe_y) ** 2) ** 0.5
                        if abs(d_toe - R) > 0.1 * H:
                            continue
                    result = analyzer.factor_of_safety(surface)
                    if result["FS"] is None:
                        continue
                    n_valid += 1
                    fs = result["FS"]
                    results.append({
                        "cx": round(cx, 2), "cy": round(cy, 2),
                        "R": round(R, 2), "FS": round(fs, 4)
                    })
                    if 0.3 < fs < best_FS:
                        best_FS = fs
                        best_surf = surface

        if best_surf is None:
            return {
                "critical_FS": None,
                "error": "No valid failure surfaces found",
                "n_surfaces_tried": self.n_cx * self.n_cy * self.n_R,
                "n_valid": n_valid,
            }

        return {
            "method": method_label,
            "critical_FS": round(best_FS, 4),
            "critical_surface": {
                "cx": round(best_surf.cx, 3),
                "cy": round(best_surf.cy, 3),
                "R":  round(best_surf.R, 3),
            },
            "n_surfaces_tried": self.n_cx * self.n_cy * self.n_R,
            "n_valid": n_valid,
            "all_results": sorted(results, key=lambda r: r["FS"])[:20],  # top 20
            "acceptable_EC7": best_FS >= 1.5,
        }

    def optimize_circle(self, initial_surface: CircularSurface | None = None) -> dict:
        """Refine critical circular surface using Nelder-Mead local optimization."""
        from scipy.optimize import minimize
        from .slices import slices_from_circle

        if initial_surface is None:
            grid_res = self.search()
            if grid_res.get("critical_FS") is None:
                return grid_res
            cs = grid_res["critical_surface"]
            x0 = np.array([cs["cx"], cs["cy"], cs["R"]], dtype=float)
        else:
            x0 = np.array([initial_surface.cx, initial_surface.cy, initial_surface.R], dtype=float)

        method_key = self.method.strip().lower()
        if method_key == "bishop":
            analyzer = BishopAnalysis(self.slope, n_slices=self.n_slices)
        elif method_key == "spencer":
            analyzer = SpencerAnalysis(self.slope, n_slices=self.n_slices)
        elif method_key in ("morgenstern_price", "morgenstern-price", "mp"):
            analyzer = MorgensternPriceAnalysis(self.slope, n_slices=self.n_slices)
        elif method_key == "janbu":
            analyzer = JanbuAnalysis(self.slope, n_slices=self.n_slices)
        else:
            analyzer = JanbuGeneralizedAnalysis(self.slope, n_slices=self.n_slices)

        toe_x, toe_y = self.slope.toe
        H = self.slope.height

        def _obj(vec: np.ndarray) -> float:
            cx, cy, R = float(vec[0]), float(vec[1]), float(vec[2])
            if R <= 0.5 * H or cy <= self.slope.ground_y(cx) + 0.1:
                return 1e6
            if (cy - R) < toe_y - self.max_depth_below_toe * H:
                return 1e6
            surf = CircularSurface(cx, cy, R)
            sl = slices_from_circle(self.slope, surf, self.n_slices)
            if len(sl) < 0.5 * self.n_slices:
                return 1e6
            res = analyzer.factor_of_safety(surf)
            fs = res.get("FS")
            if fs is None or not np.isfinite(fs) or fs <= 0.3:
                return 1e6
            return float(fs)

        opt = minimize(_obj, x0, method="Nelder-Mead", options={"maxiter": 60, "xatol": 1e-2, "fatol": 1e-4})
        best_cx, best_cy, best_R = float(opt.x[0]), float(opt.x[1]), float(opt.x[2])
        best_fs = float(opt.fun)
        return {
            "method": f"{self.method} (Nelder-Mead)",
            "critical_FS": round(best_fs, 4),
            "critical_surface": {
                "cx": round(best_cx, 3),
                "cy": round(best_cy, 3),
                "R": round(best_R, 3),
            },
            "acceptable_EC7": best_fs >= 1.5,
        }


def grid_search(slope: Slope, method: str = "bishop", **kwargs) -> dict:
    """Convenience wrapper for :class:`CriticalSurfaceSearch`."""
    return CriticalSurfaceSearch(slope=slope, method=method, **kwargs).search()


def optimize_circle(slope: Slope, method: str = "bishop", **kwargs) -> dict:
    """Convenience wrapper for local Nelder-Mead circle optimization."""
    initial_surface = kwargs.pop("initial_surface", None)
    return CriticalSurfaceSearch(slope=slope, method=method, **kwargs).optimize_circle(initial_surface=initial_surface)
