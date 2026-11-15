"""Gamma index (Low et al., 1998) for dose distributions on regular grids, 1D to 3D.

For every reference point r the gamma value is

    gamma(r) = min over r' of sqrt( |r' - r|^2 / DTA^2 + (D_eval(r') - D_ref(r))^2 / tol(r)^2 )

where tol is dose_percent % of the maximum reference dose (global) or of the local reference dose.
The evaluated distribution is linearly interpolated, and candidate positions r' are searched on a
fine grid of displacements inside a sphere of radius max_gamma * DTA.
"""
import numpy as np
from scipy.ndimage import map_coordinates


def gamma_index(reference, evaluated, spacing_mm, dose_percent=3.0, dta_mm=3.0, local=False,
                cutoff_percent=10.0, max_gamma=2.0, step_mm=None):
    """Return an array of gamma values (same shape as the inputs).

    Points where the reference dose is below cutoff_percent of its maximum are returned as NaN.
    Values above max_gamma are reported as max_gamma. Both distributions must share one grid.
    """
    ref = np.asarray(reference, dtype=float)
    ev = np.asarray(evaluated, dtype=float)
    if ref.shape != ev.shape:
        raise ValueError("reference and evaluated must have the same shape")
    if not (np.isfinite(ref).all() and np.isfinite(ev).all()):
        raise ValueError("doses must be finite")
    spacing = np.broadcast_to(np.asarray(spacing_mm, dtype=float), (ref.ndim,)).copy()
    if (spacing <= 0).any() or dose_percent <= 0 or dta_mm <= 0 or max_gamma <= 0:
        raise ValueError("spacing, dose_percent, dta_mm and max_gamma must be positive")
    dmax = ref.max()
    if dmax <= 0:
        raise ValueError("reference dose must have a positive maximum")
    if local and cutoff_percent <= 0:
        raise ValueError("local gamma needs a positive cutoff to avoid dividing by zero dose")

    tol = dose_percent / 100.0 * (ref if local else dmax)
    mask = ref >= cutoff_percent / 100.0 * dmax
    if local:
        tol = np.where(mask, tol, 1.0)             # masked points are discarded; avoid division by zero

    step = dta_mm / 10.0 if step_mm is None else float(step_mm)
    radius = max_gamma * dta_mm
    n = int(np.floor(radius / step + 1e-9))
    axis = np.arange(-n, n + 1) * step
    grids = np.meshgrid(*([axis] * ref.ndim), indexing="ij")
    disp = np.stack([g.ravel() for g in grids], axis=1)
    dist2 = (disp ** 2).sum(axis=1)
    inside = dist2 <= radius ** 2 + 1e-9
    disp, dist2 = disp[inside], dist2[inside]

    voxel = np.indices(ref.shape).astype(float)
    best = np.full(ref.shape, np.inf)
    for d, d2 in zip(disp, dist2):
        coords = voxel + (d / spacing).reshape((-1,) + (1,) * ref.ndim)
        cand = map_coordinates(ev, coords, order=1, mode="constant", cval=np.nan)
        g2 = d2 / dta_mm ** 2 + ((cand - ref) / tol) ** 2
        best = np.fmin(best, g2)              # fmin ignores candidates outside the grid (NaN)

    gamma = np.minimum(np.sqrt(best), max_gamma)
    gamma[~mask] = np.nan
    return gamma


def pass_rate(gamma, tolerance=1e-9):
    """Percentage of evaluated (non-NaN) points with gamma <= 1.

    A numerical tolerance of 1e-9 keeps points whose gamma is exactly 1 in theory (for example a
    3% dose error against a 3% tolerance) from being decided by floating-point rounding.
    """
    g = gamma[~np.isnan(gamma)]
    if g.size == 0:
        raise ValueError("no points above the dose cutoff")
    return 100.0 * float(np.mean(g <= 1.0 + tolerance))
