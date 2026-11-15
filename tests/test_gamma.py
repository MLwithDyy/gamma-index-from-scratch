"""Tests use cases where the gamma value is known analytically, plus an independent brute-force check."""
import numpy as np
import pytest
from scipy.special import erf

from gammaqa import gamma_index, pass_rate

SPACING = 1.0     # mm
DTA = 3.0         # mm


def field(x, half_width=50.0, sigma=3.0):
    """Flat field with error-function penumbra, maximum 1.0."""
    return 0.5 * (erf((x + half_width) / (sigma * np.sqrt(2))) - erf((x - half_width) / (sigma * np.sqrt(2))))


@pytest.fixture
def x():
    return np.arange(-100.0, 100.0 + SPACING, SPACING)


def test_identical_distributions_give_zero(x):
    d = field(x)
    g = gamma_index(d, d, SPACING, 3, DTA)
    assert abs(np.nanmax(g)) < 1e-12
    assert pass_rate(g) == 100.0


def test_below_cutoff_is_nan(x):
    d = field(x)
    g = gamma_index(d, d, SPACING, 3, DTA, cutoff_percent=10)
    assert np.isnan(g[d < 0.10]).all() and not np.isnan(g[d >= 0.10]).any()


@pytest.mark.parametrize("shift", [1.0, 2.0, 3.0])
def test_pure_shift_gamma_is_bounded_by_shift_over_dta(x, shift):
    # Evaluated = reference moved by `shift` mm. A search step of 1 mm puts that displacement exactly on
    # the voxel grid, so the dose difference there is zero and gamma <= shift / DTA.
    g = gamma_index(field(x), field(x - shift), SPACING, 3, DTA, step_mm=1.0)
    assert np.nanmax(g) <= shift / DTA + 1e-9


def test_global_scale_on_plateau_is_exact(x):
    d = field(x)
    g = gamma_index(d, 1.02 * d, SPACING, 3, DTA)
    plateau = np.abs(x) < 30                      # dose is 1.0 here, so gamma = 0.02 / 0.03
    assert np.allclose(g[plateau], 2.0 / 3.0, atol=1e-9)
    assert np.nanmax(g) <= 2.0 / 3.0 + 1e-9


def test_local_scale_gives_same_bound(x):
    d = field(x)
    g = gamma_index(d, 1.02 * d, SPACING, 3, DTA, local=True)
    assert np.nanmax(g) <= 2.0 / 3.0 + 1e-9


def test_scale_beyond_tolerance_fails_on_plateau(x):
    d = field(x)
    g = gamma_index(d, 1.04 * d, SPACING, 3, DTA)
    assert np.allclose(g[np.abs(x) < 30], 4.0 / 3.0, atol=1e-9)


def test_error_exactly_at_tolerance_passes(x):
    # 3% scale error against a 3% tolerance: gamma is exactly 1 on the plateau and must count as a pass.
    d = field(x)
    g = gamma_index(d, 1.03 * d, SPACING, 3, DTA)
    assert pass_rate(g) == 100.0


def test_large_shift_lowers_pass_rate(x):
    g = gamma_index(field(x), field(x - 6.0), SPACING, 3, DTA)
    assert 50.0 < pass_rate(g) < 100.0


def test_2d_diagonal_shift_bound():
    ax = np.arange(-100.0, 101.0, 1.0)
    X, Y = np.meshgrid(ax, ax, indexing="ij")
    ref = field(X) * field(Y)
    ev = field(X - 1.0) * field(Y - 1.0)
    g = gamma_index(ref, ev, SPACING, 3, DTA, step_mm=1.0)
    assert np.nanmax(g) <= np.hypot(1.0, 1.0) / DTA + 1e-9


def test_matches_brute_force_1d():
    rng = np.random.default_rng(0)
    xs = np.arange(0.0, 120.0, 1.0)
    ref = sum(a * np.exp(-0.5 * ((xs - c) / s) ** 2) for a, c, s in
              zip(rng.uniform(0.3, 1.0, 4), rng.uniform(20, 100, 4), rng.uniform(5, 15, 4)))
    ev = ref * 1.01 + 0.02 * np.sin(xs / 7.0)
    g = gamma_index(ref, ev, SPACING, 3, DTA, cutoff_percent=20)

    dmax, step = ref.max(), DTA / 10.0
    shifts = np.arange(-int(2 * DTA / step), int(2 * DTA / step) + 1) * step
    shifts = shifts[np.abs(shifts) <= 2 * DTA + 1e-9]
    expected = np.full(xs.shape, np.nan)
    for i in np.where(ref >= 0.2 * dmax)[0]:
        best = np.inf
        for s in shifts:
            xe = xs[i] + s
            if xe < xs[0] or xe > xs[-1]:
                continue
            best = min(best, np.sqrt((s / DTA) ** 2 + ((np.interp(xe, xs, ev) - ref[i]) / (0.03 * dmax)) ** 2))
        expected[i] = min(best, 2.0)
    assert np.allclose(g, expected, equal_nan=True, atol=1e-9)


def test_input_validation():
    a = np.ones(10)
    with pytest.raises(ValueError):
        gamma_index(a, np.ones(11), 1.0)
    with pytest.raises(ValueError):
        gamma_index(a, a, 0.0)
    with pytest.raises(ValueError):
        gamma_index(np.zeros(10), np.zeros(10), 1.0)
    with pytest.raises(ValueError):
        gamma_index(a, a, 1.0, local=True, cutoff_percent=0)
