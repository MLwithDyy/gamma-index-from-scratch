# Plan (fixed before running the scenarios)

Locked on: 15 November 2026

## Goal

Implement the gamma index from scratch, show that it is correct on cases with known answers, and show how pass rates behave under controlled dose-distribution errors. This is a learning and demonstration project, not clinical software.

## Scope

- Regular grids in 1D, 2D or 3D. Reference and evaluated distributions must be on the same grid (resample first otherwise).
- Gamma as defined by Low et al. (1998): dose-difference and distance-to-agreement combined, minimum over a search neighbourhood, with linear interpolation of the evaluated distribution on a fine search grid.
- Global normalisation (dose tolerance relative to the maximum reference dose) and local normalisation (relative to the local reference dose).
- Dose cutoff: points below 10% of the maximum reference dose are ignored. Gamma values are capped at 2.

Not in scope: reading DICOM files, comparing different grids, clinical validation, comparison with commercial software, 3D clinical workflows.

## Validation (unit tests)

Each test has an expected value derived on paper, not from the code:

1. Identical distributions give gamma 0 everywhere.
2. Points below the cutoff are excluded.
3. A pure shift of s mm gives gamma <= s / DTA everywhere (the shifted displacement gives zero dose difference).
4. A uniform scale error e on a flat region gives gamma = e / tolerance exactly (2% against 3% gives 0.667, 4% gives 1.333).
5. Local normalisation obeys the same bound for a scale error.
6. A 2D diagonal shift obeys gamma <= |shift| / DTA.
7. A shift larger than the DTA lowers the pass rate below 100% but not to nothing.
8. The vectorised N-D implementation agrees with an independent brute-force 1D implementation.
9. Invalid inputs raise errors.

## Scenarios (fixed)

A synthetic 2D field: 100 x 100 mm flat field with a 3 mm error-function penumbra on a 200 x 200 mm grid with 1 mm spacing. Evaluated distributions:

| Family | Values |
|---|---|
| Identical | no change |
| Shift along x | 1, 2, 4, 6 mm |
| Uniform scale | +1, +2, +4, +5 % |
| Gaussian noise | standard deviation 0.5, 1, 2 % of the maximum (seed 42) |

Criteria: 3%/3 mm, 2%/2 mm and 1%/1 mm with global normalisation, and 3%/3 mm with local normalisation. The search step is DTA/10 for global and DTA/20 for local normalisation; a convergence check compares each with half that step on a subset (see the deviation log for why these values).

## Expectations written in advance (from theory)

- Identical: 100% pass for every criterion.
- Pass rate does not increase as the perturbation grows, within a family, for a given criterion.
- Tighter criteria give lower or equal pass rates for the same perturbation.
- A shift at or below the DTA passes at 100% (bound from test 3). Shifts above the DTA fail only where the dose gradient is large (the penumbra), so the pass rate stays well above 0%.
- A uniform scale error below the tolerance passes at 100% with global normalisation; above it, flat regions fail.
- Local normalisation is at least as strict as global for the same criteria on regions below the maximum dose.
- Noise: pass rate falls as noise approaches the dose tolerance, and tighter criteria suffer first.
- Halving the search step changes pass rates by only a small amount (to be checked, not assumed; local normalisation is expected to be the most sensitive).

Anything that contradicts these expectations is reported and investigated, not hidden.

## Reporting

All scenarios in a table and figure, expectations checked against results, runtime for each grid size, and limitations: linear interpolation error on coarse grids, same-grid requirement, synthetic fields only, a single implementation compared with itself and with theory (no external reference implementation).

## Deviation log

| Date | Change | Reason |
|---|---|---|
| 15 November 2026 | Search step in the scenario script changed from DTA/5 to DTA/10 (global) and DTA/20 (local) | A development run at DTA/5 differed from DTA/10 by up to 4.8 percentage points for 2%/2 mm, and local gamma still changed from 95.5% to 99.3% between DTA/10 and DTA/20 |
| 15 November 2026 | Pass means gamma <= 1 + 1e-9 | A 3% scale error against a 3% tolerance gave gamma = 1 + 9e-16 and a pass rate of about 50% from rounding alone |
| 15 November 2026 | Disclosure | The scenario script was run during development, before this plan was committed, so the first results were seen before the expectations were checked |
