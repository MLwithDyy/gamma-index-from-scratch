# Exploratory, not in PLAN.md

import numpy as np
from scipy.special import erf
from gammaqa import gamma_index, pass_rate

f = lambda x, a=50.0, s=3.0: 0.5 * (erf((x + a) / (s * 2 ** .5)) - erf((x - a) / (s * 2 ** .5)))
ax = np.arange(-100.0, 101.0, 1.0); X, Y = np.meshgrid(ax, ax, indexing="ij")
ref = f(X) * f(Y)
ev = ref + np.random.default_rng(0).normal(0, 0.02 * ref.max(), ref.shape)    # derau 2% dari maksimum
for dta in (1.0, 0.5, 0.2, 0.1):
    g = gamma_index(ref, ev, 1.0, dose_percent=1, dta_mm=dta, step_mm=dta / 10)
    print(f"toleransi dosis 1%, DTA {dta} mm: pass rate {pass_rate(g):.2f}%")