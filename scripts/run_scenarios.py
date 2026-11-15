"""Run the scenarios fixed in PLAN.md and write results/scenarios.csv, step_convergence.csv and pass_rates.png.

Set GAMMA_QUICK=1 for a short smoke run (a few scenarios only).
"""
import os
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))   # so `gammaqa` imports without installation

import numpy as np
import pandas as pd
from scipy.special import erf

from gammaqa import gamma_index, pass_rate

QUICK = bool(os.environ.get("GAMMA_QUICK"))
OUT = Path(__file__).resolve().parent.parent / "results"
OUT.mkdir(exist_ok=True)
SPACING = 1.0
rng = np.random.default_rng(42)


def field(x, half_width=50.0, sigma=3.0):
    return 0.5 * (erf((x + half_width) / (sigma * np.sqrt(2))) - erf((x - half_width) / (sigma * np.sqrt(2))))


ax = np.arange(-100.0, 100.0 + SPACING, SPACING)
X, Y = np.meshgrid(ax, ax, indexing="ij")
ref = field(X) * field(Y)

scenarios = {"identical": ref}
for s in (1, 2, 4, 6):
    scenarios[f"shift {s} mm"] = field(X - s) * field(Y)
for e in (1, 2, 4, 5):
    scenarios[f"scale +{e}%"] = ref * (1 + e / 100)
for n in (0.5, 1, 2):
    scenarios[f"noise {n}%"] = ref + rng.normal(0.0, n / 100 * ref.max(), ref.shape)
if QUICK:
    scenarios = {k: v for k, v in scenarios.items() if k in ("identical", "shift 1 mm", "scale +4%")}

# (name, dose %, DTA mm, local, search step = DTA / divisor). Local gamma needs a finer search step
# because its tolerance is tiny where the dose is low (see the convergence table).
criteria = [("3%/3mm global", 3, 3, False, 10), ("2%/2mm global", 2, 2, False, 10),
            ("1%/1mm global", 1, 1, False, 10), ("3%/3mm local", 3, 3, True, 20)]

rows = []
for sname, ev in scenarios.items():
    for cname, dose, dta, local, div in criteria:
        t0 = time.time()
        g = gamma_index(ref, ev, SPACING, dose, dta, local=local, step_mm=dta / div)
        rows.append(dict(scenario=sname, criterion=cname, pass_rate=round(pass_rate(g), 2),
                         mean_gamma=round(float(np.nanmean(g)), 3), seconds=round(time.time() - t0, 2)))
        print(rows[-1], flush=True)

res = pd.DataFrame(rows)
res.to_csv(OUT / "scenarios.csv", index=False)
print("\n", res.pivot(index="scenario", columns="criterion", values="pass_rate").reindex(list(scenarios)).to_string())

# Convergence check: pass rate with the search step used above versus half that step
chk = []
for sname in [k for k in ("shift 1 mm", "shift 4 mm", "scale +4%", "noise 1%") if k in scenarios]:
    for cname, dose, dta, local, div in criteria:
        if cname.startswith("1%"):
            continue
        if local and sname not in ("shift 1 mm", "scale +4%"):
            continue                                     # local runs at the finest step are slow
        a = pass_rate(gamma_index(ref, scenarios[sname], SPACING, dose, dta, local=local, step_mm=dta / div))
        b = pass_rate(gamma_index(ref, scenarios[sname], SPACING, dose, dta, local=local, step_mm=dta / (2 * div)))
        chk.append(dict(scenario=sname, criterion=cname, step=f"DTA/{div}", pass_rate=round(a, 2),
                        half_step=f"DTA/{2 * div}", pass_rate_half_step=round(b, 2), difference=round(a - b, 2)))
        print(chk[-1], flush=True)
chk = pd.DataFrame(chk)
chk.to_csv(OUT / "step_convergence.csv", index=False)
print("\nStep convergence:\n", chk.to_string(index=False))

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), sharey=True)
    for a, (title, key, vals) in zip(axes, [("Shift along x (mm)", "shift {} mm", (1, 2, 4, 6)),
                                            ("Uniform scale (%)", "scale +{}%", (1, 2, 4, 5)),
                                            ("Gaussian noise (% of max)", "noise {}%", (0.5, 1, 2))]):
        for cname, *_ in criteria:
            pts = [(v, res[(res.scenario == key.format(v)) & (res.criterion == cname)].pass_rate)
                   for v in vals]
            pts = [(v, p.iloc[0]) for v, p in pts if len(p)]
            if pts:
                a.plot(*zip(*pts), marker="o", label=cname)
        a.set_xlabel(title); a.set_ylim(0, 102); a.grid(alpha=.3)
    axes[0].set_ylabel("Gamma pass rate (%)"); axes[0].legend()
    fig.tight_layout(); fig.savefig(OUT / "pass_rates.png", dpi=150)
    print("figure saved")
except ImportError:
    print("matplotlib not installed: figure skipped")
