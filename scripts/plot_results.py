"""Redraw results/pass_rates.png from results/scenarios.csv (no recomputation).

Lines that coincide (for example 3%/3mm global and local) stay visible through different
line styles and marker sizes.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

OUT = Path(__file__).resolve().parent.parent / "results"
res = pd.read_csv(OUT / "scenarios.csv")

STYLE = {"3%/3mm global": dict(ls="-", marker="o", ms=12, mfc="none", color="tab:blue"),
         "2%/2mm global": dict(ls="--", marker="s", ms=9, mfc="none", color="tab:orange"),
         "1%/1mm global": dict(ls="-.", marker="^", ms=7, color="tab:green"),
         "3%/3mm local": dict(ls=":", marker="D", ms=5, color="tab:red")}
PANELS = [("Shift along x (mm)", "shift {} mm", (1, 2, 4, 6)),
          ("Uniform scale (%)", "scale +{}%", (1, 2, 4, 5)),
          ("Gaussian noise (% of maximum dose)", "noise {}%", (0.5, 1, 2))]

fig, axes = plt.subplots(1, 3, figsize=(15, 4.8), sharey=True)
for a, (title, key, vals) in zip(axes, PANELS):
    for crit, st in STYLE.items():
        ys = []
        for v in vals:
            r = res[(res.scenario == key.format(v)) & (res.criterion == crit)]
            ys.append(r.pass_rate.iloc[0] if len(r) else float("nan"))
        a.plot(vals, ys, lw=2, label=crit, **st)
    a.set_xlabel(title)
    a.set_xticks(vals)
    a.set_ylim(0, 102)
    a.grid(alpha=0.3)
axes[0].set_ylabel("Gamma pass rate (%)")
axes[0].legend(loc="lower left")
fig.suptitle("Gamma pass rate for controlled errors in a synthetic 2D field", y=1.0)
fig.tight_layout()
fig.savefig(OUT / "pass_rates.png", dpi=150)
print("saved", OUT / "pass_rates.png")
