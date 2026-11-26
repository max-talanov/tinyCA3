"""Weight per P1 synapse class (layout B; docs/architecture_update_plan_2026-11-26.md, Table 3).

Layout B has disjoint basic patterns, so every recurrent synapse falls into exactly one class:
  modulated   X, + (inside a basic pattern), X→+, +→X (tagged only in X∪+)
  controls    =, ‖ (inside),                   =→‖, ‖→= (tagged only in =∪‖)
  none        pairs never presented together: no co-activity (postsynaptic firing alone still sets a weaker tag)
Each modulated class has an unmodulated twin with the same tagging history (X/+ ~ =/‖, X→+/+→X ~ =→‖/‖→=), so the
difference between twins is the effect of PRP availability alone.

Reads results/layout_B/seeds/pilot_s{seed}{suffix}.pkl (run_pilot.py with TINYCA3_LAYOUT=B). For every condition
found, reports W = mean 1 + Tt + z per class at 0.5, 1, 2, 4 and 8 h (mean and range over seeds), and draws
Figure B1. Writes results/layout_B/class_stats.json and results/layout_B/fig_B1_classes.png.
Usage: TINYCA3_LAYOUT=B .venv/bin/python class_stats.py
"""
import json
import pickle

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from circuit import BASIC_DISJOINT, LAYOUT, PA, PB, PC, PCOL, PD, RESULTS  # noqa: E402

SEEDS = (1, 2, 3, 4, 5)
COND = {"10 min + trigger": "_mod10_gate30", "PRP block": "_mod10_gate30_prpblock",
        "10 min + trigger, no replay": "_mod10_gate30_norep"}
T_H = (0.5, 1.0, 2.0, 4.0, 8.0)
MOD_CLASSES = (PA, PB, f"{PA}→{PB}", f"{PB}→{PA}")
CTRL_CLASSES = (PC, PD, f"{PC}→{PD}", f"{PD}→{PC}", "none")
TWIN = {PA: PC, PB: PD, f"{PA}→{PB}": f"{PC}→{PD}", f"{PB}→{PA}": f"{PD}→{PC}"}
CCOL = {PA: PCOL[PA], PB: PCOL[PB], f"{PA}→{PB}": "#e67e22", f"{PB}→{PA}": "#16a085",
        PC: PCOL[PC], PD: PCOL[PD], f"{PC}→{PD}": "#7d6608", f"{PD}→{PC}": "#5d6d7e", "none": "#000000"}
plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False})


def load(seed, suffix):
    path = RESULTS / "seeds" / f"pilot_s{seed}{suffix}.pkl"
    if not path.exists():
        return None
    with open(path, "rb") as f:
        return pickle.load(f)


def class_weight(d):
    """t (h) and mean 1 + Tt + z per class from the slow log (keys 'syn:<class>' written by run_pilot.py)."""
    t = np.asarray(d["slow_t"]) / 3600.0
    keys = [k for k in d["slow_Tt"][0] if k.startswith("syn:")]
    return t, {k[4:]: 1 + np.array([a[k] + b[k] for a, b in zip(d["slow_Tt"], d["slow_z"])]) for k in keys}


def final_weight(d):
    """Mean 1 + Tt + z per class at the end of the run, from the final state."""
    return {k[4:]: float(np.mean(1 + d["Tt_final"][m] + d["z_final"][m])) for k, m in d["masks"].items()
            if k.startswith("syn:")}


def main():
    if not BASIC_DISJOINT:
        raise SystemExit(f"layout {LAYOUT}: synapse classes need disjoint basic patterns (use TINYCA3_LAYOUT=B)")
    runs = {c: [d for d in (load(sd, suf) for sd in SEEDS) if d is not None] for c, suf in COND.items()}
    runs = {c: r for c, r in runs.items() if r}
    if not runs:
        raise SystemExit(f"no runs found in {RESULTS / 'seeds'}")
    out = {}
    for c, rs in runs.items():
        per_seed = [class_weight(d) for d in rs]
        res = {}
        for k in MOD_CLASSES + CTRL_CLASSES:
            at = np.array([np.interp(T_H, t, w[k]) for t, w in per_seed])
            fin = np.array([final_weight(d)[k] for d in rs])
            res[k] = dict(mean_at_h={h: round(float(v), 4) for h, v in zip(T_H, at.mean(0))},
                          min_8h=round(float(fin.min()), 4), max_8h=round(float(fin.max()), 4),
                          mean_8h=round(float(fin.mean()), 4))
        res["twin_difference_8h"] = {k: round(res[k]["mean_8h"] - res[v]["mean_8h"], 4) for k, v in TWIN.items()}
        out[c] = dict(seeds=[d["seed"] for d in rs], classes=res)
        print(f"\n{c}: seeds {out[c]['seeds']}")
        for k in MOD_CLASSES + CTRL_CLASSES:
            r = res[k]
            print(f"  {k:6s} W(1h)={r['mean_at_h'][1.0]:.3f}  W(8h)={r['mean_8h']:.3f} [{r['min_8h']:.3f}, {r['max_8h']:.3f}]")
        print("  modulated - twin control at 8 h:", res["twin_difference_8h"])
    with open(RESULTS / "class_stats.json", "w") as f:
        json.dump(out, f, indent=1, ensure_ascii=False)

    conds = list(runs)
    fig, axes = plt.subplots(len(conds), 2, figsize=(9, 3.0 * len(conds)), squeeze=False, sharey=True)
    for row, c in enumerate(conds):
        per_seed = [class_weight(d) for d in runs[c]]
        tg = np.linspace(per_seed[0][0][0], 8.0, 400)
        for col, (title, classes) in enumerate((("modulated (DA → X, 5-HT → +)", MOD_CLASSES),
                                                ("unmodulated controls; none = never co-presented", CTRL_CLASSES))):
            ax = axes[row, col]
            for k in classes:
                ws = np.array([np.interp(tg, t, w[k]) for t, w in per_seed])
                ls = "--" if "→" in k else (":" if k == "none" else "-")
                ax.plot(tg, ws.mean(0), color=CCOL[k], lw=1.4, ls=ls, label=k)
                if len(ws) > 1:
                    ax.fill_between(tg, ws.min(0), ws.max(0), color=CCOL[k], alpha=0.12, lw=0)
            ax.axhline(1.0, color="0.6", lw=0.6)
            ax.set_xscale("log")
            ax.set_xlim(tg[0], 8.0)
            ax.set_title(f"{c}: {title}", fontsize=8)
            ax.set_xlabel("time (h)")
            if col == 0:
                ax.set_ylabel(f"mean 1 + T̃ + z per class\n(mean of {len(per_seed)} seeds, band: range)")
            ax.legend(fontsize=6.5, frameon=False, loc="upper right")
    fig.tight_layout()
    fig.savefig(RESULTS / "fig_B1_classes.png", dpi=160)
    print("saved", RESULTS / "class_stats.json", RESULTS / "fig_B1_classes.png")


if __name__ == "__main__":
    main()
