"""Weight per P1 synapse class (layout B; docs/architecture_update_plan_2026-11-26.md, Table 3).

Layout B has disjoint basic patterns, so every recurrent synapse falls into exactly one class:
  modulated   X, + (inside a basic pattern), X→+, +→X (tagged only in X∪+)
  controls    =, ‖ (inside),                   =→‖, ‖→= (tagged only in =∪‖)
  none        pairs never presented together: no co-activity (presynaptic firing alone still sets a weaker tag)
Each modulated class has an unmodulated twin with the same sequence of presentations (X/+ ~ =/‖, X→+/+→X ~ =→‖/‖→=).
The twins differ in tag amplitude, so the effect of PRPs alone is a class against itself under --prp-block.
Layout B7 adds the mixed combination X∪‖ and two mixed classes:
  ‖→X   unmodulated presynaptic cell onto an X cell, which has PRPs: neighbour capture is possible
  X→‖   X cell onto a ‖ cell, which has no PRPs: tag only, no capture
In B7 the twins also differ in presentation count: X and ‖ are presented 3 times, + and = twice.

Reads results/layout_B/seeds/pilot_s{seed}{suffix}.pkl (run_pilot.py with TINYCA3_LAYOUT=B). For every condition
found, reports W = mean 1 + Tt + z per class at 0.5, 1, 2, 4 and 8 h (mean and range over seeds), and draws
Figure B1. Writes results/layout_B/class_stats.json and results/layout_B/fig_B1_classes.png.
Usage: TINYCA3_LAYOUT=B .venv/bin/python class_stats.py   (or TINYCA3_LAYOUT=B7: results/layout_B7)
"""
import json
import pickle

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from circuit import BASIC_DISJOINT, LAYOUT, PA, PB, PC, PD, PMIX, RESULTS  # noqa: E402

SEEDS = (1, 2, 3, 4, 5)
COND = {"10 min + trigger": "_mod10_gate30", "PRP block": "_mod10_gate30_prpblock",
        "10 min + trigger, no replay": "_mod10_gate30_norep"}
T_H = (0.5, 1.0, 2.0, 4.0, 8.0)
MOD_CLASSES = (PA, PB, f"{PA}→{PB}", f"{PB}→{PA}")
CTRL_CLASSES = (PC, PD, f"{PC}→{PD}", f"{PD}→{PC}", "none")
# Mixed classes (layout B7 only): the mixed combination PMIX = PA ∪ PD.
MIX_CLASSES = (f"{PD}→{PA}", f"{PA}→{PD}") if PMIX else ()
ALL_CLASSES = MOD_CLASSES + CTRL_CLASSES + MIX_CLASSES
TWIN = {PA: PC, PB: PD, f"{PA}→{PB}": f"{PC}→{PD}", f"{PB}→{PA}": f"{PD}→{PC}"}
# One colour per twin pair: solid = modulated class, dashed = its unmodulated twin, dotted = never co-presented.
PAIR_COL = ("#c0392b", "#2e6fd8", "#e67e22", "#16a085")
CCOL = {**{k: c for k, c in zip(TWIN, PAIR_COL)}, **{v: c for v, c in zip(TWIN.values(), PAIR_COL)}, "none": "#000000",
        f"{PD}→{PA}": "#8e44ad", f"{PA}→{PD}": "#8e44ad"}
# Line style: solid = modulated, dashed = twin, dotted = never co-presented; mixed: dash-dot onto X, loose dots onto ‖.
CLS = {**{k: "-" for k in MOD_CLASSES}, **{k: "--" for k in CTRL_CLASSES}, "none": ":",
       f"{PD}→{PA}": "-.", f"{PA}→{PD}": (0, (1, 3))}
CLABEL = {PA: f"{PA} (DA)", PB: f"{PB} (5-HT)", f"{PA}→{PB}": f"{PA}→{PB} (5-HT)", f"{PB}→{PA}": f"{PB}→{PA} (DA)",
          PC: f"{PC}, twin of {PA}", PD: f"{PD}, twin of {PB}", f"{PC}→{PD}": f"{PC}→{PD}, twin of {PA}→{PB}",
          f"{PD}→{PC}": f"{PD}→{PC}, twin of {PB}→{PA}", "none": "never co-presented",
          f"{PD}→{PA}": f"{PD}→{PA} (mixed, onto DA cells)", f"{PA}→{PD}": f"{PA}→{PD} (mixed, onto unmodulated cells)"}
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
        for k in ALL_CLASSES:
            at = np.array([np.interp(T_H, t, w[k]) for t, w in per_seed])
            fin = np.array([final_weight(d)[k] for d in rs])
            res[k] = dict(mean_at_h={h: round(float(v), 4) for h, v in zip(T_H, at.mean(0))},
                          min_8h=round(float(fin.min()), 4), max_8h=round(float(fin.max()), 4),
                          mean_8h=round(float(fin.mean()), 4))
        res["twin_difference_8h"] = {k: round(res[k]["mean_8h"] - res[v]["mean_8h"], 4) for k, v in TWIN.items()}
        out[c] = dict(seeds=[d["seed"] for d in rs], classes=res)
        print(f"\n{c}: seeds {out[c]['seeds']}")
        for k in ALL_CLASSES:
            r = res[k]
            print(f"  {k:6s} W(1h)={r['mean_at_h'][1.0]:.3f}  W(8h)={r['mean_8h']:.3f} [{r['min_8h']:.3f}, {r['max_8h']:.3f}]")
        print("  modulated - twin control at 8 h:", res["twin_difference_8h"])
    with open(RESULTS / "class_stats.json", "w") as f:
        json.dump(out, f, indent=1, ensure_ascii=False)

    conds = list(runs)
    fig, axes = plt.subplots(len(conds), 1, figsize=(7.5, 3.4 * len(conds)), squeeze=False, sharex=True, sharey=True)
    for row, c in enumerate(conds):
        per_seed = [class_weight(d) for d in runs[c]]
        tg = np.geomspace(per_seed[0][0][0], 8.0, 600)
        ax = axes[row, 0]
        for k in ALL_CLASSES:
            ws = np.array([np.interp(tg, t, w[k]) for t, w in per_seed])
            ax.plot(tg, ws.mean(0), color=CCOL[k], lw=1.6 if k in MOD_CLASSES + MIX_CLASSES else 1.2, ls=CLS[k],
                    label=CLABEL[k])
            if len(ws) > 1:
                ax.fill_between(tg, ws.min(0), ws.max(0), color=CCOL[k], alpha=0.10, lw=0)
        ax.axhline(1.0, color="0.6", lw=0.6)
        ax.set_xscale("log")
        ax.set_xlim(tg[0], 8.0)
        ax.set_title(f"{c}: solid = modulated, dashed = unmodulated twin"
                     f"{', purple = mixed' if PMIX else ''} ({len(per_seed)} seeds, band: range)",
                     fontsize=8)
        ax.set_ylabel("mean 1 + T̃ + z per synapse class")
    axes[-1, 0].set_xlabel("time from the start of encoding (h, log scale)")
    h_, l_ = axes[0, 0].get_legend_handles_labels()
    pos = {k: i for i, k in enumerate(ALL_CLASSES)}     # each modulated class next to its twin, then the rest
    order = [pos[x] for k, v in TWIN.items() for x in (k, v)] + [pos[k] for k in ("none",) + MIX_CLASSES]
    fig.legend([h_[i] for i in order], [l_[i] for i in order], fontsize=7, frameon=False, ncol=5,
               loc="lower center", bbox_to_anchor=(0.5, 0.0))
    fig.tight_layout(rect=(0, (0.07 if not PMIX else 0.10) / len(conds) * 1.5, 1, 1))
    fig.savefig(RESULTS / "fig_B1_classes.png", dpi=160)
    print("saved", RESULTS / "class_stats.json", RESULTS / "fig_B1_classes.png")


if __name__ == "__main__":
    main()
