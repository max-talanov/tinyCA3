"""Figure W1 of docs/SEROTONIN_WEAK_DRIVE_2026-11-27.md.
(a) firing of the + cells while + is encoded, against the mossy-fibre scale onto them, GIRK off and on;
(b) tag at the end of encoding for + and X→+ (synapses onto the + cells), same runs (weak_drive_sweep.json, seed 1);
(c) W = 1 + T̃ + z at 8 h for the modulated classes at scale 0.05, GIRK off and on, one line per seed (5 seeds,
    10 min + trigger; results/layout_B/seeds/*_mf0.05{,_girk}.pkl).
Writes results/layout_B/fig_weak_drive.png and results/layout_B/weak_drive_8h.json.
Usage: TINYCA3_LAYOUT=B env -u PYTHONPATH .venv/bin/python weak_drive_fig.py
"""
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from circuit import ROOT  # noqa: E402
from class_stats import SEEDS, final_weight, load  # noqa: E402

RES = ROOT / "results" / "layout_B"
GCOL = {0.0: "#7f8c8d", 6e-5: "#8e44ad"}
GLAB = {0.0: "GIRK off", 6e-5: "GIRK on (6e-5 S/cm²)"}
CLASSES = ("X", "X→+", "+", "+→X")
COND = {"control": "_mod10_gate30", "PRP block": "_mod10_gate30_prpblock", "no replay": "_mod10_gate30_norep"}
plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False})


def main():
    sweep = json.loads((RES / "weak_drive_sweep.json").read_text())
    w8 = {}
    for c, suf in COND.items():
        w8[c] = {}
        for k in CLASSES + ("=", "none"):
            off = [final_weight(load(s, suf + "_mf0.05"))[k] for s in SEEDS]
            on = [final_weight(load(s, suf + "_mf0.05_girk"))[k] for s in SEEDS]
            w8[c][k] = dict(off=[round(x, 4) for x in off], on=[round(x, 4) for x in on],
                            mean_off=round(float(np.mean(off)), 4), mean_on=round(float(np.mean(on)), 4),
                            paired_diff=[round(b - a, 4) for a, b in zip(off, on)])
    (RES / "weak_drive_8h.json").write_text(json.dumps(w8, indent=1, ensure_ascii=False))

    fig, axes = plt.subplots(1, 3, figsize=(10.2, 3.1), gridspec_kw=dict(width_ratios=[1, 1, 1.3]))
    for g in GCOL:
        rs = sorted((r for r in sweep if r["g"] == g), key=lambda r: r["scale"])
        sc = [r["scale"] for r in rs]
        axes[0].plot(sc, [r["rate_Hz"]["+:+"] for r in rs], "o-", color=GCOL[g], ms=3.5, label=GLAB[g])
        axes[1].plot(sc, [r["Tt"]["+"] for r in rs], "o-", color=GCOL[g], ms=3.5, label=f"+ ({GLAB[g]})")
        axes[1].plot(sc, [r["Tt"]["X→+"] for r in rs], "s--", color=GCOL[g], ms=3.5, label=f"X→+ ({GLAB[g]})")
    for ax in axes[:2]:
        ax.set_xscale("log")
        ax.set_xticks([0.05, 0.1, 0.2, 0.5, 1.0], ["0.05", "0.1", "0.2", "0.5", "1"])
        ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.set_xlabel("mossy-fibre weight onto + cells (× calibrated)")
        ax.axvline(0.05, color="0.8", lw=0.8, zorder=0)
    axes[0].set_ylabel("+ cell rate while + is encoded (Hz)")
    axes[0].set_title("(a) firing of the + cells (seed 1)", fontsize=8)
    axes[0].legend(frameon=False, fontsize=6.5)
    for y in (0.2, -0.2):
        axes[1].axhline(y, color="0.6", lw=0.6, ls=":")
    axes[1].axhline(0, color="0.6", lw=0.6)
    axes[1].set_ylabel("tag T̃ at the end of encoding")
    axes[1].set_title("(b) tags onto the + cells (dotted: capture threshold)", fontsize=8)
    axes[1].legend(frameon=False, fontsize=6)
    ax = axes[2]
    for i, k in enumerate(CLASSES):
        d = w8["control"][k]
        for a, b in zip(d["off"], d["on"]):
            ax.plot([i - 0.15, i + 0.15], [a, b], color="0.7", lw=0.7)
        ax.plot([i - 0.15] * len(d["off"]), d["off"], "o", color=GCOL[0.0], ms=3.5, label=GLAB[0.0] if i == 0 else None)
        ax.plot([i + 0.15] * len(d["on"]), d["on"], "o", color=GCOL[6e-5], ms=3.5, label=GLAB[6e-5] if i == 0 else None)
        ax.text(i, max(d["off"] + d["on"]) + 0.03, f"{np.mean(d['paired_diff']):+.3f}", ha="center", fontsize=7)
    ax.set_xticks(range(len(CLASSES)), CLASSES)
    ax.set_ylabel("W = 1 + T̃ + z at 8 h")
    ax.set_title("(c) W at 8 h, scale 0.05, 5 seeds (mean paired change)", fontsize=8)
    ax.set_ylim(top=ax.get_ylim()[1] + 0.05)
    ax.legend(frameon=False, fontsize=6.5, loc="center right")
    fig.tight_layout()
    fig.savefig(RES / "fig_weak_drive.png", dpi=160)
    print("saved", RES / "fig_weak_drive.png", RES / "weak_drive_8h.json")


if __name__ == "__main__":
    main()
