"""Figure D1 of docs/DA_VS_5HT_2026-11-26.md: mean W = 1 + T̃ + z inside X (DA) and + (5-HT) over time, control and
PRP block (layout B, 10 min + trigger, 5 seeds). Reads results/layout_B/class_stats.json (class_stats.py) and writes
results/layout_B/fig_da_vs_5ht.png.
Usage: env -u PYTHONPATH .venv/bin/python da_vs_5ht_fig.py
"""
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from circuit import ROOT  # noqa: E402

RES = ROOT / "results" / "layout_B"
COND = {"10 min + trigger": "-", "PRP block": "--"}
MOD = {"X": ("DA → X cells", "#2e6fd8", 2.6, 5), "+": ("5-HT → + cells", "#e67e22", 1.4, 3)}
plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False})


def main():
    d = json.loads((RES / "class_stats.json").read_text())
    fig, ax = plt.subplots(figsize=(6.4, 3.4))
    gap = 0.0
    for c, ls in COND.items():
        cls = d[c]["classes"]
        for k, (lab, col, lw, ms) in MOD.items():
            t = [float(h) for h in cls[k]["mean_at_h"]]
            w = list(cls[k]["mean_at_h"].values())
            ax.plot(t, w, ls=ls, color=col, lw=lw, marker="o", ms=ms, label=lab if c == "10 min + trigger" else None)
        a, b = cls["X"]["mean_at_h"], cls["+"]["mean_at_h"]
        gap = max(gap, max(abs(a[h] - b[h]) for h in a))
        ax.text(8.6, (cls["X"]["mean_8h"] + cls["+"]["mean_8h"]) / 2,
                f"{c}:\nX {cls['X']['mean_8h']:.3f}, + {cls['+']['mean_8h']:.3f}", va="center", fontsize=7)
    ax.set_xscale("log")
    ax.set_xticks([0.5, 1, 2, 4, 8], ["0.5", "1", "2", "4", "8"])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xlim(0.45, 8.5)
    ax.set_xlabel("time from encoding (h)")
    ax.set_ylabel("W = 1 + T̃ + z inside the pattern")
    ax.set_title(f"DA on X and 5-HT on + give the same weight curve (largest gap {gap:.3f}; "
                 "solid: control, dashed: PRP block)", fontsize=8)
    ax.legend(frameon=False, fontsize=7, loc="upper right")
    fig.tight_layout()
    fig.savefig(RES / "fig_da_vs_5ht.png", dpi=160)
    print("saved", RES / "fig_da_vs_5ht.png", "gap", round(gap, 4))


if __name__ == "__main__":
    main()
