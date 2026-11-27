"""Figure W2 of docs/SEROTONIN_WEAK_DRIVE_2026-11-27.md: the LTD test at mossy-fibre scale 0.15 onto the + cells.
(a) mean W = 1 + T̃ + z over 8 h for + (5-HT gates capture), +→X (DA gates capture) and X→+ (5-HT), with PRPs
(10 min + trigger) and under the PRP block, 5 seeds, GIRK off;
(b) mean tag T̃ of + and +→X against the capture threshold, and the captured z, control.
Writes results/layout_B/fig_ltd_test.png and results/layout_B/ltd_test.json.
Usage: TINYCA3_LAYOUT=B env -u PYTHONPATH .venv/bin/python ltd_fig.py
"""
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from circuit import ROOT  # noqa: E402
from class_stats import SEEDS, final_weight, load  # noqa: E402

RES = ROOT / "results" / "layout_B"
SUF = {"control": "_mod10_gate30_mf0.15", "PRP block": "_mod10_gate30_prpblock_mf0.15",
       "no replay": "_mod10_gate30_norep_mf0.15"}
CLS = {"+": "#e67e22", "+→X": "#16a085", "X→+": "#c0392b"}
plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False})


def series(suf, k, var):
    out = []
    for s in SEEDS:
        d = load(s, suf)
        out.append([1 + a["syn:" + k] + b["syn:" + k] if var == "W" else (a if var == "Tt" else b)["syn:" + k]
                    for a, b in zip(d["slow_Tt"], d["slow_z"])])
    return np.asarray(load(SEEDS[0], suf)["slow_t"]) / 3600.0, np.mean(out, 0)


def main():
    summary = {}
    for c, suf in SUF.items():
        for g in ("", "_girk"):
            rows = [final_weight(load(s, suf + g)) for s in SEEDS]
            summary[f"{c}{', GIRK' if g else ''}"] = {k: dict(mean=round(float(np.mean([r[k] for r in rows])), 4),
                                                       min=round(min(r[k] for r in rows), 4),
                                                       max=round(max(r[k] for r in rows), 4))
                                                   for k in ("X", "+", "X→+", "+→X", "=", "‖", "none")}
    (RES / "ltd_test.json").write_text(json.dumps(summary, indent=1, ensure_ascii=False))

    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.2))
    ax = axes[0]
    for k, col in CLS.items():
        for c, ls in (("control", "-"), ("PRP block", "--")):
            t, w = series(SUF[c], k, "W")
            ax.plot(t, w, ls=ls, color=col, lw=1.6, label=f"{k}, {c}")
    ax.axhline(1, color="0.6", lw=0.6)
    ax.set_xlabel("time from encoding (h)")
    ax.set_ylabel("W = 1 + T̃ + z")
    ax.set_title("(a) scale 0.15, 5 seeds: LTD of + and +→X lasts only with PRPs", fontsize=8)
    ax.legend(frameon=False, fontsize=6, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.17))
    ax = axes[1]
    for k in ("+", "+→X"):
        t, tt = series(SUF["control"], k, "Tt")
        _, z = series(SUF["control"], k, "z")
        ax.plot(t, tt, color=CLS[k], lw=1.6, label=f"tag T̃, {k}")
        ax.plot(t, z, color=CLS[k], lw=1.2, ls=":", label=f"captured z, {k}")
    ax.axhline(-0.2, color="k", lw=0.7, ls="--", label="capture threshold (|T̃| > 0.2)")
    ax.axhline(0, color="0.6", lw=0.6)
    ax.set_xlim(0, 2.5)
    ax.set_xlabel("time from encoding (h)")
    ax.set_title("(b) control: capture stops when |T̃| falls below 0.2", fontsize=8)
    ax.legend(frameon=False, fontsize=6)
    fig.tight_layout()
    fig.savefig(RES / "fig_ltd_test.png", dpi=160)
    print("saved", RES / "fig_ltd_test.png", RES / "ltd_test.json")


if __name__ == "__main__":
    main()
