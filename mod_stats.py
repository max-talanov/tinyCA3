"""Neuromodulator duration: DA (Λ cells) and 5-HT (V cells) on for 60, 10 or 5 min (docs/PROTOCOL.md section 3).

  60 min      pilot_s{seed}.pkl             report sections 5-7 (Lehr et al. 2022)
  10 min      pilot_s{seed}_mod10.pkl       run_pilot.py --light --mod-min 10
  5 min       pilot_s{seed}_mod5.pkl        run_pilot.py --light --mod-min 5
  PRP block   pilot_s{seed}_prpblock.pkl    the floor: no capture at all (report section 8)
Everything else is the calibrated pilot. Writes results/mod_stats.json and results/fig_pilot_13_mod_duration.png
(Figure 13, numbering continues the report of 2026-11-17).
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from circuit import PATTERN_ORDER  # noqa: E402
from prp_stats import EPOCHS_H, PCOL, SEEDS, TOL, epoch_counts, w_at, w_final  # noqa: E402
from replay_stats import ALPHA, BURST, p_chance  # noqa: E402
from stc_compare import load  # noqa: E402

RES = Path(__file__).resolve().parent / "results"
COND = {"60 min": "", "10 min": "_mod10", "5 min": "_mod5", "PRP block": "_prpblock"}
LS = {"60 min": "-", "10 min": "--", "5 min": "-.", "PRP block": ":"}
MODULATED = ("Λ", "V", "Λ∪V", "Λ∪V∪I")
plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False})


def significant(d, first=0):
    """Events from index `first` on in which a pattern is replayed with p_chance < ALPHA, per pattern."""
    n = {k: 0 for k in PATTERN_ORDER}
    for e in d["events"][first:]:
        if e["n_active"] < BURST:
            for k in e["replayed"]:
                n[k] += int(p_chance(k, e["n_active"]) < ALPHA)
    return n


def main():
    runs = {c: {sd: load(sd, suf) for sd in SEEDS} for c, suf in COND.items()}
    out = {}
    for c in COND:
        wf = {sd: w_final(runs[c][sd]) for sd in SEEDS}
        ec = {sd: epoch_counts(runs[c][sd]) for sd in SEEDS}
        z_in = {k: [float(runs[c][sd]["z_final"][runs[c][sd]["masks"][k]].mean()) for sd in SEEDS]
                for k in PATTERN_ORDER}
        tg = (0.0039, 0.5, 1, 2, 4, 8)
        out[c] = dict(
            W_8h={k: dict(mean=round(float(np.mean([wf[sd][k] for sd in SEEDS])), 4),
                          min=round(float(min(wf[sd][k] for sd in SEEDS)), 4),
                          max=round(float(max(wf[sd][k] for sd in SEEDS)), 4)) for k in PATTERN_ORDER},
            z_8h={k: round(float(np.mean(v)), 4) for k, v in z_in.items()},
            P_max=round(float(max(runs[c][sd]["slow_P"].max() for sd in SEEDS)), 4),
            W_t={k: [round(float(np.mean([w_at(runs[c][sd], k, t) for sd in SEEDS])), 4) for t in tg]
                 for k in PATTERN_ORDER},
            replays_per_epoch={k: np.sum([ec[sd][0][k] for sd in SEEDS], 0).tolist() for k in PATTERN_ORDER},
            chance_per_epoch={k: np.round(np.sum([ec[sd][1][k] for sd in SEEDS], 0), 2).tolist()
                              for k in PATTERN_ORDER},
            significant_events={k: int(sum(significant(runs[c][sd])[k] for sd in SEEDS)) for k in PATTERN_ORDER},
            significant_late={k: int(sum(significant(runs[c][sd], 40)[k] for sd in SEEDS)) for k in PATTERN_ORDER},
            above_10pct_seeds={k: int(sum(wf[sd][k] > TOL for sd in SEEDS)) for k in PATTERN_ORDER},
            bursts=[int(sum(e["n_active"] >= BURST for e in runs[c][sd]["events"])) for sd in SEEDS],
        )
    (RES / "mod_stats.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
    print(json.dumps({c: {k: out[c][k] for k in ("W_8h", "z_8h", "P_max", "above_10pct_seeds", "significant_late")}
                      for c in COND}, indent=1, ensure_ascii=False))

    # Figure 13
    fig, axes = plt.subplots(1, 3, figsize=(9.6, 3.2), gridspec_kw=dict(width_ratios=[1.3, 1.1, 1]))
    ax = axes[0]
    tg = np.linspace(0.004, 8, 500)
    for c in COND:
        for k in ("Λ", "V"):
            ws = np.array([w_at(runs[c][sd], k, tg) for sd in SEEDS]).mean(0)
            ax.plot(tg, ws, color=PCOL[k], ls=LS[c], lw=1.2, label=c if k == "Λ" else None)
    for e in EPOCHS_H:
        ax.axvspan(e, e + 0.5, color="#999", alpha=0.12, lw=0)
    ax.axhline(TOL, color="k", lw=0.6, ls=":")
    ax.text(8, TOL, "w0 + 10%", ha="right", va="bottom", fontsize=6.5)
    ax.set_xlabel("time (h)")
    ax.set_ylabel("mean 1 + T̃ + z inside the pattern")
    ax.set_title("Λ (red) and V (blue), mean of 5 seeds", fontsize=8)
    ax.legend(frameon=False, fontsize=6.5, title="DA / 5-HT on for", title_fontsize=6.5)

    ax = axes[1]
    x = np.arange(len(PATTERN_ORDER))
    wbar = 0.2
    for j, c in enumerate(COND):
        vals = np.array([[w_final(runs[c][sd])[k] for k in PATTERN_ORDER] for sd in SEEDS])
        xs = x + (j - 1.5) * wbar
        ax.bar(xs, vals.mean(0) - 1, wbar, bottom=1, color=[PCOL[k] for k in PATTERN_ORDER],
               alpha=(1.0, 0.7, 0.45, 0.2)[j], edgecolor="k", lw=0.4, label=c)
        for i in range(len(SEEDS)):
            ax.plot(xs, vals[i], "k.", ms=1.8)
    ax.axhline(TOL, color="k", lw=0.6, ls=":")
    ax.set_xticks(x, PATTERN_ORDER)
    ax.set_ylim(0.98, 1.32)
    ax.set_ylabel("W_k at 8 h")
    ax.set_title("W_k at 8 h: 60, 10, 5 min, PRP block", fontsize=8)

    ax = axes[2]
    for c in COND:
        obs = [sum(out[c]["replays_per_epoch"][k][e] for k in MODULATED) for e in range(4)]
        chn = [sum(out[c]["chance_per_epoch"][k][e] for k in MODULATED) for e in range(4)]
        ax.plot(range(4), obs, marker="o", ms=3.5, ls=LS[c], color="#8e44ad", label=c)
        ax.plot(range(4), chn, marker="o", ms=3, ls=LS[c], color="#bbbbbb", mfc="white")
    ax.plot([], [], color="#bbbbbb", marker="o", ms=3, mfc="white", label="chance (same style)")
    ax.set_xticks(range(4), [f"{e:g}–{e + 0.5:g} h" for e in EPOCHS_H])
    ax.set_ylabel("replays of Λ, V, Λ∪V, Λ∪V∪I\n(sum over 5 seeds × 20 events)")
    ax.set_title("modulated patterns, per epoch", fontsize=8)
    ax.legend(frameon=False, fontsize=6.5)
    for i, a in enumerate(axes):
        a.text(-0.02, 1.02, f"({chr(97 + i)})", transform=a.transAxes, fontsize=9, fontweight="bold",
               ha="right", va="bottom")
    fig.tight_layout()
    fig.savefig(RES / "fig_pilot_13_mod_duration.png", dpi=150)


if __name__ == "__main__":
    main()
