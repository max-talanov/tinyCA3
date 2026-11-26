"""PRP-block control (check 1 of the algorithm doc section 4; T7 (b) of the implementation plan).

The anisomycin analogue: alpha_P = 0, so no PRPs are made, and tags can never be captured. Everything else, including
encoding, DA and 5-HT and the 80 replay events, is as in the calibrated pilot.
  control     pilot_s{seed}.pkl            PRPs on (report of 2026-11-17, sections 5-7)
  PRP block   pilot_s{seed}_prpblock.pkl   run_pilot.py --light --prp-block
Pass criteria (algorithm doc section 4, check 1):
  (1) every pattern weight W_k = mean 1 + Tt + z inside pattern k is within 10% of w0 at 8 h (W_k <= 1.10), in every
      seed;
  (2) replay of every pattern stops: in the late epochs (3-3.5 h, 6-6.5 h) no pattern is read out as replayed above
      chance (summed over seeds), and no single event passes the per-event test p_chance < 0.05/6.
Chance uses the hypergeometric null of replay_stats.py. Writes results/prp_stats.json and
results/fig_pilot_12_prp_block.png (Figure 12, numbering continues the report of 2026-11-17).
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from circuit import MODULATED, P6, PA, PAB, PB, PC, PCOL, PD, RESULTS, PATTERN_ORDER  # noqa: E402
from replay_stats import ALPHA, BURST, p_chance  # noqa: E402
from stc_compare import load, pattern_weight  # noqa: E402

RES = RESULTS
SEEDS = (1, 2, 3, 4, 5)
COND = {"control": "", "PRP block": "_prpblock"}
EPOCHS_H = (0.5, 1.5, 3.0, 6.0)
TOL = 1.10
plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False})


def w_final(d):
    return {k: float(np.mean(1 + d["Tt_final"][d["masks"][k]] + d["z_final"][d["masks"][k]])) for k in PATTERN_ORDER}


def w_at(d, k, t_h):
    """Mean 1 + Tt + z inside pattern k at time(s) t_h (hours), from the slow log."""
    t, w = pattern_weight(d)
    return np.interp(t_h, t, w[k])


def epoch_counts(d):
    """Replays read out per pattern and epoch, and their chance counts (bursts excluded)."""
    obs = {k: [0] * 4 for k in PATTERN_ORDER}
    exp = {k: [0.0] * 4 for k in PATTERN_ORDER}
    for i, e in enumerate(d["events"]):
        if e["n_active"] >= BURST:
            continue
        ep = i // 20
        for k in PATTERN_ORDER:
            obs[k][ep] += int(k in e["replayed"])
            exp[k][ep] += p_chance(k, e["n_active"])
    return obs, exp


def main():
    runs = {c: {sd: load(sd, suf) for sd in SEEDS} for c, suf in COND.items()}
    out = {}
    for c in COND:
        wf = {sd: w_final(runs[c][sd]) for sd in SEEDS}
        ec = {sd: epoch_counts(runs[c][sd]) for sd in SEEDS}
        out[c] = dict(
            W_8h={k: dict(mean=float(np.mean([wf[sd][k] for sd in SEEDS])), max=float(max(wf[sd][k] for sd in SEEDS)),
                          per_seed=[round(wf[sd][k], 4) for sd in SEEDS]) for k in PATTERN_ORDER},
            z_max=float(max(np.abs(runs[c][sd]["z_final"]).max() for sd in SEEDS)),
            P_max=float(max(runs[c][sd]["slow_P"].max() for sd in SEEDS)),
            replays_per_epoch={k: np.sum([ec[sd][0][k] for sd in SEEDS], 0).tolist() for k in PATTERN_ORDER},
            chance_per_epoch={k: np.round(np.sum([ec[sd][1][k] for sd in SEEDS], 0), 2).tolist()
                              for k in PATTERN_ORDER},
            replays_total={k: [int(sum(ec[sd][0][k])) for sd in SEEDS] for k in PATTERN_ORDER},
            chance_total={k: round(float(np.mean([sum(ec[sd][1][k]) for sd in SEEDS])), 2) for k in PATTERN_ORDER},
            median_active=[float(np.median([e["n_active"] for e in runs[c][sd]["events"]])) for sd in SEEDS],
            bursts=[sum(e["n_active"] >= BURST for e in runs[c][sd]["events"]) for sd in SEEDS],
        )
        # pattern weight over time, mean over seeds, at fixed times
        tg = (0.0039, 0.5, 1, 2, 4, 8)
        out[c]["W_t"] = {k: [round(float(np.mean([w_at(runs[c][sd], k, t) for sd in SEEDS])), 4) for t in tg]
                         for k in PATTERN_ORDER}
    blk = out["PRP block"]
    late = {k: dict(obs=sum(blk["replays_per_epoch"][k][2:]), chance=round(sum(blk["chance_per_epoch"][k][2:]), 2))
            for k in PATTERN_ORDER}
    def late_significant(c):
        n = {k: 0 for k in PATTERN_ORDER}
        for sd in SEEDS:
            for i, e in enumerate(runs[c][sd]["events"]):
                if i >= 40 and e["n_active"] < BURST:
                    for k in e["replayed"]:
                        n[k] += int(p_chance(k, e["n_active"]) < ALPHA)
        return n

    out["checks"] = dict(
        late_significant_events={c: late_significant(c) for c in COND},
        weights_within_10pct=all(blk["W_8h"][k]["max"] <= TOL for k in PATTERN_ORDER),
        late_replays_vs_chance=late,
        replay_stops=all(late[k]["obs"] <= late[k]["chance"] for k in PATTERN_ORDER),
    )
    (RES / "prp_stats.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
    print(json.dumps(out, indent=1, ensure_ascii=False))

    # Figure 12
    fig, axes = plt.subplots(1, 3, figsize=(9.6, 3.2), gridspec_kw=dict(width_ratios=[1.3, 1, 1]))
    ax = axes[0]
    tg = np.linspace(0.004, 8, 500)
    for k in (PA, PB, PC, PD):
        for c, ls in (("control", "-"), ("PRP block", "--")):
            ws = np.array([w_at(runs[c][sd], k, tg) for sd in SEEDS])
            ax.plot(tg, ws.mean(0), color=PCOL[k], ls=ls, lw=1.2 if c == "control" else 1.4,
                    label=f"{k}" if c == "control" else None)
            if c == "PRP block":
                ax.fill_between(tg, ws.min(0), ws.max(0), color=PCOL[k], alpha=0.15, lw=0)
    for e in EPOCHS_H:
        ax.axvspan(e, e + 0.5, color="#999", alpha=0.12, lw=0)
    ax.axhline(TOL, color="k", lw=0.6, ls=":")
    ax.text(8, TOL, "w0 + 10%", ha="right", va="bottom", fontsize=6.5)
    ax.plot([], [], "k-", lw=1.2, label="control")
    ax.plot([], [], "k--", lw=1.4, label="PRP block")
    ax.set_xlabel("time (h)")
    ax.set_ylabel("mean 1 + T̃ + z inside the pattern")
    ax.set_title("pattern weight, 5 seeds (band: range, PRP block)", fontsize=8)
    ax.legend(frameon=False, fontsize=6.5, ncol=2)

    ax = axes[1]
    x = np.arange(len(PATTERN_ORDER))
    for j, (c, hatch) in enumerate((("control", None), ("PRP block", "///"))):
        vals = np.array([[out[c]["W_8h"][k]["per_seed"][i] for k in PATTERN_ORDER] for i in range(len(SEEDS))])
        ax.bar(x + (j - 0.5) * 0.38, vals.mean(0) - 1, 0.38, bottom=1, color=[PCOL[k] for k in PATTERN_ORDER],
               alpha=1 if c == "control" else 0.45, hatch=hatch, edgecolor="k", lw=0.4, label=c)
        for i in range(len(SEEDS)):
            ax.plot(x + (j - 0.5) * 0.38 + (i - 2) * 0.05, vals[i], "k.", ms=2.5)
    ax.axhline(TOL, color="k", lw=0.6, ls=":")
    ax.set_xticks(x, PATTERN_ORDER)
    ax.set_ylim(0.98, 1.32)
    ax.set_ylabel("W_k at 8 h")
    ax.set_title("weight at 8 h (dots: seeds)", fontsize=8)
    ax.legend(frameon=False, fontsize=6.5)

    ax = axes[2]
    for c, mk, ls in (("control", "o", "-"), ("PRP block", "s", "--")):
        cap = [sum(out[c]["replays_per_epoch"][k][e] for k in MODULATED) for e in range(4)]
        chn = [sum(out[c]["chance_per_epoch"][k][e] for k in MODULATED) for e in range(4)]
        ax.plot(range(4), cap, marker=mk, ls=ls, color="#8e44ad", label=f"{c}: observed")
        ax.plot(range(4), chn, marker=mk, ls=ls, color="#bbbbbb", mfc="white", label=f"{c}: chance")
    ax.set_xticks(range(4), [f"{e:g}–{e + 0.5:g} h" for e in EPOCHS_H])
    ax.set_ylabel(f"replays of {', '.join(MODULATED)}\n(sum over 5 seeds × 20 events)")
    ax.set_title("modulated patterns, per epoch", fontsize=8)
    ax.legend(frameon=False, fontsize=6.5)
    for i, a in enumerate(axes):
        a.text(-0.02, 1.02, f"({chr(97 + i)})", transform=a.transAxes, fontsize=9, fontweight="bold",
               ha="right", va="bottom")
    fig.tight_layout()
    fig.savefig(RES / "fig_pilot_12_prp_block.png", dpi=150)


if __name__ == "__main__":
    main()
