"""Synthesis trigger (report of 2026-11-24): DA and 5-HT released for minutes switch on a per-cell synthesis state G
that decays with tau_G = 30 min (stc.py), so PRP synthesis outlasts the modulator.

  L&T 60 min          pilot_s{seed}.pkl                 theta_pro follows a 60-min modulator (report 2026-11-17)
  10 min, no trigger  pilot_s{seed}_mod10.pkl           theta_pro follows a 10-min modulator (report section 9)
  10 min + trigger    pilot_s{seed}_mod10_gate30.pkl    run_pilot.py --light --mod-min 10 --gate 30
  5 min + trigger     pilot_s{seed}_mod5_gate30.pkl     run_pilot.py --light --mod-min 5 --gate 30
  PRP block           pilot_s{seed}_prpblock.pkl        no PRPs (report section 8)
Also an offline sweep (no replay, from the end-of-encoding tags of each seed) over tau_G and the release time.
Writes results/trigger_stats.json and results/fig_trig_14_synthesis_trigger.png (Figure 14).
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from circuit import MODULATED, PA, PATTERNS, PATTERN_ORDER, PB, PC, RESULTS, IDX  # noqa: E402
from mod_stats import MODULATED, significant  # noqa: E402
from prp_stats import EPOCHS_H, PCOL, SEEDS, TOL, epoch_counts, w_at, w_final  # noqa: E402
from replay_stats import BURST  # noqa: E402
from stc import ModPulse, STCParams, advance, modulator_at, tag_sum_per_cell  # noqa: E402
from stc_compare import load  # noqa: E402

RES = RESULTS
COND = {"L&T 60 min": "", "10 min, no trigger": "_mod10", "10 min + trigger": "_mod10_gate30",
        "5 min + trigger": "_mod5_gate30", "PRP block": "_prpblock"}
STY = {"L&T 60 min": dict(color="k", ls="-"), "10 min, no trigger": dict(color="#999999", ls="--"),
       "10 min + trigger": dict(color="#d35400", ls="-"), "5 min + trigger": dict(color="#d35400", ls="-."),
       "PRP block": dict(color="#999999", ls=":")}
LAM = [IDX[c] for c in PATTERNS[PA]]
plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False})


def offline_sweep():
    """No replay: integrate from the end-of-encoding tags of each seed for 8 h. Returns z of Λ, V, I and the end of
    synthesis (min) for each (release, tau_G)."""
    out = {}
    for rel, g in ((60, 0), (10, 0), (5, 0), (10, 15), (10, 30), (10, 45), (10, 60), (5, 30), (2, 30)):
        sp = STCParams(tau_G_min=g)
        pulses = [ModPulse("DA", 0.3, 60.0 * rel), ModPulse("5-HT", 2.6, 60.0 * rel)]
        zs, ends = {k: [] for k in (PA, PB, PC)}, []
        for sd in SEEDS:
            d = load(sd)
            t0, Tt0 = d["snapshots"]["end of encoding"]
            post = np.asarray(d["p1_post"])
            last = [0.0]

            def on(t, Tt, z, P, G, post=post, last=last, sp=sp, pulses=pulses):
                gate = G if G is not None else modulator_at(pulses, t)[0]
                if (tag_sum_per_cell(Tt, post) > sp.theta0 / (gate + 0.001)).any():
                    last[0] = t
            _, z, _, _ = advance(Tt0.copy(), np.zeros_like(Tt0), np.zeros(25), post, pulses, t0, 8 * 3600.0, sp,
                                 on_step=on)
            for k in zs:
                zs[k].append(float(z[d["masks"][k]].mean()))
            ends.append(last[0] / 60.0)
        out[f"release {rel} min, tau_G {g} min"] = dict(z={k: round(float(np.mean(v)), 3) for k, v in zs.items()},
                                                         synthesis_ends_min=round(float(np.mean(ends)), 1))
    return out


def main():
    runs = {c: {sd: load(sd, suf) for sd in SEEDS} for c, suf in COND.items()}
    out = {}
    for c in COND:
        wf = {sd: w_final(runs[c][sd]) for sd in SEEDS}
        ec = {sd: epoch_counts(runs[c][sd]) for sd in SEEDS}
        tg = (0.0039, 0.5, 1, 2, 4, 8)
        out[c] = dict(
            W_8h={k: dict(mean=round(float(np.mean([wf[sd][k] for sd in SEEDS])), 4),
                          min=round(float(min(wf[sd][k] for sd in SEEDS)), 4),
                          max=round(float(max(wf[sd][k] for sd in SEEDS)), 4)) for k in PATTERN_ORDER},
            z_8h={k: round(float(np.mean([runs[c][sd]["z_final"][runs[c][sd]["masks"][k]].mean() for sd in SEEDS])), 4)
                  for k in PATTERN_ORDER},
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
    out["offline_sweep"] = offline_sweep()
    (RES / "trigger_stats.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
    print(json.dumps({c: {k: out[c][k] for k in ("z_8h", "P_max", "above_10pct_seeds", "significant_late")}
                      for c in COND}, indent=1, ensure_ascii=False))
    print(json.dumps(out["offline_sweep"], indent=1, ensure_ascii=False))

    # Figure 14
    fig, axes = plt.subplots(2, 2, figsize=(9.6, 6.4))
    ax = axes[0, 0]
    for c in ("L&T 60 min", "10 min, no trigger", "10 min + trigger", "5 min + trigger"):
        d = runs[c][1]
        t = d["slow_t"] / 60.0
        sel = t <= 180
        ax.plot(t[sel], d["slow_P"][sel][:, LAM].mean(1), lw=1.3, label=c, **STY[c])
        if "trigger" in c and "no" not in c:
            ax.plot(t[sel], d["slow_G"][sel][:, LAM].mean(1), lw=0.9, alpha=0.6, color="#2e86c1", ls=STY[c]["ls"],
                    label=f"G, {c.split(',')[0].split(' +')[0]}")
    ax.axvspan(0, 10, color="#c0392b", alpha=0.12, lw=0)
    ax.text(5, 1.02, "DA\n10 min", ha="center", va="top", fontsize=6.5, color="#c0392b")
    ax.set_xlabel("time (min)")
    ax.set_ylabel(f"{PA} cells, mean (seed 1)")
    ax.set_ylim(0, 1.05)
    ax.set_title("PRP pool P (warm) and synthesis trigger G (blue)", fontsize=8)
    ax.legend(frameon=False, fontsize=6, loc="upper right")

    ax = axes[0, 1]
    tg = np.linspace(0.004, 8, 500)
    for c in COND:
        ws = np.array([w_at(runs[c][sd], PA, tg) for sd in SEEDS]).mean(0)
        ax.plot(tg, ws, lw=1.3, label=c, **STY[c])
    for e in EPOCHS_H:
        ax.axvspan(e, e + 0.5, color="#999", alpha=0.12, lw=0)
    ax.axhline(TOL, color="k", lw=0.6, ls=":")
    ax.text(8, TOL, "w0 + 10%", ha="right", va="bottom", fontsize=6.5)
    ax.set_xlabel("time (h)")
    ax.set_ylabel(f"mean 1 + T̃ + z inside {PA}")
    ax.set_title(f"weight of {PA}, mean of 5 seeds", fontsize=8)
    ax.legend(frameon=False, fontsize=6)

    ax = axes[1, 0]
    x = np.arange(len(PATTERN_ORDER))
    wbar = 0.16
    for j, c in enumerate(COND):
        vals = np.array([[w_final(runs[c][sd])[k] for k in PATTERN_ORDER] for sd in SEEDS])
        xs = x + (j - 2) * wbar
        ax.bar(xs, vals.mean(0) - 1, wbar, bottom=1, color=[PCOL[k] for k in PATTERN_ORDER],
               alpha=(1.0, 0.35, 0.85, 0.6, 0.15)[j], edgecolor="k", lw=0.4,
               hatch=(None, "///", None, None, None)[j])
        for i in range(len(SEEDS)):
            ax.plot(xs, vals[i], "k.", ms=1.6)
    ax.axhline(TOL, color="k", lw=0.6, ls=":")
    ax.set_xticks(x, PATTERN_ORDER)
    ax.set_ylim(0.98, 1.34)
    ax.set_ylabel("W_k at 8 h")
    ax.set_title("W_k at 8 h (dots: seeds)", fontsize=8)
    ax.text(0.99, 0.98, "bars, left to right:\n" + "\n".join(COND), transform=ax.transAxes, ha="right", va="top",
            fontsize=6)

    ax = axes[1, 1]
    for c in COND:
        obs = [sum(out[c]["replays_per_epoch"][k][e] for k in MODULATED) for e in range(4)]
        chn = [sum(out[c]["chance_per_epoch"][k][e] for k in MODULATED) for e in range(4)]
        ax.plot(range(4), obs, marker="o", ms=3.5, lw=1.3, label=c, **STY[c])
        ax.plot(range(4), chn, marker="o", ms=2.5, lw=0.8, color="#cccccc", ls=STY[c]["ls"], mfc="white")
    ax.plot([], [], color="#cccccc", marker="o", ms=2.5, mfc="white", label="chance (same line style)")
    ax.set_xticks(range(4), [f"{e:g}–{e + 0.5:g} h" for e in EPOCHS_H])
    ax.set_ylabel(f"replays of {', '.join(MODULATED)}\n(sum over 5 seeds × 20 events)")
    ax.set_title("modulated patterns, per epoch", fontsize=8)
    ax.set_ylim(0, 31)
    ax.legend(frameon=False, fontsize=6, ncol=2, loc="upper right")
    for i, a in enumerate(axes.flat):
        a.text(-0.02, 1.02, f"({chr(97 + i)})", transform=a.transAxes, fontsize=9, fontweight="bold",
               ha="right", va="bottom")
    fig.tight_layout()
    fig.savefig(RES / "fig_trig_14_synthesis_trigger.png", dpi=150)


if __name__ == "__main__":
    main()
