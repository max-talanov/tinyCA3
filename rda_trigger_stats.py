"""Replay DA with the synthesis trigger (report of 2026-11-24, section 5): sections 6 and 7 of the report of
2026-11-17 rerun on the bio-plausible protocol (--mod-min 10 --gate 30), next to the L&T runs.

  section 6 (VTA mode)            L&T suffix          trigger suffix
    no replay                     _norep              _mod10_gate30_norep
    replay                        ""                  _mod10_gate30
    replay + DA                   _rda10              _mod10_gate30_rda10
  section 7 (X4, no encoding DA for Λ)
    no DA                         _noencda            _mod10_gate30_noencda
    replay DA                     _noencda_rda10      _mod10_gate30_noencda_rda10
    encoding DA (reference)       ""                  _mod10_gate30
Writes results/rda_trigger_stats.json and results/fig_trig_15_replay_da.png (Figure 15).
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from circuit import IDX, PA, PATTERN_ORDER, RESULTS, PATTERNS  # noqa: E402
from prp_stats import PCOL, SEEDS  # noqa: E402
from replay_stats import BURST  # noqa: E402
from stc_compare import load  # noqa: E402

RES = RESULTS
MODEL = {"L&T": "", "trigger": "_mod10_gate30"}
S6 = {"no replay": "_norep", "replay": "", "replay + DA": "_rda10"}
S7 = {"no DA": "_noencda", "replay DA": "_noencda_rda10", "encoding DA": ""}
LAM = [IDX[c] for c in PATTERNS[PA]]
plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False})


def zmean(d, k):
    return float(d["z_final"][d["masks"][k]].mean())


def block(conds, pre):
    runs = {c: {sd: load(sd, pre + suf) for sd in SEEDS} for c, suf in conds.items()}
    z = {c: {k: [zmean(runs[c][sd], k) for sd in SEEDS] for k in PATTERN_ORDER} for c in conds}
    return runs, z


def main():
    out, R = {}, {}
    for m, pre in MODEL.items():
        r6, z6 = block(S6, pre)
        r7, z7 = block(S7, pre)
        R[m] = (r6, r7)
        da6 = {sd: [round(e["t_s"] / 3600, 2) for e in r6["replay + DA"][sd]["events"] if e.get("da")] for sd in SEEDS}
        da7 = {sd: [round(e["t_s"] / 3600, 2) for e in r7["replay DA"][sd]["events"] if e.get("da")] for sd in SEEDS}
        lam7 = {c: [sum(1 for e in r7[c][sd]["events"] if PA in e["replayed"] and e["n_active"] < BURST)
                    for sd in SEEDS] for c in S7}
        s6, s7 = {}, {}
        for k in PATTERN_ORDER:
            d = np.array(z6["replay + DA"][k]) - np.array(z6["replay"][k])
            s6[k] = dict(**{f"z_{c}": round(float(np.mean(z6[c][k])), 4) for c in S6},
                         d_da=round(float(d.mean()), 4), d_da_sd=round(float(d.std(ddof=1)), 4),
                         seeds_pos=int(np.sum(d > 1e-4)))
            d = np.array(z7["replay DA"][k]) - np.array(z7["no DA"][k])
            s7[k] = dict(**{f"z_{c}": round(float(np.mean(z7[c][k])), 4) for c in S7},
                         d=round(float(d.mean()), 4), d_sd=round(float(d.std(ddof=1)), 4),
                         seeds_pos=int(np.sum(d > 1e-4)))
        out[m] = dict(
            section6=dict(summary=s6, da_events_h=da6, n_da=[len(da6[sd]) for sd in SEEDS]),
            section7=dict(summary=s7, da_events_h=da7, lam_replays=lam7,
                          z_lam_per_seed={c: [round(v, 4) for v in z7[c][PA]] for c in S7},
                          P_lam_peak={c: [round(float(np.asarray(r7[c][sd]["slow_P"])[:, LAM].mean(1).max()), 3)
                                          for sd in SEEDS] for c in ("no DA", "replay DA")}))
    (RES / "rda_trigger_stats.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
    for m in MODEL:
        print("==", m)
        print(" s6 DA events per seed", out[m]["section6"]["n_da"])
        for k in PATTERN_ORDER:
            print("  s6", k, out[m]["section6"]["summary"][k])
        print(" s7 DA events", out[m]["section7"]["da_events_h"])
        print(f" s7 {PA} replays", out[m]["section7"]["lam_replays"])
        print(f" s7 z {PA} per seed", out[m]["section7"]["z_lam_per_seed"])
        print(f" s7 P peak {PA}", out[m]["section7"]["P_lam_peak"])
        for k in PATTERN_ORDER:
            print("  s7", k, out[m]["section7"]["summary"][k])

    # Figure 15
    fig, axes = plt.subplots(2, 2, figsize=(9.6, 6.2))
    x = np.arange(len(PATTERN_ORDER))
    ax = axes[0, 0]
    for j, m in enumerate(MODEL):
        s = out[m]["section6"]["summary"]
        ax.bar(x + (j - 0.5) * 0.36, [s[k]["d_da"] for k in PATTERN_ORDER], 0.36,
               yerr=[s[k]["d_da_sd"] for k in PATTERN_ORDER], capsize=2, color=[PCOL[k] for k in PATTERN_ORDER],
               alpha=(0.35, 1.0)[j], edgecolor="k", lw=0.4, label=f"{m}")
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xticks(x, PATTERN_ORDER)
    ax.set_ylabel("Δz at 8 h (replay + DA − replay)")
    ax.set_title("section 6: DA during replay, on top of encoding DA", fontsize=8)
    ax.legend(frameon=False, fontsize=6.5, title="pale: L&T 60 min; solid: 10 min + trigger", title_fontsize=6.5)

    ax = axes[0, 1]
    w = 0.13
    for j, m in enumerate(MODEL):
        for i, c in enumerate(S7):
            v = out[m]["section7"]["z_lam_per_seed"][c]
            xi = (j * 3 + i) * 1.0
            ax.bar(xi, np.mean(v), 0.8, color=("#999999", "#c0392b", "#e9a3a3")[i], alpha=(0.45, 1.0)[j],
                   edgecolor="k", lw=0.4)
            ax.plot(xi + np.linspace(-0.25, 0.25, len(v)), v, "k.", ms=3)
    ax.set_xticks(range(6), [f"{c}\n{m}" for m in MODEL for c in S7], fontsize=6.5)
    ax.set_ylabel(f"z of {PA} at 8 h")
    ax.set_title(f"section 7 (X4): {PA} with DA only from replay", fontsize=8)

    ax = axes[1, 0]
    r7 = R["trigger"][1]
    for c, col in (("no DA", "#999999"), ("replay DA", "#c0392b")):
        for sd in SEEDS:
            d = r7[c][sd]
            ax.plot(np.asarray(d["slow_t"]) / 3600, [v[PA] for v in d["slow_z"]], color=col, lw=0.8,
                    label=c if sd == 1 else None)
    for sd in SEEDS:
        for t in out["trigger"]["section7"]["da_events_h"][sd]:
            ax.axvline(t, color="#c0392b", lw=0.4, alpha=0.4)
    ax.set_xlabel("time after encoding (h)")
    ax.set_ylabel(f"mean z of {PA}")
    ax.set_title(f"X4 with trigger: capture of {PA} (lines: replay DA events)", fontsize=8)
    ax.legend(frameon=False, fontsize=6.5)

    ax = axes[1, 1]
    for sd in SEEDS:
        d = r7["replay DA"][sd]
        t = np.asarray(d["slow_t"]) / 3600
        n = len(out["trigger"]["section7"]["da_events_h"][sd])
        ln, = ax.plot(t, np.asarray(d["slow_P"])[:, LAM].mean(1), lw=0.9, label=f"seed {sd} ({n} DA)")
        ax.plot(t, np.asarray(d["slow_G"])[:, LAM].mean(1), lw=0.6, ls=":", color=ln.get_color())
    ax.set_xlabel("time after encoding (h)")
    ax.set_ylabel(f"{PA} cells, mean")
    ax.set_title("X4 with trigger: PRP pool P (solid) and trigger G (dotted)", fontsize=8)
    ax.legend(frameon=False, fontsize=6)
    for i, a in enumerate(axes.flat):
        a.text(-0.02, 1.02, f"({chr(97 + i)})", transform=a.transAxes, fontsize=9, fontweight="bold",
               ha="right", va="bottom")
    fig.tight_layout()
    fig.savefig(RES / "fig_trig_15_replay_da.png", dpi=150)


if __name__ == "__main__":
    main()
