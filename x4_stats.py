"""Experiment X4: Λ encoded without DA, so that its only DA comes from replay.

Conditions per seed (results/seeds/):
  no DA        pilot_s{seed}_noencda.pkl         no encoding DA, replay, no replay DA (5-HT on V unchanged)
  replay DA    pilot_s{seed}_noencda_rda10.pkl   no encoding DA, every event that replays Λ releases DA onto the
                                                 Λ cells for 10 min
  reference    pilot_s{seed}.pkl                 encoding DA (section 5)
Writes results/x4_stats.json and results/fig_pilot_9_x4.png (Figure 9).
"""
import json
import pickle
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from circuit import IDX, PATTERN_ORDER, PATTERNS  # noqa: E402

RES = Path(__file__).resolve().parent / "results"
SEEDS = (1, 2, 3, 4, 5)
COND = {"no DA": "_noencda", "replay DA": "_noencda_rda10", "encoding DA": ""}
CCOL = {"no DA": "#7f7f7f", "replay DA": "#c0392b", "encoding DA": "#e9a3a3"}


def load(seed, suffix):
    with open(RES / "seeds" / f"pilot_s{seed}{suffix}.pkl", "rb") as f:
        return pickle.load(f)


def main():
    data = {c: {sd: load(sd, suf) for sd in SEEDS} for c, suf in COND.items()}
    lam = [IDX[c] for c in PATTERNS["Λ"]]
    z = {c: {sd: {k: float(np.mean(d["z_final"][d["masks"][k]])) for k in PATTERN_ORDER}
             for sd, d in dd.items()} for c, dd in data.items()}
    da_events = {sd: [round(e["t_s"] / 3600, 2) for e in data["replay DA"][sd]["events"] if e.get("da")] for sd in SEEDS}
    lam_rep = {c: {sd: sum(1 for e in data[c][sd]["events"] if "Λ" in e["replayed"] and e["n_active"] < 23)
                   for sd in SEEDS} for c in COND}
    summary = {}
    for k in PATTERN_ORDER:
        a = {c: np.array([z[c][sd][k] for sd in SEEDS]) for c in COND}
        d = a["replay DA"] - a["no DA"]
        summary[k] = dict(**{f"z_{c}": float(a[c].mean()) for c in COND}, d_mean=float(d.mean()),
                          d_sd=float(d.std(ddof=1)), seeds_d_pos=int(np.sum(d > 1e-6)))
    (RES / "x4_stats.json").write_text(json.dumps(dict(z=z, da_events=da_events, lam_replays=lam_rep,
                                                       summary=summary), indent=1, ensure_ascii=False))
    for k in PATTERN_ORDER:
        s = summary[k]
        print(f"{k:6s} z noDA {s['z_no DA']:.3f}  replayDA {s['z_replay DA']:.3f}  encDA {s['z_encoding DA']:.3f}  "
              f"Δ {s['d_mean']:+.3f}±{s['d_sd']:.3f} ({s['seeds_d_pos']}/5 >0)")
    print("DA events (h):", da_events)
    print("Λ replays per run:", lam_rep)

    plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 3, figsize=(9.6, 3.0), gridspec_kw=dict(width_ratios=[1.2, 1.2, 1]))
    ax = axes[0]
    x = np.arange(len(PATTERN_ORDER))
    for i, c in enumerate(COND):
        m = [summary[k][f"z_{c}"] for k in PATTERN_ORDER]
        ax.bar(x + (i - 1) * 0.27, m, 0.27, color=CCOL[c], label=c)
        for sd in SEEDS:
            ax.plot(x + (i - 1) * 0.27, [z[c][sd][k] for k in PATTERN_ORDER], "k.", ms=2)
    ax.set_xticks(x, PATTERN_ORDER)
    ax.set_ylabel("mean z at 8 h")
    ax.legend(frameon=False, fontsize=7)
    ax.set_title("captured weight per condition", fontsize=8)
    ax = axes[1]
    for c in ("no DA", "replay DA"):
        for sd in SEEDS:
            d = data[c][sd]
            ax.plot(np.array(d["slow_t"]) / 3600, [v["Λ"] for v in d["slow_z"]], color=CCOL[c], lw=0.8,
                    label=c if sd == 1 else None)
    for sd in SEEDS:
        for t in da_events[sd]:
            ax.axvline(t, color=CCOL["replay DA"], lw=0.4, alpha=0.4)
    ax.set_xlabel("time after encoding (h)")
    ax.set_ylabel("mean z of Λ")
    ax.legend(frameon=False, fontsize=7)
    ax.set_title("capture of Λ (vertical lines: replay DA events)", fontsize=8)
    ax = axes[2]
    for sd in SEEDS:
        d = data["replay DA"][sd]
        t = np.array(d["slow_t"]) / 3600
        ax.plot(t, np.array(d["slow_P"])[:, lam].mean(axis=1), lw=0.8, label=f"seed {sd} ({len(da_events[sd])} DA)")
    ax.set_xlabel("time after encoding (h)")
    ax.set_ylabel("PRP pool P of Λ cells")
    ax.legend(frameon=False, fontsize=6)
    ax.set_title("PRPs made by replay DA", fontsize=8)
    for i, a_ in enumerate(axes):
        a_.text(-0.02, 1.02, f"({chr(97 + i)})", transform=a_.transAxes, fontsize=9, fontweight="bold",
                ha="right", va="bottom")
    fig.tight_layout()
    fig.savefig(RES / "fig_pilot_9_x4.png", dpi=150)


if __name__ == "__main__":
    main()
