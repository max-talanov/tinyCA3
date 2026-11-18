"""DA during replay (VTA mode): compare three conditions per seed.

  control : no SWR events               (results/seeds/pilot_s{seed}_norep.pkl)
  replay  : 80 SWR events               (results/seeds/pilot_s{seed}.pkl)
  replay+DA: the same, and every event that replays Λ (non-burst) releases DA onto the Λ cells for 10 min
            (results/seeds/pilot_s{seed}_rda10.pkl)
Prints mean z at 8 h per pattern and the paired differences; writes results/da_stats.json and
results/fig_pilot_8_replay_da.png (Figure 8).
"""
import json
import pickle
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from circuit import PATTERN_ORDER  # noqa: E402

RES = Path(__file__).resolve().parent / "results"
SEEDS = (1, 2, 3, 4, 5)
COND = {"control": "_norep", "replay": "", "replay+DA": "_rda10"}
PCOL = {"Λ": "#c0392b", "V": "#2e6fd8", "I": "#2e9e4f", "O": "#7f7f7f", "Λ∪V": "#8e44ad", "Λ∪V∪I": "#a0522d"}


def load(seed, suffix):
    with open(RES / "seeds" / f"pilot_s{seed}{suffix}.pkl", "rb") as f:
        return pickle.load(f)


def main():
    z = {c: {} for c in COND}
    n_da, slow = {}, {}
    for sd in SEEDS:
        for c, suf in COND.items():
            d = load(sd, suf)
            z[c][sd] = {k: float(np.mean(d["z_final"][d["masks"][k]])) for k in PATTERN_ORDER}
            if c == "replay+DA":
                n_da[sd] = sum(1 for e in d["events"] if e.get("da"))
                slow[sd] = (np.array(d["slow_t"]), np.array([x["Λ"] for x in d["slow_z"]]))
            if c == "replay" and sd == 1:
                slow_rep = (np.array(d["slow_t"]), np.array([x["Λ"] for x in d["slow_z"]]))
    summary = {}
    for k in PATTERN_ORDER:
        a = {c: np.array([z[c][sd][k] for sd in SEEDS]) for c in COND}
        d_da = a["replay+DA"] - a["replay"]
        d_ctl = a["replay+DA"] - a["control"]
        summary[k] = dict(**{f"z_{c}": float(a[c].mean()) for c in COND},
                          d_da_mean=float(d_da.mean()), d_da_sd=float(d_da.std(ddof=1)), seeds_d_da_pos=int(np.sum(d_da > 1e-6)),
                          d_ctl_mean=float(d_ctl.mean()), d_ctl_sd=float(d_ctl.std(ddof=1)),
                          seeds_d_ctl_pos=int(np.sum(d_ctl > 1e-6)))
    (RES / "da_stats.json").write_text(json.dumps(dict(z=z, n_da_events=n_da, summary=summary), indent=1,
                                                  ensure_ascii=False))
    for k in PATTERN_ORDER:
        s = summary[k]
        print(f"{k:6s} z ctl {s['z_control']:.3f} rep {s['z_replay']:.3f} rep+DA {s['z_replay+DA']:.3f} | "
              f"DA effect {s['d_da_mean']:+.3f}±{s['d_da_sd']:.3f} ({s['seeds_d_da_pos']}/5 >0) | "
              f"vs control {s['d_ctl_mean']:+.3f}±{s['d_ctl_sd']:.3f} ({s['seeds_d_ctl_pos']}/5 >0)")
    print("DA events per seed:", n_da)

    plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 3, figsize=(9.6, 3.0), gridspec_kw=dict(width_ratios=[1.4, 1, 1]))
    ax = axes[0]
    x = np.arange(len(PATTERN_ORDER))
    for j, k in enumerate(PATTERN_ORDER):
        for sd in SEEDS:
            ax.plot([j - 0.25, j, j + 0.25], [z["control"][sd][k], z["replay"][sd][k], z["replay+DA"][sd][k]],
                    color=PCOL[k], lw=0.7, marker=".", ms=4)
    ax.set_xticks(x, PATTERN_ORDER)
    ax.set_ylabel("mean z at 8 h")
    ax.set_title("control → replay → replay + DA (one line per seed)", fontsize=8)
    ax = axes[1]
    for j, k in enumerate(PATTERN_ORDER):
        dd = [z["replay+DA"][sd][k] - z["replay"][sd][k] for sd in SEEDS]
        ax.bar(j, np.mean(dd), 0.6, yerr=np.std(dd, ddof=1), color=PCOL[k], capsize=2)
        ax.plot(j + np.random.default_rng(j).uniform(-0.15, 0.15, len(SEEDS)), dd, "k.", ms=3)
    ax.axhline(0, color="k", lw=0.7)
    ax.set_xticks(x, PATTERN_ORDER)
    ax.set_ylabel("Δz (replay + DA − replay)")
    ax.set_title("effect of DA during replay", fontsize=8)
    ax = axes[2]
    for sd in SEEDS:
        t, zz = slow[sd]
        ax.plot(t / 3600, zz, lw=0.8, label=f"seed {sd} ({n_da[sd]} DA events)")
    t, zz = slow_rep
    ax.plot(t / 3600, zz, "k--", lw=0.8, label="seed 1, replay without DA")
    ax.set_xlabel("time after encoding (h)")
    ax.set_ylabel("mean z of Λ")
    ax.set_title("capture of Λ over time", fontsize=8)
    ax.legend(frameon=False, fontsize=6)
    for i, a_ in enumerate(axes):
        a_.text(-0.02, 1.02, f"({chr(97 + i)})", transform=a_.transAxes, fontsize=9, fontweight="bold",
                ha="right", va="bottom")
    fig.tight_layout()
    fig.savefig(RES / "fig_pilot_8_replay_da.png", dpi=150)


if __name__ == "__main__":
    main()
