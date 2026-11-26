"""Null readout and seed statistics for the pilot (results/seeds/pilot_s{seed}[_norep].pkl).

Null readout. An event "replays" pattern k when >= 80% of k's cells spike. With n_active random cells out of 25, the
chance that this happens for a pattern of n_k cells is hypergeometric:
    p_chance(k, n_active) = P(X >= ceil(0.8 n_k)),  X ~ Hypergeom(population 25, successes n_k, draws n_active)
(the same as shuffling the pattern over random cells of equal size). Per seed and pattern:
    observed  = events read out as replay of k (bursts with >= 23 active cells excluded)
    expected  = sum of p_chance over the same events (chance level)
    significant = replays with p_chance < 0.05 / 6 (Bonferroni over the six patterns)
Paired control. The same seed run without SWR events (--no-replay) gives the captured weight z without replay; the
difference replay - control is the effect of replay on capture.
Writes results/replay_stats.json and results/fig_pilot_7_seeds.png (Figure 7).
--suffix S reads pilot_s{seed}S.pkl and pilot_s{seed}S_norep.pkl (e.g. --suffix _mod10_gate30) and adds S to the
output names.
"""
import argparse
import json
import math
import pickle
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.stats import hypergeom  # noqa: E402

from circuit import PATTERN_ORDER, PCOL, RESULTS, PATTERNS  # noqa: E402

RES = RESULTS
SUFFIX = ""
SEEDS = (1, 2, 3, 4, 5)
BURST = 23
ALPHA = 0.05 / len(PATTERN_ORDER)


def p_chance(k, n_active, n_cells=25, frac=0.8):
    n_k = len(PATTERNS[k])
    need = math.ceil(frac * n_k)
    return float(hypergeom.sf(need - 1, n_cells, n_k, n_active))


def load(seed, norep=False):
    with open(RES / "seeds" / f"pilot_s{seed}{SUFFIX}{'_norep' if norep else ''}.pkl", "rb") as f:
        return pickle.load(f)


def seed_stats(d):
    evs = [e for e in d["events"] if e["n_active"] < BURST]
    out = dict(n_events=len(d["events"]), n_bursts=len(d["events"]) - len(evs),
               median_active=float(np.median([e["n_active"] for e in d["events"]])))
    for k in PATTERN_ORDER:
        p = [p_chance(k, e["n_active"]) for e in evs]
        rep = [k in e["replayed"] for e in evs]
        out[k] = dict(observed=int(sum(rep)), expected=float(sum(p)),
                      significant=int(sum(r and pp < ALPHA for r, pp in zip(rep, p))))
    return out


def z_by_pattern(d):
    return {k: float(np.mean(d["z_final"][d["masks"][k]])) for k in PATTERN_ORDER}


def main():
    global SUFFIX
    ap = argparse.ArgumentParser()
    ap.add_argument("--suffix", default="")
    SUFFIX = ap.parse_args().suffix
    stats, z_rep, z_ctl = {}, {}, {}
    for sd in SEEDS:
        d, c = load(sd), load(sd, norep=True)
        stats[sd] = seed_stats(d)
        z_rep[sd], z_ctl[sd] = z_by_pattern(d), z_by_pattern(c)
    summary = {}
    for k in PATTERN_ORDER:
        obs = np.array([stats[sd][k]["observed"] for sd in SEEDS])
        exp = np.array([stats[sd][k]["expected"] for sd in SEEDS])
        sig = np.array([stats[sd][k]["significant"] for sd in SEEDS])
        dz = np.array([z_rep[sd][k] - z_ctl[sd][k] for sd in SEEDS])
        summary[k] = dict(observed_mean=float(obs.mean()), observed_sd=float(obs.std(ddof=1)),
                          expected_mean=float(exp.mean()), expected_sd=float(exp.std(ddof=1)),
                          significant_mean=float(sig.mean()), significant_sd=float(sig.std(ddof=1)),
                          seeds_obs_gt_exp=int(np.sum(obs > exp)),
                          z_replay_mean=float(np.mean([z_rep[sd][k] for sd in SEEDS])),
                          z_control_mean=float(np.mean([z_ctl[sd][k] for sd in SEEDS])),
                          dz_mean=float(dz.mean()), dz_sd=float(dz.std(ddof=1)), seeds_dz_negative=int(np.sum(dz < 0)))
    other = dict(bursts=[stats[sd]["n_bursts"] for sd in SEEDS],
                 median_active=[stats[sd]["median_active"] for sd in SEEDS])
    (RES / f"replay_stats{SUFFIX}.json").write_text(json.dumps(dict(per_seed=stats, summary=summary, other=other,
                                                           z_replay=z_rep, z_control=z_ctl), indent=1,
                                                      ensure_ascii=False))

    plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 3, figsize=(9.6, 3.0))
    x = np.arange(len(PATTERN_ORDER))
    w = 0.38
    ax = axes[0]
    m = [summary[k]["observed_mean"] for k in PATTERN_ORDER]
    sd = [summary[k]["observed_sd"] for k in PATTERN_ORDER]
    ax.bar(x - w / 2, m, w, yerr=sd, color=[PCOL[k] for k in PATTERN_ORDER], capsize=2, label="observed")
    me = [summary[k]["expected_mean"] for k in PATTERN_ORDER]
    se = [summary[k]["expected_sd"] for k in PATTERN_ORDER]
    ax.bar(x + w / 2, me, w, yerr=se, color="white", edgecolor="k", hatch="///", capsize=2, label="chance")
    for sdn in SEEDS:
        ax.plot(x - w / 2 + np.random.default_rng(sdn).uniform(-0.1, 0.1, len(x)),
                [stats[sdn][k]["observed"] for k in PATTERN_ORDER], "k.", ms=3)
    ax.set_xticks(x, PATTERN_ORDER)
    ax.set_ylabel("events per 8-h run (of 80)")
    ax.legend(frameon=False, fontsize=7)
    ax.set_title("replay read out (≥ 80% of cells)", fontsize=8)
    ax = axes[1]
    ex = np.array([[stats[sdn][k]["observed"] - stats[sdn][k]["expected"] for k in PATTERN_ORDER] for sdn in SEEDS])
    ax.bar(x, ex.mean(0), 0.6, yerr=ex.std(0, ddof=1), color=[PCOL[k] for k in PATTERN_ORDER], capsize=2)
    for i, sdn in enumerate(SEEDS):
        ax.plot(x + np.random.default_rng(sdn).uniform(-0.15, 0.15, len(x)), ex[i], "k.", ms=3)
    ax.axhline(0, color="k", lw=0.7)
    ax.set_xticks(x, PATTERN_ORDER)
    ax.set_ylabel("observed − chance (events per run)")
    ax.set_title("excess over chance", fontsize=8)
    ax = axes[2]
    for j, k in enumerate(PATTERN_ORDER):
        for sdn in SEEDS:
            ax.plot([j - 0.15, j + 0.15], [z_ctl[sdn][k], z_rep[sdn][k]], color=PCOL[k], lw=0.8, marker=".", ms=4)
    ax.set_xticks(x, PATTERN_ORDER)
    ax.set_ylabel("mean z at 8 h")
    ax.set_title("z at 8 h: no replay (left) → replay (right)", fontsize=8, pad=10)
    for i, ax in enumerate(axes):
        ax.text(-0.02, 1.02, f"({chr(97 + i)})", transform=ax.transAxes, fontsize=9, fontweight="bold",
                ha="right", va="bottom")
    fig.tight_layout()
    fig.savefig(RES / f"fig_pilot_7_seeds{SUFFIX}.png", dpi=150)
    np.set_printoptions(precision=2)
    for k in PATTERN_ORDER:
        s = summary[k]
        print(f"{k:6s} obs {s['observed_mean']:.1f}±{s['observed_sd']:.1f}  chance {s['expected_mean']:.2f}±"
              f"{s['expected_sd']:.2f}  sig {s['significant_mean']:.1f}±{s['significant_sd']:.1f}  "
              f"obs>chance {s['seeds_obs_gt_exp']}/5  z rep {s['z_replay_mean']:.3f} ctl {s['z_control_mean']:.3f}"
              f"  dz {s['dz_mean']:+.3f}±{s['dz_sd']:.3f} ({s['seeds_dz_negative']}/5 negative)")
    print("bursts", other["bursts"], "median active", other["median_active"])


if __name__ == "__main__":
    main()
