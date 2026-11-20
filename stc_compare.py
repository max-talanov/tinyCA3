"""Synapses with and without STC: STP against STC over time, and the weights during replay (report of 2026-11-20).

Runs (results/seeds/, all made by run_pilot.py --light, which records the probe synapses of run_pilot.PROBES):
  with STC     pilot_s{seed}.pkl          the calibrated pilot of the 2026-11-17 report, sections 5-7
  without STC  pilot_s{seed}_nostc.pkl    --no-stc: P1 synapses with plastic = 0, so Tt = z = 0 and only the
                                          Tsodyks-Markram release (STP) changes the transmitted weight
Seed 1 of the "with STC" run is rerun with probes (Figure 10); seeds 2-5 are the existing runs (Figure 11 uses only
their pattern means and events).

Weight of one P1 synapse at presynaptic spike n, relative to the baseline NetCon weight w0 of the encoding state:
    W_n = s_state * (rel_n / U) * (1 + Tt(t_n) + z(t_n))
          s_state: 1 in encoding, replay_p1_scale (2.5) in replay blocks (brain-state switch, report section 4)
          rel_n / U: STP, the release of spike n relative to the first spike of a rested synapse (ms to s)
          1 + Tt + z: STC, tag plus captured weight (minutes to hours); 1 without STC
Writes results/fig_cmp_10_stp_stc.png, results/fig_cmp_11_replay.png and results/stc_compare.json.
"""
import json
import pickle
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from circuit import IDX, PATTERN_ORDER, CircuitParams  # noqa: E402
from reference.stc_numpy import advance_slow  # noqa: E402
from replay_stats import BURST, p_chance  # noqa: E402
from stc import STCParams  # noqa: E402

RES = Path(__file__).resolve().parent / "results"
SEEDS = (1, 2, 3, 4, 5)
PCOL = {"Λ": "#c0392b", "V": "#2e6fd8", "I": "#2e9e4f", "O": "#7f7f7f", "Λ∪V": "#8e44ad", "Λ∪V∪I": "#a0522d"}
STC_COL, STP_COL = "#c0392b", "#555555"
EPOCHS_H = (0.5, 1.5, 3.0, 6.0)
ENC_ON_MS = {k: 300.0 + i * 2300.0 for i, k in enumerate(PATTERN_ORDER)}
REPLAY_SCALE = CircuitParams().replay_p1_scale
plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False})


def load(seed, suffix=""):
    with open(RES / "seeds" / f"pilot_s{seed}{suffix}.pkl", "rb") as f:
        return pickle.load(f)


def spikes_of(block, name, U):
    """Per-spike values of one probe synapse in one online block: time (ms in block), rel/U, 1 + Tt + z."""
    p = block[name]
    return p["spk_t"], p["spk_rel"] / U, 1.0 + p["spk_Tt"] + p["spk_z"]


def probe_slow(d, name, sp=STCParams()):
    """Tt and z of one probe synapse over the whole run: recorded inside online blocks, integrated offline between
    them with the logged PRP pool of its postsynaptic cell (logged every ~60 s, held between log points)."""
    post = IDX[d["probe_pairs"][name][1]]
    slow_t, slow_P = d["slow_t"], d["slow_P"][:, post]
    blocks = d["probes"]
    t_out, w_out = [], []
    for i, b in enumerate(blocks):
        t0 = b["t0_s"]
        tb = t0 + b["t"] / 1000.0
        Tt, z = b[name]["Tt"], b[name]["z"]
        t_out.append(tb)
        w_out.append(1 + Tt + z)
        t_next = blocks[i + 1]["t0_s"] if i + 1 < len(blocks) else 8 * 3600.0
        t, Tt_, z_ = tb[-1], Tt[-1], z[-1]
        ts, ws = [], []
        while t < t_next - 1e-6:
            dt = min(sp.dt_slow_s, t_next - t)
            P = slow_P[max(np.searchsorted(slow_t, t, side="right") - 1, 0)]
            Tt_, z_ = advance_slow(Tt_, z_, P, dt * 1000.0, k_rel=sp.k_rel, tau_h=sp.tau_h_ms,
                                   tau_z=sp.tau_z_min * 60000.0)
            t += dt
            ts.append(t)
            ws.append(1 + Tt_ + z_)
        t_out.append(np.array(ts))
        w_out.append(np.array(ws))
    return np.concatenate(t_out), np.concatenate(w_out)


def pattern_weight(d):
    """Mean 1 + Tt + z of the synapses inside each pattern, from the slow log (every ~60 s)."""
    t = np.asarray(d["slow_t"]) / 3600.0
    return t, {k: np.array([a[k] + b[k] for a, b in zip(d["slow_Tt"], d["slow_z"])]) + 1 for k in PATTERN_ORDER}


def at_events(d):
    """Mean 1 + Tt + z per pattern at the start of each replay event (slow-log value just before the event)."""
    t, w = pattern_weight(d)
    ev_t = np.array([e["t_s"] for e in d["events"]]) / 3600.0
    i = np.clip(np.searchsorted(t, ev_t, side="right") - 1, 0, len(t) - 1)
    return ev_t, {k: w[k][i] for k in PATTERN_ORDER}


def replay_counts(d):
    evs = [e for e in d["events"] if e["n_active"] < BURST]
    out = {k: dict(observed=int(sum(k in e["replayed"] for e in evs)),
                   expected=float(sum(p_chance(k, e["n_active"]) for e in evs))) for k in PATTERN_ORDER}
    out["bursts"] = len(d["events"]) - len(evs)
    out["median_active"] = float(np.median([e["n_active"] for e in d["events"]]))
    return out


def letters(axes):
    for i, ax in enumerate(axes):
        ax.text(-0.02, 1.02, f"({chr(97 + i)})", transform=ax.transAxes, fontsize=9, fontweight="bold",
                ha="right", va="bottom")


# ------------------------------------------------------------------------------------------------------------------
def figure_10(stc, nostc, summary):
    U = stc["U"]
    enc, enc0 = stc["probes"][0], nostc["probes"][0]
    fig = plt.figure(figsize=(9.6, 6.6))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.05], hspace=0.42, wspace=0.22)
    axes = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, 0]),
            fig.add_subplot(gs[1, 1])]

    # (a) one Λ encoding (3 bursts of 10 spikes at 100 Hz): STP inside each burst, STC between bursts
    ax = axes[0]
    for blk, lab, col, mk in ((enc0, "without STC (STP only)", STP_COL, "s"), (enc, "with STC", STC_COL, "o")):
        t, stp, stcf = spikes_of(blk, "Λ", U)
        m = (t > 250) & (t < 1500)
        ax.plot(t[m], (stp * stcf)[m], "-", color=col, lw=0.6, marker=mk, ms=2.5, label=lab)
    t, stp, stcf = spikes_of(enc, "Λ", U)
    m = (t > 250) & (t < 1500)
    ax.step(t[m], stcf[m], where="post", color=STC_COL, lw=1.2, ls="--", label="STC factor 1 + T̃ + z")
    ax.axhline(1, color="k", lw=0.5, ls=":")
    ax.set_xlabel("time in encoding block (ms)")
    ax.set_ylabel("weight per spike W / w0")
    ax.set_title("Λ synapse ex_22 → ex_24, first Λ encoding", fontsize=8)
    ax.legend(frameon=False, fontsize=6.5, loc="upper right")

    # (b) whole encoding block: the STC factor steps up at each encoding that contains the synapse
    ax = axes[1]
    for k, t_on in ENC_ON_MS.items():
        ax.axvspan(t_on / 1000, (t_on + 1100) / 1000, color=PCOL[k], alpha=0.10, lw=0)
        ax.text((t_on + 550) / 1000, 0.02, k, transform=ax.get_xaxis_transform(), ha="center", fontsize=7,
                color=PCOL[k])
    for blk, col, mk in ((enc0, STP_COL, "s"), (enc, STC_COL, "o")):
        t, stp, stcf = spikes_of(blk, "Λ", U)
        ax.plot(t / 1000, stp * stcf, ".", color=col, ms=2)
    for name, ls in (("Λ", "-"), ("I", "-.")):
        b = enc[name]
        ax.plot(enc["t"] / 1000, 1 + b["Tt"] + b["z"], color=PCOL[name], ls=ls, lw=1.2,
                label=f"STC factor, {name} synapse")
    ax.plot(enc0["t"] / 1000, 1 + enc0["Λ"]["Tt"] + enc0["Λ"]["z"], color=STP_COL, lw=1.0, ls="--",
            label="STC factor without STC (= 1)")
    ax.set_xlabel("time in encoding block (s)")
    ax.set_ylabel("W / w0 (dots: spikes of the Λ synapse)")
    ax.set_title("encoding block, 0–14.1 s", fontsize=8)
    ax.legend(frameon=False, fontsize=6.5, loc="upper left")

    # (c) one log time axis: STP is the spread of the spike dots, STC the slow line
    ax = axes[2]
    for d, col, mk, lab in ((nostc, STP_COL, "s", "without STC"), (stc, STC_COL, "o", "with STC")):
        ts, ws = [], []
        for blk in d["probes"]:
            t, stp, stcf = spikes_of(blk, "Λ", U)
            sc = 1.0 if blk["kind"] == "encoding" else REPLAY_SCALE
            ts.append(blk["t0_s"] + t / 1000)
            ws.append(sc * stp * stcf)
        ax.plot(np.concatenate(ts), np.concatenate(ws), mk, color=col, ms=2, alpha=0.6, mew=0,
                label=f"spikes, {lab}")
    tt, wt = probe_slow(stc, "Λ")
    ax.plot(tt[tt > 0.05], wt[tt > 0.05], color=STC_COL, lw=1.4, label="STC factor 1 + T̃ + z")
    ax.axhline(1, color=STP_COL, lw=1.0, ls="--", label="STC factor without STC")
    for e in EPOCHS_H:
        ax.axvspan(e * 3600, e * 3600 + 1800, color="#999", alpha=0.12, lw=0)
    ax.set_xscale("log")
    ax.set_xlim(0.2, 8 * 3600)
    ax.set_xticks([1, 10, 60, 600, 3600, 4 * 3600], ["1 s", "10 s", "1 min", "10 min", "1 h", "4 h"])
    ax.set_xlabel("time from start of encoding (log)")
    ax.set_ylabel("W / w0, Λ synapse")
    ax.set_title("Λ synapse, 0.2 s to 8 h (grey: replay epochs, P1 × 2.5)", fontsize=8)
    ax.legend(frameon=False, fontsize=6.5, loc="upper left", ncol=2)

    # (d) 8 h, pattern means, all five seeds
    ax = axes[3]
    for k in ("Λ", "V", "I", "O"):
        ws = []
        for sd in SEEDS:
            t, w = pattern_weight(summary["runs"][sd]["stc"])
            ws.append(np.interp(np.linspace(0.004, 8, 500), t, w[k]))
        ws = np.array(ws)
        tg = np.linspace(0.004, 8, 500)
        ax.fill_between(tg, ws.min(0), ws.max(0), color=PCOL[k], alpha=0.2, lw=0)
        ax.plot(tg, ws.mean(0), color=PCOL[k], lw=1.3, label=f"{k}, with STC")
    for sd in SEEDS:
        t, w = pattern_weight(summary["runs"][sd]["nostc"])
        ax.plot(t, w["Λ"], color=STP_COL, lw=1.0, ls="--", label="every pattern, without STC" if sd == 1 else None)
    for e in EPOCHS_H:
        ax.axvspan(e, e + 0.5, color="#999", alpha=0.12, lw=0)
    ax.axvspan(0, 1, ymin=0.96, ymax=1, color="#e9a3a3", lw=0)
    ax.text(0.5, 0.985, "DA, 5-HT", transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=6.5)
    ax.set_xlabel("time (h)")
    ax.set_ylabel("mean 1 + T̃ + z inside the pattern")
    ax.set_title("8 h, five seeds (line: mean, band: range)", fontsize=8)
    ax.legend(frameon=False, fontsize=6.5, loc="upper right")
    letters(axes)
    fig.savefig(RES / "fig_cmp_10_stp_stc.png", dpi=150, bbox_inches="tight")


def figure_11(summary):
    runs = summary["runs"]
    fig, axes = plt.subplots(1, 3, figsize=(9.6, 3.2), gridspec_kw=dict(width_ratios=[1.35, 1, 1]))

    # (a) pattern weight at each of the 80 events
    ax = axes[0]
    for k in ("Λ", "V", "I", "O"):
        w = np.array([at_events(runs[sd]["stc"])[1][k] for sd in SEEDS])
        x = np.arange(w.shape[1])
        ax.fill_between(x, w.min(0), w.max(0), color=PCOL[k], alpha=0.2, lw=0)
        ax.plot(x, w.mean(0), color=PCOL[k], lw=1.2, marker=".", ms=2.5, label=k)
    w0 = np.array([at_events(runs[sd]["nostc"])[1]["Λ"] for sd in SEEDS])
    ax.plot(np.arange(w0.shape[1]), w0.mean(0), color=STP_COL, ls="--", lw=1, label="all, without STC")
    for i, e in enumerate(EPOCHS_H):
        ax.axvline(20 * i - 0.5, color="k", lw=0.4)
        ax.text(20 * i + 10, 1.0, f"{e:g}–{e + 0.5:g} h", transform=ax.get_xaxis_transform(), ha="center",
                va="bottom", fontsize=6.5)
    ax.set_xlabel("replay event (20 per epoch)")
    ax.set_ylabel("mean 1 + T̃ + z at the event")
    ax.set_title("pattern weight at each event, 5 seeds", fontsize=8, pad=12)
    ax.legend(frameon=False, fontsize=6.5, ncol=2, loc="upper right")

    # (b) replays read out per pattern, with and without STC, against chance
    ax = axes[1]
    x = np.arange(len(PATTERN_ORDER))
    wbar = 0.27
    for j, (cond, col, hatch, lab) in enumerate((("stc", None, None, "with STC"),
                                                 ("nostc", "#bbbbbb", None, "without STC"),
                                                 ("chance_stc", "white", "///", "chance (with STC)"))):
        if cond == "chance_stc":
            vals = np.array([[summary["counts"][sd]["stc"][k]["expected"] for k in PATTERN_ORDER] for sd in SEEDS])
        else:
            vals = np.array([[summary["counts"][sd][cond][k]["observed"] for k in PATTERN_ORDER] for sd in SEEDS])
        colors = [PCOL[k] for k in PATTERN_ORDER] if col is None else col
        ax.bar(x + (j - 1) * wbar, vals.mean(0), wbar, color=colors, edgecolor="k", lw=0.4, hatch=hatch, label=lab)
        if cond != "chance_stc":
            for i in range(len(SEEDS)):
                ax.plot(x + (j - 1) * wbar + (i - 2) * 0.04, vals[i], "k.", ms=2.5)
    ax.set_xticks(x, PATTERN_ORDER)
    ax.set_ylabel("events read out as replay (of 80)")
    ax.set_title("replay content (bars: mean, dots: seeds)", fontsize=8)
    ax.legend(frameon=False, fontsize=6.5)

    # (c) cells recruited per event
    ax = axes[2]
    bins = np.arange(0, 27) - 0.5
    for cond, col, lab in (("nostc", STP_COL, "without STC"), ("stc", STC_COL, "with STC")):
        n = np.concatenate([[e["n_active"] for e in runs[sd][cond]["events"]] for sd in SEEDS])
        ax.hist(n, bins, color=col, alpha=0.5, label=f"{lab} (median {np.median(n):g})")
    ax.set_xlabel("PYR active in the event")
    ax.set_ylabel("events (5 seeds × 80)")
    ax.set_title("size of the replay events", fontsize=8)
    ax.legend(frameon=False, fontsize=6.5)
    letters(axes)
    fig.tight_layout()
    fig.savefig(RES / "fig_cmp_11_replay.png", dpi=150)


def main():
    runs = {sd: dict(stc=load(sd), nostc=load(sd, "_nostc")) for sd in SEEDS}
    counts = {sd: {c: replay_counts(runs[sd][c]) for c in ("stc", "nostc")} for sd in SEEDS}
    summary = dict(runs=runs, counts=counts)
    stc1, nostc1 = runs[1]["stc"], runs[1]["nostc"]
    U = stc1["U"]

    # numbers for the report
    out = dict(U=U, replay_scale=REPLAY_SCALE)
    for lab, d in (("stc", stc1), ("nostc", nostc1)):
        t, stp, stcf = spikes_of(d["probes"][0], "Λ", U)
        m = (t > 250) & (t < 1500)
        out[f"seed1_{lab}_first_lambda"] = dict(stp_min=float(stp[m].min()), stp_max=float(stp[m].max()),
                                                stc_first=float(stcf[m][0]), stc_last=float(stcf[m][-1]),
                                                W_max=float((stp * stcf)[m].max()), n_spikes=int(m.sum()))
        out[f"seed1_{lab}_encoding_end"] = {n: float(1 + d["probes"][0][n]["Tt"][-1] + d["probes"][0][n]["z"][-1])
                                            for n in d["probe_pairs"]}
        stp_all = np.concatenate([spikes_of(b, "Λ", U)[1] for b in d["probes"]])
        out[f"seed1_{lab}_stp_range_all_spikes"] = [float(stp_all.min()), float(stp_all.max())]
        out[f"seed1_{lab}_n_spikes_Λ_probe"] = int(stp_all.size)
    tt, wt = probe_slow(stc1, "Λ")
    out["seed1_stc_probe_Λ_at_h"] = {h: float(np.interp(h * 3600, tt, wt)) for h in (0.5, 1, 2, 4, 8)}
    for c in ("stc", "nostc"):
        pw = {}
        for k in PATTERN_ORDER:
            vals = []
            for sd in SEEDS:
                t, w = pattern_weight(runs[sd][c])
                vals.append([float(np.interp(h, t, w[k])) for h in (0.0039, 1, 2, 4, 8)])
            pw[k] = np.mean(vals, 0).round(4).tolist()
        out[f"pattern_weight_{c}_at_0_1_2_4_8h"] = pw
        ev = {k: [] for k in PATTERN_ORDER}
        for sd in SEEDS:
            _, w = at_events(runs[sd][c])
            for k in PATTERN_ORDER:
                ev[k].append([float(w[k][20 * i:20 * i + 20].mean()) for i in range(4)])
        out[f"event_weight_{c}_per_epoch"] = {k: np.mean(v, 0).round(4).tolist() for k, v in ev.items()}
        out[f"replays_{c}"] = {k: dict(mean=float(np.mean([counts[sd][c][k]["observed"] for sd in SEEDS])),
                                       sd=float(np.std([counts[sd][c][k]["observed"] for sd in SEEDS], ddof=1)),
                                       chance=float(np.mean([counts[sd][c][k]["expected"] for sd in SEEDS])),
                                       per_seed=[counts[sd][c][k]["observed"] for sd in SEEDS])
                               for k in PATTERN_ORDER}
        out[f"bursts_{c}"] = [counts[sd][c]["bursts"] for sd in SEEDS]
        out[f"median_active_{c}"] = [counts[sd][c]["median_active"] for sd in SEEDS]
        out[f"z_final_max_{c}"] = float(max(np.abs(runs[sd][c]["z_final"]).max() for sd in SEEDS))
        out[f"Tt_final_max_{c}"] = float(max(np.abs(runs[sd][c]["Tt_final"]).max() for sd in SEEDS))
    (RES / "stc_compare.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
    print(json.dumps(out, indent=1, ensure_ascii=False))
    figure_10(stc1, nostc1, summary)
    figure_11(summary)


if __name__ == "__main__":
    main()
