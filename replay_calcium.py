"""Why replay only depresses: the calcium of the tag rule (neuron/TMExp2SynNMDA_STC.mod) in replay events and in
encoding, reconstructed from the recorded spikes of the seed-1 trace run (run_pilot.py --trace), for
docs/REPLAY_DEPRESSION_2026-11-27.md.

Rule: c jumps c_pre (1.0, CircuitParams) at a presynaptic spike after d_rc + tc_delay (1.5 + 18.8 ms) and c_post
(0.2758) at a postsynaptic spike, decays with tau_c = 48.8 ms. dTt/dt contains +gamma_p (Tt_max - Tt) H(c - theta_p)
and -gamma_d (1 + Tt) H(c - theta_d), with theta_d = 1.2, theta_p = 3, all over tau_h.
For every synapse and event: number of pre and post spikes, peak c, time above theta_d and theta_p, and the actual
change of W in the event (W_after - W_before). The same for the 18 presentation bursts of encoding.
Writes results/layout_B/replay_calcium.json and results/layout_B/fig_replay_calcium.png.
Usage: TINYCA3_LAYOUT=B env -u PYTHONPATH .venv/bin/python replay_calcium.py
"""
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from circuit import CELL_IDS, ROOT, CircuitParams  # noqa: E402
from class_stats import CCOL, load  # noqa: E402

RES = ROOT / "results" / "layout_B"
CP = CircuitParams()
C_PRE, C_POST, TAU_C = CP.c_pre, CP.c_post, 48.8
DELAY = CP.d_rc + 18.8
TH_D, TH_P = 1.2, 3.0
G_D, G_P, TAU_H = 313.1, 1645.6, 688400.0
DT = 0.1
plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False})


def calcium(sp_pre, sp_post, t):
    c = np.zeros_like(t)
    for s in sp_pre:
        m = t >= s + DELAY
        c[m] += C_PRE * np.exp(-(t[m] - s - DELAY) / TAU_C)
    for s in sp_post:
        m = t >= s
        c[m] += C_POST * np.exp(-(t[m] - s) / TAU_C)
    return c


def exposure(sp_pre, sp_post, t0, t1):
    t = np.arange(t0, t1, DT)
    c = calcium(sp_pre, sp_post, t)
    return float(c.max()), float((c > TH_D).sum() * DT), float((c > TH_P).sum() * DT)


def main():
    d = load(1, "_mod10_gate30_trace")
    pre, post = d["p1_pre"], d["p1_post"]
    cls = np.empty(len(pre), dtype=object)
    for key, m in d["masks"].items():
        if key.startswith("syn:"):
            cls[m] = key[4:]
    rows = []                                    # event, synapse, n_pre, n_post, c_max, t_d, t_p, dW
    for ei, e in enumerate(d["events"]):
        sp = {c: np.asarray(e["spikes"][("PYR", c)]) for c in CELL_IDS}
        dW = e["W_after"] - e["W_before"]
        for i in range(len(pre)):
            a, b = sp[CELL_IDS[pre[i]]], sp[CELL_IDS[post[i]]]
            if len(a) or len(b):
                rows.append((ei, i, len(a), len(b), *exposure(a, b, 0.0, 450.0), dW[i]))
    r = np.array(rows)
    enc = d["enc"]
    sp = {c: np.asarray(enc["spikes"][("PYR", c)]) for c in CELL_IDS}
    erows = []                                   # per presentation and synapse with both cells active
    for k, ton in enc["enc_on"].items():
        for i in range(len(pre)):
            a = sp[CELL_IDS[pre[i]]]
            b = sp[CELL_IDS[post[i]]]
            a, b = a[(a >= ton) & (a < ton + 1300)], b[(b >= ton) & (b < ton + 1300)]
            if len(a) > 3 and len(b) > 3:
                erows.append((len(a), len(b), *exposure(a, b, ton, ton + 1400)))
    er = np.array(erows)

    # spikes per active cell in an event and in an encoding burst
    rep_n = [int(np.sum((s > 195) & (s < 440))) for e in d["events"] for s in
             (np.asarray(e["spikes"][("PYR", c)]) for c in CELL_IDS)]
    rep_n = [n for n in rep_n if n]
    on = r[:, 5] > 0
    slope, icpt = np.polyfit(r[on, 5], r[on, 7], 1)
    out = dict(
        params=dict(c_pre=C_PRE, c_post=C_POST, tau_c=TAU_C, theta_d=TH_D, theta_p=TH_P, delay_ms=DELAY),
        replay=dict(synapse_events=int(len(r)), c_max=round(float(r[:, 4].max()), 3),
                    above_theta_p=int((r[:, 6] > 0).sum()), spikes_per_active_cell=np.bincount(rep_n).tolist(),
                    by_n_pre={int(n): dict(n=int((r[:, 2] == n).sum()),
                                           c_max_median=round(float(np.median(r[r[:, 2] == n, 4])), 3),
                                           frac_above_theta_d=round(float(np.mean(r[r[:, 2] == n, 5] > 0)), 3),
                                           mean_dW=round(float(r[r[:, 2] == n, 7].mean()), 5))
                              for n in np.unique(r[:, 2]).astype(int)},
                    dW_per_ms_above_theta_d=round(float(slope), 6),
                    corr=round(float(np.corrcoef(r[on, 5], r[on, 7])[0, 1]), 3),
                    rule_range_per_ms=[round(-G_D / TAU_H * 1.3, 6), round(-G_D / TAU_H * 1.7, 6)],
                    max_abs_dW_below_theta_d=float(np.abs(r[~on, 7]).max()),
                    two_or_more_pre_by_class={c: round(float(np.mean(r[cls[r[:, 1].astype(int)] == c, 2] >= 2)), 3)
                                              for c in ("X", "+", "X→+", "+→X", "=", "‖", "none")}),
        encoding=dict(synapse_presentations=int(len(er)), c_max_median=round(float(np.median(er[:, 2])), 3),
                      ms_above_theta_p_median=round(float(np.median(er[:, 4])), 1),
                      ms_above_theta_d_median=round(float(np.median(er[:, 3])), 1)))

    # minimal burst for potentiation: n presynaptic spikes at ISI, with the postsynaptic cell firing with them
    isis = np.arange(2, 31, 1.0)
    ns = np.arange(1, 9)
    grid = np.zeros((len(ns), len(isis)))
    tg = np.arange(0, 500, DT)
    for a, n in enumerate(ns):
        for b, isi in enumerate(isis):
            s = np.arange(n) * isi + 10
            grid[a, b] = calcium(s, s, tg).max()
    out["burst_c_max"] = dict(isi_ms=isis.tolist(), n=ns.tolist(), c_max=np.round(grid, 3).tolist())
    (RES / "replay_calcium.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))

    fig, axes = plt.subplots(2, 2, figsize=(9.4, 6.8))
    ax = axes[0, 0]
    t = np.arange(0, 140, DT)
    for n, col in ((1, "#7f8c8d"), (2, "#e67e22"), (3, "#c0392b")):
        s = np.arange(n) * 6.0
        ax.plot(t, calcium(s, s, t), color=col, lw=1.5, label=f"replay: {n} spike{'s' if n > 1 else ''}, pre and post, 6 ms apart")
    s = np.arange(9) * 10.5
    ax.plot(t, calcium(s, s, t), color="#2e6fd8", lw=1.5, label="encoding: 9-spike burst at 95 Hz")
    ax.axhspan(TH_D, TH_P, color="#f5b7b1", alpha=0.35, lw=0)
    ax.axhline(TH_P, color="#27ae60", lw=1, ls="--")
    ax.text(138, TH_P + 0.08, "θ_p = 3: potentiation", ha="right", color="#27ae60", fontsize=7)
    ax.text(138, TH_D + 0.08, "θ_d = 1.2: depression", ha="right", color="#c0392b", fontsize=7)
    ax.set_xlabel("ms after the first spike")
    ax.set_ylabel("synaptic calcium c")
    ax.set_title("(a) calcium of one synapse for the spike patterns of replay and encoding", fontsize=8)
    ax.legend(frameon=False, fontsize=6, loc="upper left")
    ax = axes[0, 1]
    bins = np.arange(0.5, 12.5, 1)
    enc_n = [int(np.sum((s >= ton) & (s < ton + 100))) for ton0 in enc["enc_on"].values() for c in CELL_IDS
             for s in [np.asarray(enc["spikes"][("PYR", c)])] for ton in (ton0, ton0 + 500, ton0 + 1000)]
    enc_n = [n for n in enc_n if n > 2]
    ax.hist(rep_n, bins=bins, color="#e67e22", alpha=0.8, label=f"replay event ({len(rep_n)} cell firings)")
    ax.hist(enc_n, bins=bins, color="#2e6fd8", alpha=0.6, label=f"encoding burst, active cells ({len(enc_n)})")
    ax.set_yscale("log")
    ax.set_xlabel("spikes per cell")
    ax.set_ylabel("count")
    ax.set_title("(b) spikes a cell fires in one replay event and in one 100-ms encoding burst", fontsize=8)
    ax.legend(frameon=False, fontsize=6.5)
    ax = axes[1, 0]
    for n, col in ((1, "#7f8c8d"), (2, "#e67e22"), (3, "#c0392b"), (4, "#8e44ad")):
        m = r[:, 2] == n
        ax.plot(r[m, 5], r[m, 7] * 1000, ".", ms=2.5, color=col, alpha=0.5, label=f"{n} presynaptic spike{'s' if n > 1 else ''}")
    xx = np.linspace(0, r[:, 5].max(), 10)
    ax.plot(xx, (-G_D / TAU_H * 1.3 * xx) * 1000, "k--", lw=0.8, label="rule, 1 + T̃ = 1.3")
    ax.plot(xx, (-G_D / TAU_H * 1.7 * xx) * 1000, "k:", lw=0.8, label="rule, 1 + T̃ = 1.7")
    ax.set_xlabel("time with c > θ_d in the event (ms)")
    ax.set_ylabel("ΔW of the synapse in the event (× 10⁻³)")
    ax.set_title(f"(c) replay, seed 1: every synapse-event (r = {out['replay']['corr']})", fontsize=8)
    ax.legend(frameon=False, fontsize=6)
    ax = axes[1, 1]
    im = ax.pcolormesh(isis, ns, grid, cmap="viridis", shading="nearest")
    cs = ax.contour(isis, ns, grid, levels=[TH_D, TH_P], colors=["#f5b7b1", "w"], linewidths=1.2)
    ax.clabel(cs, fmt={TH_D: "θ_d", TH_P: "θ_p"}, fontsize=7)
    ax.plot([6], [1], "o", color="#e67e22", ms=5)
    ax.plot([6], [2], "o", color="#e67e22", ms=5)
    ax.plot([10.5], [8], "s", color="#2e6fd8", ms=5)
    ax.set_xlabel("interspike interval (ms)")
    ax.set_ylabel("spikes (pre and post together)")
    ax.set_title("(d) peak calcium of a pre+post burst: θ_p needs 3 spikes ≤ 8 ms or 4 ≤ 17 ms apart", fontsize=8)
    fig.colorbar(im, ax=ax, label="peak c")
    fig.tight_layout()
    fig.savefig(RES / "fig_replay_calcium.png", dpi=160)
    print(json.dumps({k: v for k, v in out.items() if k != "burst_c_max"}, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
