"""Figure P1 of docs/PROTOCOL.md: the pilot protocol at five time scales.

(a) the whole 8-h run; (b) the first 90 min (modulator windows, PRP build-up); (c) the encoding block (14.1 s in layout B, 16.4 s in B7);
(d) one replay epoch; (e) one replay event. Timings are read from run_pilot.py and circuit.py, so the figure
follows the code. Writes results/fig_protocol.png.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

from circuit import PA, PATTERN_ORDER, PB, PCOL, RESULTS, CircuitParams  # noqa: E402
from run_pilot import ENC_START_MS, ENC_STEP_MS, EPOCHS_H, EVENT_GAP_S, EVENT_MS, N_EVENTS, SWR_AT_MS  # noqa: E402

RES = RESULTS
C_DA, C_HT, C_REP, C_OFF, C_ENC = "#c0392b", "#2e6fd8", "#e67e22", "#d9d9d9", "#555555"
MOD_MIN = 10.0          # bio-plausible modulator duration of the protocol (alternatives drawn: 5 and 60 min)
TAU_P_MIN = 60.0
plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.spines.left": False})


def lane(ax, y, x0, x1, color, label=None, h=0.6, **kw):
    ax.add_patch(Rectangle((x0, y - h / 2), x1 - x0, h, color=color, lw=0, **kw))
    if label:
        ax.text(x0, y + h / 2 + 0.05, label, fontsize=6.5, va="bottom")


def main():
    cp = CircuitParams()
    t_enc_s = (ENC_START_MS + len(PATTERN_ORDER) * ENC_STEP_MS) / 1000.0
    on = {k: (ENC_START_MS + i * ENC_STEP_MS) / 1000.0 for i, k in enumerate(PATTERN_ORDER)}
    fig = plt.figure(figsize=(10, 9.2))
    gs = fig.add_gridspec(4, 2, height_ratios=[1.05, 1.25, 1.25, 1.1], hspace=0.62, wspace=0.18)

    # (a) whole run --------------------------------------------------------------------------------------------
    ax = fig.add_subplot(gs[0, :])
    rows = {"encoding": 3, f"DA → {PA} cells": 2, f"5-HT → {PB} cells": 1.4, "replay epochs": 0.6, "offline STC": -0.2}
    lane(ax, rows["offline STC"], 0, 8, C_OFF)
    ax.plot([t_enc_s / 3600] * 2, [2.7, 3.3], color=C_ENC, lw=3)
    ax.text(0.05, 3, f"{(ENC_START_MS + len(PATTERN_ORDER) * ENC_STEP_MS) / 1000:.1f} s, "
                     f"{len(PATTERN_ORDER)} patterns (c)", va="center", fontsize=7)
    for d, ls, a in ((MOD_MIN, "-", 1), (60, ":", 0.35)):
        lane(ax, rows[f"DA → {PA} cells"], on[PA] / 3600, on[PA] / 3600 + d / 60, C_DA, alpha=a, h=0.4)
        lane(ax, rows[f"5-HT → {PB} cells"], on[PB] / 3600, on[PB] / 3600 + d / 60, C_HT, alpha=a, h=0.4)
    ax.text(1.05, rows[f"DA → {PA} cells"], f"released {MOD_MIN:g} min (solid) + synthesis trigger; 60 min L&T form (pale)",
            va="center", fontsize=6.5)
    for i, e in enumerate(EPOCHS_H):
        lane(ax, rows["replay epochs"], e, e + N_EVENTS * EVENT_GAP_S / 3600, C_REP)
        ax.text(e + 0.25, rows["replay epochs"], f"E{i + 1}", ha="center", va="center", fontsize=7, color="white",
                fontweight="bold")
    ax.text(6.6, rows["replay epochs"], f"{N_EVENTS} SWR events per epoch (d)", va="center", fontsize=6.5)
    ax.annotate("readout\nat 8 h", (8, -0.2), (7.6, 2.4), fontsize=6.5, ha="center",
                arrowprops=dict(arrowstyle="->", lw=0.6))
    ax.set_yticks(list(rows.values()), list(rows.keys()))
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(-0.05, 8.05)
    ax.set_ylim(-0.7, 3.5)
    ax.set_xlabel("time (h)")
    ax.set_title("(a) whole run: 8 h of real time, no compression", fontsize=8, loc="left")

    # (b) first 90 min: modulator windows and PRP pool -----------------------------------------------------------
    ax = fig.add_subplot(gs[1, :])
    t = np.linspace(0, 90, 901)
    for d, ls in ((5, "--"), (MOD_MIN, "-"), (60, ":")):
        P = np.where(t < d, 1 - np.exp(-t / TAU_P_MIN), (1 - np.exp(-d / TAU_P_MIN)) * np.exp(-(t - d) / TAU_P_MIN))
        ax.plot(t, P, color="k", ls=ls, lw=1.1, label=f"{d:g} min, L&T (θ_pro follows M)")
    # synthesis trigger: 10-min release, G = dG/dt = M (1 - G) / 2 min - G / 30 min (stc.py)
    G, Gs = 0.0, []
    for ti in t:
        Gs.append(G)
        M = 1.0 if ti < MOD_MIN else 0.0
        G += 0.1 * (M * (1 - G) / 2.0 - G / 30.0)
    ax.plot(t, Gs, color="#2e86c1", lw=1.0, alpha=0.8, label="trigger G, 10-min release, τ_G = 30 min")
    end = 68.0                                     # synthesis ends ~68 min (mean of 5 seeds, report of 2026-11-24)
    P = np.where(t < end, 1 - np.exp(-t / TAU_P_MIN), (1 - np.exp(-end / TAU_P_MIN)) * np.exp(-(t - end) / TAU_P_MIN))
    ax.plot(t, P, color="#d35400", lw=1.6, label="P, 10-min release + trigger (bio-plausible)")
    ax.axvspan(30, 60, color=C_REP, alpha=0.15, lw=0)
    ax.text(45, 0.88, "replay epoch E1 (0.5–1 h)", ha="center", fontsize=7, color="#a0521d")
    ax.set_xlim(0, 90)
    ax.set_ylim(0, 1.0)
    ax.set_xlabel("time (min)")
    ax.set_ylabel("PRP pool P / α_P\n(cell with tags above θ)")
    ax.legend(frameon=False, fontsize=6.5, title="DA or 5-HT released for", title_fontsize=6.5,
              loc="upper left", bbox_to_anchor=(1.0, 1.0))
    ax.spines["left"].set_visible(True)
    ax.set_title("(b) first 90 min: the modulator (or the trigger G it starts) lowers the PRP threshold; P builds and decays with τ_P = 60 min",
                 fontsize=8, loc="left")

    # (c) encoding block ------------------------------------------------------------------------------------------
    ax = fig.add_subplot(gs[2, 0])
    for i, k in enumerate(PATTERN_ORDER):
        y = len(PATTERN_ORDER) - i
        for b in range(3):
            t0 = on[k] + b * 0.5
            lane(ax, y, t0, t0 + 0.1, PCOL[k], h=0.55)
        ax.text(on[k] + 1.25, y, f"EC_{k}", va="center", fontsize=6.5, color=PCOL[k])
    ax.annotate("DA on", (on[PA], 7.0), (on[PA], 7.7), color=C_DA, fontsize=6.5, ha="center",
                arrowprops=dict(arrowstyle="->", color=C_DA, lw=0.8))
    ax.annotate("5-HT on", (on[PB], 7.0), (on[PB], 7.7), color=C_HT, fontsize=6.5, ha="center",
                arrowprops=dict(arrowstyle="->", color=C_HT, lw=0.8))
    ax.set_yticks([])
    ax.set_xlim(0, t_enc_s)
    ax.set_ylim(0.4, 8.2)
    ax.set_xlabel("time (s)")
    ax.set_title("(c) encoding block, 14.1 s: 3 × 100 ms at 100 Hz into EC_k,\n0.5 s apart; next pattern 2.3 s later",
                 fontsize=8, loc="left")

    # (d) one replay epoch ------------------------------------------------------------------------------------------
    ax = fig.add_subplot(gs[2, 1])
    span = N_EVENTS * EVENT_GAP_S / 60
    lane(ax, 0, 0, span, C_OFF, h=0.5)
    for i in range(N_EVENTS):
        ax.plot([i * EVENT_GAP_S / 60] * 2, [-0.45, 0.45], color=C_REP, lw=2)
    ax.text(0, 0.62, "SWR event (e), 450 ms", color="#a0521d", fontsize=6.5)
    ax.text(EVENT_GAP_S / 60 / 2, -0.62, "90 s offline", ha="center", va="top", fontsize=6.5)
    ax.annotate("", (0, -0.5), (EVENT_GAP_S / 60, -0.5), arrowprops=dict(arrowstyle="<->", lw=0.6))
    ax.set_yticks([])
    ax.set_xlim(-0.5, span + 0.5)
    ax.set_ylim(-1.1, 1.1)
    ax.set_xlabel("time from epoch start (min)")
    ax.set_title(f"(d) one replay epoch: {N_EVENTS} events, {EVENT_GAP_S:g} s apart (30 min)", fontsize=8, loc="left")

    # (e) one replay event -----------------------------------------------------------------------------------------
    ax = fig.add_subplot(gs[3, 0])
    lane(ax, 2, 0, EVENT_MS, "#f3d9b1", h=0.5)
    ax.text(5, 2.38, f"replay state: P1 × {cp.replay_p1_scale:g}, BC → PYR × {cp.replay_inh_scale:g}",
            fontsize=6.5)
    ax.plot([SWR_AT_MS] * 2, [0.65, 1.35], color="k", lw=2)
    ax.text(SWR_AT_MS + 6, 1.0, f"SWT trigger: {cp.swt_active} of {cp.n_swt} cells fire once,\n"
            f"seed ~3–8 random PYR (jitter ≤ {cp.swt_jitter:g} ms)", va="center", fontsize=6.5)
    lane(ax, 0, SWR_AT_MS - 5, EVENT_MS - 10, C_REP, h=0.4, alpha=0.6)
    ax.text(SWR_AT_MS - 5, -0.4, "readout window 195–440 ms", fontsize=6.5, va="top")
    ax.set_yticks([2, 1, 0], ["state", "trigger", "readout"])
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(0, EVENT_MS)
    ax.set_ylim(-0.9, 2.8)
    ax.set_xlabel("time in event (ms)")
    ax.set_title("(e) one replay event, 450 ms of NEURON simulation", fontsize=8, loc="left")

    # (f) readout rules -----------------------------------------------------------------------------------------
    ax = fig.add_subplot(gs[3, 1])
    ax.axis("off")
    ax.set_title("(f) readout of an event", fontsize=8, loc="left")
    txt = ("active cells = PYR that spike in the readout window\n"
           "pattern k replayed  ⇔  ≥ 80% of its cells active\n"
           "burst (not replay)  ⇔  ≥ 23 of 25 PYR active\n"
           "chance: hypergeometric draw of the same number of cells\n"
           f"per-event significance: p_chance < 0.05 / {len(PATTERN_ORDER)}\n\n"
           "weights: W_k = mean 1 + T̃ + z of the P1 synapses inside k")
    ax.text(0, 0.95, txt, va="top", fontsize=7, linespacing=1.6, family="DejaVu Sans")

    fig.savefig(RES / "fig_protocol.png", dpi=150, bbox_inches="tight")


if __name__ == "__main__":
    main()
