"""Figures of the pilot run (results/pilot.npz, results/pilot_events.pkl) -> results/fig_pilot_*.png.

Numbering continues docs/RESULTS_2026-11-17.md (Figure 1 is the T1 f-I curve):
  Figure 2  encoding block: patterns and modulators, spikes and population potentials by cell type, tags
  Figure 3  replay events: four SWR events at close range
  Figure 4  the 8-h dynamics: DA/5-HT, SWR events and replays; tags, captured weights, total weights, PRP pools
  Figure 5  P1 weight matrices (1 + Tt + z) at three times
  Figure 6  replay share per pattern and epoch
The "population potential" is the mean somatic membrane potential of a cell type, a proxy for the LFP, not a computed
extracellular potential.
"""
import pickle
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from circuit import CELL_IDS, IDX, PATTERN_ORDER, PATTERNS  # noqa: E402

RES = Path(__file__).resolve().parent / "results"
BURST = 23   # events with >= 23 of 25 PYR active are whole-network bursts, not pattern replays
PCOL = {"Λ": "#c0392b", "V": "#2e6fd8", "I": "#2e9e4f", "O": "#7f7f7f", "Λ∪V": "#8e44ad", "Λ∪V∪I": "#a0522d"}
TCOL = {"PYR": "#d35400", "BC": "#1f5fbf", "EC": "#16a085", "DG": "#8e44ad", "SWT": "#555555"}
DA_COL, HT_COL = "#c0392b", "#2e6fd8"
plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False})


def _letters(axes, x=-0.02, y=1.0):
    """Panel indexes (a), (b), ... at the top-left of each axes, in reading order."""
    for i, ax in enumerate(axes):
        ax.text(x, y, f"({chr(97 + i)})", transform=ax.transAxes, fontsize=9, fontweight="bold",
                ha="right", va="bottom")


def _load():
    d = np.load(RES / "pilot.npz")
    with open(RES / "pilot_events.pkl", "rb") as f:
        e = pickle.load(f)
    return d, e


def _types():
    return [("PYR", list(CELL_IDS)), ("BC", list(CELL_IDS)), ("EC", list(PATTERN_ORDER)), ("DG", list(PATTERN_ORDER)),
            ("SWT", list(range(10)))]


def _raster(ax, spikes, t_scale=1.0, xlim=None):
    row, ticks, labels = 0, [], []
    for typ, keys in _types():
        start = row
        for k in keys:
            st = spikes.get((typ, k), np.array([]))
            if xlim is not None:
                st = st[(st >= xlim[0]) & (st <= xlim[1])]
            col = TCOL[typ]
            if typ == "PYR":
                member = [p for p in ("Λ", "V", "I", "O") if k in PATTERNS[p]]
                col = PCOL[member[0]] if member else TCOL["PYR"]
            big = typ in ("EC", "DG", "SWT")
            ax.plot(st * t_scale, np.full(len(st), row), "|", color=col, ms=5.0 if big else 2.5, mew=1.4 if big else 0.6)
            row += 1
        ticks.append((start + row - 1) / 2)
        labels.append(typ)
        ax.axhline(row - 0.5, color="#dddddd", lw=0.5)
    ax.set_yticks(ticks, labels)
    ax.set_ylim(-1, row)
    ax.invert_yaxis()


def _pop_v(ax_list, t, v, t_scale=1.0):
    for ax, typ in zip(ax_list, ("PYR", "BC", "EC", "DG")):
        keys = [k for k in v if k[0] == typ]
        m = np.mean([v[k] for k in keys], axis=0)
        ax.plot(t * t_scale, m, color=TCOL[typ], lw=0.6)
        ax.set_ylabel(f"{typ}\n(mV)", rotation=0, ha="right", va="center")


# ---------------------------------------------------------------------------------------------------------------------
def fig_encoding(e):
    enc, enc_on = e["enc"], e["enc_on"]
    t = enc["t"] / 1000.0
    fig, axes = plt.subplots(7, 1, figsize=(7.2, 8.4), sharex=True,
                             gridspec_kw=dict(height_ratios=[0.6, 3.0, 0.6, 0.6, 0.6, 0.6, 1.6]))
    ax = axes[0]
    for k, t_on in enc_on.items():
        a = t_on / 1000.0
        ax.axvspan(a, a + 1.3, color=PCOL[k], alpha=0.15, lw=0)
        ax.text(a + 0.65, 0.65, k, color=PCOL[k], ha="center", va="center", fontsize=11, fontweight="bold")
    t_end = t[-1]
    for p in e["pulses"]:
        y, col, glyph = (0.15, DA_COL, "Λ") if p["source"] == "DA" else (0.0, HT_COL, "V")
        ax.plot([p["t_on_s"], t_end], [y, y], color=col, lw=3, solid_capstyle="butt")
        ax.text(p["t_on_s"], y + 0.12, f"{p['source']} {glyph}", color=col, fontsize=7, va="bottom")
    ax.set_ylim(-0.1, 1.0)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    _raster(axes[1], enc["spikes"], t_scale=1e-3)
    axes[1].set_ylabel("spikes")
    _pop_v(axes[2:6], t, enc["v"])
    for k in PATTERN_ORDER:
        m = e["masks"][k]
        tr = enc["Tt_trace"][m].mean(axis=0)
        axes[6].plot(np.arange(len(tr)) * 0.005, tr, color=PCOL[k], lw=1, label=k)
    axes[6].axhline(0.2, color="k", ls=":", lw=0.7)
    axes[6].axhline(-0.2, color="k", ls=":", lw=0.7)
    axes[6].set_ylabel("mean T̃\n(tag)")
    axes[6].legend(ncol=6, frameon=False, fontsize=7, loc="lower right")
    axes[6].set_xlabel("time (s)")
    _letters(axes)
    fig.align_ylabels(axes)
    fig.tight_layout()
    fig.savefig(RES / "fig_pilot_2_encoding.png", dpi=150)
    plt.close(fig)


def _pick_events(events, epochs_h=(0.5, 1.5, 3.0, 6.0)):
    picked = []
    for ep in epochs_h:
        evs = [x for x in events if ep * 3600 <= x["t_s"] < ep * 3600 + 1800]
        rep = [x for x in evs if x["replayed"]]
        picked.append(rep[0] if rep else evs[0])
    return picked


def fig_replay_events(e):
    evs = _pick_events(e["events"])
    fig, axes = plt.subplots(5, len(evs), figsize=(9.0, 6.6), sharex=True,
                             gridspec_kw=dict(height_ratios=[3.0, 0.7, 0.7, 0.7, 0.7]))
    for j, ev in enumerate(evs):
        t = ev["t"]
        _raster(axes[0, j], ev["spikes"])
        _pop_v(axes[1:, j], t, ev["v"])
        for ax in axes[:, j]:
            ax.axvline(200.0, color="k", lw=0.6, ls="--")
            ax.set_xlim(100, 450)
        glyphs = " + ".join(ev["replayed"]) if ev["replayed"] else "no replay"
        col = PCOL[ev["replayed"][0]] if ev["replayed"] else "k"
        axes[0, j].set_title(f"t = {ev['t_s'] / 3600:.2f} h\n{glyphs}", color=col, fontsize=9)
        axes[-1, j].set_xlabel("time in event (ms)")
        if j:
            for ax in axes[:, j]:
                ax.set_ylabel("")
                ax.set_yticklabels([])
    _letters(axes[0, :], x=-0.02, y=1.18)
    fig.tight_layout()
    fig.savefig(RES / "fig_pilot_3_replay_events.png", dpi=150)
    plt.close(fig)


def fig_overview(d, e):
    th = d["slow_t"] / 3600.0
    fig, axes = plt.subplots(5, 1, figsize=(7.2, 8.0), sharex=True,
                             gridspec_kw=dict(height_ratios=[1.3, 1, 1, 1, 1]))
    ax = axes[0]
    for p in e["pulses"]:
        y, col, glyph = (2.2, DA_COL, "Λ") if p["source"] == "DA" else (1.7, HT_COL, "V")
        a, b = p["t_on_s"] / 3600, (p["t_on_s"] + p["dur_s"]) / 3600
        ax.plot([a, b], [y, y], color=col, lw=4, solid_capstyle="butt")
        ax.text(a - 0.05, y, glyph, color=col, va="center", ha="right", fontsize=12, fontweight="bold")
        ax.text(7.95, y, f"{glyph}: {p['source']} release onto {glyph} cells, {a:.0f}–{b:.0f} h", color=col,
                va="center", ha="right", fontsize=7)
    ax.text(0.0, 2.85, "encoding at 0 h: " + "  ".join(PATTERN_ORDER), fontsize=7, va="center")
    for ep in (0.5, 1.5, 3.0, 6.0):
        for a_ in axes:
            a_.axvspan(ep, ep + 0.5, color="#f0e6c8", alpha=0.6, lw=0)
    for ev in e["events"]:
        x = ev["t_s"] / 3600
        if ev["n_active"] >= BURST:
            ax.plot([x, x], [0.2, 1.1], color="k", lw=0.9)
        elif ev["replayed"]:
            k = max(ev["replayed"], key=lambda q: len(PATTERNS[q]))     # largest pattern read out
            ax.plot([x, x], [0.2, 1.1], color=PCOL[k], lw=1.2)
        else:
            ax.plot([x, x], [0.2, 0.6], color="#999999", lw=0.6)
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], color=PCOL[k], lw=2, label=k) for k in PATTERN_ORDER]
    handles.append(Line2D([], [], color="k", lw=2, label="all 25 (burst)"))
    ax.legend(handles=handles, ncol=7, frameon=False, fontsize=7, loc="lower right", bbox_to_anchor=(1.0, 0.18),
              title="replay read out as (largest pattern ≥ 80% active):", title_fontsize=7, handlelength=1.0,
              columnspacing=0.8)
    ax.set_ylim(0, 3.1)
    ax.set_yticks([0.6, 1.95], ["SWR events", "modulators"])
    for k in PATTERN_ORDER:
        tt = np.array([x[k] for x in e["slow_Tt"]])
        zz = np.array([x[k] for x in e["slow_z"]])
        axes[1].plot(th, tt, color=PCOL[k], lw=1, label=k)
        axes[2].plot(th, zz, color=PCOL[k], lw=1)
        axes[3].plot(th, 1 + tt + zz, color=PCOL[k], lw=1)
    axes[1].axhline(0.2, color="k", ls=":", lw=0.7)
    axes[1].set_ylabel("mean T̃\n(tag)")
    axes[1].legend(ncol=6, frameon=False, fontsize=7, loc="upper right")
    axes[2].set_ylabel("mean z\n(captured)")
    axes[3].set_ylabel("weight\n1 + T̃ + z")
    P = d["slow_P"]
    groups = {"Λ cells (DA)": set(PATTERNS["Λ"]), "V cells (5-HT)": set(PATTERNS["V"]),
              "I only": set(PATTERNS["I"]) - set(PATTERNS["Λ∪V"]),
              "O only": set(PATTERNS["O"]) - set(PATTERNS["Λ∪V"])}
    gcol = {"Λ cells (DA)": PCOL["Λ"], "V cells (5-HT)": PCOL["V"], "I only": PCOL["I"], "O only": PCOL["O"]}
    gsty = {"Λ cells (DA)": dict(lw=3.0, ls="-", alpha=0.5), "V cells (5-HT)": dict(lw=1.2, ls="--"),
            "I only": dict(lw=3.0, ls="-", alpha=0.5), "O only": dict(lw=1.2, ls="--")}
    for g, cells in groups.items():
        axes[4].plot(th, P[:, [IDX[c] for c in cells]].mean(axis=1), color=gcol[g], label=g, **gsty[g])
    axes[4].set_ylabel("PRP pool P")
    axes[4].legend(frameon=False, fontsize=7)
    axes[4].set_xlabel("time after encoding (h)")
    axes[4].set_xlim(-0.1, 8.0)
    _letters(axes)
    fig.align_ylabels(axes)
    fig.tight_layout()
    fig.savefig(RES / "fig_pilot_4_overview.png", dpi=150)
    plt.close(fig)


def fig_weights(d, e):
    snaps = e["snapshots"]
    fig, axes = plt.subplots(1, len(snaps), figsize=(9.0, 3.3))
    order = [cid for k in ("Λ", "V", "I", "O") for cid in PATTERNS[k]]
    order = list(dict.fromkeys(order)) + [c for c in CELL_IDS if c not in set(order)]
    pos = {cid: i for i, cid in enumerate(order)}
    for ax, (label, (t_s, w)) in zip(axes, snaps.items()):
        W = np.full((25, 25), np.nan)
        for s, (a, b) in enumerate(zip(d["p1_pre"], d["p1_post"])):
            W[pos[CELL_IDS[b]], pos[CELL_IDS[a]]] = 1 + w[s]
        im = ax.imshow(W, cmap="magma", vmin=0.9, vmax=1.7)
        ax.set_title(f"t = {t_s / 3600:.2f} h", fontsize=8, pad=16)
        ax.set_xlabel("presynaptic cell")
        ax.set_ylabel("postsynaptic cell")
        start = 0
        for k in ("Λ", "V", "I", "O"):
            n_new = len([c for c in PATTERNS[k] if pos[c] >= start])
            ax.add_patch(plt.Rectangle((start - 0.5, start - 0.5), n_new, n_new, fill=False, ec=PCOL[k], lw=1.2))
            blk = W[start:start + n_new, start:start + n_new]
            ax.text(start + n_new - 1.0, start + 0.3, f"{np.nanmean(blk):.2f}", color="w",
                    ha="right", va="top", fontsize=7, fontweight="bold",
                    bbox=dict(boxstyle="round,pad=0.1", fc="k", ec="none", alpha=0.5))
            ax.text(start + n_new / 2, -1.2, k, color=PCOL[k], ha="center", fontsize=9, fontweight="bold")
            start += n_new
        ax.set_xticks([])
        ax.set_yticks([])
    _letters(axes, x=-0.08, y=1.12)
    fig.colorbar(im, ax=axes, shrink=0.8, label="weight multiplier 1 + T̃ + z")
    fig.savefig(RES / "fig_pilot_5_weights.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def fig_replay_share(e):
    epochs = (0.5, 1.5, 3.0, 6.0)
    share = np.zeros((len(epochs), len(PATTERN_ORDER)))
    bursts = []
    for i, ep in enumerate(epochs):
        evs = [x for x in e["events"] if ep * 3600 <= x["t_s"] < ep * 3600 + 1800]
        for j, k in enumerate(PATTERN_ORDER):
            share[i, j] = np.mean([k in x["replayed"] and x["n_active"] < BURST for x in evs]) if evs else 0
        bursts.append(np.mean([x["n_active"] >= BURST for x in evs]) if evs else 0)
    fig, ax = plt.subplots(figsize=(6.0, 2.8))
    w = 0.13
    for j, k in enumerate(PATTERN_ORDER):
        ax.bar(np.arange(len(epochs)) + (j - 2.5) * w, share[:, j], w, color=PCOL[k], label=k)
    ax.bar(np.arange(len(epochs)) + 3.5 * w, bursts, w, color="k", label="all (burst)")
    ax.set_xticks(range(len(epochs)), [f"{ep:g}–{ep + 0.5:g} h" for ep in epochs])
    ax.set_ylabel("fraction of SWR events\nreplaying the pattern")
    ax.set_xlabel("replay epoch")
    ax.legend(ncol=7, frameon=False, fontsize=7, loc="upper center", bbox_to_anchor=(0.5, 1.18))
    fig.tight_layout()
    fig.savefig(RES / "fig_pilot_6_replay_share.png", dpi=150)
    plt.close(fig)
    return share


def main():
    d, e = _load()
    fig_encoding(e)
    fig_replay_events(e)
    fig_overview(d, e)
    fig_weights(d, e)
    share = fig_replay_share(e)
    np.set_printoptions(precision=2, suppress=True)
    print("replay share (epochs x patterns", PATTERN_ORDER, ")\n", share)
    last = e["slow_Tt"][-1], e["slow_z"][-1]
    print("final mean Tt:", {k: round(v, 3) for k, v in last[0].items()})
    print("final mean z :", {k: round(v, 3) for k, v in last[1].items()})


if __name__ == "__main__":
    main()
