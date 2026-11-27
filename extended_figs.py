"""Figures E1-E6 of docs/DA_5HT_REPLAY_EXTENDED_2026-11-27.md: weight dynamics of the DA- and 5-HT-gated synapse classes
with replay, and the electrical activity of every cell group (layout B, 10 min + trigger, calibrated drive).

Reads the trace runs of seed 1 (run_pilot.py --seed 1 --light --trace --mod-min 10 --gate 30 [--prp-block |
--no-replay]: slow state every 10 s, spikes and voltages of encoding and of every replay event, every synapse's
W = 1 + T̃ + z before and after each event) and the 5-seed light runs (class_stats.COND).
  E1  whole run: W per class, replay and no replay, PRP block; modulators, trigger G and PRP pool P
  E2  first phase: encoding (5 ms tags) to the end of replay epoch 1 (1.0 h), with the change of W in every event
  E3  middle: replay epoch 3 (3.0-3.5 h) and around it
  E4  end: replay epoch 4 (6.0-6.5 h) to 8 h
  E5  up/down: per-synapse W distributions, potentiated and depressed fractions, per-event up/down changes, 5 seeds
  E6  activity: rasters and rates of PYR (X, +, =, ‖), BC, EC, DG and SWT in encoding and in replay events
Writes results/layout_B/extended/fig_E{1..6}_*.png and extended_numbers.json.
Usage: TINYCA3_LAYOUT=B env -u PYTHONPATH .venv/bin/python extended_figs.py
"""
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from circuit import CELL_IDS, PATTERN_ORDER, PATTERNS, PCOL, ROOT  # noqa: E402
from class_stats import CCOL, CLS, COND, SEEDS, load  # noqa: E402

RES = ROOT / "results" / "layout_B"
OUT = RES / "extended"
TR = {"replay": "_mod10_gate30_trace", "no replay": "_mod10_gate30_norep_trace",
      "PRP block": "_mod10_gate30_prpblock_trace"}
CLASSES = ("X", "+→X", "+", "X→+", "=", "‖", "=→‖", "‖→=", "none")
LAB = {"X": "X (DA)", "+→X": "+→X (DA)", "+": "+ (5-HT)", "X→+": "X→+ (5-HT)", "=": "= (none)", "‖": "‖ (none)",
       "=→‖": "=→‖ (none)", "‖→=": "‖→= (none)", "none": "never co-presented"}
GATE = {"DA": ("X", "+→X"), "5-HT": ("+", "X→+"), "none": ("=", "‖", "=→‖", "‖→=")}
EPOCHS_H = (0.5, 1.5, 3.0, 6.0)
EPOCH_LEN_H = 19 * 90 / 3600.0
GROUPS = [("PYR", g) for g in ("X", "+", "=", "‖")] + [("BC", None), ("EC", None), ("DG", None), ("SWT", None)]
GCOL = {"PYR X": PCOL["X"], "PYR +": PCOL["+"], "PYR =": PCOL["="], "PYR ‖": PCOL["‖"], "BC": "#555555",
        "EC": "#d4a017", "DG": "#16a085", "SWT": "#8e44ad"}
plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False})


def W_series(d, k):
    t = np.asarray(d["slow_t"]) / 3600.0
    return t, np.array([1 + a["syn:" + k] + b["syn:" + k] for a, b in zip(d["slow_Tt"], d["slow_z"])])


def shade_epochs(ax, t0=0, t1=8):
    for e in EPOCHS_H:
        if e < t1 and e + EPOCH_LEN_H > t0:
            ax.axvspan(e, e + EPOCH_LEN_H, color="#f2e6c9", lw=0, zorder=0)


def ev_dW(d, k):
    m = d["masks"]["syn:" + k]
    return np.array([float(np.mean(e["W_after"][m] - e["W_before"][m])) for e in d["events"]])


def group_keys(typ, grp):
    if typ == "PYR":
        return [("PYR", c) for c in PATTERNS[grp]]
    if typ == "BC":
        return [("BC", c) for c in CELL_IDS]
    if typ in ("EC", "DG"):
        return [(typ, k) for k in PATTERN_ORDER]
    return [("SWT", i) for i in range(10)]


def gname(typ, grp):
    return f"PYR {grp}" if grp else typ


# ------------------------------------------------------------------------------------------------------------------
def fig_whole(D, nums):
    d = D["replay"]
    fig, axes = plt.subplots(3, 1, figsize=(9, 8), sharex=True, gridspec_kw=dict(height_ratios=[1, 3, 1.4]))
    ax = axes[0]
    t = np.asarray(d["slow_t"]) / 3600.0
    for m, col, idx in (("DA", "#c0392b", "X"), ("5-HT", "#2e6fd8", "+")):
        cells = [CELL_IDS.index(c) for c in PATTERNS[idx]]
        ax.plot(t, d["slow_G"][:, cells].mean(1), color=col, lw=1.4, label=f"trigger G in {idx} cells ({m})")
        ax.plot(t, d["slow_P"][:, cells].mean(1), color=col, lw=1.0, ls="--", label=f"PRP pool P in {idx} cells")
    other = [CELL_IDS.index(c) for c in PATTERNS["="] + PATTERNS["‖"]]
    ax.plot(t, d["slow_P"][:, other].mean(1), color="0.5", lw=1.0, ls="--", label="P in = and ‖ cells")
    for p in d["pulses"]:
        ax.axvspan(p["t_on_s"] / 3600, (p["t_on_s"] + p["dur_s"]) / 3600, color="#c0392b" if p["source"] == "DA"
                   else "#2e6fd8", alpha=0.25, lw=0)
    shade_epochs(ax)
    ax.set_ylabel("G, P")
    ax.set_title("(a) modulators (bands: 10-min DA and 5-HT release), synthesis trigger and PRPs, seed 1", fontsize=8)
    ax.legend(frameon=False, fontsize=6, ncol=3)
    ax = axes[1]
    shade_epochs(ax)
    for k in CLASSES:
        tt, w = W_series(d, k)
        ax.plot(tt, w, color=CCOL[k], ls=CLS[k], lw=1.5, label=LAB[k])
        tt, w = W_series(D["no replay"], k)
        if k in GATE["DA"] + GATE["5-HT"]:
            ax.plot(tt, w, color=CCOL[k], ls=(0, (1, 1)), lw=0.9)
    ax.set_ylabel("W = 1 + T̃ + z (class mean)")
    ax.set_title("(b) seed 1: with replay (lines), without replay (fine dots, modulated classes); shaded: replay epochs",
                 fontsize=8)
    ax.legend(frameon=False, fontsize=6, ncol=3, loc="upper right")
    ax = axes[2]
    shade_epochs(ax)
    for gate, ks in (("DA", GATE["DA"]), ("5-HT", GATE["5-HT"]), ("none", GATE["none"])):
        for c, ls in (("replay", "-"), ("PRP block", "--")):
            ws = []
            for s in SEEDS:
                dd = load(s, COND["10 min + trigger" if c == "replay" else "PRP block"])
                ws.append(np.mean([W_series(dd, k)[1] for k in ks], 0))
            tt = np.asarray(load(SEEDS[0], COND["PRP block"])["slow_t"]) / 3600.0
            ax.plot(tt, np.mean(ws, 0), ls=ls, lw=1.5,
                    color={"DA": "#c0392b", "5-HT": "#2e6fd8", "none": "0.5"}[gate],
                    label=f"{gate}-gated, {'PRPs' if c == 'replay' else 'PRP block'}")
            nums.setdefault("E1_gate_means_8h", {})[f"{gate}, {c}"] = round(float(np.mean(ws, 0)[-1]), 4)
    ax.set_xlabel("time from encoding (h)")
    ax.set_ylabel("W, 5 seeds")
    ax.set_title("(c) mean over the classes each modulator gates, 5 seeds: PRPs (solid) and PRP block (dashed)",
                 fontsize=8)
    ax.legend(frameon=False, fontsize=6, ncol=3)
    ax.set_xlim(0, 8.05)
    fig.tight_layout()
    fig.savefig(OUT / "fig_E1_whole.png", dpi=160)
    plt.close(fig)


def zoom(D, t0, t1, title, fname, enc_panel=False):
    d = D["replay"]
    ncol = 3 if enc_panel else 2
    fig, axes = plt.subplots(1, ncol, figsize=(4.1 * ncol, 3.6), gridspec_kw=dict(width_ratios=[1.2] * (ncol - 1) + [1]))
    ai = 0
    if enc_panel:
        ax = axes[0]
        e = d["enc"]
        tt = np.arange(e["Tt_trace"].shape[1]) * 5.0 / 1000.0
        for k in CLASSES:
            ax.plot(tt, e["Tt_trace"][d["masks"]["syn:" + k]].mean(0) + 1, color=CCOL[k], ls=CLS[k], lw=1.2,
                    label=LAB[k])
        for k, ton in e["enc_on"].items():
            ax.axvline(ton / 1000, color=PCOL[k], lw=0.6, alpha=0.6)
            ax.text(ton / 1000, 1.72, k, color=PCOL[k], fontsize=7)
        for p in d["pulses"]:
            ax.annotate(f"{p['source']} on", xy=(p["t_on_s"], 0.98), fontsize=6,
                        color="#c0392b" if p["source"] == "DA" else "#2e6fd8")
        ax.set_xlabel("time (s)")
        ax.set_ylabel("W = 1 + T̃ (z = 0 during encoding)")
        ax.set_title("(a) encoding, tags every 5 ms", fontsize=8)
        ax.set_ylim(0.95, 1.78)
        ax.legend(frameon=False, fontsize=5.5, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.17))
        ai = 1
    ax = axes[ai]
    shade_epochs(ax, t0, t1)
    for k in CLASSES:
        tt, w = W_series(d, k)
        sel = (tt >= t0) & (tt <= t1)
        ax.plot(tt[sel], w[sel], color=CCOL[k], ls=CLS[k], lw=1.4, label=LAB[k])
        if k in GATE["DA"] + GATE["5-HT"]:
            tt2, w2 = W_series(D["no replay"], k)
            sel = (tt2 >= t0) & (tt2 <= t1)
            ax.plot(tt2[sel], w2[sel], color=CCOL[k], ls=(0, (1, 1)), lw=0.9)
    evt = np.array([e["t_s"] for e in d["events"]]) / 3600.0
    insel = (evt >= t0) & (evt <= t1)
    for x, e in zip(evt[insel], np.array(d["events"], dtype=object)[insel]):
        ax.axvline(x, ymax=0.04, color="#8e44ad" if e["replayed"] else "0.7", lw=1.0)
    ax.set_xlabel("time from encoding (h)")
    ax.set_ylabel("W = 1 + T̃ + z")
    ax.set_title(f"({'b' if enc_panel else 'a'}) {title}; ticks: events (purple: a pattern replayed)", fontsize=8)
    if not enc_panel:
        ax.legend(frameon=False, fontsize=5.5, ncol=2)
    ax = axes[ai + 1]
    idx = np.nonzero(insel)[0]
    for k in ("X", "+→X", "+", "X→+", "=", "‖"):
        ax.plot(evt[idx], ev_dW(d, k)[idx] * 1000, marker="o", ms=2.5, lw=0.8, color=CCOL[k], ls=CLS[k], label=LAB[k])
    ax.axhline(0, color="0.6", lw=0.6)
    ax.set_xlabel("event time (h)")
    ax.set_ylabel("ΔW in the event (× 10⁻³, class mean)")
    ax.set_title(f"({'c' if enc_panel else 'b'}) change of W in each 450-ms event", fontsize=8)
    ax.legend(frameon=False, fontsize=5.5, ncol=2)
    fig.tight_layout()
    fig.savefig(OUT / fname, dpi=160)
    plt.close(fig)


def fig_updown(D, nums):
    fig, axes = plt.subplots(2, 2, figsize=(9.4, 6.6))
    d = D["replay"]
    ax = axes[0, 0]
    bins = np.linspace(1.40, 1.75, 71)
    for (name, (ts, W)), col in zip(d["snapshots"].items(), ("#95a5a6", "#2980b9", "#2c3e50")):
        for gate, ls in (("DA", "-"), ("5-HT", "--")):
            m = np.any([d["masks"]["syn:" + k] for k in GATE[gate]], 0)
            ax.hist(1 + W[m], bins=bins, histtype="step", color=col, ls=ls, lw=1.2,
                    label=f"{gate}-gated, {name}")
    ax.set_xlabel("W = 1 + T̃ + z per synapse")
    ax.set_ylabel("synapses")
    ax.set_title("(a) seed 1: W of DA-gated (solid) and 5-HT-gated (dashed) synapses at three times", fontsize=8)
    ax.legend(frameon=False, fontsize=5.5)
    ax = axes[0, 1]
    bounds = (0.0, 0.5, 1.5, 3.0, 6.0, 8.0)
    evt_h = np.array([e["t_s"] for e in d["events"]]) / 3600.0
    ks = ("X", "+→X", "+", "X→+", "=", "‖")
    decomp = {}
    width = 0.13
    for j, k in enumerate(ks):
        tt, w = W_series(d, k)
        dW = ev_dW(d, k)
        rp, off = [], []
        for a, b in zip(bounds[:-1], bounds[1:]):
            total = float(np.interp(b, tt, w) - np.interp(a, tt, w))
            r = float(dW[(evt_h >= a) & (evt_h < b)].sum())
            rp.append(r)
            off.append(total - r)
        decomp[k] = dict(offline=[round(v, 4) for v in off], replay=[round(v, 4) for v in rp])
        x = np.arange(len(bounds) - 1) + (j - 2.5) * width
        ax.bar(x, off, width=width, color=CCOL[k], alpha=0.45 if CLS[k] != "-" else 0.9,
               hatch=None if CLS[k] == "-" else "///", edgecolor=CCOL[k], lw=0.5, label=LAB[k])
        ax.bar(x, rp, width=width, bottom=[min(o, 0) for o in off], color="k", alpha=0.85, lw=0)
    ax.axhline(0, color="0.6", lw=0.6)
    ax.set_xticks(range(len(bounds) - 1), [f"{a:g}-{b:g} h" for a, b in zip(bounds[:-1], bounds[1:])])
    ax.set_ylabel("change of class-mean W in the interval")
    ax.set_title("(b) seed 1: offline change (colour) and change in replay events (black)", fontsize=8)
    ax.legend(frameon=False, fontsize=5.5, ncol=2)
    nums["E5_decomposition"] = dict(intervals_h=bounds, classes=decomp)
    ax = axes[1, 0]
    up_n, dn_n, up_s, dn_s = [], [], [], []
    for e in d["events"]:
        dW = e["W_after"] - e["W_before"]
        up_n.append(int(np.sum(dW > 1e-3)))
        dn_n.append(int(np.sum(dW < -1e-3)))
        up_s.append(float(dW[dW > 0].sum()))
        dn_s.append(float(dW[dW < 0].sum()))
    x = np.arange(len(d["events"]))
    ax.bar(x, up_n, color="#27ae60", width=0.8, label="synapses up (ΔW > 0.001)")
    ax.bar(x, [-v for v in dn_n], color="#c0392b", width=0.8, label="synapses down (ΔW < −0.001)")
    for i, e in enumerate(d["events"]):
        if e["replayed"]:
            ax.text(i, max(up_n) * 1.02, "·", color="#8e44ad", ha="center", fontsize=10)
    for b in (20, 40, 60):
        ax.axvline(b - 0.5, color="0.7", lw=0.6)
    ax.set_xlabel("event number (20 per epoch: 0.5, 1.5, 3, 6 h)")
    ax.set_ylabel("synapses (of 600)")
    ax.set_title("(c) seed 1: synapses moved up or down in each event (dots: a pattern replayed)", fontsize=8)
    ax.legend(frameon=False, fontsize=6)
    nums["E5_event_updown"] = dict(up_n=up_n, down_n=dn_n, up_sum=[round(v, 4) for v in up_s],
                                   down_sum=[round(v, 4) for v in dn_s])
    ax = axes[1, 1]
    tt = np.asarray(load(SEEDS[0], COND["10 min + trigger"])["slow_t"]) / 3600.0
    for c, ls in (("10 min + trigger", "-"), ("10 min + trigger, no replay", ":"), ("PRP block", "--")):
        r = []
        for s in SEEDS:
            dd = load(s, COND[c])
            r.append(np.mean([W_series(dd, k)[1] for k in GATE["DA"]], 0) /
                     np.mean([W_series(dd, k)[1] for k in GATE["5-HT"]], 0))
        r = np.array(r)
        tt = np.asarray(load(SEEDS[0], COND[c])["slow_t"]) / 3600.0
        ax.plot(tt, r.mean(0), ls=ls, color="k", lw=1.3, label=c)
        if c == "10 min + trigger":
            ax.fill_between(tt, r.min(0), r.max(0), color="0.85", alpha=0.7, lw=0)
        nums.setdefault("E5_DA_over_5HT", {})[c] = dict(max=round(float(r.mean(0).max()), 4),
                                                       min=round(float(r.mean(0).min()), 4),
                                                       at_8h=round(float(r.mean(0)[-1]), 4))
    shade_epochs(ax)
    ax.axhline(1, color="0.6", lw=0.6)
    ax.set_xlabel("time from encoding (h)")
    ax.set_ylabel("W(DA-gated) / W(5-HT-gated)")
    ax.set_title("(d) 5 seeds: W(DA-gated) / W(5-HT-gated); band: seed range, replay", fontsize=8)
    ax.legend(frameon=False, fontsize=6)
    fig.tight_layout()
    fig.savefig(OUT / "fig_E5_updown.png", dpi=160)
    plt.close(fig)


def raster(ax, spikes, t0, t1, scale=1.0):
    y, ticks = 0, []
    for typ, grp in GROUPS:
        keys = group_keys(typ, grp)
        for i, key in enumerate(keys):
            st = np.asarray(spikes.get(key, []))
            st = st[(st >= t0) & (st <= t1)]
            ax.plot(st * scale, np.full(st.size, y + i), "|", ms=2.2, mew=0.6, color=GCOL[gname(typ, grp)])
        ticks.append((y + len(keys) / 2, gname(typ, grp)))
        y += len(keys) + 2
    ax.set_yticks([a for a, _ in ticks], [b for _, b in ticks], fontsize=6)
    ax.set_ylim(-1, y)
    ax.invert_yaxis()


def fig_activity(D, nums):
    d = D["replay"]
    e = d["enc"]
    fig = plt.figure(figsize=(10, 9.5))
    gs = fig.add_gridspec(3, 3, height_ratios=[1.5, 1.2, 1.1])
    ax = fig.add_subplot(gs[0, :])
    raster(ax, e["spikes"], 0, e["t_enc_ms"], 1 / 1000)
    for k, ton in e["enc_on"].items():
        ax.axvline(ton / 1000, color=PCOL[k], lw=0.6, alpha=0.5)
        ax.text(ton / 1000, -1.5, k, color=PCOL[k], fontsize=7)
    ax.set_xlabel("time (s)")
    ax.set_title("(a) encoding block, seed 1: spikes of every cell group", fontsize=8)
    # rates per group in each presentation (3 bursts in 1.3 s after onset)
    ax = fig.add_subplot(gs[1, :2])
    rates = {}
    for typ, grp in GROUPS:
        keys = group_keys(typ, grp)
        rates[gname(typ, grp)] = [float(np.mean([np.sum((np.asarray(e["spikes"][key]) >= ton) &
                                                        (np.asarray(e["spikes"][key]) < ton + 1300)) for key in keys]) / 0.3)
                                  for ton in e["enc_on"].values()]
    w = 0.1
    for i, (g, r) in enumerate(rates.items()):
        ax.bar(np.arange(len(PATTERN_ORDER)) + (i - 4) * w, r, width=w, color=GCOL[g], label=g)
    ax.set_xticks(range(len(PATTERN_ORDER)), [f"{k} presented" for k in PATTERN_ORDER])
    ax.set_ylabel("mean rate per cell in the 3 bursts (Hz)")
    ax.set_title("(b) encoding: mean rate of each group while each pattern is presented", fontsize=8)
    ax.legend(frameon=False, fontsize=5.5, ncol=4)
    nums["E6_encoding_rates_Hz"] = {g: dict(zip(PATTERN_ORDER, [round(v, 1) for v in r])) for g, r in rates.items()}
    # participation in replay events
    ax = fig.add_subplot(gs[1, 2])
    evt = np.array([ev["t_s"] for ev in d["events"]]) / 3600.0
    part = {}
    for typ, grp in GROUPS:
        keys = group_keys(typ, grp)
        part[gname(typ, grp)] = [float(np.mean([np.any((np.asarray(ev["spikes"][key]) > 195) &
                                                       (np.asarray(ev["spikes"][key]) < 440)) for key in keys]))
                                 for ev in d["events"]]
    names = [g for g in part if g not in ("EC", "DG")]
    im = ax.imshow(np.array([part[g] for g in names]), aspect="auto", cmap="Greys", vmin=0, vmax=1,
                   interpolation="nearest")
    for i, ev in enumerate(d["events"]):
        if ev["replayed"]:
            ax.text(i, -0.9, "▾", ha="center", fontsize=5, color="#8e44ad")
    ax.set_yticks(range(len(names)), names, fontsize=6)
    for b in (20, 40, 60):
        ax.axvline(b - 0.5, color="#e67e22", lw=0.6)
    ax.set_ylim(len(names) - 0.5, -1.3)
    ax.set_xlabel("event number (epochs of 20)")
    ax.set_title("(c) fraction of each group firing per event (▾: replay)", fontsize=8)
    fig.colorbar(im, ax=ax, fraction=0.05, pad=0.02)
    nums["E6_replay_participation_mean"] = {g: round(float(np.mean(p)), 3) for g, p in part.items()}
    # one event per epoch: raster and a voltage trace
    picks = []
    for ep in range(4):
        idx = [i for i in range(ep * 20, ep * 20 + 20) if d["events"][i]["replayed"]]
        picks.append(idx[0] if idx else ep * 20)
    sub = gs[2, :].subgridspec(1, 4)
    for j, i in enumerate(picks):
        ev = d["events"][i]
        ax = fig.add_subplot(sub[0, j])
        raster(ax, ev["spikes"], 195, 260)
        ax.set_xlim(195, 260)
        if j:
            ax.set_yticklabels([])
        ax.set_xlabel("ms in event")
        ax.set_title(f"({'defg'[j]}) event {i + 1}, {ev['t_s'] / 3600:.2f} h: "
                     f"{', '.join(ev['replayed']) or 'no pattern'}", fontsize=7)
    fig.tight_layout()
    fig.savefig(OUT / "fig_E6_activity.png", dpi=160)
    plt.close(fig)
    # voltage examples
    fig, axes = plt.subplots(2, 1, figsize=(9, 4.6))
    ax = axes[0]
    t = e["t"] / 1000.0
    VC = {"PYR 22 (X)": PCOL["X"], "PYR 23 (+)": PCOL["+"], "PYR 11 (=)": PCOL["="], "BC 22": "#555555"}
    for key, lab in ((("PYR", 22), "PYR 22 (X)"), (("PYR", 23), "PYR 23 (+)"), (("PYR", 11), "PYR 11 (=)"),
                     (("BC", 22), "BC 22")):
        n = min(len(t), len(e["v"][key]))
        ax.plot(t[:n], e["v"][key][:n], lw=0.5, label=lab, color=VC[lab], alpha=0.85)
    ax.set_xlim(0, e["t_enc_ms"] / 1000)
    ax.set_ylabel("V (mV)")
    ax.set_title("(a) membrane potential during encoding, seed 1", fontsize=8)
    ax.legend(frameon=False, fontsize=6, ncol=5)
    ax = axes[1]
    ev = d["events"][picks[0]]
    swt = next(i for i in range(10) if len(ev["spikes"][("SWT", i)]))
    VC[f"SWT {swt}"] = GCOL["SWT"]
    for key, lab in ((("PYR", 22), "PYR 22 (X)"), (("PYR", 23), "PYR 23 (+)"), (("PYR", 11), "PYR 11 (=)"),
                     (("BC", 22), "BC 22"), (("SWT", swt), f"SWT {swt}")):
        n = min(len(ev["t"]), len(ev["v"][key]))
        ax.plot(ev["t"][:n], ev["v"][key][:n], lw=0.8, label=lab, color=VC[lab])
    ax.set_xlabel("ms")
    ax.set_ylabel("V (mV)")
    ax.set_title(f"(b) membrane potential in event {picks[0] + 1} ({ev['t_s'] / 3600:.2f} h)", fontsize=8)
    ax.legend(frameon=False, fontsize=6, ncol=5)
    fig.tight_layout()
    fig.savefig(OUT / "fig_E6b_voltage.png", dpi=160)
    plt.close(fig)


def main():
    OUT.mkdir(exist_ok=True)
    D = {c: load(1, s) for c, s in TR.items()}
    nums = {}
    fig_whole(D, nums)
    zoom(D, 0.0, 1.05, "encoding to the end of epoch 1", "fig_E2_first_phase.png", enc_panel=True)
    zoom(D, 2.0, 4.0, "middle: 2-4 h, epoch 3 at 3.0-3.5 h", "fig_E3_middle.png")
    zoom(D, 5.5, 8.0, "end: 5.5-8 h, epoch 4 at 6.0-6.5 h", "fig_E4_end.png")
    fig_updown(D, nums)
    fig_activity(D, nums)
    d = D["replay"]
    nums["W_seed1"] = {c: {k: {f"{h:g} h": round(float(np.interp(h, *W_series(D[c], k))), 4)
                               for h in (0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 3.5, 4.0, 6.0, 6.5, 8.0)}
                           for k in CLASSES} for c in TR}
    nums["event_dW_by_epoch_x1e3"] = {k: [round(float(ev_dW(d, k)[i * 20:(i + 1) * 20].sum() * 1000), 3)
                                          for i in range(4)] for k in CLASSES}
    nums["replayed_by_epoch"] = [[ev["replayed"] for ev in d["events"][i * 20:(i + 1) * 20]] for i in range(4)]
    (OUT / "extended_numbers.json").write_text(json.dumps(nums, indent=1, ensure_ascii=False))
    print("saved", OUT)


if __name__ == "__main__":
    main()
