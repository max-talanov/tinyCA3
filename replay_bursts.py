"""Option B of docs/REPLAY_DEPRESSION_2026-11-27.md: can replay be made bursty enough to potentiate, without runaway?
Reruns only the 80 replay events of the seed-1 trace run (pilot_s1_mod10_gate30_trace.pkl) with the weights they
started from, for one setting of the replay state (CircuitParams.replay_inh_scale, replay_p1_scale).

Each event starts from the stored W_before = 1 + T̃ + z of every synapse, restored as T̃ = W_before - 1, z = 0, P = 0:
transmission is exactly as in the full run (each quantum is scaled by 1 + T̃ + z), the sign of the tag change is
exact, its size is approximate (the depression rate scales with 1 + T̃). Events are independent: weights are not
carried from one event to the next, so this measures what each setting does to the same starting state.
Per event: active cells, spikes per active cell, replayed patterns, and the change of W of every synapse.
Usage: TINYCA3_LAYOUT=B env -u PYTHONPATH .venv/bin/python replay_bursts.py INH P1     (one setting, writes JSON)
       TINYCA3_LAYOUT=B env -u PYTHONPATH .venv/bin/python replay_bursts.py --summary  (table and figure)
"""
import json
import sys

import numpy as np

from circuit import CELL_IDS, PATTERN_ORDER, PATTERNS, ROOT

RES = ROOT / "results" / "layout_B" / "replay_bursts"
SRC = "_mod10_gate30_trace"
EVENT_MS, SWR_AT_MS = 450.0, 200.0
REPLAY_FRACTION = 0.8
CLASSES = ("X", "+", "X→+", "+→X", "=", "‖", "=→‖", "‖→=", "none")


def run(inh, p1):
    from circuit import Circuit, CircuitParams
    from class_stats import load
    d = load(1, SRC)
    circ = Circuit(CircuitParams(seed=1, replay_inh_scale=inh, replay_p1_scale=p1))
    n = len(circ.p1)
    circ.set_state(replay=True)
    masks = {k: d["masks"]["syn:" + k] for k in CLASSES}
    out = []
    for block, e in enumerate(d["events"], start=1):
        Tt0 = e["W_before"].astype(float) - 1.0
        circ.schedule_swr(SWR_AT_MS, event=block)
        r = circ.run_block(EVENT_MS, Tt0, np.zeros(n), np.zeros(len(CELL_IDS)), block=block)
        circ.clear_schedules()
        dW = (r["Tt"] + r["z"]) - Tt0
        nsp = {c: int(np.sum((s > SWR_AT_MS - 5) & (s < EVENT_MS - 10)))
               for c, s in ((c, r["spikes"][("PYR", c)]) for c in CELL_IDS)}
        active = {c for c, k in nsp.items() if k}
        frac = {k: len(active & set(PATTERNS[k])) / len(PATTERNS[k]) for k in PATTERN_ORDER}
        bc = int(sum(np.sum((s > SWR_AT_MS - 5) & (s < EVENT_MS - 10))
                     for s in (r["spikes"][("BC", c)] for c in CELL_IDS)))
        out.append(dict(t_s=e["t_s"], n_active=len(active), active=sorted(active), spikes=nsp, bc_spikes=bc,
                        replayed=[k for k in PATTERN_ORDER if frac[k] >= REPLAY_FRACTION],
                        n_up=int(np.sum(dW > 1e-3)), n_down=int(np.sum(dW < -1e-3)),
                        dW_class={k: float(dW[m].mean()) for k, m in masks.items()},
                        dW_up_class={k: int(np.sum(dW[m] > 1e-3)) for k, m in masks.items()},
                        same_active_as_full_run=sorted(active) == list(e["active"])))
        print(f"inh {inh:g} p1 {p1:g} event {block:2d}: active {len(active):2d} up {out[-1]['n_up']:3d} "
              f"down {out[-1]['n_down']:3d} {out[-1]['replayed']}", flush=True)
    RES.mkdir(parents=True, exist_ok=True)
    (RES / f"inh{inh:g}_p1{p1:g}.json").write_text(json.dumps(out, ensure_ascii=False))


def summary():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rows = []
    for f in sorted(RES.glob("inh*_p1*.json")):
        inh, p1 = (float(x) for x in f.stem[3:].split("_p1"))
        ev = json.loads(f.read_text())
        sp = [k for e in ev for k in e["spikes"].values() if k]
        rows.append(dict(inh=inh, p1=p1,
                         spikes_per_active_cell=round(float(np.mean(sp)), 2) if sp else 0.0,
                         frac_cells_3plus=round(float(np.mean(np.array(sp) >= 3)), 3) if sp else 0.0,
                         mean_active=round(float(np.mean([e["n_active"] for e in ev])), 1),
                         runaway_events=sum(e["n_active"] >= 23 for e in ev),
                         replay_events=sum(bool(e["replayed"]) for e in ev),
                         single_pattern_events=sum(len([k for k in e["replayed"] if k in ("X", "+", "=", "‖")]) == 1
                                                   and e["n_active"] < 23 for e in ev),
                         syn_up=sum(e["n_up"] for e in ev), syn_down=sum(e["n_down"] for e in ev),
                         dW_sum={k: round(float(sum(e["dW_class"][k] for e in ev)), 4) for k in CLASSES},
                         same_as_full_run=sum(e["same_active_as_full_run"] for e in ev)))
    (RES / "summary.json").write_text(json.dumps(rows, indent=1, ensure_ascii=False))
    inhs = sorted({r["inh"] for r in rows})
    p1s = sorted({r["p1"] for r in rows})
    keys = [("spikes_per_active_cell", "spikes per active cell"), ("runaway_events", "events recruiting ≥ 23 of 25 cells"),
            ("syn_up", "synapse-events potentiated"), ("single_pattern_events", "events replaying one basic pattern")]
    fig, axes = plt.subplots(1, 4, figsize=(13, 3.2))
    for ax, (k, lab) in zip(axes, keys):
        g = np.full((len(inhs), len(p1s)), np.nan)
        for r in rows:
            g[inhs.index(r["inh"]), p1s.index(r["p1"])] = r[k]
        im = ax.imshow(g, origin="lower", cmap="viridis", aspect="auto")
        for i in range(len(inhs)):
            for j in range(len(p1s)):
                ax.text(j, i, f"{g[i, j]:g}", ha="center", va="center", fontsize=7,
                        color="w" if g[i, j] < np.nanmax(g) * 0.6 else "k")
        ax.set_xticks(range(len(p1s)), [f"{p:g}" for p in p1s])
        ax.set_yticks(range(len(inhs)), [f"{v:g}" for v in inhs])
        ax.set_xlabel("replay_p1_scale (recurrent)")
        ax.set_ylabel("replay_inh_scale (basket → PYR)")
        ax.set_title(lab, fontsize=8)
        ax.plot(p1s.index(2.5), inhs.index(2.0), "s", mfc="none", mec="r", ms=16)
    fig.suptitle("Replay of seed 1 (80 events, same starting weights) under different replay states; red: current",
                 fontsize=8)
    fig.tight_layout()
    fig.savefig(RES / "fig_replay_bursts.png", dpi=160)
    for r in rows:
        print({k: v for k, v in r.items() if k != "dW_sum"}, r["dW_sum"]["X"], r["dW_sum"]["+"], r["dW_sum"]["="])


if __name__ == "__main__":
    if sys.argv[1] == "--summary":
        summary()
    else:
        run(float(sys.argv[1]), float(sys.argv[2]))
