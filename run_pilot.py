"""Pilot: encoding of the six patterns, DA/5-HT, then 8 h of slow STC dynamics with four replay epochs.

Single seed, hand-calibrated weights (circuit.CircuitParams). A pilot: its numbers are not accepted results.
Saves results/pilot.npz (+ results/pilot_events.pkl) for make_figures.py.

Timeline (real time, no compression):
  0 s       online encoding block: Λ, V, I, O, Λ∪V, Λ∪V∪I in turn, each 3 x 100 ms at 100 Hz into EC_k, 0.4 s apart,
            next pattern 1 s later (L&T 2021 learning protocol); DA (Λ cells) and 5-HT (V cells) on for 60 min
            from the first encoding that contains Λ or V
  0.5, 1.5, 3, 6 h  replay epochs: 20 triggered SWR events each, 90 s apart (30 min)
  between   offline integration of Tt, z, P every 10 s
Usage: .venv/bin/python run_pilot.py
"""
import pickle
import time
from pathlib import Path

import numpy as np

from circuit import CELL_IDS, PATTERN_ORDER, PATTERNS, Circuit, CircuitParams, pattern_synapse_mask
from stc import ModPulse, STCParams, advance, modulator_at

OUT = Path(__file__).resolve().parent / "results"
T_END_S = 8 * 3600.0
EPOCHS_H = (0.5, 1.5, 3.0, 6.0)
N_EVENTS, EVENT_GAP_S = 20, 90.0
ENC_START_MS, ENC_STEP_MS = 300.0, 2300.0
EVENT_MS, SWR_AT_MS = 450.0, 200.0
REPLAY_FRACTION = 0.8          # pattern k counts as replayed if >= 80% of its cells spike in the event window


def main():
    OUT.mkdir(exist_ok=True)
    cp, sp = CircuitParams(), STCParams()
    circ = Circuit(cp)
    n = len(circ.p1)
    post_idx = circ.p1_post
    masks = {k: pattern_synapse_mask(circ, k) for k in PATTERN_ORDER}
    Tt, z, P = np.zeros(n), np.zeros(n), np.zeros(len(CELL_IDS))

    # --- encoding (one online block) ---------------------------------------------------------------------------
    enc_on = {k: ENC_START_MS + i * ENC_STEP_MS for i, k in enumerate(PATTERN_ORDER)}
    for k, t_on in enc_on.items():
        circ.schedule_encoding(k, t_on)
    t_enc_ms = ENC_START_MS + len(PATTERN_ORDER) * ENC_STEP_MS
    from cells import h
    tt_rec = [h.Vector() for _ in range(n)]
    for vec, (_, _, s) in zip(tt_rec, circ.p1):
        vec.record(s._ref_Tt, 5.0)
    t0 = time.time()
    circ.set_state(replay=False)                        # encoding state
    enc = circ.run_block(t_enc_ms, Tt, z, P, block=0)
    circ.clear_schedules()
    print(f"encoding block {t_enc_ms / 1000:.1f} s simulated in {time.time() - t0:.0f} s", flush=True)
    enc["Tt_trace"] = np.array([np.array(v) for v in tt_rec])
    for v in tt_rec:
        v.play_remove()
    Tt, z = enc["Tt"], enc["z"]
    t_s = t_enc_ms / 1000.0
    pulses = [ModPulse("DA", enc_on["Λ"] / 1000.0, 3600.0), ModPulse("5-HT", enc_on["V"] / 1000.0, 3600.0)]

    # --- slow dynamics and replay ------------------------------------------------------------------------------
    slow = dict(t=[], Tt=[], z=[], P=[], M_DA=[], M_5HT=[])

    def log(t, Tt_, z_, P_):
        if not slow["t"] or t - slow["t"][-1] >= 59.9:
            _, src = modulator_at(pulses, t)
            slow["t"].append(t)
            slow["Tt"].append({k: float(Tt_[m].mean()) for k, m in masks.items()})
            slow["z"].append({k: float(z_[m].mean()) for k, m in masks.items()})
            slow["P"].append(P_.copy())
            slow["M_DA"].append(src["DA"])
            slow["M_5HT"].append(src["5-HT"])

    log(t_s, Tt, z, P)
    snapshots = {"end of encoding": (t_s, (Tt + z).copy())}
    events = []
    block = 1
    event_times = [h_ * 3600.0 + i * EVENT_GAP_S for h_ in EPOCHS_H for i in range(N_EVENTS)]
    for ev_t in event_times + [T_END_S]:
        Tt, z, P = advance(Tt, z, P, post_idx, pulses, t_s, ev_t, sp, on_step=log)
        t_s = ev_t
        if ev_t >= T_END_S:
            break
        circ.set_state(replay=True)                     # replay state
        circ.schedule_swr(SWR_AT_MS, event=block)
        r = circ.run_block(EVENT_MS, Tt, z, P, block=block)
        circ.clear_schedules()
        block += 1
        Tt, z = r["Tt"], r["z"]
        win = lambda st: np.any((st > SWR_AT_MS - 5) & (st < EVENT_MS - 10))
        active = {cid for cid in CELL_IDS if win(r["spikes"][("PYR", cid)])}
        frac = {k: len(active & set(PATTERNS[k])) / len(PATTERNS[k]) for k in PATTERN_ORDER}
        replayed = [k for k in PATTERN_ORDER if frac[k] >= REPLAY_FRACTION]
        events.append(dict(t_s=t_s, frac=frac, replayed=replayed, n_active=len(active),
                           spikes={k: v for k, v in r["spikes"].items()}, v=r["v"], t=r["t"]))
        print(f"t={t_s / 3600:.2f} h  active {len(active):2d}  replayed {replayed}", flush=True)
        if abs(t_s - EPOCHS_H[1] * 3600.0 - (N_EVENTS - 1) * EVENT_GAP_S) < 1:
            snapshots["after epoch 2 (~1.98 h)"] = (t_s, (Tt + z).copy())
    log(t_s, Tt, z, P)
    snapshots["8 h"] = (t_s, (Tt + z).copy())

    np.savez_compressed(OUT / "pilot.npz", slow_t=np.array(slow["t"]), slow_P=np.array(slow["P"]),
                        slow_M_DA=np.array(slow["M_DA"]), slow_M_5HT=np.array(slow["M_5HT"]),
                        p1_pre=circ.p1_pre, p1_post=circ.p1_post)
    with open(OUT / "pilot_events.pkl", "wb") as f:
        pickle.dump(dict(enc=dict(t=enc["t"], spikes=enc["spikes"], v=enc["v"], Tt_trace=enc["Tt_trace"]),
                         enc_on=enc_on, t_enc_ms=t_enc_ms, slow_Tt=slow["Tt"], slow_z=slow["z"], events=events,
                         snapshots=snapshots, pulses=[p.__dict__ for p in pulses], params=cp.__dict__,
                         stc_params=sp.__dict__, masks={k: m for k, m in masks.items()}), f)
    print("saved", OUT / "pilot.npz", OUT / "pilot_events.pkl")


if __name__ == "__main__":
    main()
