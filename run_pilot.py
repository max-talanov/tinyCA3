"""Pilot: encoding of the six patterns, DA/5-HT, then 8 h of slow STC dynamics with four replay epochs.

Single seed, hand-calibrated weights (circuit.CircuitParams). A pilot: its numbers are not accepted results.
Saves results/pilot.npz (+ results/pilot_events.pkl) for make_figures.py.

Timeline (real time, no compression):
  0 s       online encoding block: the six patterns in turn (layout B: X, +, =, ‖, X∪+, =∪‖; v1: Λ, V, I, O,
            Λ∪V, Λ∪V∪I), each 3 x 100 ms at 100 Hz into EC_k, 0.4 s apart,
            next pattern 1 s later (L&T 2021 learning protocol); DA (X cells) and 5-HT (+ cells) on for 60 min
            from the encoding of X and +
  0.5, 1.5, 3, 6 h  replay epochs: 20 triggered SWR events each, 90 s apart (30 min)
  between   offline integration of Tt, z, P every 10 s
Usage: .venv/bin/python run_pilot.py                      (seed 1, full output for make_figures.py)
       .venv/bin/python run_pilot.py --seed 3 --light     (results/seeds/pilot_s3.pkl, no voltages or spikes)
       .venv/bin/python run_pilot.py --seed 3 --light --no-replay   (paired control: no SWR events, same slow dynamics)
       .venv/bin/python run_pilot.py --seed 3 --light --no-stc   (STP only: P1 plastic = 0, no tags, no capture;
            results/seeds/pilot_s3_nostc.pkl)
       .venv/bin/python run_pilot.py --seed 3 --light --prp-block   (anisomycin analogue: alpha_P = 0, no PRPs, so
            no capture; check 1 of the algorithm doc section 4, T7 (b) of the plan; results/seeds/pilot_s3_prpblock.pkl)
       .venv/bin/python run_pilot.py --seed 3 --light --mod-min 10   (DA and 5-HT on for 10 min instead of 60;
            results/seeds/pilot_s3_mod10.pkl)
       .venv/bin/python run_pilot.py --seed 3 --light --mod-min 10 --gate 30   (DA and 5-HT for 10 min switch on a
            synthesis trigger that decays with 30 min; results/seeds/pilot_s3_mod10_gate30.pkl)
       .venv/bin/python run_pilot.py --seed 3 --light --replay-da 10   (VTA mode: an event that replays the DA pattern
            releases DA onto its cells for 10 min; algorithm doc section 2.6, McNamara et al. 2014 via S&A Ch. 14)
The pattern layout is set by TINYCA3_LAYOUT (circuit.py, default B); output goes to circuit.RESULTS.
The seed sets the background noise, the sparse trigger wiring (P8) and the trigger subsets of every event.
"""
import argparse
import pickle
import time
from pathlib import Path

import numpy as np

from circuit import (CELL_IDS, LAYOUT, MOD_PATTERN, PATTERN_ORDER, PATTERNS, PROBES, RESULTS, Circuit, CircuitParams,
                     class_synapse_masks, pattern_synapse_mask)
from stc import ModPulse, STCParams, advance, modulator_at

OUT = RESULTS
T_END_S = 8 * 3600.0
EPOCHS_H = (0.5, 1.5, 3.0, 6.0)
N_EVENTS, EVENT_GAP_S = 20, 90.0
ENC_START_MS, ENC_STEP_MS = 300.0, 2300.0
EVENT_MS, SWR_AT_MS = 450.0, 200.0
# Probe synapses (pre, post): release, tag and captured weight recorded every 0.5 ms in every online block, so the
# per-spike efficacy (release / U) * (1 + Tt + z) can be read out (comparison report of 2026-11-20). PROBES comes
# from the layout in circuit.py: one probe per synapse class.
REPLAY_FRACTION = 0.8          # pattern k counts as replayed if >= 80% of its cells spike in the event window


def compact_probes(t0_s, kind, t, rec, every=20):
    """Probe traces of one online block, compacted: the values at each presynaptic spike (spk_*) and Tt, z every
    `every` samples (10 ms). Sample 0 is taken inside finitialize, before run_block restores Tt and z, so it is
    dropped."""
    out = dict(t0_s=t0_s, kind=kind, t=t[1::every].astype(np.float32))
    for name, r in rec.items():
        idx = np.nonzero(np.diff(r["n_pre"]) > 0)[0] + 1
        out[name] = dict(spk_t=t[idx], spk_rel=r["last_rel"][idx], spk_Tt=r["Tt"][idx], spk_z=r["z"][idx],
                         Tt=r["Tt"][1::every].astype(np.float32), z=r["z"][1::every].astype(np.float32))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--light", action="store_true", help="store events without spikes and voltages")
    ap.add_argument("--no-replay", action="store_true", help="skip the SWR events (paired control)")
    ap.add_argument("--replay-da", type=float, default=0.0,
                    help="minutes of DA onto the DA-pattern cells after an event that replays that pattern (0 = off)")
    ap.add_argument("--replay-da-amp", type=float, default=1.0)
    ap.add_argument("--no-encoding-da", action="store_true",
                    help="no DA pulse at encoding (5-HT unchanged): the DA pattern gets DA only from replay (experiment X4)")
    ap.add_argument("--mod-min", type=float, default=60.0,
                    help="minutes of DA and 5-HT from the encoding of their patterns (default 60, Lehr 2022)")
    ap.add_argument("--gate", type=float, default=0.0,
                    help="synthesis trigger: the modulator switches on a per-cell synthesis state that decays with "
                         "this time constant (min); 0 = off, theta_pro follows the modulator (L&T 2021)")
    ap.add_argument("--prp-block", action="store_true",
                    help="PRP synthesis blocked (alpha_P = 0): tags are set but never captured")
    ap.add_argument("--no-stc", action="store_true",
                    help="STP only: P1 synapses with plastic = 0 (no tag, no capture), everything else unchanged")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    cp, sp = CircuitParams(seed=args.seed), STCParams(alpha_P=0.0 if args.prp_block else 1.0,
                                                         tau_G_min=args.gate)
    circ = Circuit(cp)
    n = len(circ.p1)
    post_idx = circ.p1_post
    masks = {k: pattern_synapse_mask(circ, k) for k in PATTERN_ORDER}
    masks.update(class_synapse_masks(circ))            # layout B: per synapse class, keys 'syn:<class>'
    if args.no_stc:
        for _, _, s in circ.p1:
            s.plastic = 0
    from cells import h
    syn_of = {(a, b): s for a, b, s in circ.p1}
    probe_rec = {}
    for name, ab in PROBES.items():
        s = syn_of[ab]
        probe_rec[name] = {}
        for var in ("last_rel", "Tt", "z", "n_pre"):
            v = h.Vector()
            v.record(getattr(s, f"_ref_{var}"), 0.5)
            probe_rec[name][var] = v
    probes = []          # one dict per online block: t (ms in block), t0_s (block start), and the probe traces

    def grab_probes(t0_s, kind):
        probes.append(compact_probes(t0_s, kind, np.array(circ.trec),
                                     {name: {var: np.array(v) for var, v in rec.items()}
                                      for name, rec in probe_rec.items()}))
    Tt, z, P = np.zeros(n), np.zeros(n), np.zeros(len(CELL_IDS))

    # --- encoding (one online block) ---------------------------------------------------------------------------
    enc_on = {k: ENC_START_MS + i * ENC_STEP_MS for i, k in enumerate(PATTERN_ORDER)}
    for k, t_on in enc_on.items():
        circ.schedule_encoding(k, t_on)
    t_enc_ms = ENC_START_MS + len(PATTERN_ORDER) * ENC_STEP_MS
    tt_rec = [h.Vector() for _ in range(n)]
    for vec, (_, _, s) in zip(tt_rec, circ.p1):
        vec.record(s._ref_Tt, 5.0)
    t0 = time.time()
    circ.set_state(replay=False)                        # encoding state
    enc = circ.run_block(t_enc_ms, Tt, z, P, block=0)
    circ.clear_schedules()
    print(f"encoding block {t_enc_ms / 1000:.1f} s simulated in {time.time() - t0:.0f} s", flush=True)
    grab_probes(0.0, "encoding")
    enc["Tt_trace"] = np.array([np.array(v) for v in tt_rec])
    for v in tt_rec:
        v.play_remove()
    Tt, z = enc["Tt"], enc["z"]
    t_s = t_enc_ms / 1000.0
    pulses = [ModPulse(m, enc_on[k] / 1000.0, 60.0 * args.mod_min) for m, k in MOD_PATTERN.items()]
    if args.no_encoding_da:
        pulses = [p for p in pulses if p.source != "DA"]

    # --- slow dynamics and replay ------------------------------------------------------------------------------
    slow = dict(t=[], Tt=[], z=[], P=[], G=[], M_DA=[], M_5HT=[])
    G = None

    def log(t, Tt_, z_, P_, G_=None):
        if not slow["t"] or t - slow["t"][-1] >= 59.9:
            _, src = modulator_at(pulses, t)
            slow["t"].append(t)
            slow["Tt"].append({k: float(Tt_[m].mean()) for k, m in masks.items()})
            slow["z"].append({k: float(z_[m].mean()) for k, m in masks.items()})
            slow["P"].append(P_.copy())
            slow["G"].append(np.zeros_like(P_) if G_ is None else G_.copy())
            slow["M_DA"].append(src["DA"])
            slow["M_5HT"].append(src["5-HT"])

    log(t_s, Tt, z, P)
    snapshots = {"end of encoding": (t_s, (Tt + z).copy())}
    events = []
    block = 1
    event_times = [h_ * 3600.0 + i * EVENT_GAP_S for h_ in EPOCHS_H for i in range(N_EVENTS)]
    for ev_t in event_times + [T_END_S]:
        Tt, z, P, G = advance(Tt, z, P, post_idx, pulses, t_s, ev_t, sp, on_step=log, G=G)
        t_s = ev_t
        if ev_t >= T_END_S:
            break
        if args.no_replay:
            continue
        circ.set_state(replay=True)                     # replay state
        circ.schedule_swr(SWR_AT_MS, event=block)
        r = circ.run_block(EVENT_MS, Tt, z, P, block=block)
        grab_probes(t_s, "replay")
        circ.clear_schedules()
        block += 1
        Tt, z = r["Tt"], r["z"]
        win = lambda st: np.any((st > SWR_AT_MS - 5) & (st < EVENT_MS - 10))
        active = {cid for cid in CELL_IDS if win(r["spikes"][("PYR", cid)])}
        frac = {k: len(active & set(PATTERNS[k])) / len(PATTERNS[k]) for k in PATTERN_ORDER}
        replayed = [k for k in PATTERN_ORDER if frac[k] >= REPLAY_FRACTION]
        ev = dict(t_s=t_s, frac=frac, replayed=replayed, n_active=len(active), active=sorted(active))
        if not args.light:
            ev.update(spikes={k: v for k, v in r["spikes"].items()}, v=r["v"], t=r["t"])
        events.append(ev)
        if args.replay_da > 0 and MOD_PATTERN["DA"] in replayed and len(active) < 23:
            pulses.append(ModPulse("DA", t_s, 60.0 * args.replay_da, args.replay_da_amp))
            ev["da"] = True
        print(f"t={t_s / 3600:.2f} h  active {len(active):2d}  replayed {replayed}", flush=True)
        if abs(t_s - EPOCHS_H[1] * 3600.0 - (N_EVENTS - 1) * EVENT_GAP_S) < 1:
            snapshots["after epoch 2 (~1.98 h)"] = (t_s, (Tt + z).copy())
    log(t_s, Tt, z, P, G)
    snapshots["8 h"] = (t_s, (Tt + z).copy())

    if args.light:
        (OUT / "seeds").mkdir(exist_ok=True)
        suffix = "_norep" if args.no_replay else (f"_rda{args.replay_da:g}" if args.replay_da > 0 else "")
        if args.no_encoding_da:
            suffix = "_noencda" + suffix
        if args.no_stc:
            suffix = "_nostc" + suffix
        if args.prp_block:
            suffix = "_prpblock" + suffix
        if args.gate > 0:
            suffix = f"_gate{args.gate:g}" + suffix
        if args.mod_min != 60.0:
            suffix = f"_mod{args.mod_min:g}" + suffix
        name = OUT / "seeds" / f"pilot_s{args.seed}{suffix}.pkl"
        with open(name, "wb") as f:
            pickle.dump(dict(layout=LAYOUT, seed=args.seed, no_replay=args.no_replay, replay_da=args.replay_da,
                             pulses=[p.__dict__ for p in pulses], slow_t=np.array(slow["t"]),
                             slow_Tt=slow["Tt"], slow_z=slow["z"], slow_P=np.array(slow["P"]),
                             slow_G=np.array(slow["G"]), gate=args.gate, events=events,
                             z_final=z, Tt_final=Tt, snapshots=snapshots, masks=masks,
                             p1_pre=circ.p1_pre, p1_post=circ.p1_post, params=cp.__dict__, no_stc=args.no_stc, prp_block=args.prp_block, mod_min=args.mod_min,
                             probes=probes, probe_pairs=PROBES, U=float(circ.p1[0][2].U)), f)
        print("saved", name)
        return
    np.savez_compressed(OUT / "pilot.npz", slow_t=np.array(slow["t"]), slow_P=np.array(slow["P"]),
                        slow_M_DA=np.array(slow["M_DA"]), slow_M_5HT=np.array(slow["M_5HT"]),
                        p1_pre=circ.p1_pre, p1_post=circ.p1_post)
    with open(OUT / "pilot_events.pkl", "wb") as f:
        pickle.dump(dict(layout=LAYOUT, enc=dict(t=enc["t"], spikes=enc["spikes"], v=enc["v"], Tt_trace=enc["Tt_trace"]),
                         enc_on=enc_on, t_enc_ms=t_enc_ms, slow_Tt=slow["Tt"], slow_z=slow["z"], events=events,
                         snapshots=snapshots, pulses=[p.__dict__ for p in pulses], params=cp.__dict__,
                         stc_params=sp.__dict__, masks={k: m for k, m in masks.items()}), f)
    print("saved", OUT / "pilot.npz", OUT / "pilot_events.pkl")


if __name__ == "__main__":
    main()
