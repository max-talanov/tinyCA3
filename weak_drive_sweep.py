"""Encoding-only sweep for docs/SEROTONIN_WEAK_DRIVE_2026-11-27.md: the mossy-fibre weight onto the + (5-HT) cells is
scaled by s (CircuitParams.mf_scale_5ht), with the 5-HT1A -> Kir3 (GIRK) current off and on (6e-5 S/cm2 in the + cells
from the 5-HT onset). Seed 1, layout B. For each run: firing of the + cells during the bursts of + and X∪+, and the mean
tag per synapse class at the end of encoding.
Writes results/layout_B/weak_drive_sweep.json.
Usage: TINYCA3_LAYOUT=B env -u PYTHONPATH .venv/bin/python weak_drive_sweep.py            (all runs, 8 in parallel)
       TINYCA3_LAYOUT=B env -u PYTHONPATH .venv/bin/python weak_drive_sweep.py s g        (one run, prints JSON)
"""
import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

import numpy as np

from circuit import ROOT

SCALES = (1.0, 0.5, 0.3, 0.2, 0.15, 0.1, 0.05)
G = (0.0, 6e-5)
OUT = ROOT / "results" / "layout_B" / "weak_drive_sweep.json"


def one(scale, g):
    from cells import h
    from circuit import (CELL_IDS, MOD_PATTERN, PATTERN_ORDER, PATTERNS, Circuit,
                         CircuitParams, class_synapse_masks)
    from run_pilot import ENC_START_MS, ENC_STEP_MS
    from stc import MOD_CELLS
    circ = Circuit(CircuitParams(seed=1, mf_scale_5ht=scale))
    n = len(circ.p1)
    enc_on = {k: ENC_START_MS + i * ENC_STEP_MS for i, k in enumerate(PATTERN_ORDER)}
    for k, t in enc_on.items():
        circ.schedule_encoding(k, t)
    keep = []
    if g > 0:
        for cid in MOD_CELLS["5-HT"]:
            for key in ("soma", "dend"):
                for seg in circ.pyr[cid][key]:
                    tv, gv = h.Vector([0.0, enc_on[MOD_PATTERN["5-HT"]]]), h.Vector([0.0, g])
                    gv.play(seg._ref_gmax_ch_GIRK, tv, 0)
                    keep.append((tv, gv))
    circ.set_state(replay=False)
    r = circ.run_block(ENC_START_MS + len(PATTERN_ORDER) * ENC_STEP_MS, np.zeros(n), np.zeros(n),
                       np.zeros(len(CELL_IDS)), block=0)
    rates = {}
    for k in PATTERN_ORDER:                 # 3 bursts of 100 ms in the 1.3 s after the onset
        for grp in ("X", "+"):
            cells = [c for c in PATTERNS[grp] if c in PATTERNS[k]]
            if cells:
                sp = [r["spikes"][("PYR", c)] for c in cells]
                rates[f"{k}:{grp}"] = float(np.mean([np.sum((s >= enc_on[k]) & (s < enc_on[k] + 1300.0)) for s in sp]) / 0.3)
    return dict(scale=scale, g=g, rate_Hz=rates,
                Tt={k[4:]: float(r["Tt"][m].mean()) for k, m in class_synapse_masks(circ).items()})


def main():
    if len(sys.argv) == 3:
        print(json.dumps(one(float(sys.argv[1]), float(sys.argv[2]))))
        return

    def run(sg):
        p = subprocess.run([sys.executable, __file__, str(sg[0]), str(sg[1])], capture_output=True, text=True, check=True)
        return json.loads(p.stdout.strip().splitlines()[-1])
    with ThreadPoolExecutor(8) as ex:
        res = list(ex.map(run, [(s, g) for s in SCALES for g in G]))
    OUT.write_text(json.dumps(res, indent=1, ensure_ascii=False))
    print("saved", OUT)


if __name__ == "__main__":
    main()
