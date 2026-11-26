"""5-HT1A -> Kir3 (GIRK) current of the pyramidal cells (neuron/ch_GIRK.mod; Johnston et al. 2014)."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from cells import PyrParams, build_bc, build_pyr, h  # noqa: E402

G_CAL = 1.5e-4      # S/cm2, the value used with run_pilot.py --ht-girk


def _rest_and_spikes(cell, g, amp):
    for key in ("soma", "dend"):
        for seg in cell[key]:
            seg.ch_GIRK.gmax = g
    ic = h.IClamp(cell["soma"](0.5))
    ic.delay, ic.dur, ic.amp = 500.0, 1000.0, amp
    v, t = h.Vector().record(cell["soma"](0.5)._ref_v), h.Vector().record(h._ref_t)
    h.dt = 0.025
    h.finitialize(-66.0)
    h.continuerun(1500.0)
    v, t = np.array(v), np.array(t)
    spikes = int(np.sum((v[1:] > 0) & (v[:-1] <= 0) & (t[1:] > 500)))
    return float(v[(t > 400) & (t < 500)].mean()), spikes


def test_girk_off_by_default_and_only_in_pyramidal_cells():
    pyr = build_pyr(PyrParams(k_scale=2.0), name="tg_pyr")
    assert all(seg.ch_GIRK.gmax == 0.0 for key in ("soma", "dend") for seg in pyr[key])
    bc = build_bc(name="tg_bc")
    sec = next(iter(bc.values())) if isinstance(bc, dict) else bc
    assert not h.ismembrane("ch_GIRK", sec=sec)


def test_girk_hyperpolarizes_and_lowers_firing():
    pyr = build_pyr(PyrParams(k_scale=2.0), name="tg_pyr2")
    r0, s0 = _rest_and_spikes(pyr, 0.0, 0.1)
    r1, s1 = _rest_and_spikes(pyr, G_CAL, 0.1)
    assert 4.0 < r0 - r1 < 8.0                      # calibrated to ~6 mV
    assert s1 < s0
