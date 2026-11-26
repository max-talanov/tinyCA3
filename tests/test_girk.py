"""5-HT1A -> Kir3 (GIRK) current of the pyramidal cells (neuron/ch_GIRK.mod; Johnston et al. 2014)."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from cells import PyrParams, build_bc, build_pyr, h  # noqa: E402

G_CAL = 6e-5        # S/cm2, the value used with run_pilot.py --ht-girk (4.3 mV, Johnston et al. 2014)


def _rest_and_spikes(cell, g, amp):
    for key in ("soma", "dend"):
        for seg in cell[key]:
            seg.ch_GIRK.gmax = g
    ic = h.IClamp(cell["soma"](0.5))  # noqa: F841 (kept alive for the run)
    ic.delay, ic.dur, ic.amp = 0.0, 1500.0, amp
    v, t = h.Vector().record(cell["soma"](0.5)._ref_v), h.Vector().record(h._ref_t)
    h.dt = 0.025
    h.finitialize(-66.0)
    h.continuerun(1500.0)
    v, t = np.array(v), np.array(t)
    spikes = int(np.sum((v[1:] > 0) & (v[:-1] <= 0) & (t[1:] > 500)))
    return float(v[(t > 1000) & (t < 1500)].mean()), spikes


def test_girk_off_by_default_and_only_in_pyramidal_cells():
    pyr = build_pyr(PyrParams(k_scale=2.0), name="tg_pyr")
    assert all(seg.ch_GIRK.gmax == 0.0 for key in ("soma", "dend") for seg in pyr[key])
    bc = build_bc(name="tg_bc")
    sec = next(iter(bc.values())) if isinstance(bc, dict) else bc
    assert not h.ismembrane("ch_GIRK", sec=sec)


def test_girk_matches_measured_hyperpolarization_and_lowers_firing():
    pyr = build_pyr(PyrParams(k_scale=2.0), name="tg_pyr2")
    v0, _ = _rest_and_spikes(pyr, 0.0, 0.08)        # held near -45 mV, subthreshold
    v1, _ = _rest_and_spikes(pyr, G_CAL, 0.08)
    assert 3.5 < v0 - v1 < 5.5                      # measured 4.3 +/- 1.6 mV (Johnston et al. 2014, Fig. 4C)
    _, s0 = _rest_and_spikes(pyr, 0.0, 0.1)
    _, s1 = _rest_and_spikes(pyr, G_CAL, 0.1)
    assert s1 < s0
