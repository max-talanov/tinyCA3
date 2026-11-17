"""T1: the three Hodgkin-Huxley-type cell types.

1. The basket-cell port (cells.build_bc) reproduces the original Bezaire class_pvbasketcell.hoc: soma voltage within
   0.01 mV and equal spike counts for a hyperpolarising and a suprathreshold step (34 C, dt 0.025 ms).
2. Every cell type has a DC rheobase below 2 nA and fires repetitively above it.
3. The basket cell fires > 100 Hz with little accommodation (last/first ISI < 1.2) at 5x rheobase.
4. An RS unit (EC, DG, SWT) fires exactly one spike to a single AMPA-only synaptic event over a weight window of at
   least 3x (here 0.003-0.016 uS-equivalent NetCon weight); with NMDA it fires trains, so these inputs are AMPA-only.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from cells import (build_bc, build_pyr, build_rs, f_i_curve, h, rheobase,  # noqa: E402
                   spike_times)


def _clear():
    for sec in list(h.allsec()):
        h.delete_section(sec=sec)


def _step(sec, amp, t_stop=400.0):
    st = h.IClamp(sec(0.5))
    st.delay, st.dur, st.amp = 100.0, 200.0, amp
    v = h.Vector()
    v.record(sec(0.5)._ref_v)
    h.cvode_active(0)
    h.dt = 0.025
    h.celsius = 34.0
    h.finitialize(-65.0)
    h.continuerun(t_stop)
    return np.array(v)


@pytest.mark.parametrize("amp_na", [-0.05, 0.5])
def test_bc_port_matches_template(amp_na):
    if not hasattr(h, "numCellTypes"):
        h("numCellTypes = 0")
        h("objref cellType[1]")
    if not hasattr(h, "pvbasketcell"):
        h.load_file(str(ROOT / "reference" / "bezaire2016" / "class_pvbasketcell.hoc"))
    ref = h.pvbasketcell(0, 0, 0)
    v_ref = _step(ref.soma[0], amp_na)
    del ref
    _clear()
    port = build_bc()
    v_port = _step(port["soma"], amp_na)
    del port
    _clear()
    n = lambda v: int(np.sum((v[1:] > 0) & (v[:-1] <= 0)))
    assert n(v_port) == n(v_ref)
    assert np.max(np.abs(v_port - v_ref)) < 0.01


@pytest.mark.parametrize("builder", [build_pyr, build_bc, build_rs])
def test_rheobase_and_repetitive_firing(builder):
    c = builder()
    rb = rheobase(c["soma"], tol=1e-3)
    assert rb is not None and 0 < rb < 2.0
    rows = f_i_curve(c["soma"], [2 * rb])
    assert rows[0]["n_spikes"] >= 5
    del c
    _clear()


def test_bc_fast_spiking():
    c = build_bc()
    rb = rheobase(c["soma"], tol=1e-3)
    row = f_i_curve(c["soma"], [5 * rb])[0]
    assert row["rate_hz"] > 100.0
    assert row["adaptation"] < 1.2
    del c
    _clear()


def _spikes_per_event(cell, weight, nmda=False):
    s = cell["soma"]
    syn = h.TMExp2SynNMDACA1_SCExc(s(0.5))
    if not nmda:
        syn.gmax_nmda = 0
    ns = h.NetStim()
    ns.start, ns.number, ns.noise = 200.0, 1, 0
    nc = h.NetCon(ns, syn, 0, 0, weight)
    tv, vv = h.Vector(), h.Vector()
    tv.record(h._ref_t)
    vv.record(s(0.5)._ref_v)
    h.dt = 0.025
    h.finitialize(-66.0)
    h.continuerun(400.0)
    n = len(spike_times(np.array(tv), np.array(vv)))
    del nc, ns, syn
    return n


@pytest.mark.parametrize("weight", [0.003, 0.006, 0.016])
def test_rs_single_spike_window(weight):
    c = build_rs()
    assert _spikes_per_event(c, weight) == 1
    del c
    _clear()
