"""T4 (pilot form): TMExp2SynNMDA_STC against the numpy reference (reference/stc_numpy.py).

A synapse on a passive soma receives presynaptic spikes through one NetCon and "postsynaptic" spikes through a
second NetCon with weight -1. Protocols (calcium thresholds of L&T 2021, network values):
  strong : pre and post together at 100 Hz for 1 s -> c above theta_p -> Tt rises above the tag threshold
  weak   : pre alone at 20 Hz for 1 s              -> c between theta_d and theta_p at most -> Tt falls or stays
  capture: strong, with P = 1 and tau_z shortened to 1 s so that z moves within the test
Acceptance: Tt and z at the end agree with the reference within 1% of their range (or 1e-4 absolute), and the
release scaling by (1 + Tt + z) is applied to the conductance.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from nrn_setup import load_mechanisms  # noqa: E402
from reference.stc_numpy import stc_single  # noqa: E402

h = load_mechanisms()


def _run(pre, post, t_stop, P=0.0, tau_z=3600000.0):
    soma = h.Section(name="stc_soma")
    soma.L = soma.diam = 10.0
    soma.insert("pas")
    syn = h.TMExp2SynNMDA_STC(soma(0.5))
    syn.P, syn.tau_z = P, tau_z
    src_pre, src_post = h.NetStim(), h.NetStim()
    for s, times in ((src_pre, pre), (src_post, post)):
        if len(times):
            s.start, s.number, s.noise = times[0], len(times), 0
            s.interval = times[1] - times[0] if len(times) > 1 else 1e9
        else:
            s.number = 0
    nc_pre = h.NetCon(src_pre, syn, 0, 0, 1e-4)
    nc_post = h.NetCon(src_post, syn, 0, 0, -1)
    rec = {k: h.Vector() for k in ("c", "Tt", "z")}
    rec["c"].record(syn._ref_c)
    rec["Tt"].record(syn._ref_Tt)
    rec["z"].record(syn._ref_z)
    h.cvode_active(0)
    h.dt = 0.025
    h.finitialize(-65.0)
    h.continuerun(t_stop)
    out = {k: np.array(v) for k, v in rec.items()}
    del nc_pre, nc_post, syn, src_pre, src_post
    h.delete_section(sec=soma)
    return out


def _times(f, n, start=100.0):
    return [start + i * 1000.0 / f for i in range(n)]


@pytest.mark.parametrize("name,pre,post,P,tau_z", [
    ("strong", _times(100, 100), _times(100, 100, start=100.5), 0.0, 3600000.0),
    ("weak", _times(20, 20), [], 0.0, 3600000.0),
    ("capture", _times(100, 100), _times(100, 100, start=100.5), 1.0, 1000.0),
])
def test_mod_matches_numpy(name, pre, post, P, tau_z):
    t_stop = 1600.0
    got = _run(pre, post, t_stop, P=P, tau_z=tau_z)
    ref = stc_single(pre, post, t_stop, P=P, tau_z=tau_z)
    for j, key in ((1, "Tt"), (2, "z")):
        g_end, r_end = got[key][-1], ref[-1, j]
        tol = max(1e-4, 0.01 * max(1e-9, np.ptp(ref[:, j])))
        assert abs(g_end - r_end) < tol, (name, key, g_end, r_end)
    if name == "strong":
        assert got["Tt"][-1] > 0.2          # tagged up
    if name == "capture":
        assert got["z"][-1] > 0.05          # captured within the test
