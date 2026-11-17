"""T0: environment and the existing STP mechanism.

Acceptance: NEURON imports from the venv, every mechanism in neuron/ loads, and the release per spike of
TMExp2SynNMDACA1_SCExc at 10, 20, 40 and 80 Hz matches the numpy reference (reference/tm_numpy.py) to 1e-9.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from nrn_setup import load_mechanisms  # noqa: E402
from reference.tm_numpy import tm_release  # noqa: E402

h = load_mechanisms()

MECHS = ["TMExp2SynNMDACA1_SCExc", "ch_Navp", "ch_Navaxonp", "ch_Kdrp", "ch_KvAproxp", "ch_KvAdistp", "ch_HCNp",
         "ch_Kdrfast", "ch_KvA", "ch_CavN", "ch_CavL", "ch_KCaS", "ch_KvCaB", "ch_leak"]


def test_mechanisms_loaded():
    mt_point, mt_density = h.MechanismType(1), h.MechanismType(0)
    names = set()
    for mt in (mt_point, mt_density):
        s = h.ref("")
        for i in range(int(mt.count())):
            mt.select(i)
            mt.selected(s)
            names.add(s[0])
    missing = [m for m in MECHS if m not in names]
    assert not missing, missing


@pytest.mark.parametrize("freq_hz", [10.0, 20.0, 40.0, 80.0])
def test_tm_release_matches_numpy(freq_hz):
    soma = h.Section(name="tm_soma")
    soma.L = soma.diam = 5.0
    soma.insert("pas")
    syn = h.TMExp2SynNMDACA1_SCExc(soma(0.5))
    isi = 1000.0 / freq_hz
    times = [200.0 + i * isi for i in range(10)]
    vs = h.NetStim()
    vs.start, vs.interval, vs.number, vs.noise = times[0], isi, len(times), 0
    nc = h.NetCon(vs, syn, 0, 0, 1e-4)
    h.cvode_active(0)
    h.dt = 0.025
    h.celsius = 34.0
    h.finitialize(-65.0)
    got = []
    for t in times:
        h.continuerun(t + 0.05)
        got.append(syn.last_rel)
    ref = tm_release(times, syn.U, syn.tau_f, syn.tau_d)
    assert np.max(np.abs(np.array(got) - ref)) < 1e-9
    del nc, syn, vs
    h.delete_section(sec=soma)
