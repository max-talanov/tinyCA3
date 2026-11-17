"""T1: Hodgkin-Huxley-type cells of tinyCA3. No integrate-and-fire or artificial cells.

- PYR  (CA3 pyramidal cell): soma + one apical dendrite + axon initial segment, Bezaire et al. 2016 pyramidal channels
       (ch_Navp, ch_Kdrp, ch_KvAproxp, ch_HCNp; ch_Navaxonp in the axon). Densities are the template's somatic values
       (params.BEZ_PYR_SOMA). Status: CA1 model values used as placeholders for CA3.
- BC   (CA3 basket cell): the published Bezaire et al. 2016 PV basket cell (morphology and channels), ported from
       class_pvbasketcell.hoc and checked against it (tests/test_t1_cells.py). Status: CA1 model, placeholder for CA3.
- RS   (regular-spiking unit for EC_k, DG_k and SWT): one compartment with the pyramidal somatic channel set.
       Status: placeholder; these are input/readout units, not models of entorhinal or granule cells.

Rheobase is measured with injected DC current (`rheobase`), never with synaptic drive.
"""
import math
from dataclasses import dataclass

import numpy as np

from nrn_setup import load_mechanisms
from params import BEZ_PYR_SOMA, CELSIUS, Geometry

h = load_mechanisms()

_G = Geometry()


def _dlambda_nseg(L, diam, Ra, cm, f=100.0, d_lambda=0.1):
    lam = 1e5 * math.sqrt(diam / (4 * math.pi * f * Ra * cm))
    return int((L / (d_lambda * lam) + 0.9) / 2) * 2 + 1


# ------------------------------------------------------------------------------------------------------------------
# Pyramidal cell
# ------------------------------------------------------------------------------------------------------------------
@dataclass(frozen=True)
class PyrParams:
    """Scale factors on the Bezaire somatic densities (1 = template value) and the reduced geometry."""
    na_scale: float = 1.0
    k_scale: float = 1.0
    ka_scale: float = 1.0
    h_scale: float = 1.0
    leak_scale: float = 1.0
    e_pas: float = BEZ_PYR_SOMA["v_rest"]
    soma_L: float = _G.pyr_soma_L
    soma_diam: float = _G.pyr_soma_diam
    dend_L: float = _G.pyr_dend_L
    dend_diam: float = _G.pyr_dend_diam
    axon_L: float = _G.pyr_axon_L
    axon_diam: float = _G.pyr_axon_diam


def build_pyr(p: PyrParams = PyrParams(), name="pyr"):
    """Soma, apical dendrite (from soma(1)) and axon initial segment (from soma(0)).

    Channels follow the Bezaire template's rules: soma and dendrite carry Navp, Kdrp, KvAproxp and HCNp (dendritic
    Navp with ar2 = 0.8 as in its apical tree; HCN and KvA at their somatic densities, i.e. without the template's
    distance gradients); the axon carries Navaxonp, Kdrp and KvAproxp at 0.2x; dendrites have Rm/2 and 2*Cm.
    Returns {"soma", "dend", "axon"}. Keep the dict: NEURON deletes unreferenced sections."""
    b = BEZ_PYR_SOMA
    h.celsius = CELSIUS
    soma = h.Section(name=f"{name}_soma")
    dend = h.Section(name=f"{name}_dend")
    axon = h.Section(name=f"{name}_axon")
    soma.L, soma.diam = p.soma_L, p.soma_diam
    dend.L, dend.diam = p.dend_L, p.dend_diam
    axon.L, axon.diam = p.axon_L, p.axon_diam
    dend.connect(soma(1))
    axon.connect(soma(0))
    soma.Ra, soma.cm = b["Ra_soma"], b["Cm"]
    dend.Ra, dend.cm = b["Ra_soma"], 2 * b["Cm"]
    axon.Ra, axon.cm = b["Ra_axon"], b["Cm"]
    for sec in (soma, dend, axon):
        sec.nseg = _dlambda_nseg(sec.L, sec.diam, sec.Ra, sec.cm)
        sec.insert("pas")
    for sec, rm in ((soma, b["Rm"]), (dend, b["Rm"] / 2), (axon, b["Rm"])):
        for seg in sec:
            seg.pas.g, seg.pas.e = p.leak_scale / rm, p.e_pas
    for sec in (soma, dend):
        for mech in ("ch_Navp", "ch_Kdrp", "ch_KvAproxp", "ch_HCNp"):
            sec.insert(mech)
        for seg in sec:
            seg.ch_Navp.gmax = b["gNav"] * p.na_scale
            seg.ch_Kdrp.gmax = b["gKdr"] * p.k_scale
            seg.ch_KvAproxp.gmax = b["gKvAprox"] * p.ka_scale
            seg.ch_HCNp.gmax = b["gHCN"] * p.h_scale
            seg.ch_HCNp.e = b["e_hcn"]
    for seg in dend:
        seg.ch_Navp.ar2 = 0.8
    for mech in ("ch_Navaxonp", "ch_Kdrp", "ch_KvAproxp"):
        axon.insert(mech)
    for seg in axon:
        seg.ch_Navaxonp.gmax = b["gNavaxon"] * p.na_scale
        seg.ch_Kdrp.gmax = b["gKdr"] * p.k_scale
        seg.ch_KvAproxp.gmax = b["gKvAprox"] * 0.2 * p.ka_scale
    for sec in (soma, dend, axon):
        sec.ena, sec.ek = b["ena"], b["ek"]
    return {"soma": soma, "dend": dend, "axon": axon}


# ------------------------------------------------------------------------------------------------------------------
# Basket cell: published Bezaire PV basket cell
# ------------------------------------------------------------------------------------------------------------------
# Morphology of class_pvbasketcell.hoc: soma 20 x 10 um; two apical trees of five sections (4 -> 1 um) from soma(1),
# two basal trees of three sections (2 -> 1 um) from soma(0); section lengths are the distances between the
# template's 3-D points (z = 0). cm 1.4 uF/cm2, Ra 100 ohm cm. Same port as tinyHippo-bio-plasticity cells.py.
_TREES = {
    "apA": ((0, 20), (38.9418, 112.106), (77.8837, 204.212), (155.767, 388.424), (194.709, 480.531), (233.651, 572.637)),
    "apB": ((0, 20), (-38.9418, 112.106), (-77.8837, 204.212), (-155.767, 388.424), (-194.709, 480.531),
            (-233.651, 572.637)),
    "baA": ((0, 0), (-38.9418, -92.1061), (-77.8837, -184.212), (-116.826, -276.318)),
    "baB": ((0, 0), (38.9419, -92.1061), (77.8837, -184.212), (116.826, -276.318)),
}
_DIAMS = {"ap": (4.0, 3.0, 2.0, 1.5, 1.0), "ba": (2.0, 1.5, 1.0)}
_PV_CH = dict(Na=("ch_Navaxonp", 0.15), Kdr=("ch_Kdrfast", 0.013), KvA=("ch_KvA", 0.00015),
              CavN=("ch_CavN", 0.0008), CavL=("ch_CavL", 0.005), KvCaB=("ch_KvCaB", 2e-7), KCaS=("ch_KCaS", 2e-6))
_PV = dict(rm=5555.0, e_leak=-65.0, ena=55.0, ek=-90.0, cm=1.4, Ra=100.0, soma=(20.0, 10.0))
_CA_OUT, _CA_IN = 2.0, 5e-6  # mM; template's fixed eca (about 171 mV at 34 C)


def _tree(prefix, root_end):
    pts, diams = _TREES[prefix], _DIAMS[prefix[:2]]
    return tuple((f"{prefix}{i}", math.dist(pts[i], pts[i + 1]), d, None if i == 0 else f"{prefix}{i - 1}",
                  root_end if i == 0 else 1.0) for i, d in enumerate(diams))


_PV_BRANCHES = _tree("apA", 1.0) + _tree("apB", 1.0) + _tree("baA", 0.0) + _tree("baB", 0.0)


@dataclass(frozen=True)
class BCParams:
    """Scale factors on the published PV basket cell (1 = published)."""
    na_scale: float = 1.0
    k_scale: float = 1.0
    leak_scale: float = 1.0
    compact: bool = False   # True: one segment per section (network use, ~3x faster); False: Bezaire's d_lambda nseg


def _bez_nseg(L, diam):
    """Bezaire's set_nseg, evaluated before Ra and cm are set (NEURON defaults 35.4 ohm cm, 1 uF/cm2)."""
    lam = 1e5 * math.sqrt(diam / (4 * math.pi * 100.0 * 35.4 * 1.0))
    return int((L / (0.1 * lam) + 0.9) / 2) * 2 + 1


def build_bc(p: BCParams = BCParams(), name="bc"):
    """Bezaire 2016 PV basket cell. Returns {"soma", "dend"}."""
    h.celsius = CELSIUS
    soma = h.Section(name=f"{name}_soma")
    soma.L, soma.diam = _PV["soma"]
    secs, dend = {"soma": soma}, {}
    for bname, L, d, parent, end in _PV_BRANCHES:
        sec = h.Section(name=f"{name}_{bname}")
        sec.L, sec.diam = L, d
        sec.connect(soma(end) if parent is None else dend[parent](1.0))
        dend[bname] = secs[bname] = sec
    scale = dict(Na=p.na_scale, Kdr=p.k_scale)
    eca = 8.314 * (273.15 + CELSIUS) / (2 * 9.649e4) * math.log(_CA_OUT / _CA_IN) * 1000
    for sec in secs.values():
        sec.nseg = 1 if p.compact else _bez_nseg(sec.L, sec.diam)
        sec.Ra, sec.cm = _PV["Ra"], _PV["cm"]
        for key, (mech, g) in _PV_CH.items():
            sec.insert(mech)
            for seg in sec:
                getattr(seg, mech).gmax = g * scale.get(key, 1.0)
        sec.insert("ch_leak")
        for seg in sec:
            seg.ch_leak.gmax = p.leak_scale / _PV["rm"]
            seg.ch_leak.e = _PV["e_leak"]
        sec.ena, sec.ek, sec.eca = _PV["ena"], _PV["ek"], eca
    return {"soma": soma, "dend": dend}


# ------------------------------------------------------------------------------------------------------------------
# Regular-spiking unit for EC, DG and SWT
# ------------------------------------------------------------------------------------------------------------------
@dataclass(frozen=True)
class RSParams:
    na_scale: float = 1.0
    k_scale: float = 1.0
    ka_scale: float = 1.0
    leak_scale: float = 1.0
    e_pas: float = BEZ_PYR_SOMA["v_rest"]
    L: float = _G.rs_soma_L
    diam: float = _G.rs_soma_diam


def build_rs(p: RSParams = RSParams(), name="rs"):
    """One compartment with the pyramidal somatic channels (Navp, Kdrp, KvAproxp, leak). Returns {"soma"}."""
    b = BEZ_PYR_SOMA
    h.celsius = CELSIUS
    soma = h.Section(name=f"{name}_soma")
    soma.L, soma.diam = p.L, p.diam
    soma.Ra, soma.cm = b["Ra_soma"], b["Cm"]
    for mech in ("pas", "ch_Navp", "ch_Kdrp", "ch_KvAproxp"):
        soma.insert(mech)
    for seg in soma:
        seg.pas.g, seg.pas.e = p.leak_scale / b["Rm"], p.e_pas
        seg.ch_Navp.gmax = b["gNav"] * p.na_scale
        seg.ch_Kdrp.gmax = b["gKdr"] * p.k_scale
        seg.ch_KvAproxp.gmax = b["gKvAprox"] * p.ka_scale
    soma.ena, soma.ek = b["ena"], b["ek"]
    return {"soma": soma}


# ------------------------------------------------------------------------------------------------------------------
# Measurements (DC current only)
# ------------------------------------------------------------------------------------------------------------------
def _run_step(sec, amp, delay=100.0, dur=500.0, tail=100.0, v_init=None):
    stim = h.IClamp(sec(0.5))
    stim.delay, stim.dur, stim.amp = delay, dur, amp
    tv, vv = h.Vector(), h.Vector()
    tv.record(h._ref_t)
    vv.record(sec(0.5)._ref_v)
    h.cvode_active(0)
    h.dt = 0.025
    h.celsius = CELSIUS
    h.finitialize(v_init if v_init is not None else -66.0)
    h.continuerun(delay + dur + tail)
    t, v = np.array(tv), np.array(vv)
    del stim
    return t, v


def spike_times(t, v, thr=0.0):
    idx = np.where((v[1:] > thr) & (v[:-1] <= thr))[0]
    return t[idx + 1]


def resting_potential(sec, settle=1000.0):
    t, v = _run_step(sec, 0.0, delay=settle, dur=0.0, tail=0.0)
    return float(v[-1])


def input_resistance(sec, amp=-0.02, delay=1000.0, dur=500.0):
    """Steady-state Rin (MOhm) from a small hyperpolarising DC step, after a settling period."""
    t, v = _run_step(sec, amp, delay=delay, dur=dur, tail=0.0)
    v0 = v[(t > delay - 50) & (t < delay)].mean()
    v1 = v[(t > delay + dur - 50) & (t <= delay + dur)].mean()
    return float((v1 - v0) / amp)


def rheobase(sec, i_max=2.0, tol=1e-4, dur=500.0, settle=500.0):
    """Smallest 500-ms DC current (nA) that evokes at least one spike, by bisection."""
    def spikes(i):
        t, v = _run_step(sec, i, delay=settle, dur=dur)
        return len(spike_times(t, v)) > 0
    if spikes(0.0):
        return 0.0
    lo, hi = 0.0, i_max
    if not spikes(hi):
        return None
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        lo, hi = (lo, mid) if spikes(mid) else (mid, hi)
    return hi


def f_i_curve(sec, amps, dur=500.0, settle=500.0):
    """For each DC amplitude: spike count, mean rate (Hz), first and last ISI (ms) and adaptation index
    (last ISI / first ISI) over the 500-ms step."""
    rows = []
    for a in amps:
        t, v = _run_step(sec, a, delay=settle, dur=dur)
        st = spike_times(t, v)
        st = st[(st >= settle) & (st <= settle + dur)]
        isi = np.diff(st)
        rows.append(dict(amp_nA=float(a), n_spikes=int(len(st)), rate_hz=len(st) / (dur / 1000.0),
                         first_isi_ms=float(isi[0]) if len(isi) else None,
                         last_isi_ms=float(isi[-1]) if len(isi) else None,
                         adaptation=float(isi[-1] / isi[0]) if len(isi) > 1 else None,
                         peak_mV=float(v.max())))
    return rows
