"""T3 (pilot form): the tinyCA3 circuit of docs/ca3_basic_circuit_architecture.md, built with Hodgkin-Huxley cells.

Projections (implementation plan Table 3):
  P1 PYR -> PYR    all-to-all, no autapses, TMExp2SynNMDA_STC on the postsynaptic dendrite (STP + STC)
  P2 PYR -> BC     1:1, AMPA only
  P3 BC  -> PYR    1:all except self (self_inhibition off), GABA_A (TM synapse, gmax_nmda = 0, e = E_GABA)
  P4 EC_k -> DG_k  1:1, AMPA only
  P5 DG_k -> PYR   mossy fibre to every cell of pattern k, AMPA + NMDA, proximal dendrite
  P6 PYR -> EC_k   cells of pattern k to EC_k, AMPA only, strongly depressing (limits the loop to one reinstatement)
  P8 SWT -> PYR    every trigger cell to every PYR, AMPA only, jittered delays
Every weight is a pilot value (status free) chosen by calibrate_pilot.py; none is a CA3 measurement.
"""
from dataclasses import dataclass, field

import numpy as np

from cells import BCParams, PyrParams, build_bc, build_pyr, build_rs, h

# Patterns read from the architecture doc Table 3 (cell = 10*row + column).
PATTERNS = {
    "Λ": (13, 22, 24, 32, 34, 41, 45, 51, 55),
    "V": (11, 15, 21, 25, 32, 34, 42, 44, 53),
    "I": (13, 23, 33, 43, 53),
    "O": tuple(sorted({10 * r + c for r in range(1, 6) for c in range(1, 6) if r in (1, 5) or c in (1, 5)})),
}
PATTERNS["Λ∪V"] = tuple(sorted(set(PATTERNS["Λ"]) | set(PATTERNS["V"])))
PATTERNS["Λ∪V∪I"] = tuple(sorted(set(PATTERNS["Λ∪V"]) | set(PATTERNS["I"])))
PATTERN_ORDER = ("Λ", "V", "I", "O", "Λ∪V", "Λ∪V∪I")
MODULATOR_OF = {"Λ": "DA", "V": "5-HT", "Λ∪V": "DA+5-HT", "Λ∪V∪I": "DA+5-HT"}
CELL_IDS = tuple(10 * r + c for r in range(1, 6) for c in range(1, 6))
IDX = {cid: i for i, cid in enumerate(CELL_IDS)}


@dataclass
class CircuitParams:
    w0: float = 0.0006            # P1 baseline NetCon weight (uS-equivalent), calibrated: rule of plan section 2.4
    w_pb: float = 0.03            # P2
    w_bp: float = 0.01            # P3, per synapse
    e_gaba: float = -75.0         # mV, placeholder
    gaba_tau1: float = 0.3        # ms, placeholder
    gaba_tau2: float = 8.0        # ms, placeholder
    w_ed: float = 0.006           # P4
    w_mf: float = 0.35            # P5, AMPA only: keeps members near 90 Hz in the large patterns despite inhibition
    mf_nmda: bool = False         # P5 NMDA component (True gives seconds of reverberation after encoding)
    p1_gmax_nmda: float = 0.3     # P1 NMDA/AMPA conductance ratio, a quarter of the CA1 value 1.23 (limits reverberation)
    w_ce_total: float = 0.00268   # P6, divided over the cells of the pattern: EC fires at ~90% (threshold 0.00236)
    w_swt: float = 0.0016         # P8, per synapse; 10 jittered inputs fire a quiet PYR at 0.00175
    n_swt: int = 10
    swt_jitter: float = 15.0      # ms, delays uniform in [1, 1 + jitter]
    d_rc: float = 1.5             # ms, all delays are placeholders
    d_pb: float = 1.0
    d_bp: float = 1.0
    d_ed: float = 3.0
    d_mf: float = 3.0
    d_ce: float = 3.0
    w_bg: float = 0.004           # background Poisson input onto PYR dendrites
    bg_rate: float = 4.0          # Hz per PYR (gives ~0.5-1 Hz firing)
    w_drive: float = 0.008        # encoding drive NetStim -> EC_k
    k_rel: float = 0.1            # tag decay of L&T 2021 (~1.9 h); 0.265 = Frey & Morris 1998 half-life ~30 min
    c_pre: float = 1.0            # L&T 2021 Table 2 slice values (Graupner & Brunel 2012 fit); the 0.6x network
    c_post: float = 0.2758        # values never reach theta_p with these HH cells (pilot calibration, RESULTS 3)
    pyr_k_scale: float = 2.0      # Kdr x2: lets PYR follow 100 Hz mossy-fibre input (x1 follows 2 of 10 spikes)
    seed: int = 1


class Circuit:
    def __init__(self, p: CircuitParams = CircuitParams()):
        self.p = p
        self.rng = np.random.default_rng(p.seed)
        self.pyr = {cid: build_pyr(PyrParams(k_scale=p.pyr_k_scale), name=f"ex_{cid}") for cid in CELL_IDS}
        self.bc = {cid: build_bc(BCParams(compact=True), name=f"in_{cid}") for cid in CELL_IDS}
        self.ec = {k: build_rs(name=f"EC_{i + 1}") for i, k in enumerate(PATTERN_ORDER)}
        self.dg = {k: build_rs(name=f"DG_{i + 1}") for i, k in enumerate(PATTERN_ORDER)}
        self.swt = [build_rs(name=f"SWT_{i}") for i in range(p.n_swt)]
        self.keep = []
        self.enc_stims = []
        self._wire()
        self._record()

    # -------------------------------------------------------------------------------------------------------------
    def _syn(self, sec, x, nmda=True, stp=None, gaba=False):
        s = h.TMExp2SynNMDACA1_SCExc(sec(x))
        if not nmda or gaba:
            s.gmax_nmda = 0
        if gaba:
            s.e, s.tau1_ampa, s.tau2_ampa = self.p.e_gaba, self.p.gaba_tau1, self.p.gaba_tau2
        if stp:
            s.U, s.tau_f, s.tau_d = stp
        self.keep.append(s)
        return s

    def _nc(self, src_sec, target, weight, delay):
        nc = h.NetCon(src_sec(0.5)._ref_v, target, sec=src_sec)
        nc.threshold, nc.delay, nc.weight[0] = -10.0, delay, weight
        self.keep.append(nc)
        return nc

    def _wire(self):
        p = self.p
        # P1 recurrent collaterals with STC, plus the postsynaptic-spike NetCons that feed calcium
        self.p1 = []            # list of (pre_id, post_id, synapse)
        for post in CELL_IDS:
            dend = self.pyr[post]["dend"]
            for pre in CELL_IDS:
                if pre == post:
                    continue
                s = h.TMExp2SynNMDA_STC(dend(0.5))
                s.k_rel, s.c_pre, s.c_post, s.gmax_nmda = p.k_rel, p.c_pre, p.c_post, p.p1_gmax_nmda
                self._nc(self.pyr[pre]["soma"], s, p.w0, p.d_rc)
                self._nc(self.pyr[post]["soma"], s, -1.0, 0.0)      # postsynaptic spike -> calcium
                self.p1.append((pre, post, s))
        self.p1_pre = np.array([IDX[a] for a, _, _ in self.p1])
        self.p1_post = np.array([IDX[b] for _, b, _ in self.p1])
        # P2, P3
        for cid in CELL_IDS:
            self._nc(self.pyr[cid]["soma"], self._syn(self.bc[cid]["soma"], 0.5, nmda=False), p.w_pb, p.d_pb)
            for tgt in CELL_IDS:
                if tgt != cid:
                    self._nc(self.bc[cid]["soma"], self._syn(self.pyr[tgt]["soma"], 0.5, gaba=True), p.w_bp, p.d_bp)
        # P4, P5, P6
        for k, cells in PATTERNS.items():
            self._nc(self.ec[k]["soma"], self._syn(self.dg[k]["soma"], 0.5, nmda=False), p.w_ed, p.d_ed)
            ec_syn = self._syn(self.ec[k]["soma"], 0.5, nmda=False, stp=(0.9, 1.0, 1500.0))
            for cid in cells:
                self._nc(self.dg[k]["soma"], self._syn(self.pyr[cid]["dend"], 0.05, nmda=p.mf_nmda), p.w_mf, p.d_mf)
                self._nc(self.pyr[cid]["soma"], ec_syn, p.w_ce_total / len(cells), p.d_ce)
        # P8 trigger: one synapse per PYR, each trigger cell contacts it with its own jittered delay
        for cid in CELL_IDS:
            s = self._syn(self.pyr[cid]["dend"], 0.5, nmda=False)
            for sw in self.swt:
                self._nc(sw["soma"], s, p.w_swt, 1.0 + self.rng.uniform(0, p.swt_jitter))
        # drives (NetStims are inputs only): EC encoding drive, SWT burst, PYR background
        self.ec_drive = {}
        for k in PATTERN_ORDER:
            ns = h.NetStim()
            ns.number, ns.noise = 0, 0
            nc = h.NetCon(ns, self._syn(self.ec[k]["soma"], 0.5, nmda=False, stp=(1.0, 1.0, 1.0)))
            nc.weight[0] = p.w_drive
            self.keep += [ns, nc]
            self.ec_drive[k] = ns
        self.swt_drive = []
        for sw in self.swt:
            ns = h.NetStim()
            ns.number, ns.noise = 0, 0
            nc = h.NetCon(ns, self._syn(sw["soma"], 0.5, nmda=False, stp=(1.0, 1.0, 1.0)))
            nc.weight[0] = 0.008
            self.keep += [ns, nc]
            self.swt_drive.append(ns)
        self.bg = []
        for i, cid in enumerate(CELL_IDS):
            ns = h.NetStim()
            ns.interval, ns.noise, ns.number, ns.start = 1000.0 / p.bg_rate, 1.0, 1e9, 0.0
            ns.noiseFromRandom123(i, 0, p.seed)
            nc = h.NetCon(ns, self._syn(self.pyr[cid]["dend"], 0.5, nmda=False, stp=(1.0, 1.0, 1.0)))
            nc.weight[0] = p.w_bg
            self.keep += [ns, nc]
            self.bg.append(ns)

    def _record(self, dt_rec=0.5):
        self.spk = {}
        groups = [("PYR", cid, self.pyr[cid]) for cid in CELL_IDS] + [("BC", cid, self.bc[cid]) for cid in CELL_IDS]
        groups += [("EC", k, self.ec[k]) for k in PATTERN_ORDER] + [("DG", k, self.dg[k]) for k in PATTERN_ORDER]
        groups += [("SWT", i, s) for i, s in enumerate(self.swt)]
        self.vrec = {}
        for typ, key, cell in groups:
            v = h.Vector()
            nc = h.NetCon(cell["soma"](0.5)._ref_v, None, sec=cell["soma"])
            nc.threshold = -10.0
            nc.record(v)
            self.keep.append(nc)
            self.spk[(typ, key)] = v
            vv = h.Vector()
            vv.record(cell["soma"](0.5)._ref_v, dt_rec)
            self.vrec[(typ, key)] = vv
        self.trec = h.Vector()
        self.trec.record(h._ref_t, dt_rec)

    # -------------------------------------------------------------------------------------------------------------
    def get_stc(self):
        return (np.array([s.Tt for _, _, s in self.p1]), np.array([s.z for _, _, s in self.p1]))

    def set_stc(self, Tt, z, P_cell):
        for i, (_, _, s) in enumerate(self.p1):
            s.Tt, s.z, s.P = Tt[i], z[i], P_cell[self.p1_post[i]]

    def set_tt(self, value_by_syn):
        for i, (_, _, s) in enumerate(self.p1):
            s.Tt = value_by_syn[i]

    def schedule_encoding(self, k, t_on, n_pulses=3, pulse=100.0, gap=400.0, rate=100.0):
        """Drive EC_k with n_pulses bursts of `pulse` ms at `rate` Hz, `gap` ms apart (L&T 2021 learning protocol).
        One NetStim plays one burst, so extra bursts use extra NetStims."""
        stims = []
        for i in range(n_pulses):
            ns = h.NetStim()
            ns.start, ns.interval, ns.number, ns.noise = t_on + i * (pulse + gap), 1000.0 / rate, pulse * rate / 1000.0, 0
            nc = h.NetCon(ns, self._syn(self.ec[k]["soma"], 0.5, nmda=False, stp=(1.0, 1.0, 1.0)))
            nc.weight[0] = self.p.w_drive
            stims += [ns, nc]
            self.enc_stims.append(ns)
        self.keep += stims

    def schedule_swr(self, t_on, n_spikes=1):
        for ns in self.swt_drive:
            ns.start, ns.number, ns.interval, ns.noise = t_on, n_spikes, 5.0, 0

    def clear_schedules(self):
        for ns in self.swt_drive + self.enc_stims:
            ns.number = 0

    def run_block(self, t_stop, Tt, z, P_cell, v_init=-66.0, block=0):
        """finitialize, restore the slow state (finitialize resets it), run, and return spikes, voltages, Tt and z.
        `block` selects an independent background-noise stream (Random123 id2), so successive blocks differ."""
        for i, ns in enumerate(self.bg):
            ns.noiseFromRandom123(i, block, self.p.seed)
        h.cvode_active(0)
        h.dt = 0.025
        for v in self.spk.values():
            v.resize(0)
        h.finitialize(v_init)
        self.set_stc(Tt, z, P_cell)
        h.fcurrent()
        h.continuerun(t_stop)
        spikes = {key: np.array(v) for key, v in self.spk.items()}
        volts = {key: np.array(v) for key, v in self.vrec.items()}
        Tt2, z2 = self.get_stc()
        return dict(t=np.array(self.trec), spikes=spikes, v=volts, Tt=Tt2, z=z2)


def pattern_synapse_mask(circ: Circuit, k):
    cells = set(PATTERNS[k])
    return np.array([(a in cells) and (b in cells) for a, b, _ in circ.p1])
