"""T5 (pilot form): PRP pools, neuromodulator schedules and the offline (slow) integration.

Pool granularity "cell": one PRP pool P_j per postsynaptic PYR (implementation plan Table 2 of the algorithm doc).
    tau_P dP_j/dt = -P_j + alpha_P * Theta[ sum_{s -> j} |Tt_s| - theta_pro(M_j) ],  theta_pro(M) = theta0 / (M + 0.001)
(Luboeinski & Tetzlaff 2021 Eq. 15; Lehr et al. 2022 Eq. 10). M_j is the modulator reaching cell j: DA reaches the cells
of X (v1: Λ), 5-HT the cells of + (v1: V) (architecture doc, "Mapping of the modulators"), each as rectangular pulses.
Synthesis trigger (tau_G_min > 0, report of 2026-11-24): the modulator does not set theta_pro directly. It switches on
a per-cell synthesis state G_j (the D1/D5 / 5-HT receptor -> cAMP/PKA -> translation step), which then outlasts it:
    dG_j/dt = M_j (1 - G_j) / tau_on - G_j / tau_G,   theta_pro = theta0 / (G_j + 0.001)
so a modulator of minutes opens a synthesis window of about tau_G. tau_G_min = 0 gives the L&T form (G = M).
Between online blocks there is no calcium, so Tt relaxes exponentially and z follows the capture equation with the
current P (reference/stc_numpy.advance_slow). No time compression: every step is real time.
"""
from dataclasses import dataclass, field

import numpy as np

from circuit import CELL_IDS, IDX, MOD_PATTERN, PATTERNS
from reference.stc_numpy import advance_slow


@dataclass
class STCParams:
    tau_P_min: float = 60.0        # L&T 2021 (Frey & Morris 1998 give a PRP half-life of 1-2 h)
    alpha_P: float = 1.0           # L&T 2021
    theta0: float = 0.5 * 1.001    # theta_pro(M = 1) = 0.5 (L&T 2021 theta_pro = 0.5 h0)
    k_rel: float = 0.1             # L&T 2021; 0.265 = tag half-life ~30 min (Frey & Morris 1998, S&A p. 380)
    tau_h_ms: float = 688400.0
    tau_z_min: float = 60.0
    dt_slow_s: float = 10.0
    tau_G_min: float = 0.0         # synthesis-trigger decay; 0 = off (theta_pro follows M directly, L&T 2021)
    tau_on_min: float = 2.0        # synthesis-trigger rise while the modulator is on (free)


@dataclass
class ModPulse:
    source: str      # "DA" or "5-HT"
    t_on_s: float
    dur_s: float
    amp: float = 1.0


MOD_CELLS = {m: set(PATTERNS[k]) for m, k in MOD_PATTERN.items()}   # v1: Λ and V; layout B: X and +


def modulator_at(pulses, t_s):
    """Per-cell modulator level M_j(t) and per-source levels."""
    M = np.zeros(len(CELL_IDS))
    src = {"DA": 0.0, "5-HT": 0.0}
    for p in pulses:
        if p.t_on_s <= t_s < p.t_on_s + p.dur_s:
            src[p.source] += p.amp
            for cid in MOD_CELLS[p.source]:
                M[IDX[cid]] += p.amp
    return M, src


def tag_sum_per_cell(Tt, post_idx):
    return np.bincount(post_idx, weights=np.abs(Tt), minlength=len(CELL_IDS))


def advance(Tt, z, P, post_idx, pulses, t0_s, t1_s, sp: STCParams, on_step=None, G=None):
    """Integrate Tt, z, P (and the synthesis trigger G) from t0 to t1 (s) offline, in steps of dt_slow_s.
    Returns the new state Tt, z, P, G (G is None when the trigger is off)."""
    t = t0_s
    if sp.tau_G_min > 0 and G is None:
        G = np.zeros_like(P)
    while t < t1_s - 1e-9:
        dt = min(sp.dt_slow_s, t1_s - t)
        M, _ = modulator_at(pulses, t)
        if sp.tau_G_min > 0:
            a = M / (60.0 * sp.tau_on_min)                 # exact step of dG/dt = a - b G with M fixed over dt
            b = a + 1.0 / (60.0 * sp.tau_G_min)
            G = a / b + (G - a / b) * np.exp(-b * dt)
            gate = G
        else:
            gate = M
        theta = sp.theta0 / (gate + 0.001)
        drive = sp.alpha_P * (tag_sum_per_cell(Tt, post_idx) > theta)
        P = drive + (P - drive) * np.exp(-dt / (60.0 * sp.tau_P_min))
        Tt, z = advance_slow(Tt, z, P[post_idx], dt * 1000.0, k_rel=sp.k_rel, tau_h=sp.tau_h_ms,
                             tau_z=sp.tau_z_min * 60000.0)
        t += dt
        if on_step:
            on_step(t, Tt, z, P, G)
    return Tt, z, P, G
