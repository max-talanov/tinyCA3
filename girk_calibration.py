"""Calibration of the 5-HT1A -> Kir3 (GIRK) current (neuron/ch_GIRK.mod) to Johnston, McBain & Fisahn 2014,
J Physiol 592:4187 (doi:10.1113/jphysiol.2014.279083).

(a) rectification: least-squares fit of vhalf and vslope to the normalized I-V of the 5-HT-induced current
    (their Fig. 5B, points read from the plot, E_K = -85 mV there);
(b) size: hyperpolarization of one model pyramidal cell held near -45 mV (subthreshold) against gmax, with the measured
    4.3 +/- 1.6 mV (their Fig. 4C, active CA3 pyramidal cells near -49 mV, n = 14).
Writes results/layout_B/girk_calibration.json and results/layout_B/fig_girk_calibration.png.
Usage: env -u PYTHONPATH .venv/bin/python girk_calibration.py
"""
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.optimize import least_squares  # noqa: E402

from cells import PyrParams, build_pyr, h  # noqa: E402
from circuit import ROOT  # noqa: E402

OUT = ROOT / "results" / "layout_B"
V_FIG5B = np.array([-110, -100, -90, -80, -70, -60, -50, -40, -30.0])
I_FIG5B = np.array([-1.0, -0.55, -0.15, 0.05, 0.25, 0.35, 0.5, 0.5, 0.3])   # normalized, read from Fig. 5B
EK_SLICE = -85.0
DV_MEAS, DV_SEM = 4.3, 1.6
G_USED = 6e-5
HOLD_NA = 0.08
plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False})


def iv(v, vh, k, ek=EK_SLICE):
    g = 1.0 / (1.0 + np.exp((v - ek - vh) / k))
    return g * (v - ek)


def fit_rectification():
    def res(p):
        return iv(V_FIG5B, *p) / abs(iv(-110.0, *p)) - I_FIG5B
    r = least_squares(res, [15.0, 10.0], bounds=([0.0, 2.0], [100.0, 60.0]))
    vh, k = r.x
    ri = iv(-60.0, vh, k) / iv(-110.0, vh, k)
    return float(vh), float(k), float(ri)


def hyperpolarization(gs):
    cell = build_pyr(PyrParams(k_scale=2.0), name="cal")
    ic = h.IClamp(cell["soma"](0.5))
    ic.delay, ic.dur, ic.amp = 0.0, 2500.0, HOLD_NA

    def vm(g):
        for key in ("soma", "dend"):
            for seg in cell[key]:
                seg.ch_GIRK.gmax = g
        v, t = h.Vector().record(cell["soma"](0.5)._ref_v), h.Vector().record(h._ref_t)
        h.dt = 0.025
        h.finitialize(-66.0)
        h.continuerun(2500.0)
        v, t = np.array(v), np.array(t)
        return float(v[t > 2000].mean())
    v0 = vm(0.0)
    return v0, [v0 - vm(g) for g in gs]


def main():
    vh, k, ri = fit_rectification()
    gs = np.array([0, 1e-5, 2e-5, 4e-5, 6e-5, 8e-5, 1e-4])
    v0, dv = hyperpolarization(gs)
    dv_used = hyperpolarization([G_USED])[1][0]
    out = dict(vhalf_mV=round(vh, 2), vslope_mV=round(k, 2), rectification_index_fit=round(ri, 2),
               rectification_index_measured_median=-0.62, held_Vm_mV=round(v0, 1), gmax_used=G_USED,
               dV_used_mV=round(dv_used, 2), dV_measured_mV=[DV_MEAS, DV_SEM],
               dV_by_gmax={f"{g:g}": round(d, 2) for g, d in zip(gs, dv)})
    (OUT / "girk_calibration.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))

    fig, axes = plt.subplots(1, 2, figsize=(7.6, 2.9))
    ax = axes[0]
    vv = np.linspace(-115, -25, 200)
    ax.plot(V_FIG5B, I_FIG5B, "ko", ms=4, label="Johnston et al. 2014, Fig. 5B")
    ax.plot(vv, iv(vv, vh, k) / abs(iv(-110.0, vh, k)), color="#8e44ad", lw=1.6,
            label=f"model: vhalf {vh:.0f} mV, vslope {k:.1f} mV")
    ax.axhline(0, color="0.6", lw=0.6)
    ax.set_xlabel("membrane potential (mV)")
    ax.set_ylabel("5-HT current (normalized to -110 mV)")
    ax.set_title("(a) inward rectification", fontsize=8)
    ax.legend(frameon=False, fontsize=6.5)
    ax = axes[1]
    ax.axhspan(DV_MEAS - DV_SEM, DV_MEAS + DV_SEM, color="0.85", lw=0)
    ax.axhline(DV_MEAS, color="k", lw=0.8, ls="--", label="measured 4.3 ± 1.6 mV (Fig. 4C)")
    ax.plot(gs * 1e5, dv, "o-", color="#8e44ad", ms=3.5, label=f"model cell held at {v0:.0f} mV")
    ax.plot([G_USED * 1e5], [dv_used], "s", color="#8e44ad", ms=7, mfc="white", label=f"used: {dv_used:.1f} mV")
    ax.set_xlabel("gmax (10⁻⁵ S/cm²)")
    ax.set_ylabel("hyperpolarization (mV)")
    ax.set_title("(b) size of the effect", fontsize=8)
    ax.legend(frameon=False, fontsize=6.5)
    fig.tight_layout()
    fig.savefig(OUT / "fig_girk_calibration.png", dpi=160)


if __name__ == "__main__":
    main()
