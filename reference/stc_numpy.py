"""Numpy reference for the STC part of neuron/TMExp2SynNMDA_STC.mod (forward Euler, same dt).

c, Tt and z of one synapse, given presynaptic and postsynaptic spike times (ms). Presynaptic calcium arrives
tc_delay after the spike. Also `advance_slow`: the offline integration used between online blocks (no spikes, c = 0),
which the pilot uses for minutes to hours.
"""
import numpy as np

STC_DEFAULTS = dict(c_pre=0.6, c_post=0.1655, tc_delay=18.8, tau_c=48.8, theta_p=3.0, theta_d=1.2, gamma_p=1645.6,
                    gamma_d=313.1, tau_h=688400.0, k_rel=0.1, Tt_max=1.3805, theta_tag=0.2, tau_z=3600000.0)


def stc_single(pre_ms, post_ms, t_stop, P=0.0, dt=0.025, **kw):
    p = dict(STC_DEFAULTS, **kw)
    n = int(round(t_stop / dt))
    pre_idx = {int(round((t + p["tc_delay"]) / dt)) for t in pre_ms}
    post_idx = {int(round(t / dt)) for t in post_ms}
    c = Tt = z = 0.0
    out = np.empty((n, 3))
    for k in range(n):
        dc = -c / p["tau_c"]
        dT = (-p["k_rel"] * Tt + p["gamma_p"] * (p["Tt_max"] - Tt) * (c > p["theta_p"])
              - p["gamma_d"] * (1 + Tt) * (c > p["theta_d"])) / p["tau_h"]
        dz = P * ((1 - z) * (Tt > p["theta_tag"]) - (z + 0.5) * (Tt < -p["theta_tag"])) / p["tau_z"]
        c, Tt, z = c + dt * dc, Tt + dt * dT, z + dt * dz
        if k + 1 in pre_idx:
            c += p["c_pre"]
        if k + 1 in post_idx:
            c += p["c_post"]
        out[k] = c, Tt, z
    return out


def advance_slow(Tt, z, P, dt_ms, k_rel=0.1, tau_h=688400.0, theta_tag=0.2, tau_z=3600000.0):
    """One offline step (no calcium): Tt relaxes exactly; z by forward Euler with the current P (array or scalar)."""
    Tt_new = Tt * np.exp(-k_rel * dt_ms / tau_h)
    dz = P * ((1 - z) * (Tt > theta_tag) - (z + 0.5) * (Tt < -theta_tag)) / tau_z
    return Tt_new, z + dt_ms * dz
