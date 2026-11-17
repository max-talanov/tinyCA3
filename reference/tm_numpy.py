"""Numpy reference for the Tsodyks-Markram release of neuron/TMExp2SynNMDACA1_SCExc.mod (NET_RECEIVE block).

Per spike n after an interval dt_n (first spike: u = 0, R = 1):
    u <- u exp(-dt/tau_f);  R <- 1 - (1 - R) exp(-dt/tau_d);  u <- u + U (1 - u)
    rel = u R;  R <- R - rel
The conductance increment is proportional to rel / U.
"""
import numpy as np


def tm_release(spike_times_ms, U, tau_f, tau_d):
    """Release fraction rel at each spike time (ms)."""
    u, R, last = 0.0, 1.0, None
    out = []
    for t in spike_times_ms:
        if last is not None:
            dt = t - last
            u = u * np.exp(-dt / tau_f)
            R = 1.0 - (1.0 - R) * np.exp(-dt / tau_d)
        u = u + U * (1.0 - u)
        rel = min(max(u * R, 0.0), 1.0)
        R -= rel
        last = t
        out.append(rel)
    return np.array(out)
