COMMENT
Tsodyks-Markram AMPA + NMDA synapse (same release and conductance equations as TMExp2SynNMDACA1_SCExc.mod) with
synaptic tagging and capture on top (tinyCA3 implementation plan section 3.2-3.3; Luboeinski & Tetzlaff 2021
Eqs. 11-14, Lehr et al. 2022 Eq. 10, in units of the baseline weight w0 = NetCon weight[0]).

State per synapse:
  c   postsynaptic calcium (dimensionless), jumps c_pre after a presynaptic spike (delay tc_delay) and c_post at a
      postsynaptic spike
  Tt  early-phase weight change in units of w0 (the tag): tagged if |Tt| > theta_tag
  z   late-phase (captured) weight in units of w0, in [-0.5, 1]
  P   PRP amount of the postsynaptic pool, set from Python (RANGE), held constant within a block
Each released quantum is scaled by (1 + Tt + z).

Events (NET_RECEIVE(weight)):
  flag 0, weight >= 0 : presynaptic spike -> TM release, conductance increment, schedules flag 2 after tc_delay
  flag 0, weight < 0  : postsynaptic spike (NetCon from the postsynaptic soma with weight -1) -> c += c_post
  flag 2              : delayed presynaptic calcium -> c += c_pre
Integration: forward Euler (METHOD euler) so that reference/stc_numpy.py can reproduce it step by step.
No plasticity noise (L&T's xi term is omitted in this pilot).
ENDCOMMENT

NEURON {
    POINT_PROCESS TMExp2SynNMDA_STC
    RANGE e, U, tau_f, tau_d
    RANGE tau1_ampa, tau2_ampa, gmax_ampa, tau1_nmda, tau2_nmda, gmax_nmda
    RANGE mg, alpha, beta
    RANGE i, i_ampa, i_nmda, g, g_ampa, g_nmda, g_nmda_eff, mggate
    RANGE R, u, last_rel
    RANGE c_pre, c_post, tc_delay, tau_c, theta_p, theta_d, gamma_p, gamma_d, tau_h, k_rel, Tt_max
    RANGE theta_tag, tau_z, P, plastic, n_pre, n_post
    NONSPECIFIC_CURRENT i
}

UNITS {
    (nA) = (nanoamp)
    (uS) = (microsiemens)
    (mV) = (millivolt)
    (mM) = (milli/liter)
}

PARAMETER {
    e = 0 (mV)
    U     = 0.1528 <1e-9,1>
    tau_f = 138.4  (ms) <1e-9,1e9>
    tau_d = 101.0  (ms) <1e-9,1e9>
    tau1_ampa = 0.5 (ms) <1e-9,1e9>
    tau2_ampa = 3.0 (ms) <1e-9,1e9>
    tau1_nmda = 2.93  (ms) <1e-9,1e9>
    tau2_nmda = 148.5 (ms) <1e-9,1e9>
    gmax_ampa = 1 (uS) <0,1e9>
    gmax_nmda = 1.23 (uS) <0,1e9>
    mg    = 1.0  (mM) <0,10>
    alpha = 0.062 (/mV)
    beta  = 3.57 (1)

    : --- STC (Luboeinski & Tetzlaff 2021, Table 2, network values; divided by h0) ---
    c_pre    = 0.6
    c_post   = 0.1655
    tc_delay = 18.8   (ms)
    tau_c    = 48.8   (ms)
    theta_p  = 3
    theta_d  = 1.2
    gamma_p  = 1645.6
    gamma_d  = 313.1
    tau_h    = 688400 (ms)
    k_rel    = 0.1          : 0.1 = L&T (tag ~1.9 h); 0.265 = Frey & Morris 1998 half-life ~30 min
    Tt_max   = 1.3805       : 1 nC / h0 - 1
    theta_tag = 0.2
    tau_z    = 3600000 (ms) : 60 min
    P        = 0
    plastic  = 1            : 0 freezes Tt and z (STP-only behaviour)
}

ASSIGNED {
    v (mV)
    i (nA)
    i_ampa (nA)
    i_nmda (nA)
    g (uS)
    g_ampa (uS)
    g_nmda (uS)
    g_nmda_eff (uS)
    factor_ampa
    factor_nmda
    last (ms)
    deltat (ms)
    rel
    scale
    R
    u
    last_rel
    mggate
    n_pre
    n_post
}

STATE {
    A_ampa (uS)
    B_ampa (uS)
    A_nmda (uS)
    B_nmda (uS)
    c
    Tt
    z
}

INITIAL {
    LOCAL tp
    if (tau1_ampa/tau2_ampa > 0.9999) { tau1_ampa = 0.9999*tau2_ampa }
    if (tau1_nmda/tau2_nmda > 0.9999) { tau1_nmda = 0.9999*tau2_nmda }
    A_ampa = 0
    B_ampa = 0
    A_nmda = 0
    B_nmda = 0
    c = 0
    : Tt and z keep their values across finitialize only if set after it; Python restores them each block.
    Tt = 0
    z = 0
    R = 1
    u = 0
    last = -1e9
    last_rel = 0
    n_pre = 0
    n_post = 0
    tp = (tau1_ampa*tau2_ampa)/(tau2_ampa - tau1_ampa) * log(tau2_ampa/tau1_ampa)
    factor_ampa = 1 / (-exp(-tp/tau1_ampa) + exp(-tp/tau2_ampa))
    tp = (tau1_nmda*tau2_nmda)/(tau2_nmda - tau1_nmda) * log(tau2_nmda/tau1_nmda)
    factor_nmda = 1 / (-exp(-tp/tau1_nmda) + exp(-tp/tau2_nmda))
}

BREAKPOINT {
    SOLVE state METHOD euler
    g_ampa = B_ampa - A_ampa
    g_nmda = B_nmda - A_nmda
    mggate     = 1 / (1 + (mg/beta)*exp(-alpha*v))
    g_nmda_eff = g_nmda * mggate
    g = g_ampa + g_nmda_eff
    i_ampa = g_ampa     * (v - e)
    i_nmda = g_nmda_eff * (v - e)
    i      = i_ampa + i_nmda
}

FUNCTION heav(x) {
    if (x > 0) { heav = 1 } else { heav = 0 }
}

DERIVATIVE state {
    A_ampa' = -A_ampa/tau1_ampa
    B_ampa' = -B_ampa/tau2_ampa
    A_nmda' = -A_nmda/tau1_nmda
    B_nmda' = -B_nmda/tau2_nmda
    c' = -c/tau_c
    Tt' = plastic * (-k_rel*Tt + gamma_p*(Tt_max - Tt)*heav(c - theta_p) - gamma_d*(1 + Tt)*heav(c - theta_d)) / tau_h
    z' = plastic * P * ((1 - z)*heav(Tt - theta_tag) - (z + 0.5)*heav(-Tt - theta_tag)) / tau_z
}

NET_RECEIVE(weight) {
    LOCAL uu, RR, w
    if (flag == 2) {
        c = c + c_pre
    } else if (weight < 0) {
        c = c + c_post
        n_post = n_post + 1
    } else {
        deltat = t - last
        if (deltat < 0) { deltat = 0 }
        uu = u * exp(-deltat/tau_f)
        RR = 1 - (1 - R) * exp(-deltat/tau_d)
        uu = uu + U*(1 - uu)
        rel = uu * RR
        if (rel < 0) { rel = 0 }
        if (rel > 1) { rel = 1 }
        scale = rel / U
        RR = RR - rel
        u = uu
        R = RR
        last = t
        last_rel = rel
        n_pre = n_pre + 1
        w = weight * (1 + Tt + z)
        if (w < 0) { w = 0 }
        A_ampa = A_ampa + w * gmax_ampa * scale * factor_ampa
        B_ampa = B_ampa + w * gmax_ampa * scale * factor_ampa
        A_nmda = A_nmda + w * gmax_nmda * scale * factor_nmda
        B_nmda = B_nmda + w * gmax_nmda * scale * factor_nmda
        net_send(tc_delay, 2)
    }
}
