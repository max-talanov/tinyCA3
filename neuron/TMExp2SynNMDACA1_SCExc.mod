COMMENT
Tsodyks-Markram synaptic model with AMPA + NMDA components.
Calibrated for Schaffer Collateral -> CA1 pyramidal cell (SC->Exc).

Parameters from hippocampushub.eu (SC->Exc, CA1):

  Short Term Plasticity - optimized fit to Wierenga & Wadman 2003 Fig. 4A+4B:
    U     = 0.1528  (hippocampushub: 0.14 +/- 0.08)
    tau_f = 138.4 ms (hippocampushub: 129 +/- 68)
    tau_d = 101.0 ms (hippocampushub: 186 +/- 71)
    Optimization: differential_evolution on plateau levels (Fig.4A)
    and paired-pulse ratios (Fig.4B), weighted by SEM
    Cost reduction: 49.9 -> 5.8 (~9x improvement)
    Species: Wistar rat, 14-28 d, CA1

  NMDA kinetics (Andrasfalvy & Magee 2001, J Neurosci):
    tau1_nmda = 2.93  ms  (rise)
    tau2_nmda = 148.5 ms  (decay)
    Species: Sprague Dawley rat, 6-12 w, CA1, n=52

  NMDA/AMPA conductance ratio = 1.23 (Le Roux et al. 2013, at +40 mV, Mg-free)
    => gmax_nmda = 1.23 * gmax_ampa
    Species: Mus musculus, 17-23 d, CA1, n=12

  AMPA kinetics (standard CA1 Schaffer ~30 degC):
    tau1_ampa = 0.5 ms, tau2_ampa = 3.0 ms

  Mg block: Jahr & Stevens (1990)
    mg=1.0 mM, alpha=0.062 /mV, beta=3.57

Filosofia gmax (identica al MOD originale Tsodyks-Markram):
  gmax_ampa = 1 uS  (placeholder; weight=1 nella NetCon)
  gmax_nmda = 1.23 uS  (= NMDA/AMPA ratio * gmax_ampa)
  La scala assoluta e' interamente controllata dal weight della NetCon.

  Per ottenere il primo picco EPSC desiderato (es. I_target pA a -70 mV):
    g_target  = I_target [pA] / 70000 [pA/uS]   (driving force 70 mV)
    weight    = g_target [uS] / gmax_ampa [uS]
              = g_target [uS]     (perche' gmax_ampa=1)

  Esempi:
    I_target = 135 pA  (Table 1 Wierenga 2003)
      => g = 135/70000 = 0.00193 uS  => weight = 0.00193
    I_target = 200 pA
      => weight = 200/70000 = 0.00286

Note: i parametri TM (U, tau_d, tau_f) descrivono la dinamica AMPA
(registrazioni a -70 mV, NMDA bloccata da Mg). La stessa dinamica
viene applicata anche alla componente NMDA (stesso terminale presinaptico).
A -70 mV la componente NMDA e' ~90% bloccata dal Mg; diventa rilevante
quando il neurone si depolarizza durante l'attivita'.
ENDCOMMENT

NEURON {
    POINT_PROCESS TMExp2SynNMDACA1_SCExc
    RANGE e, U, tau_f, tau_d
    RANGE tau1_ampa, tau2_ampa, gmax_ampa
    RANGE tau1_nmda, tau2_nmda, gmax_nmda
    RANGE mg, alpha, beta
    RANGE i, i_ampa, i_nmda, g, g_ampa, g_nmda, g_nmda_eff, mggate
    RANGE R, u, last_rel
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

    : --- Tsodyks-Markram STP (Wierenga & Wadman 2003, SC->Exc CA1) ---
    U     = 0.1528 <1e-9,1>           : optimized fit Fig.4A+4B (hippocampushub: 0.14)
    tau_f = 138.4  (ms) <1e-9,1e9>    : optimized fit Fig.4A+4B (hippocampushub: 129ms)
    tau_d = 101.0  (ms) <1e-9,1e9>    : optimized fit Fig.4A+4B (hippocampushub: 186ms)

    : --- AMPA kinetics (~30 degC, CA1 Schaffer) ---
    tau1_ampa = 0.5 (ms) <1e-9,1e9>
    tau2_ampa = 3.0 (ms) <1e-9,1e9>

    : --- NMDA kinetics (Andrasfalvy & Magee 2001) ---
    tau1_nmda = 2.93  (ms) <1e-9,1e9>
    tau2_nmda = 148.5 (ms) <1e-9,1e9>

    : --- Conductances: placeholder, scala con weight NetCon ---
    : gmax_nmda = 1.23 * gmax_ampa  (Le Roux et al. 2013, NMDA/AMPA ratio)
    : weight NetCon = I_target[pA] / 70000  per Vh=-70mV, e=0mV
    gmax_ampa = 1 (uS) <0,1e9>
    gmax_nmda = 1.23 (uS) <0,1e9>

    : --- Mg2+ block (Jahr & Stevens 1990) ---
    mg    = 1.0  (mM) <0,10>
    alpha = 0.062 (/mV)
    beta  = 3.57 (1)
}

ASSIGNED {
    v (mV)
    i      (nA)
    i_ampa (nA)
    i_nmda (nA)
    g          (uS)
    g_ampa     (uS)
    g_nmda     (uS)
    g_nmda_eff (uS)
    factor_ampa
    factor_nmda
    last   (ms)
    deltat (ms)
    rel
    scale
    R
    u
    last_rel
    mggate
}

STATE {
    A_ampa (uS)
    B_ampa (uS)
    A_nmda (uS)
    B_nmda (uS)
}

INITIAL {
    LOCAL tp

    if (tau1_ampa/tau2_ampa > 0.9999) { tau1_ampa = 0.9999*tau2_ampa }
    if (tau1_nmda/tau2_nmda > 0.9999) { tau1_nmda = 0.9999*tau2_nmda }

    A_ampa = 0
    B_ampa = 0
    A_nmda = 0
    B_nmda = 0

    R        = 1
    u        = 0
    last     = -1e9
    last_rel = 0

    tp = (tau1_ampa*tau2_ampa)/(tau2_ampa - tau1_ampa) * log(tau2_ampa/tau1_ampa)
    factor_ampa = 1 / (-exp(-tp/tau1_ampa) + exp(-tp/tau2_ampa))

    tp = (tau1_nmda*tau2_nmda)/(tau2_nmda - tau1_nmda) * log(tau2_nmda/tau1_nmda)
    factor_nmda = 1 / (-exp(-tp/tau1_nmda) + exp(-tp/tau2_nmda))
}

BREAKPOINT {
    SOLVE state METHOD cnexp

    g_ampa = B_ampa - A_ampa
    g_nmda = B_nmda - A_nmda

    mggate     = 1 / (1 + (mg/beta)*exp(-alpha*v))
    g_nmda_eff = g_nmda * mggate

    g = g_ampa + g_nmda_eff

    i_ampa = g_ampa     * (v - e)
    i_nmda = g_nmda_eff * (v - e)
    i      = i_ampa + i_nmda
}

DERIVATIVE state {
    A_ampa' = -A_ampa/tau1_ampa
    B_ampa' = -B_ampa/tau2_ampa
    A_nmda' = -A_nmda/tau1_nmda
    B_nmda' = -B_nmda/tau2_nmda
}

NET_RECEIVE(weight (1)) {
    LOCAL uu, RR

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

    u        = uu
    R        = RR
    last     = t
    last_rel = rel

    A_ampa = A_ampa + weight * gmax_ampa * scale * factor_ampa
    B_ampa = B_ampa + weight * gmax_ampa * scale * factor_ampa

    A_nmda = A_nmda + weight * gmax_nmda * scale * factor_nmda
    B_nmda = B_nmda + weight * gmax_nmda * scale * factor_nmda
}
