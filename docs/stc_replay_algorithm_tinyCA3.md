# Synaptic tagging and capture with replay in tinyCA3: the algorithm

Nov 12, 2026

*Copied from [tinyHippo-bio-plasticity](https://github.com/max-talanov/tinyHippo-bio-plasticity/blob/main/docs/stc_replay_algorithm_tinyCA3.md). References to CLAUDE.md, `params.py`, the S1 rule and "the plasticity model" refer to that repository.*

## The main idea

This document turns the tag/capture literature
([synaptic_tagging_literature_review.md](synaptic_tagging_literature_review.md)) into one concrete
algorithm for the tinyCA3 circuit
([ca3_basic_circuit_architecture.md](ca3_basic_circuit_architecture.md),
"the CA3 doc"):

- 25 pyramidal cells (PYR) with all-to-all recurrent collaterals (RC, projection P1, 600 synapses);
- 25 basket cells giving global feedback inhibition;
- EC and DG units that impose six patterns;
- dopamine (DA) on the Λ pattern and serotonin (5-HT) on the V pattern;
- a CA3 → EC loop (P6) that replays stored patterns.

STC and replay act on each other in two directions, and the algorithm makes both explicit:

1. **STC decides what is replayed.** Replay is pattern completion through the P1 weights. A pattern
   whose tags were captured keeps its raised weights and is still replayed hours later. A pattern
   whose tags decayed stops being replayed once its early weights fall back (~1–2 h).
2. **Replay decides what is captured.** Each replay event co-activates the cells of a pattern.
   That re-sets (refreshes) its tags, can deliver a modulator, and can supply PRPs. A refreshed tag
   can capture PRPs that arrive later than the original tag would have survived. Replay of one
   pattern also lowers the chance that another pattern is captured: through a limited PRP pool, and
   through SWR-triggered depression of what is not replayed.

The synapse model is the published recurrent-network STC model of Luboeinski & Tetzlaff (2021) and
Lehr et al. (2022). Of the models reviewed in Sajikumar & Abel (2024, Ch. 23), only these have
complete equations and parameters for a recurrent spiking network. It is written in the notation
of the CA3 doc (tag T, PRP pool P, stored weight w_cap), and switches select the alternatives the
literature leaves open.

Nothing below is built. Every value is marked with its source, or as a free parameter.

## 1. State variables

**Table 1.** State of the circuit. *s* indexes a P1 synapse i → j, *k* a PRP pool, *m* a
neuromodulator.

| Symbol | Where | Meaning | Source of the form |
|---|---|---|---|
| u_s, x_s | Every synapse | Tsodyks–Markram facilitation and resources | Plasticity model Eq. 1 |
| c_s | P1 | Postsynaptic calcium at synapse s | L&T 2021 Eq. 12 |
| T_s | P1 | Early-phase weight change, signed: T > 0 is E-LTP, T < 0 is E-LTD. It is also the tag | L&T 2021 Eq. 11 with T = h − h0 |
| z_s | P1 | Late-phase (captured) weight, z ∈ [−0.5, 1] | L&T 2021 Eqs. 13–14 |
| P_k | One per pool | PRP amount | L&T 2021 Eq. 15 |
| M_m(t) | Per source, m ∈ {DA, 5-HT} | Neuromodulator signal | Lehr 2022 Eq. 10; CA3 doc Eq. 2 |

The weight that a spike delivers through P1 is CA3 doc Eq. 1, with α = 1 and w_cap = w0·z:

```
Δg_n = [ w0 + T_s(t_n) + w0·z_s(t_n) ] · A_n                                  (1)
```

Equivalence: L&T 2021 write w = h + h0·z. With h = w0 + T and h0 = w0 this is the same weight.

## 2. Synapse equations

### 2.1 Induction: calcium sets the sign (fast, during spiking)

```
dc_s/dt = −c_s/τ_c + c_pre Σ δ(t − t_i − t_c,delay) + c_post Σ δ(t − t_j)       (2)

τ_h dT_s/dt = −k_rel·T_s
              + γ_p (T_max − T_s) Θ[c_s − θ_p]
              − γ_d (w0 + T_s)    Θ[c_s − θ_d]
              + ξ_s(t)                                                         (3)
```

- **Origin of Eq. 3:** L&T 2021 Eq. 11, rewritten with T = h − h0 and T_max = 1 nC − h0. ξ is
  their calcium-gated noise, σ_pl = 0.290436 nC s^½.
- **Calcium decides the sign:**
  - c ≥ θ_p: potentiation.
  - θ_d ≤ c < θ_p: depression.
  - c < θ_d: T only relaxes.

  This is the two-threshold form proposed as Eq. 8 in the plasticity model.
- **Option:** the plasticity model's segment signal g_Ca(V_seg) (its Eq. 2) can stand in for
  Eq. 2 here (switch `induction = "calcium" | "vseg"`). Calcium is the default, because its
  parameters are fitted to slice data (Graupner & Brunel 2012) and published for a recurrent
  network.
- **k_rel sets the tag lifetime.** It is the one place where the literature disagrees, so it is a
  parameter (§5):
  - k_rel = 0.1 is L&T: T relaxes with τ_h/0.1 = 1.9 h.
  - k_rel = 0.265 gives τ = 43 min, the ~30 min tag half-life of Frey & Morris 1998 (S&A p. 380).

### 2.2 Tag

```
tag⁺_s ⇔ T_s >  θ_tag,      tag⁻_s ⇔ T_s < −θ_tag,      θ_tag = 0.2·w0        (4)
```

- The tag is a threshold on the early-phase change, as in L&T 2021 and Clopath et al. 2008.
- It decays only through Eq. 3. There is no separate tag variable.
- Two optional extensions:
  - `tag_reset`: low-frequency input within t_w of induction multiplies T by ρ (plasticity
    model Eq. 10).
  - A separate tag variable (Ziegler et al. 2015), needed only if the dissociation between tag and
    early LTP becomes a target.

### 2.3 PRP pool (slow)

```
τ_P dP_k/dt = −P_k + α_P · Θ[ Σ_{s∈k} |T_s| − θ_pro(M_k) ] − γ_c · C_k        (5)

θ_pro(M) = θ_pro,0 / (M + ε)                                                 (6)
```

- **Trigger:** PRPs are made when enough early change has built up in the pool (Frey & Morris
  1997; Clopath 2008; L&T 2021 Eq. 15). The neuromodulator lowers the threshold (Lehr 2022 Eq. 10:
  θ_pro = 1/(NM + 0.001)).
- **ε = 0.001** (Lehr 2022). With M = 0 the threshold is then ~1000 θ_pro,0, so in practice no
  PRPs are made without a modulator. This reproduces the CA3 doc's rule that capture needs a
  modulator (its Eq. 2) without a separate gate.
- **θ_pro,0 is a free scale.** L&T 2021 used a fixed θ_pro = 0.5 h0 (0.210037 nC). Lehr 2022
  replaced it with 1/(NM + 0.001). Choose θ_pro,0 so that θ_pro(M⁰) = 0.5 w0 at the modulator
  amplitude M⁰ used during encoding. This makes the encoding condition match L&T 2021.
- **C_k** is the consumption by capture, nonzero only when `pool = "cell_limited"` (§2.5).
- **Modulator reach** depends on the pool granularity (§2.5):
  `M_k(t) = Σ_m κ_m · M_m(t) · 𝟙[pool k is reached by m]`.

### 2.4 Capture (slow)

```
τ_z dz_s/dt =   P_k(s) · (1 − z_s)   · Θ[ T_s − θ_tag]
              − P_k(s) · (z_s + 0.5) · Θ[−T_s − θ_tag]                       (7)
```

- **Origin:** L&T 2021 Eqs. 13–14.
  - Late LTP needs a positive tag; late LTD needs a negative tag.
  - z saturates at 1 and −0.5.
  - The saturation fixes the plasticity model's open limit 6 (no upper bound on w_cap).
- **Time constants:** τ_z = 60 min (L&T 2021). z does not decay. This is acceptable over hours to
  days. Turnover (plasticity model Eq. 16) is a later extension.

### 2.5 Pool granularity: the switch that pattern 6 tests

**Table 2.** Which synapses share a PRP pool, and which modulator reaches it.

| `pool` | Pool k | Modulator m reaches pool k if | Literature |
|---|---|---|---|
| `"synapse"` | One per synapse s | s ∈ S_m (both cells in a pattern modulated by m; CA3 doc Eq. 2) | Purely local control; no STC between synapses |
| `"cell"` (default) | One per postsynaptic PYR j | j belongs to a pattern modulated by m (volume transmission reaches the cell) | L&T 2021; Clopath 2008; Ziegler 2015; Frey–Morris design |
| `"cell_limited"` | As `"cell"`, plus consumption C_k = Σ_{s∈k} \|dz_s/dt\| · τ_P / z_cost | As `"cell"` | Fonseca et al. 2004 (S&A p. 169); plasticity model Eq. 11 (D9) |

- **Dendrites:** the tinyCA3 PYR have no separate apical and basal compartments yet. A `"cell"` pool
  stands in for the dendritic compartment that, in CA1, limits capture (Alarcon et al. 2006,
  S&A p. 144).
- **Split pools:** `prp_pools = "split"` gives separate LTP and LTD pools (P⁺ fills z > 0, P⁻
  fills z < 0), as Sajikumar et al. 2005/2007 report (S&A p. 496).

### 2.6 Neuromodulators

M_m(t) is a sum of rectangular pulses of amplitude M_m⁰ and duration D_m. Two modes follow the two
dopamine systems in S&A Ch. 14 (pp. 272–273):

- **`da_mode = "VTA"`** (default for the CA3 doc):
  - DA is released while its pattern (Λ, and patterns 5 and 6) is presented.
  - With `replay_modulation = True`, DA is also released during SWR events that replay Λ.
  - Reach is pattern-specific (Table 2).
  - Source: McNamara et al. 2014, where VTA stimulation boosts later reactivation only of patterns
    present during the stimulation.
- **`da_mode = "LC"`** (the grace period):
  - DA starts at a novelty event and lasts D_DA. It reaches *all* PYR pools.
  - Any pattern encoded within the window can capture.
  - Source: Takeuchi et al. 2016, where LC-TH⁺ activation makes unrelated memories encoded close in
    time persist.

**5-HT** follows the same equations with its own κ_5HT. Benoy et al. 2022 (S&A p. 242) found, in
rat CA2, that 5-HT turns weak LTD into late LTD and allows cross-tagging. With
`prp_pools = "split"`, 5-HT therefore feeds P⁻ (and P⁺ through cross-tagging, κ⁺_5HT ≤ κ⁻_5HT).
This is our reading of one CA2 study, not a CA3 result.

## 3. The two-phase simulation loop

The simulation alternates between **online** phases (spiking, ms resolution) and **offline**
phases (slow variables only, between replay events). No timescale is compressed.

- **Why the offline phase is valid:** L&T 2021 show that early-phase plasticity does not occur
  below ~3 Hz firing (S&A p. 487). Between replay events the slow variables can therefore be
  integrated without spikes.
- **Check:** the background firing of tinyCA3 must be verified to stay below that rate. If the
  circuit is instead run with any compression factor, the factor must be stated (CLAUDE.md).

```
ALGORITHM  STC_with_replay(circuit, protocol, params)

  init:  T_s = 0, z_s = 0, P_k = 0, u_s = U, x_s = 1 for all P1 synapses;
         M_m(t) from protocol.modulator_schedule
         t = 0

  for each block b in protocol.blocks (ordered by start time):

    # ---- offline: advance slow variables to the start of the block ----
    advance_slow(t → b.t_start)                      # §3.1

    if b.type == ENCODE(k, n_pulses=3, pulse=0.1 s, gap=0.4 s):
        # online: EC_k → DG_k → mossy fibres impose pattern k
        run_online(duration = b.duration, drive = EC_k)   # §3.2

    if b.type == REPLAY_EPOCH(duration, rate r_SWR):
        for each SWR event e at times from Poisson(r_SWR) or schedule:
            advance_slow(t → t_e)
            R_e = run_swr_event(e)                     # §3.3
            apply_selection(R_e)                       # §3.4 (optional D3)
            log(e, R_e, per-pattern activation a_k)

    if b.type == RECALL_TEST(k, cue_fraction = 0.5):
        run_online(duration = 0.1 s, drive = half of pattern k's PYR directly)
        log Q_k = (ν_ans − ν_ctrl) / ν_as                # L&T 2021

    t = b.t_end

  advance_slow(t → protocol.t_end)
```

### 3.1 `advance_slow(t0 → t1)`: offline integration

With no spikes, c relaxes to 0 within ~0.2 s. Eq. 3 reduces to linear relaxation of T, and Eqs. 5
and 7 are ODEs in P and z. Integrate with dt_slow = 1 s, or exactly piecewise:

```
T_s(t)  = T_s(t0) · exp(−k_rel (t − t0) / τ_h)
P_k     : τ_P dP/dt = −P + α_P Θ[Σ|T| − θ_pro(M_k(t))] − γ_c C_k
z_s     : Eq. 7 with the current P_k(s) and the tag sign of T_s(t)
```

Event times inside [t0, t1] must be split at:

- the moments where Σ|T| crosses θ_pro(M_k) (PRP synthesis switches on or off);
- the moments where |T_s| crosses θ_tag (the tag is lost);
- the edges of the modulator pulses.

### 3.2 `run_online`: spiking phase (NEURON)

- Every P1 synapse runs Eqs. 1–4 per spike (TM release, calcium, T).
- P and z are held fixed during an online block (they change on a scale of minutes). Eq. 5 is
  evaluated at the end of the block.
- Background input: the Ornstein–Uhlenbeck current of L&T 2021 Eq. 9, with I_0 tuned so PYR fire
  at ~0.2–1.5 Hz (as in L&T; S&A p. 480). This rate is a target to check, not a CA3 value.

### 3.3 `run_swr_event(e)`: one replay event

```
1. Initiation   : for t_init (≈10 ms), add a brief depolarising current to every PYR
                  with amplitude drawn per cell (noise), or to a random seed set of n_seed
                  cells. No pattern is favoured by the drive itself.
2. Completion   : run_online(t_SWR ≈ 100 ms). Strong P1 weights recruit the rest of a
                  pattern; basket-cell feedback (P2/P3) limits the active set (soft k-WTA).
3. Loop         : P6 is live: if ≥ θ_EC of pattern k's cells fire, EC_k → DG_k fires and
                  re-imposes pattern k within the same event (CA3 doc §P6).
4. Readout      : a_k = fraction of pattern k's cells that spiked;
                  R_e = { k : a_k ≥ θ_EC }  (the replayed set; can be empty or multiple).
5. Plasticity   : already applied inside steps 2–3 by Eqs. 2–4. Co-firing in replayed
                  patterns raises c and resets or raises T (tag refresh).
6. Modulation   : if replay_modulation and a modulated pattern is in R_e, add a pulse
                  M_m += M_replay for D_replay (VTA-like; Kim & Park 2025 pp. 428–430).
return R_e
```

How each part affects replay:

- **Which patterns replay** depends on the stored weight within each pattern,
  W_k = mean over s ∈ k of (w0 + T_s + w0·z_s). W_k contains both the uncaptured early part (T)
  and the captured part (z), which is how STC enters replay.
- **Tag refresh** comes from steps 2–3. A replayed pattern's synapses reach c ≥ θ_p again if
  replay firing is dense enough. Reactivated memories carry tags again (S&A pp. 246–248), and
  reactivated firing with ripple input induces LTP in slices (Sadowski et al. 2016, cited in the
  plasticity model).
- **PRP supply by replay** happens when replay of a strong pattern drives Σ|T| in its cells above
  θ_pro(M). Overlapping weak patterns on the same cells can then capture. This is the circuit form
  of Cassini et al. 2013, where a reactivated strong memory rescues a weak one (S&A p. 246).

### 3.4 `apply_selection(R_e)`: optional depression of what is not replayed

With `swr_depression = True`, after each event and for every P1 synapse s that does not lie within
a pattern in R_e:

```
T_s ← T_s · (1 − δ_SWR),      z_s ← z_s · (1 − δ_SWR · f_z)                 (8)
```

- This is plasticity model Eq. 12 (D3; Norimoto et al. 2018).
- f_z ∈ [0, 1] sets how much ripple depression reaches captured weight. f_z = 0 protects captured
  synapses. That is the reading of Kim & Park 2025 (Fig. 2C): reactivated assemblies keep their
  strength while the rest is downscaled.

## 4. How the coupling can fail, and the checks that catch it

Replay that re-tags can sustain early LTP on its own, with no PRPs. That would contradict the
protein-synthesis dependence of consolidation. Two checks are mandatory before any result is
reported:

1. **PRP block** (anisomycin analogue, α_P = 0).
   - Run the full protocol, including replay.
   - **Pass:** every W_k returns to within 10% of w0 by 8 h, and replay of every pattern stops.
   - **Fail:** replay alone keeps a pattern alive. Lower the replay drive or the replay density.
     Alternatively, add the L&T 2022 condition that LTP needs pre- and postsynaptic rates above
     ν_th, which uses ripple-like bursts but not sparse replay firing.
   - This is the dashed trace of Fig. 2 in [hippocampal_plasticity_model.md](https://github.com/max-talanov/tinyHippo-bio-plasticity/blob/main/hippocampal_plasticity_model.md).
2. **Runaway.**
   - Capture → stronger W_k → more replay → more capture is a positive loop. It is bounded only by
     z ≤ 1 (Eq. 7), basket-cell inhibition and θ_EC.
   - **Pass:** no pattern takes over (> 90% of SWR events) unless it was the only one captured.
     The P6 loop gain must stay below 1 (CA3 doc).

## 5. Parameters

**Table 3.** Parameters. "L&T" values were fitted to rat CA1 slices at in vitro [Ca²⁺]ₒ and
rescaled by 0.6 for in vivo use. None is a CA3 value.

| Parameter | Value | Source | Status |
|---|---|---|---|
| τ_c, t_c,delay | 48.8 ms, 18.8 ms | L&T 2021 Table 2 (Graupner & Brunel 2012) | Published; CA1 slice fit |
| c_pre, c_post | 0.6, 0.1655 (network) | L&T 2021 Table 2 (×0.6, Higgins et al. 2014) | Published |
| θ_p, θ_d | 3, 1.2 | L&T 2021 | Published |
| γ_p, γ_d | 1645.6, 313.1 | L&T 2021 | Published |
| τ_h | 688.4 s | L&T 2021 | Published |
| k_rel | 0.1 (L&T, τ ≈ 1.9 h) or 0.265 (τ ≈ 43 min) | L&T 2021; Frey & Morris 1998 via S&A p. 380 | **Choice:** run both |
| θ_tag | 0.2 w0 | L&T 2021 | Published |
| θ_pro,0, ε | 0.5 w0 (scale), 0.001 | L&T 2021; Lehr 2022 | Scale free |
| τ_P | 60 min (L&T) or 1.4–2.9 h | L&T 2021; Frey & Morris 1998 via S&A p. 380 | **Choice:** run both |
| τ_z, α_P | 60 min, 1 | L&T 2021 | Published |
| w0 | Baseline P1 weight of the CA3 doc | — | Set so that a 50% cue does not complete a pattern before learning |
| M⁰_DA, D_DA | Free; Lehr 2022 used 30 and 60 min pulses | Lehr 2022 | Free |
| κ_DA, κ_5HT | Free | — | Free; 5-HT role from CA2 only |
| r_SWR, t_SWR, t_init, n_seed | Free | — | **Unsourced:** take from a primary ripple study before use |
| θ_EC | Free | CA3 doc | Proposal |
| δ_SWR, f_z | Free | Plasticity model Eq. 12 | Candidate |
| γ_c, z_cost | Free | Plasticity model Eq. 11 | Candidate |

## 6. Experiments and predictions

Each experiment runs ≥ 5 seeds (CLAUDE.md) and reports each pattern separately.

**Table 4.** Protocols. Defaults are `pool = "cell"`, `da_mode = "VTA"` and
`swr_depression = False`, unless stated.

| # | Protocol | Prediction | What it decides |
|---|---|---|---|
| X0 | Single synapse: STET, WTET, SLFS, WLFS, and the two-pathway weak-before/after-strong protocols, without the network | Reproduce Fig. 2 of L&T 2021 (late LTP/LTD only after strong stimulation; weak captured within the window) | S1 acceptance: verify against a numpy reference of Eqs. 2–7 |
| X1 | Encode patterns 1–6 once (the CA3 doc Table 4 design), no replay; recall at 10 s and 8 h | Λ, V captured; I, O decay to w0; pattern 6 depends on `pool` | Baseline without replay |
| X2 | As X1, plus replay epochs (e.g. 30 min every 1 h for 8 h) | Captured patterns dominate replay after ~2 h. I and O are replayed in the first ~1 h, then drop out | **STC → replay:** content of replay follows capture |
| X3 | PRP block (α_P = 0) with replay | No pattern survives at 8 h | Check 1 (§4) |
| X4 | Delayed modulation: encode Λ *without* DA. DA arrives only during replay that starts Δ after encoding, Δ ∈ {0, 15, 30, 60, 120, 240 min} | Capture for Δ below the tag lifetime. Replay that refreshes the tag widens the window compared with the same DA without replay | **Replay → STC:** does replay extend the tag window? |
| X5 | Behavioural tagging in the circuit: encode I (weak, no modulator); deliver DA at Δ ∈ ±{15, 30, 60, 120, 240} min, in `da_mode = "VTA"` vs `"LC"` | LC: I captured within the window, both before and after (symmetric; S&A p. 368). VTA: I captured only through cells shared with Λ (ex_13) | LC vs VTA reading of DA |
| X6 | Pattern 6 under the three `pool` settings | `synapse`: I never captured. `cell`: I→Λ/V-cell synapses captured, I→I-only (23, 33, 43) not. `cell_limited`: I captures less as Λ and V consume P | CA3 doc's pool test |
| X7 | `cell_limited`; encode Λ and O (they overlap in 5 cells); replay only Λ | O loses capture where it shares pools with Λ | Competition through replay (Fonseca 2004; Sajikumar 2014) |
| X8 | `swr_depression = True`, f_z ∈ {0, 0.5} | Non-replayed patterns are depressed. f_z = 0 keeps captured but non-replayed patterns | Selection by replay (D3) |
| X9 | Encode Λ, V, O in sequence 3 s apart vs 8 h apart (consolidate each before the next) | The 3-s spacing shows recency in spontaneous replay. The 8-h spacing shows more replay and less recency | L&T 2022 recency result in a 25-cell circuit |
| X10 | Priming: one 0.1-s cue to a captured pattern 8 h after encoding | Its replay share rises for ~1–2 h (raised T), then returns | L&T 2022 priming |

**Measures**, logged per pattern k:

- W_k(t), split into its early part (mean T) and its captured part (mean w0·z);
- the fraction of P1 synapses tagged;
- P_k(t);
- the replay share (the fraction of SWR events with k ∈ R_e);
- recall quality Q_k.

## 7. Implementation outline

Following CLAUDE.md, each mechanism is a dataclass plus a `build_*` function:

- **`STCSynParams` + `build_stc_synapse`**
  - An NMODL point process extending `TMSyn` with c, T and z (Eqs. 1–4, 7).
  - P_k(s) is read through a POINTER, so z can be integrated inside NEURON during online blocks.
- **`PRPPoolParams` + `build_prp_pools(circuit, pool)`**
  - One pool per synapse or per PYR (Eqs. 5–6).
  - Held in Python, written to the synapses' POINTER variable between blocks.
- **`ModulatorSchedule` + `build_modulators(da_mode, replay_modulation)`**: M_m(t) and the reach
  masks of Table 2.
- **`ReplayParams` + `run_swr_event` and `advance_slow`**: §3.1–3.4. `advance_slow` is pure numpy
  and is the reference against which the NMODL slow terms are verified (S1 rule).

## 8. Limits

1. **The synapse fit comes from another region and condition.** All plasticity values come from
   rat CA1 slices (Graupner & Brunel 2012), rescaled for in vivo calcium. There is no CA3 STC fit.
2. **The tag lifetime is not settled.** It ranges from a 30-min half-life (slices) to 2–3 h
   (behaviour). k_rel and τ_P are run at both ends.
3. **25 cells are not a ripple.** An SWR event here is a forced, brief reactivation. Ripple
   frequency, rate and duration are not modelled and have no source in this document yet.
4. **5-HT** rests on one CA2 study (Benoy et al. 2022) and is read as an LTD-side gate.
5. **One compartment per cell.** The apical–basal limit on capture (Alarcon et al. 2006) and
   clustering on single dendritic branches (Govindarajan et al. 2011) cannot appear.
6. **No reconsolidation or late maintenance tag.** The ~12 h BDNF-dependent maintenance phase
   (S&A p. 381) and the lability of reactivated captured weights (S&A pp. 246–248) are not in the
   equations. A captured weight cannot be destabilised by replay.

## Sources

- Luboeinski J, Tetzlaff C (2021) *Commun Biol* 4:275. doi:10.1038/s42003-021-01778-y
- Lehr AB, Luboeinski J, Tetzlaff C (2022) *Sci Rep* 12:17772. doi:10.1038/s41598-022-22430-7
- Sajikumar S, Abel T (eds.) (2024) *Synaptic Tagging and Capture*, 2nd ed. Springer.
  doi:10.1007/978-3-031-54864-2. Chapters 1, 7, 8, 12, 14, 15, 19, 23.
- Primary studies cited through the book or the papers above (not read here):
  - Frey & Morris 1997, 1998
  - Clopath et al. 2008
  - Ziegler et al. 2015
  - Graupner & Brunel 2012
  - Higgins et al. 2014
  - Moncada & Viola 2007
  - Takeuchi et al. 2016
  - McNamara et al. 2014
  - Cassini et al. 2013
  - Fonseca et al. 2004
  - Alarcon et al. 2006
  - Benoy et al. 2022
  - Sadowski et al. 2016
  - Norimoto et al. 2018
- Circuit: [tinyCA3 basic circuit architecture](ca3_basic_circuit_architecture.md)
- Plasticity model: [plasticity_up_down_model_2026-11-06.md](https://github.com/max-talanov/tinyHippo-bio-plasticity/blob/main/docs/plasticity_up_down_model_2026-11-06.md)
