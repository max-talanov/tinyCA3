# Synaptic tagging and capture (STC): a short literature review

*Copied from [tinyHippo-bio-plasticity](https://github.com/max-talanov/tinyHippo-bio-plasticity/blob/main/docs/synaptic_tagging_literature_review.md) and restricted to mammalian data: the invertebrate
studies (Aplysia, Drosophila, crayfish) discussed in the source reviews are left out.*

**Sources read:** the full texts of the sources below. Page numbers refer to the journal or book
pagination. Most are reviews, so the numbers quoted here are second-hand: each is attributed
to the primary study the review cites, and must be checked against that primary paper (with its
recording conditions: species, age, [Ca²⁺]ₒ, temperature) before it goes into `params.py`.

1. Martin KC, Kosik KS (2002). Synaptic tagging — who's it? *Nat Rev Neurosci* 3:813–820.
   doi:10.1038/nrn942. PMID 12360325.
2. Rogerson T, Cai DJ, Frank A, Sano Y, Shobe J, Lopez-Aranda MF, Silva AJ (2014). Synaptic tagging
   during memory allocation. *Nat Rev Neurosci* 15:157–169. doi:10.1038/nrn3667. PMID 24496410.
   PMCID PMC3992944.
3. Kim J, Park M (2025). Systems memory consolidation during sleep: oscillations, neuromodulators,
   and synaptic remodeling. *BMB Rep* 58(10):425–436. doi:10.5483/BMBRep.2025-0033. PMID 40962324.
   PMCID PMC12576410.
4. Sajikumar S, Abel T (eds.) (2024). *Synaptic Tagging and Capture: From Synapses to Behavior*,
   2nd ed. Springer. doi:10.1007/978-3-031-54864-2. Chapters 1, 7, 8, 12, 14, 15, 19 and 23 read.
   Cited as "S&A p. n".
5. Luboeinski J, Tetzlaff C (2021), *Commun Biol* 4:275; Lehr AB, Luboeinski J, Tetzlaff C (2022),
   *Sci Rep* 12:17772. Methods and parameter tables read. Cited as "L&T 2021" and "Lehr 2022".

## The STC hypothesis in one paragraph

Synapses activated by a weak stimulus get a local, transient, protein-synthesis-independent
**tag**. Strong stimulation anywhere on the same neuron triggers transcription and translation of
**plasticity-related products (PRPs)**, which are delivered throughout the cell but strengthen
only tagged synapses. So a weak input that alone gives early LTP (E-LTP) gets late LTP (L-LTP) if
a strong input on the same neuron overlaps it in time. The synapse, not the nucleus, becomes the
unit of lasting plasticity (Martin & Kosik, pp. 813–815, 819; Rogerson et al., p. 161 and Fig. 3).

## 1. Martin & Kosik (2002): what is the tag?

### Evidence for tagging in rat hippocampus (pp. 814–815)

- **Rat hippocampal slice, CA1, two independent inputs S1/S2** (Frey & Morris 1997, 1998):
  - 1 tetanic train gives LTP that decays after ~1.5 h.
  - 3 trains give LTP that lasts ≥8 h and depends on transcription.
  - A weak input expresses L-LTP if it is delivered 1–2 h before, or less than 2.5–3 h after, the
    strong stimulation.
  - Anisomycin during the S1 tetanus does not block persistent LTP in S2. The tag itself is
    protein-synthesis independent (Fig. 2c).
- **Two temporal features of the tag** (p. 815):
  - After strong stimulation, the cell can turn transient changes into persistent ones for hours
    (in rat CA1, a weak input up to 2.5–3 h later).
  - The tag set by a weak stimulus persists for 1–2 h.

### Criteria for a tag (p. 815)

A tag must:

1. be spatially restricted (local);
2. be time-limited and reversible;
3. interact with the cell-wide events that follow strong stimulation (usually gene expression) to
   give lasting, synapse-specific strengthening.

### Candidate tags (pp. 815–819)

| Class | Candidates and evidence |
|---|---|
| Fast synaptic changes | AMPA receptor clustering; lateral mobility of NMDA receptors; Homer-driven insertion of mGluR5; palmitate cycling on PSD-95 |
| Kinases | CaMKII, made autonomous by autophosphorylation. PKMζ, whose persistent activity needs protein synthesis, is necessary and sufficient for L-LTP maintenance. PKA, activated through D1/D5 receptors and cAMP in rodent hippocampus (S&A p. 294) |
| Adhesion molecules | NCAM, cadherins, β-catenin, δ-catenin |
| Actin network | Remodelling after repeated tetanus lasts ≥18 h (Colicos et al. 2001) |
| Local translation and degradation | mRNA released from RNA granules into polysomes on depolarisation. More spines contain polyribosomes after LTP. Dendritic targeting of *BDNF*/*TrkB* mRNA. *Arc* mRNA accumulates at stimulated sites. Local ubiquitin-mediated protein degradation |

### Key qualifications

- **Local translation is not a general tag** (pp. 818–819). It is not required for synaptic tagging
  at rodent hippocampal synapses (Frey & Morris 1997, 1998).
- **Captured products** (p. 819):
  - Largely unknown. MicroRNAs are a speculative candidate.
  - Mice expressing constitutively active CREB (VP16-CREB) in CA1 have a lower threshold for L-LTP,
    and their L-LTP becomes independent of transcription (Barco et al. 2002). This suggests the PRPs
    are CREB targets.
- **Conclusion** (p. 819):
  - There are probably several molecular tags. Different stimuli recruit different ones, and they act
    over different time domains, from hours to days.
  - A single tag molecule would be "extremely surprising". The common feature of all candidates is
    that they are spatially restricted and persistent.

**For the model:** a single abstract tag variable per synapse, decaying over ~1–2 h and
protein-synthesis independent, matches the rodent CA1 data this review summarises. The review
argues against identifying it with any one molecule.

## 2. Rogerson et al. (2014): tagging as part of memory allocation

### Neuronal allocation (pp. 157–160)

- **CREB raises excitability, and that biases which neurons store a memory.**
  - Lateral amygdala: neurons with virally raised CREB are ~3× more likely to be recruited into an
    auditory fear memory trace (Han et al. 2007).
  - Only ~⅓ of recorded lateral amygdala neurons show synaptic changes after fear conditioning
    (Rumpel et al. 2005).
  - CA1: virus transduced ~25% of neurons. Raising CREB before training, but not after, enhances
    the memory.
  - CA1 with constitutively active CREB: smaller afterhyperpolarisation, higher excitability, lower
    LTP threshold (Lopez de Armentia et al. 2007).
  - CA1 neurons recruited as place cells have lower spike thresholds (Epsztein et al. 2011).
    Preplay is consistent with excitability deciding allocation in advance.
- **Linking memories in time** (p. 160, Fig. 2):
  - A first memory raises excitability in the neurons that store it.
  - So, shortly afterwards, a second memory is stored largely in the same neurons. Recalling one then
    recalls the other.
  - As excitability decays, later memories no longer overlap.
  - Supported by Garner et al. 2012 (*Fos*-driven DREADDs).

### Synaptic allocation and STC (pp. 161–163, Figs 3–4)

- **Classic slice data** (Frey & Morris 1997, 1998), as reported by the review:
  - A weak S2 tetanus preceded 1 h earlier by L-LTP induction in S1 gives L-LTP in S2.
  - Weak tetanization of S2 that gives only E-LTP can come up to 1 h before the strong S1 tetanus
    and still be converted. The tag lasts up to ~1 h.
- **Tag and products work for LTD too** (Sajikumar & Frey 2004): short-lived LTD at one set of
  synapses becomes late LTD when L-LTP is induced at another set on the same neurons.
- **Molecules:**
  - Setting the tag needs CaMKII and actin remodelling.
  - Making the PRPs needs the CaMKII–CaMKIV–CREB pathway.
  - Candidate PRPs: ARC, GluR1, HOMER1A, PKMζ.
  - L-LTP can also be induced by low-frequency stimulation combined with dopamine D1/D5 agonists,
    BDNF or carbachol.
- **Crosstalk between spines (Fig. 4):**
  - LTP at one spine lowers the LTP threshold at neighbouring spines. The effect falls with both
    distance and time (Harvey et al. 2008).
  - Active RAS persists ~5 min and spreads ~10 µm. Active RHOA persists ~5 min and spreads ~5 µm
    (Harvey et al. 2008; Murakoshi et al. 2011).
  - In the Fig. 4 scheme, a weakly firing input must arrive within ~10 min.
- **Single-spine STC** (Govindarajan et al. 2011):
  - E-LTP at one spine is converted to L-LTP when L-LTP was induced at a nearby spine. The
    conversion needs protein synthesis and is temporally asymmetric.
  - There is little STC if the spines are >70 µm apart, and less across dendritic branches.
  - Simultaneous L-LTP at two nearby spines reduces spine growth. This suggests competition for a
    limited pool of PRPs.

### Metaplasticity and behavioural tagging (pp. 164–165)

- **Metaplasticity** links bouts of plasticity separated by minutes to days. Excitability changes
  driven by CREB are one proposed mechanism for heterosynaptic metaplasticity.
- **Behavioural tagging, inhibitory avoidance** (Moncada & Viola 2007):
  - A weak memory, lasting hours, becomes long-term, lasting days, when the rat explores a novel
    environment 1 h before training.
  - This depends on dopamine and on protein synthesis.
- **Behavioural tagging, food location** (Wang et al. 2010):
  - Novelty 30 min after a weak encoding episode makes the memory last.
  - Needs hippocampal D1/D5 receptors and protein synthesis.
  - The slice equivalent: strong tetanization converts E-LTP to L-LTP on a convergent pathway.
- **Counter-example** (Izquierdo et al. 1999):
  - A novel open field 1 h *after* inhibitory avoidance training impairs the memory.
  - There is no effect at 5 min or at 6 h, or when the field is familiar.
  - The exact parameters (strength, interval, novelty) matter.

### Clustering, synergy, open questions (pp. 165–167)

- **Spine clustering:**
  - Spines added during motor learning cluster on dendrites, and clustered spines are more stable
    (Fu et al. 2012).
  - Molecules involved: RHOA, CDC42, CaMKII, TIAM1, β-PIX.
- **Synaptic and neuronal allocation work together** (Fig. 3):
  - For ≥1 h after strong training, neuronal ensemble 1 has PRPs that can turn a weak memory into a
    strong one, so STC also decides which *neurons* store a memory.
  - Conversely, local excitability can bias which dendrites take part.
- **The authors' own caveat:** there is still no direct evidence that the behavioural interactions
  actually use STC. They might instead reflect protein-synthesis-dependent rises in neuromodulators
  such as dopamine.

**For the model:**

- STC decides where changes are stored, not only whether they last.
- The paper suggests these model features:
  - tags and PRP pools that are local to a branch or spine neighbourhood, with distance limits of
    ~10 µm (fast crosstalk) and ~70 µm (STC);
  - a limited PRP pool that synapses compete for;
  - dopamine (D1/D5) as a trigger for PRP synthesis;
  - neuron-level excitability as an allocation variable separate from the tag.

## 3. Kim & Park (2025): sleep, neuromodulators and synaptic remodelling

**Not an STC paper.** It never discusses synapse-level tagging. "Emotional tagging" (pp. 425, 427)
means REM sleep prioritising emotional memories. What it adds is the setting in which tagged
changes would be consolidated or downscaled.

- **Oscillations** (p. 426):
  - slow oscillations, 0.1–4 Hz;
  - spindles, 10–15 Hz, lasting 0.5–2 s;
  - sharp-wave ripples, 150–250 Hz, from CA3 to CA1;
  - REM theta, 5–8 Hz.
- **NREM:** coupling of slow oscillations, spindles and ripples drives replay and transfer from
  hippocampus to cortex (p. 427).
- **Slow oscillations vs delta waves:** they play opposing roles. Slow oscillations promote
  consolidation and delta waves promote forgetting (Kim et al. 2019).
- **Neuromodulators** (pp. 427–429):
  - In NREM, noradrenaline oscillates at ~0.02 Hz. Bursts of the locus coeruleus (LC) end individual
    spindles, while sustained low noradrenaline lets spindles cluster (Fig. 2B caption).
  - Noradrenaline is lowest in REM.
  - Suppressing noradrenaline too early, or raising it during NREM, impairs consolidation.
  - Dopamine surges in the nucleus accumbens and basolateral amygdala at the NREM→REM transition.
  - Dopamine from the VTA during NREM aligns with memory-related activity and may decide which
    memories are prioritised.
- **Synaptic regulation** (pp. 429–430):
  - Synaptic Homeostasis Hypothesis (Tononi & Cirelli): synapses potentiate overall during waking and
    are downscaled during sleep. Axon–spine contact area and synaptic AMPA receptors fall.
  - Homer1a drives this downscaling (Diering et al. 2017). This is the same molecule as the candidate
    PRP HOMER1A in Rogerson et al.
  - Sleep also *forms* and *preserves* synapses. Sleep promotes branch-specific spine formation after
    learning (Yang et al. 2014). REM sleep selectively prunes and maintains new spines (Li et al.
    2017).
  - Fig. 2C caption: NREM selectively downscales highly active neurons, REM weakens the network more
    broadly, and assemblies reactivated in NREM keep their strength.
  - The main text instead says NREM mainly potentiates and REM mainly prunes (p. 430). The review is
    not fully consistent on how NREM and REM divide the work.

**For the model:** relevant to the downscaling pathways (D1–D13) and to a later consolidation stage,
not to S4–S6 directly. It gives a reason to gate PRP synthesis or capture by neuromodulators
(dopamine, noradrenaline), and a reason to protect tagged or captured synapses from global
downscaling.

## 4. Sajikumar & Abel (eds., 2024): *Synaptic Tagging and Capture*, 2nd ed.

Springer. doi:10.1007/978-3-031-54864-2. A 23-chapter book; the chapters relevant to the model
were read in full text. Pages are the printed book pages.

### Timing of tag and PRPs

- **Ch. 19** (Viola, Ballarini & Moncada, p. 380):
  - Frey & Morris (1998) measured a tag **half-life of ~30 min** and a PRP half-life of
    **1 to 2 h**.
  - These are the first explicit half-lives in this review. They correspond to τ_tag ≈ 43 min and
    τ_P ≈ 1.4–2.9 h (τ = t½ / ln 2).
- **Ch. 1** (Prodan & Morris, p. 30): STC studies broadly agree that tags last approximately 1 h
  at 32 °C. This is the only temperature given for a tag lifetime.
- **Behavioural tags last longer:**
  - From the BT experiments as a whole, learning tags have a half-life of 30 min to 2–3 h,
    depending on the task (Ch. 19, p. 380).
  - Inhibitory avoidance (IA): a 5-min novel open field promotes long-term memory about 1 h before
    or after weak training. It does not between 30 min before and the moment of training (Moncada &
    Viola 2007; Ch. 19, p. 368). The window is symmetric, as in slices.
  - Spatial object recognition: the window is 1 h before to 2 h after.
  - Conditioned taste aversion: 1 h before or 2.5 h after.
  - Contextual fear: 3 h before to 1 h after (Chen et al. 2020).
- **A second, late window** (Ch. 19, p. 381): memory persistence at 7 days needs a BDNF- and
  protein-synthesis-dependent phase about 12 h after training (Bekinschtein et al. 2007). The
  authors propose a late "maintenance tag".
- **Sleep:** sleep deprivation impairs late-LTP, disrupts STC, and blocks behavioural tagging of
  weak IA by novelty (Wong et al. 2020; Ch. 19, p. 377).

### Where capture can happen, and competition

- **Compartments (Ch. 7,** Alarcon, p. 144**):**
  - Capture works between two inputs on apical dendrites, and between two on basal dendrites.
  - It fails between an apical and a basal input (Alarcon et al. 2006).
  - LTD is compartmentalised too.
  - So the tag/PRP unit is a dendritic compartment, not necessarily the whole cell.
- **Competition (Ch. 8,** Fonseca, p. 169**):**
  - When PRPs are limited, tagged groups compete. This was shown with protein-synthesis inhibitors
    and also with natural patterns that only limit the initial pool (Fonseca et al. 2004).
  - Reactivating one previously activated group adds tags and puts the other groups under
    pressure.
  - Competition is proportional to the potentiation at the winning input. Tag "activity" scales
    with activation, and stronger tags capture more.

### Dopamine and serotonin as the PRP gate

- **Dopamine (Ch. 15,** Navakkode, p. 294**):** D1/D5 receptors enable cellular consolidation by
  inducing protein synthesis through cAMP/PKA. They act as a gate.
  - L-LTP is impaired in D1-receptor knockout mice.
  - A D1/D5 agonist alone gives slow-onset potentiation that lasts ~6 h *in vitro* and takes part
    in associativity.
- **Two dopamine systems (Ch. 14,** Takeuchi, pp. 272–273**):**
  - **LC → hippocampus, the "grace period":** optogenetic activation of LC-TH⁺ (not VTA-TH⁺)
    neurons mimics novelty. It makes spatial memories encoded close in time persist, and the effect
    is sensitive to D1/D5 (Takeuchi et al. 2016). This works over time, across memories: a
    cell-wide PRP gate.
  - **VTA → hippocampus, pattern-specific replay:** stimulating VTA-TH⁺ axons while the animal
    explores a novel environment enhances later SWR reactivation, and memory persistence. Only the
    activity patterns present *during* the stimulation benefit, not those immediately before
    (McNamara et al. 2014). This is a pattern-specific gate that acts through replay.
- **Serotonin (Ch. 12,** Bin Ibrahim & Sajikumar, p. 242**):**
  - In rat CA2, weak low-frequency stimulation that alone gives transient LTD becomes long-lasting
    LTD when paired with serotonin.
  - This serotonin-facilitated LTD takes part in tagging and capture with a second input, and in
    *cross*-tagging: weak LTP in another input becomes lasting (Benoy et al. 2022).
  - This is the first source in this review for a role of 5-HT in STC. It points to LTD rather
    than LTP, and comes from CA2, not CA3.

### Reactivation and reconsolidation (Ch. 12, pp. 246–248)

- **A reactivated strong memory supplies PRPs:** briefly re-exposing an animal to a fear-conditioned
  context right after weak spatial object recognition training makes the weak memory persist
  (Cassini et al. 2013).
  - It needs the strong memory to become labile: blocking L-type calcium channels, or shortening
    re-exposure to 90 s, abolishes the effect.
  - It needs NMDA receptors during the weak task, to set the tag.
- **Reactivated memories carry tags again** and need PRPs to reconsolidate. Blocking protein
  synthesis during reactivation causes loss, unless PRPs were supplied by an earlier strong event
  (Rabinovich Orlandi et al. 2020).
- **Proposed:** synapses of a new memory and of a reactivated one compete, winner-take-all
  (Sajikumar et al. 2014).

### Computational models (Ch. 23, Luboeinski & Tetzlaff, pp. 471–498)

**Common ground** (p. 474):

- All models reproduce the slice data: tag and early LTP/LTD decay over ~1–2 h, and late changes
  last many hours.
- All of them, except the compartment models, assume an unlimited PRP pool for the whole neuron.

| Model | Tag | PRP trigger | Notes |
|---|---|---|---|
| Clopath et al. 2008 | Three Markov states (neutral/up/down); a tag is the early state | Enough early change, with a dopamine-dependent threshold; PRPs decay over hours | Three needed parts: early weight/tag, PRP, late weight |
| Barrett et al. 2009 | Six states; early change can occur without a tag | — | Depotentiation and tag reset |
| Ziegler et al. 2015 | Separate tag variable; weight, tag and scaffold are all bistable | Dopamine only (one value per neuron) | Slow-onset LTP; behavioural tagging in feed-forward fear model, including no consolidation when novelty coincides with the shock |
| O'Donnell & Sejnowski 2014 | Exponential tag | Alpha-function PRP, local to dendrite | CA3→CA1 feed-forward; weak inputs must share segments with strong ones; clustering |
| Kastellakis et al. 2016 | Calcium-based | Alpha-function, somatic and/or dendritic | Behavioural tagging; homeostasis |
| Päpper et al. 2011 | Abstract | — | First recurrent STC network; multi-timescale STC raises storage capacity |
| Li et al. 2016 → Luboeinski & Tetzlaff 2021, 2022; Lehr et al. 2022 | Early weight beyond a threshold | Summed early change above θ_pro; θ_pro depends on a neuromodulator | Recurrent spiking network; the basis of the algorithm in [stc_replay_algorithm_tinyCA3.md](stc_replay_algorithm_tinyCA3.md) |

**Results relevant to replay:**

- **Consolidation and improvement** (pp. 488–490):
  - STC consolidates cell assemblies, and recall 8 h later can be *better* than at 10 s.
  - "Passive" improvement comes from the late weights. "Active" improvement comes from early
    plasticity during recall.
- **Neuromodulation acts retroactively** (pp. 490–491):
  - Low neuromodulator: no consolidation.
  - Moderate: good recall of the patterns as they were presented (input-defined).
  - High: the assembly grows to recruit new cells, and recall of self-organised patterns and
    temporal traces is good.
  - Neuromodulation works only if it falls within the tag window. Pattern completion peaks with
    neuromodulation ~60–90 min after learning, as in behavioural tagging.
- **Interaction between memories** (pp. 491–493):
  - LTD makes newer assemblies weaken older ones, so the newer ones reactivate spontaneously more
    often (recency).
  - Consolidating each assembly before learning the next gives more reactivation and less recency.
  - A brief priming stimulus raises the early weights and spontaneous reactivation of one assembly
    for a few hours.
- **Open points named by the authors** (pp. 494–496):
  - Tag reset and the partial dissociation between tag and early LTP (as in Ziegler 2015).
  - Limited, consumed PRP pools and competition.
  - Separate PRP pools for LTP and LTD (Sajikumar et al. 2005, 2007).
  - A mix of somatic and dendritic protein synthesis.

## 5. Luboeinski & Tetzlaff (2021); Lehr, Luboeinski & Tetzlaff (2022)

- Luboeinski J, Tetzlaff C (2021). Memory consolidation and improvement by synaptic tagging and
  capture in recurrent neural networks. *Commun Biol* 4:275. doi:10.1038/s42003-021-01778-y.
  PMCID PMC7977149.
- Lehr AB, Luboeinski J, Tetzlaff C (2022). Neuromodulator-dependent synaptic tagging and capture
  retroactively controls neural coding in spiking neural networks. *Sci Rep* 12:17772.
  doi:10.1038/s41598-022-22430-7. PMCID PMC9588040.

Read in full text (Methods and Tables 1–2 of the 2021 paper; Methods of the 2022 paper). These
papers give the only complete, published equation set for STC in a recurrent spiking network.

**Equations** (2021 Eqs. 10–15; 2022 Eqs. 9–10):

```
w = h + h0·z                                                              (E→E only)
τ_h dh/dt = 0.1(h0 − h) + γ_p(1 nC − h)Θ[c − θ_p] − γ_d h Θ[c − θ_d] + ξ(t)
dc/dt     = −c/τ_c + c_pre Σ δ(t − t_pre − t_c,delay) + c_post Σ δ(t − t_post)
tagged    ⇔ |h − h0| > θ_tag
τ_z dz/dt = p (1 − z) Θ[(h − h0) − θ_tag] − p (z + 0.5) Θ[(h0 − h) − θ_tag]
τ_p dp/dt = −p + α Θ[ Σ_j |h_j − h0| − θ_pro ]                          (one p per neuron)
θ_pro(NM) = 1 / (NM + 0.001)                                              (2022)
```

**Parameters** (Table 2 of the 2021 paper):

| Parameter | Value | Note |
|---|---|---|
| τ_c | 48.8 ms | Graupner & Brunel 2012 fit to slice data |
| c_pre, c_post | 1, 0.2758 (network: 0.6, 0.1655) | Network values scaled by 0.6 for in vivo [Ca²⁺]ₒ (Higgins et al. 2014) |
| t_c,delay | 18.8 ms | — |
| θ_p, θ_d | 3, 1.2 | — |
| γ_p, γ_d | 1645.6, 313.1 | — |
| τ_h | 688.4 s | With the 0.1 factor, h relaxes with τ_h/0.1 = 6884 s ≈ 1.9 h |
| σ_pl | 0.290436 nC s^½ | Noise in plasticity |
| h0 | 0.420075 nC | — |
| θ_tag | 0.0840149 nC | = 0.2 h0 |
| θ_pro | 0.210037 nC | = 0.5 h0 |
| τ_p, τ_z | 60 min, 60 min | — |
| α | 1 | — |

**Network and protocol:**

- 1600 excitatory and 400 inhibitory LIF neurons, connection probability 0.1, Δt = 0.2 ms.
- Learning: 3 pulses of 0.1 s, 0.4 s apart. Recall: one 0.1-s pulse to half of the assembly.
- Pattern completion is measured as Q = (ν_ans − ν_ctrl)/ν_as.
- The slow phase is computed without spiking. This is valid because the early weights do not
  change below ~3 Hz firing (Ch. 23, p. 487).

## The tagging mechanism, step by step

Assembled from the three reviews. "M&K" is Martin & Kosik 2002, "R" is Rogerson et al. 2014, "K&P" is
Kim & Park 2025.

1. **Tag setting (weak input).**
   - Activity that is too weak for lasting plasticity still produces E-LTP (or early LTD) and sets
     a tag at the activated synapses.
   - Setting the tag needs no new protein (anisomycin does not block it; M&K Fig. 2c).
   - It involves CaMKII and actin remodelling (R p. 163). Other candidates are PKMζ, PKA, adhesion
     molecules, ion channels, and activated local translation (M&K pp. 815–818).
   - The tag is local, transient and reversible (M&K p. 815). It is probably a synaptic state made of
     several components rather than one molecule (M&K p. 819).
2. **Tag decay.**
   - Without capture the tag decays within ~1–2 h (M&K p. 815) or ~1 h (R p. 162).
   - The E-LTP it supports also decays: 1 tetanic train gives LTP that is gone after ~1.5 h
     (M&K Fig. 2b).
3. **PRP synthesis (strong input or neuromodulator).**
   - Strong stimulation of any input on the neuron (e.g. 3 tetanic trains) activates the
     CaMKII–CaMKIV–CREB pathway and transcription (R p. 163; M&K p. 814).
   - The same can be triggered by weak stimulation combined with dopamine D1/D5 agonists, BDNF or
     carbachol (R p. 163).
   - Products (candidates: ARC, GluR1, HOMER1A, PKMζ; R p. 163) are delivered throughout the cell
     and remain available for hours (M&K p. 814).
4. **Capture.**
   - A synapse that still carries its tag when PRPs are present converts E-LTP into
     protein-synthesis-dependent L-LTP, lasting ≥8 h in slices (M&K Fig. 2).
   - The order can go either way, weak before strong or strong before weak, within the windows in
     the synthesis table below.
   - Untagged synapses ignore the PRPs. That keeps the change input-specific even though the PRPs
     spread across the whole cell (R p. 161).
5. **Spatial and resource limits.**
   - Fast crosstalk: RAS and RHOA spread ~5–10 µm for ~5 min after LTP at one spine (R p. 163).
   - STC itself: weak beyond ~70 µm, and weaker across dendritic branches.
   - When two nearby spines are potentiated together they compete for a limited PRP pool
     (Govindarajan et al. 2011; R p. 164).
6. **Bidirectionality.** Early LTD is captured into late LTD in the same way, including by PRPs from
   L-LTP at other synapses (cross-capture; R p. 162).

```
weak input ──► E-LTP + TAG (local, no new protein, decays ~1–2 h)
                    │
                    ├── no PRPs before the tag decays ──► weight returns to baseline
                    │
strong input / ──► CaMKII–CaMKIV–CREB ──► transcription ──► PRPs (cell-wide, limited)
DA (D1/D5),                                                    │
BDNF, ACh                                                      ▼
                    └── tag still present when PRPs arrive ──► CAPTURE ──► L-LTP (stable)
```

## Influence on memory consolidation

The two meanings of "consolidation" follow
[hippocampal_plasticity_model.md](https://github.com/max-talanov/tinyHippo-bio-plasticity/blob/main/hippocampal_plasticity_model.md) (Terminology section).

### Cellular consolidation

- **STC is the proposed synaptic correlate.** It turns a transient, protein-synthesis-independent
  change into a lasting, protein-synthesis-dependent one, the same transition that defines cellular
  consolidation. That is a hypothesis, not an identity.
- **Behavioural tagging:** the evidence that the same rules work at the level of behaviour (R pp. 164–165).
  - A weak inhibitory-avoidance memory, lasting hours, became long-term, lasting days, when the rat
    explored a *novel* environment 1 h before training. A familiar environment did not do this. The
    effect needed dopamine and protein synthesis (Moncada & Viola 2007).
  - The same held for spatial object recognition, contextual fear conditioning and conditioned taste
    aversion (Ballarini et al. 2009).
  - A weakly encoded food location became long-lasting when novelty followed 30 min later. This
    needed hippocampal D1/D5 receptors and protein synthesis. In slices, the matching strong tetanus
    converted E-LTP to L-LTP on a convergent pathway (Wang et al. 2010).
- **What the model takes from this:**
  - A strong or novel event within the window supplies the PRPs that rescue a weak memory.
  - Dopamine is the gate that links novelty to PRP synthesis.
- **Limits:**
  - Novelty 1 h *after* inhibitory-avoidance training *impaired* the memory. There was no effect at
    5 min or at 6 h (Izquierdo et al. 1999; R p. 165).
  - The authors say there is no direct evidence that behavioural tagging actually runs through STC,
    rather than through protein-synthesis-dependent changes in neuromodulators (R p. 167).

### Memory allocation and linking

- **STC decides where a memory is stored** (R pp. 161–166, Fig. 3):
  - For ≥1 h after a strong memory, the neurons that stored it hold PRPs.
  - A weak memory arriving then gets L-LTP in those neurons but not in others.
  - So STC selects the storing *neurons*, not just the synapses.
- **Excitability works alongside it.** CREB-driven excitability makes recently active neurons more
  likely to store the next memory. Memories formed close together in time therefore overlap and are
  linked, so recalling one recalls the other (R Fig. 2).
- **Clustering:** crosstalk over short distances and the ~70 µm STC limit favour clusters of
  potentiated spines on one branch. Clustered new spines are more stable (Fu et al. 2012; R p. 165).
- **Metaplasticity:** one memory changes the plasticity thresholds that apply to the next, over
  minutes to days (R p. 164).

### System consolidation and sleep

None of the three reviews links STC directly to sleep. Kim & Park describe the setting in which
captured changes are then consolidated or pruned:

- **Replay:**
  - In NREM sleep, sharp-wave ripples, coupled to spindles and slow oscillations, replay recent
    experience and drive transfer from hippocampus to cortex (K&P p. 427).
  - The cortical side then needs its own cellular consolidation.
- **Neuromodulator gating:**
  - VTA dopamine during NREM lines up with memory-related activity and may decide which memories are
    prioritised (K&P pp. 428–430).
  - Noradrenaline from the locus coeruleus oscillates at ~0.02 Hz and controls when spindles occur.
  - Dopamine is the same signal that gates PRP synthesis in behavioural tagging. So it is a plausible
    common gate for capture during waking and for prioritisation during sleep. This is our inference,
    not a claim of the reviews.
- **Downscaling vs preservation:**
  - Sleep weakens synapses overall (Synaptic Homeostasis Hypothesis; Homer1a-driven scaling,
    Diering et al. 2017).
  - Assemblies reactivated in NREM keep their strength (K&P Fig. 2C).
  - Sleep promotes branch-specific spine formation after learning (Yang et al. 2014), and REM
    selectively prunes or keeps new spines (Li et al. 2017).
  - **For the model:** captured, replayed synapses should be protected from global downscaling;
    uncaptured ones should not be.
- **Timescale:**
  - Capture needs PRPs and tags within hours of each other. Only replay in the first hours after
    learning can drive capture.
  - Replay over days is the proposed driver of system consolidation, not of STC.

## Synthesis for S4–S6 (tag/capture)

| Element | What the three reviews report | Source |
|---|---|---|
| Tag | Local, transient, needs no new protein; CaMKII and actin; probably several tags, not one molecule | M&K pp. 815, 819; R p. 163 |
| Tag lifetime (CA1) | ~1–2 h (M&K); "up to 1 h" (R); ~1 h at 32 °C (Ch. 1); **half-life ~30 min** (Frey & Morris 1998) | M&K p. 815; R p. 162; S&A pp. 30, 380 |
| PRP lifetime | **Half-life 1–2 h** (Frey & Morris 1998); τ_p = 60 min in the L&T model | S&A p. 380; L&T 2021 Table 2 |
| Behavioural tag | Half-life 30 min to 2–3 h, depending on the task; symmetric windows of ~1 h (IA) to 1–2.5 h | S&A pp. 368, 380 |
| Compartment | Capture within apical or within basal dendrites, not between them | S&A p. 144 |
| Dopamine systems | LC: grace period across time (cell-wide). VTA: boosts later replay of the patterns present during release (pattern-specific) | S&A pp. 272–273 |
| Serotonin | Converts weak LTD to late LTD in CA2 and enables cross-tagging | S&A p. 242 (Benoy et al. 2022) |
| Reactivation | A reactivated strong memory supplies PRPs to a weak one; reactivated memories need PRPs again | S&A pp. 246–248 |
| Weak → strong window (CA1) | Weak 1–2 h before strong (M&K); up to 1 h before (R) | M&K p. 814; R p. 162 |
| Strong → weak window (CA1) | Weak <2.5–3 h after strong (M&K); ≥1 h after (R, Fig. 3) | M&K p. 814; R p. 162 |
| E-LTP / L-LTP decay | 1 train decays after ~1.5 h; 3 trains persist ≥8 h | M&K Fig. 2 |
| PRPs | Need transcription and translation (CaMKII–CaMKIV–CREB); candidates ARC, GluR1, HOMER1A, PKMζ; limited pool, competition | R pp. 163–164 |
| Spatial range | Fast crosstalk ~5–10 µm for ~5 min (RAS/RHOA); STC weak beyond ~70 µm and across branches | R pp. 163–164 |
| LTD | Tags and products also convert early LTD into late LTD (cross-capture) | R p. 162 |
| Neuromodulation | Dopamine D1/D5 needed for behavioural tagging; dopamine and noradrenaline in sleep gate which memories are prioritised | R pp. 163–165; K&P pp. 428–430 |
| Allocation | CREB-driven excitability chooses neurons; STC chooses synapses (and neurons) | R pp. 158–166 |
| Sleep | Homeostatic downscaling (Homer1a) and targeted preservation and pruning of synapses | K&P pp. 429–430 |

### Caveats

- **The two reviews give different time windows.** Martin & Kosik give a 1–2 h tag and accept a
  weak input up to 3 h after the strong one. Rogerson et al. cite 1 h in both directions. Both cite
  Frey & Morris 1997/1998, so take the window from those primary papers, together with their
  recording conditions.
- **No recording conditions are given** (temperature, [Ca²⁺]ₒ, age) for any quoted number. All
  values are second-hand from reviews.
- **Behavioural tagging may not use STC.** Rogerson et al. say directly that it has not been shown
  to.
