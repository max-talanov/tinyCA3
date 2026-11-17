# tinyCA3 results

One numbered section per result. Stages and acceptance criteria are in
[docs/implementation_plan.md](docs/implementation_plan.md) §4. All runs so far are deterministic (no noise), so
repeating them adds nothing; seeds start to matter at T6.

## 1. T0: environment and the STP mechanism (2026-11-17)

**Setup.**

- Python 3.13.3 venv with pip `neuron==9.0.2`, numpy, scipy, matplotlib and pytest (`requirements.txt`).
- `PYTHONPATH` must be unset. The NEURON 8.2.7 app bundle cannot be imported under Homebrew's Python 3.14, and it
  shadows pip's NEURON when it is on `PYTHONPATH`.
- `nrn_setup.load_mechanisms` compiles `neuron/` on first use. That folder holds the existing
  `TMExp2SynNMDACA1_SCExc.mod` and 13 Bezaire et al. 2016 channel files, copied unchanged from
  tinyHippo-bio-plasticity (CC0; provenance in `neuron/BEZAIRE_CHANNELS.md`).

**Acceptance** (`tests/test_t0_environment.py`, 5 tests, all pass):

- All 14 mechanisms load.
- The Tsodyks–Markram release per spike at 10, 20, 40 and 80 Hz (10 spikes each) matches the numpy reference
  `reference/tm_numpy.py` to < 1e-9.
- At 80 Hz, release relative to spike 1 is 1.000, 1.535, 1.598, 1.406, 1.165 for spikes 1–5. The same values were
  obtained earlier under NEURON 8.2.7.

**Status:** T0 accepted.

## 2. T1: Hodgkin–Huxley-type cells (2026-11-17)

**Cells** (`cells.py`; all at 34 °C, dt 0.025 ms):

| Type | Used for | Build |
|---|---|---|
| PYR | 25 CA3 pyramidal cells | Soma (20 × 20 µm), one apical dendrite (300 × 3 µm), axon initial segment (50 × 1 µm). Bezaire pyramidal channels at the template's somatic densities (`params.BEZ_PYR_SOMA`) |
| BC | 25 basket cells | The published Bezaire PV basket cell: 17 sections, seven channel types |
| RS | EC_k, DG_k, SWT | One compartment (20 × 20 µm) with the pyramidal somatic channels |

**Measurements** (`run_t1_cells.py` → `results/t1_cells.json`, `results/t1_fI.png`). All are from DC current
injection into the soma: 500-ms steps after 500 ms of settling. Rin comes from a −20 pA step.

| Type | V_rest (mV) | R_in (MΩ) | Rheobase (pA) | Rate at 1 / 1.5 / 2 / 3 / 5 × rheobase (Hz) | Last/first ISI at 5× |
|---|---|---|---|---|---|
| PYR | −63.1 | 168.3 | 75.6 | 2 / 22 / 30 / 42 / 60 | 1.04 |
| BC | −65.0 | 52.2 | 288.2 | 2 / 38 / 56 / 86 / 130 | 1.03 |
| RS | −68.4 | 2107.6 | 16.8 | 2 / 24 / 32 / 46 / 62 | 1.03 |

**Acceptance** (`tests/test_t1_cells.py`, 9 tests, all pass):

1. **BC port against the original template.** `build_bc` agrees with `reference/bezaire2016/class_pvbasketcell.hoc`
   within 0.01 mV, with equal spike counts, for −50 pA and +500 pA steps.
2. **Rheobase and repetition.** Every type has a DC rheobase below 2 nA and fires ≥ 5 spikes at 2× rheobase.
3. **BC fast spiking.** BC fires 130 Hz at 5× rheobase, with last/first ISI 1.03. The criterion was > 100 Hz with
   little accommodation.
4. **RS single spikes.** An RS unit fires exactly one spike to one AMPA-only synaptic event for NetCon weights
   0.003–0.016.
   - With the NMDA component on, the same event fires a train: 3 spikes at weight 0.004 and 7 at 0.008. The 148-ms NMDA
     decay keeps depolarising the cell after the first spike.
   - So every input onto EC, DG and SWT units is AMPA-only (gmax_nmda = 0). The plan's Tables 2–3 are updated.

**What these cells are not.**

- **No CA3 values.** All densities are Bezaire CA1 model values, used as placeholders. Not one value is a CA3
  measurement.
- **The PYR does not adapt or burst.** Last/first ISI is ≈ 1.0, while real CA3 pyramidal cells adapt and many burst.
  The plan lists I_CaL/I_KCa as a later extension.
- **The PYR's HCN and A-type currents are flat.** They sit at their somatic densities along the dendrite, without the
  template's distance gradients.
- **The RS unit is extremely compact.** Rin is 2.1 GΩ and rheobase 17 pA, so it is very sensitive to background
  input. That is acceptable for detonator-type input and readout units, but background drive onto RS units must be
  calibrated in T3, not copied from the PYR.

**Status:** T1 accepted against the plan's criteria. The CA3-specific intrinsic properties above are open.
