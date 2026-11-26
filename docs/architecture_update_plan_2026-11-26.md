# tinyCA3: plan of the architecture update (disjoint patterns, layout B)

Nov 26, 2026 · plan, nothing built yet

## 1. The problem

### 1.1 Pattern contamination

The current circuit ([ca3_basic_circuit_architecture.md](ca3_basic_circuit_architecture.md), Table 3; `circuit.PATTERNS`)
uses four basic patterns, Λ, V, I and O, and these share cells. The PRP pool belongs to the postsynaptic cell (`stc.py`:
`P[post_idx]`). The modulator reaches whole cells by pattern membership (`MOD_CELLS = {DA: Λ, 5-HT: V}`). So any
control synapse that ends on a Λ or V cell gets PRPs and captures.

**Table 1.** Shared cells and contaminated control synapses in the current layout (v1).

| Pair | Shared cells |
| --- | --- |
| Λ, V | 32, 34 |
| Λ, I / V, I | 13 / 53 |
| Λ, O / V, O | 13, 41, 45, 51, 55 / 11, 15, 21, 25, 53 |

| Control | Internal synapses | Ending on a Λ or V cell |
| --- | --- | --- |
| I | 20 | 8 (40%) |
| O | 240 | 150 (62%) |

So W_I and W_O mix synapses that are only tagged with synapses that captured. The clean "tag decays to w_0" result is
shown only by hand-picked probes (`run_pilot.PROBES`). Synapses onto 32 and 34 get both DA and 5-HT, so a DA-only or
5-HT-only effect cannot be read out cleanly.

### 1.2 The v1 patterns do not match the figure

Compared cell by cell with `figs/CA3_basic_circuit.png`, the Λ, V, I and O shapes in Table 3 of the architecture doc
are not the patterns drawn there. The figure has four disjoint basic patterns, and its bottom row holds two
combinations of two basic patterns each (Table 2). All results so far (reports from 2026-11-17 to 2026-11-24) were
obtained with the misread v1 patterns. They stay valid as experiments on v1, but they are not experiments on the
circuit in the figure.

## 2. The new layout (B): the patterns of the figure

The grid stays 5 × 5, with 25 PYR and 25 BC. The populations and projections P1–P8 are unchanged. Only the pattern
sets change.

```text
 1  X  +DA      2  +  +5-HT    3  =          4  ‖          5  X∪+  +DA+5-HT   6  =∪‖
 . . . . .      . . . . .      # # # # #     . . . . .     . . . . .          # # # # #
 . # . # .      . . # . .      . . . . .     # . . . #     . # # # .          # . . . #
 . . . . .      . # # # .      . . . . .     # . . . #     . # # # .          # . . . #
 . # . # .      . . # . .      . . . . .     # . . . #     . # # # .          # . . . #
 . . . . .      . . . . .      # # # # #     . . . . .     . . . . .          # # # # #
```

**Figure 1.** The six patterns of layout B, rows 1 to 5 from the top. EC/DG panel order as in the figure:
EC1 = X, EC2 = +, EC3 = =, EC4 = ‖, EC5 = X∪+, EC6 = =∪‖. The bottom-right panel is labelled EC5 in the figure;
it is read as EC6.

**Table 2.** Patterns of layout B.

| k | Name | Cells | n | Synapses n(n − 1) | Modulator | Role |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | X | 22, 24, 42, 44 | 4 | 12 | DA | Capture with DA |
| 2 | + | 23, 32, 33, 34, 43 | 5 | 20 | 5-HT | Capture with 5-HT |
| 3 | = | rows 1 and 5 (11–15, 51–55) | 10 | 90 | None | Tag, no PRP |
| 4 | ‖ | 21, 31, 41, 25, 35, 45 | 6 | 30 | None | Tag, no PRP |
| 5 | X∪+ | patterns 1 and 2 (the inner 3 × 3 block) | 9 | 72 | DA + 5-HT | Combination of the two modulated patterns |
| 6 | =∪‖ | patterns 3 and 4 (border) | 16 | 240 | None | Combination of the two unmodulated patterns |

Properties:

- The four basic patterns are pairwise disjoint (checked from the cell lists). No cell of = or ‖ ever receives DA or
  5-HT, so the controls have no access to PRPs.
- Each combination joins two basic patterns of the same kind: modulated with modulated (X∪+), unmodulated with
  unmodulated (=∪‖). No combination mixes a control with a modulated pattern, so the controls stay clean in every
  phase of encoding.
- No synapse receives both modulators. A cross synapse in X∪+ gets the modulator of its postsynaptic cell only.
- Every cell belongs to exactly one basic pattern (4 + 5 + 10 + 6 = 25), so there are no naive cells. In the first
  draft of the figure, + had no centre; the corrected figure (`figs/CA3_basic_circuit.png`, 2026-11-26) adds cell 33
  to +, which makes X∪+ the full inner 3 × 3 block.

### 2.1 Synapse classes and the expected outcome

**Table 3.** Synapse classes of P1 and the expected weight at 8 h (tagging by co-activity, capture gated by the
modulator, PRPs per postsynaptic cell).

| Class | Count | Tagged by | PRPs | Expected at 8 h |
| --- | --- | --- | --- | --- |
| X internal | 12 | Patterns 1, 5 | DA | Captured |
| + internal | 20 | Patterns 2, 5 | 5-HT | Captured |
| X → + cell, + → X cell | 20 + 20 | Pattern 5 only | 5-HT / DA (postsynaptic cell) | Captured; tagged once, so weaker than internal |
| = internal | 90 | Patterns 3, 6 | None | Back to w_0 with τ_tag |
| ‖ internal | 30 | Patterns 4, 6 | None | Back to w_0 |
| = ↔ ‖ cross | 60 + 60 | Pattern 6 only | None | Back to w_0; tagged once |
| Never co-presented: (X ∪ +) ↔ (= ∪ ‖) | 144 + 144 | No co-activity; a weaker tag from presynaptic firing alone (calcium rule, c_pre = 1) | Onto X/+ cells only | w_0: the no-co-activity baseline |

Total: 12 + 20 + 40 + 90 + 30 + 120 + 288 = 600 ✓.

The design is symmetric. The upper row of Table 3 (modulated) and the middle rows (controls) have the same structure:
two basic patterns, each presented on its own and again in a combination, plus cross synapses tagged only in the
combination. So every modulated class has a matched unmodulated class. The difference between them is the effect of
PRP availability (capture) on a tag of the same history.

### 2.2 Points to check

1. **Small modulated patterns.** X has 4 cells (12 synapses) and + has 5 (20). Check that:
   (a) the tag sum per cell (3 or 4 incoming internal synapses after pattern 1 or 2, 8 after pattern 5) passes θ_pro
   while the modulator is on;
   (b) replay of a 4- or 5-cell pattern is not too easy. The 80% criterion means all 4 cells of X (4 of 5 for +), and
   the chance of that among 25 cells is higher than for 9 cells.
   If (b) fails, read replay for X∪+ (9 cells) rather than for X and + alone.
2. **Size mismatch.** The modulated basic patterns have 4 and 5 cells. The controls have 10 and 6. Compare per
   synapse (mean 1 + T̃ + z per class), not per pattern total. A size-matched pair is X∪+ (9) against = (10).
3. **Encoding drive.** w_mf = 0.35 was calibrated on patterns of 5 to 19 cells. The 4-cell pattern X gets less feedback
   inhibition, so check that its rate does not go far above ~90 Hz.
4. **P6 replay feedback.** `w_ce_total` is divided over the pattern cells. For 4 cells, each synapse is larger. Check
   that EC_k still needs ~90% of its cells.

## 3. Implementation plan

Keep the old patterns selectable (`LAYOUT = "v1" | "B"`, default "B"), so that the v1 reports stay reproducible.

### Step 1. Patterns in one place (`circuit.py`)
- `PATTERN_SETS = {"v1": {...}, "B": {...}}`, with `PATTERNS`, `PATTERN_ORDER` and `MODULATOR_OF` taken from the
  chosen layout.
- Code keys (ASCII): `X`, `PLUS`, `H` (=), `VV` (‖), `X+PLUS`, `H+VV`. Display labels X, +, =, ‖, X∪+, =∪‖.
- Move the pattern colours (`PCOL`, now copied in 6 scripts) and the list of modulated patterns (`MODULATED` in
  `mod_stats.py`, literals in `prp_stats.py`, `trigger_stats.py`) into `circuit.py`.
- Add `synapse_class(circ)`: labels each P1 synapse with its class from Table 3.
- Assertions at import: the basic patterns are pairwise disjoint, the controls are disjoint from `MOD_CELLS`, and each
  combination is exactly the union of its two parts.

### Step 2. Modulation (`stc.py`)
- `MOD_CELLS = {"DA": X, "5-HT": PLUS}`, taken from `MODULATOR_OF` instead of hard-coded names.

### Step 3. Protocol and probes (`run_pilot.py`)
- Encoding order is pattern 1 to 6, as before (X, +, =, ‖, X∪+, =∪‖). DA starts with X and 5-HT with +, as with Λ and V.
- `PROBES`: one synapse per class, e.g. X 22 → 44, + 23 → 43, X → + 22 → 23, = 12 → 14, ‖ 31 → 35, = → ‖ 13 → 31,
  never co-presented 22 → 12.

### Step 4. Calibration (re-check)
- One seed: rates during encoding per pattern, EC_k reinstatement, resting rate, replay event size. Change w_mf or
  `w_ce_total` only if a check of 2.2 fails.

### Step 5. Readouts and statistics
- `stc_compare.py`, `prp_stats.py`, `da_stats.py`, `mod_stats.py`, `trigger_stats.py`, `x4_stats.py`,
  `rda_trigger_stats.py`, `replay_stats.py`: replace the Λ/V/I/O literals with the pattern sets from `circuit.py`.
- New main readout: W per synapse class (Table 3) over time and at 8 h.
- `replay_stats.py`: population 25 and burst threshold 23 stay. Pattern size comes from `PATTERNS`.

### Step 6. Figures
- `make_figures.py`, `make_protocol_fig.py`: grid plots coloured by the new patterns. New main figure: the weight of
  each synapse class against time, modulated classes next to their matched controls.
- The circuit figure is the corrected one of 2026-11-26. Its bottom-right panel is still labelled EC5/DG5 and is read
  as EC6/DG6.

### Step 7. Tests (`tests/test_patterns.py`)
- Disjointness, the basic patterns covering all 25 cells, unions, `MOD_CELLS` ∩ (= ∪ ‖) = ∅, class counts summing to 600.
- Offline STC: with no replay, after 8 h every control class has z = 0 and P = 0 on its postsynaptic cells.

### Step 8. Runs
- Pilot, one seed.
- 5 seeds with the bio-plausible protocol (10-min DA/5-HT with the synthesis trigger), as in the report of
  2026-11-24.
- PRP-block control (5 seeds), as in commit db0ea60.

### Step 9. Documents
- `ca3_basic_circuit_architecture.md`: replace the pattern section, Tables 3 and 4 and the overlap and pattern-6
  text with section 2 of this plan. Note that v1 was a misreading of the figure.
- `PROTOCOL.md`: pattern names in the phases and the readouts.
- Add a note at the top of the reports of 2026-11-17, 11-20 and 11-24 that they use the v1 patterns.
- New report comparing v1 and B.

## 4. Open choices

1. Whether the "pattern 6 neighbour-capture" question of v1 (can an unmodulated input capture from modulated
   neighbours?) is still wanted. Layout B does not test it, because no combination mixes the two kinds. A mixed
   combination such as X∪‖ could be added later as a seventh pattern.
