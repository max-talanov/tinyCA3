# 5-HT1A → Kir3 with a weaker encoding drive onto the + cells (report of 2026-11-27)

Follow-up to [the 5-HT1A → Kir3 report](SEROTONIN_GIRK_2026-11-26.md). With the calibrated drive, the GIRK current
changed nothing: the + cells fire at 90 Hz either way. That report's first next step was to weaken the encoding drive
onto the + cells, so that a 4.3 mV hyperpolarization can change their firing. This report does that.

## Summary

| # | Question | Answer |
| --- | --- | --- |
| 1 | What changes? | The mossy-fibre weight onto the + cells is scaled by a factor s (`--mf-5ht s`, `CircuitParams.mf_scale_5ht`). It applies wherever the + cells are presented, alone and in X∪+. Every other cell keeps the calibrated drive |
| 2 | Where can the GIRK current act? | Only at the weakest drive tried. From s = 1 down to 0.1, the + cells fire the same number of spikes with the current off and on, and the tags agree to 0.012. At s = 0.05 the + cells fire 31 Hz instead of 34, and the X→+ tag falls from 0.457 to 0.428 |
| 3 | What does the weaker drive do by itself? | It turns the tag of the + synapses from LTP into LTD. Between s = 0.2 and 0.3 the + tag drops from 0.57 to about 0; at s = 0.15–0.1 it is about −0.25, past the −0.2 capture threshold |
| 4 | Effect at 8 h, s = 0.05 (30 runs) | **Yes, small and consistent.** GIRK lowers X→+ by 0.044 (1.437 → 1.393) and X by 0.013 (1.508 → 1.496), in all 5 seeds. The effect is as large without replay and gone under the PRP block. + and +→X stay at w0 with or without the current (tags below the capture threshold) |
| 5 | Where does it come from? | Encoding. The tags differ by the same amounts at the end of encoding, and 390 of 400 replay events recruit the same cells. 5-HT then gates PRP synthesis in the + cells, so the smaller X→+ tags are captured as a smaller z |

## 1. Protocol

Layout B, `--mod-min 10 --gate 30` as in the earlier reports. The only change is the mossy-fibre (P5) weight onto the
5 + cells (23, 32, 33, 34, 43): w_mf × s. The EC → CA3 input (P6) and the drive onto every other pattern are
unchanged.

- **Encoding sweep** (`weak_drive_sweep.py`, seed 1, encoding block only): s ∈ {1, 0.5, 0.3, 0.2, 0.15, 0.1, 0.05},
  with GIRK off and at 6 × 10⁻⁵ S/cm², 14 runs. Results in `results/layout_B/weak_drive_sweep.json`.
- **Full runs** at s = 0.05: 5 seeds × {control, PRP block, no replay} × {GIRK off, on}, 30 runs of 8 h.
  `run_pilot.py --seed S --light --mod-min 10 --gate 30 --mf-5ht 0.05 [--ht-girk 6e-5] [--prp-block | --no-replay]`.
  Files `results/layout_B/seeds/pilot_s{1..5}_mod10_gate30{,_prpblock,_norep}_mf0.05{,_girk}.pkl`.
- **Readout:** `class_stats.py --mf-5ht 0.05 [--girk]` (`class_stats_mf0.05{,_girk}.json`,
  `fig_B1_classes_mf0.05{,_girk}.png`) and `weak_drive_fig.py` (Figure W1, `weak_drive_8h.json`).

s = 0.05 was chosen because it is the only scale in the sweep where the current changed the firing of the + cells.

## 2. Encoding sweep

Seed 1. Rates are spikes in the three 100-ms bursts of a presentation, per + cell, divided by 0.3 s.

| s | + rate, + presented (off → on, Hz) | + rate, X∪+ presented (off → on) | X rate, X∪+ presented (off → on) | Tag + (off → on) | Tag X→+ (off → on) | Tag +→X (off → on) |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 90 → 90 | 90 → 90 | 90.8 → 90.8 | 0.638 → 0.638 | 0.553 → 0.553 | 0.557 → 0.557 |
| 0.5 | 90 → 90 | 90 → 90 | 90.8 → 90.8 | 0.637 → 0.637 | 0.554 → 0.554 | 0.554 → 0.555 |
| 0.3 | 90 → 90 | 80.7 → 80.0 | 91.7 → 90.8 | 0.573 → 0.573 | 0.535 → 0.536 | 0.518 → 0.519 |
| 0.2 | 70.7 → 70.0 | 61.3 → 60.7 | 90 → 90 | 0.029 → 0.017 | 0.507 → 0.507 | −0.098 → −0.104 |
| 0.15 | 60 → 60 | 40.7 → 40.7 | 90 → 90 | −0.257 → −0.259 | 0.490 → 0.490 | −0.284 → −0.284 |
| 0.1 | 50 → 50 | 34.0 → 34.0 | 90 → 90 | −0.226 → −0.229 | 0.489 → 0.490 | −0.264 → −0.265 |
| 0.05 | 34.0 → 31.3 | 17.3 → 20.7 | 90 → 86.7 | −0.109 → −0.096 | 0.457 → 0.428 | −0.151 → −0.147 |

![Three panels: (a) firing of the + cells against the mossy-fibre scale, GIRK off and on, overlapping except at 0.05; (b) tags of the + and X→+ synapses at the end of encoding against the scale; (c) W at 8 h per seed for X, X→+, + and +→X at scale 0.05, GIRK off and on.](../results/layout_B/fig_weak_drive.png)

**Figure W1.** `weak_drive_fig.py`. (a) and (b): the encoding sweep, seed 1. (c): W at 8 h at s = 0.05, control,
5 seeds; grey lines pair the runs of one seed, numbers are the mean paired change.

- **The current is still nearly silent.** Each mossy-fibre pulse either fires a + cell or not. A few millivolts
  changes that only for the occasional pulse, and only once the drive is close to threshold. Down to s = 0.1 the
  firing is identical; the largest tag change is 0.012 (the + tag at s = 0.2, where it crosses zero).
- **At s = 0.05 the effect reaches the tags.** The + cells lose spikes when + is presented (34 → 31 Hz), and the
  X→+ tag falls by 0.029. The X cells, which have no GIRK, fire less in X∪+ (90 → 86.7 Hz), and the X tag falls by
  0.023. This goes through the network. The + cells fire slightly more in X∪+ (17.3 → 20.7 Hz), so it is not simply
  fewer + spikes. A likely contributor for X→+ is the dendritic GIRK current, which hyperpolarizes the + dendrites and
  so admits less NMDA calcium per spike. That was not tested separately.
- **A weaker drive alone gives LTD tags.** Between s = 0.3 and 0.2 the + tag drops from 0.57 to near 0, and the +→X
  tag turns negative. At s = 0.15–0.1 both are about −0.25, past the capture threshold of 0.2 in absolute value. With
  10 min of 5-HT on the + cells, that LTD tag would be captured. This is the CA2 situation (5-HT makes weak LTD last;
  Benoy et al. 2022), produced here by the drive and not by 5-HT. At s = 0.05 the tags are −0.11 and −0.15, below
  the threshold.

## 3. Results at 8 h, s = 0.05

Mean W = 1 + T̃ + z per class at 8 h, 5 seeds (`weak_drive_8h.json`); the calibrated drive without GIRK
(`class_stats.json`) for reference.

| Class | Calibrated drive | s = 0.05, GIRK off | s = 0.05, GIRK on | Paired change (5 seeds) |
| --- | --- | --- | --- | --- |
| X (DA) | 1.542 | 1.508 | 1.496 | −0.013 (−0.008 to −0.021) |
| X→+ (5-HT) | 1.530 | 1.437 | 1.393 | −0.044 (−0.024 to −0.091) |
| + (5-HT) | 1.546 | 0.998 | 0.999 | +0.0004 |
| +→X (DA) | 1.529 | 0.997 | 0.997 | +0.0002 |
| = (control) | 1.006 | 1.006 | 1.006 | 0.0000 |
| never co-presented | 1.001 | 1.000 | 1.000 | −0.0001 |

| Condition | X: off → on | X→+: off → on |
| --- | --- | --- |
| Control (10 min + trigger) | 1.508 → 1.496 | 1.437 → 1.393 |
| No replay | 1.537 → 1.526 | 1.453 → 1.412 |
| PRP block | 1.0046 → 1.0043 | 1.0044 → 1.0038 |

- **The effect is set at encoding.** Across the 5 seeds, the GIRK current changes the tags at the end of encoding by
  −0.023 (X) and −0.041 (X→+), about the same as the 8-h changes. Without replay the 8-h changes are almost the same
  (−0.011 and −0.041).
- **It is kept by capture.** Under the PRP block every class decays to w0 and the difference disappears (below
  0.001). The captured z at 8 h is 0.506 → 0.494 for X and 0.434 → 0.391 for X→+.
- **Replay barely changes.** 390 of 400 events recruit the same cells with the current off and on. Replay counts per
  pattern agree within 1. The weaker drive itself reduces replay of the patterns containing +: + 68 → 38 and X∪+
  16 → 6 events over the 5 seeds, against the calibrated drive.
- **+ and +→X are not captured.** Their tags (−0.11 and −0.15) are below the capture threshold, so they decay to w0
  with or without 5-HT and GIRK.

## 4. What this means, and next steps

5-HT1A → Kir3 can change what is stored, but only once the + cells are driven close to threshold. Then it lowers
the LTP tags onto the + cells by about 9% and the captured weights by about 3% of W. It also lowers the X weights
through the network. Its effect is subtractive, as in the slice (fewer spikes, lower membrane potential). With the
calibrated drive it is invisible.

The larger finding is the sign change of the + tags with a weaker drive. Near s = 0.1–0.15 the + pattern sets LTD
tags strong enough to be captured, and 5-HT gates their capture. That gives the LTD protocol for + from the earlier
reports without a new plasticity rule.

- [ ] Rerun at s = 0.15 (LTD tags of −0.26, captured under 5-HT) with GIRK off and on and with the PRP block: does
  5-HT make the + LTD last, as in CA2 (Benoy et al. 2022)?
- [ ] Separate the dendritic and somatic GIRK (insert only in the soma) to test whether NMDA calcium explains the X→+
  change.
- [ ] Slow GIRK onset (about 5 min, Johnston et al. 2014 Fig. 1A): at s = 0.05 the current is on during all of
  encoding; with a 5-min onset the current should miss encoding, and the effect should disappear.
- [ ] Spiking background between events, where the slice effect (5 → 1.9 Hz) was measured.
