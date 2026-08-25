# sk2bGrow — manuscript outline

Structure mirrors Pilea (Microbiome 2026, doi:10.1186/s40168-026-02374-0) so a
reviewer can compare section-for-section. Deviations are marked **[NEW]**.

Status legend: ✅ data in hand · 🟡 partial · ❌ not started

---

## Title (working)

*Coordinate-aware peak-to-trough estimation from deterministic 2bRAD anchors*

The title deliberately foregrounds **coordinate-aware estimation**, not the
deterministic sketch. See "What we can and cannot claim" below — the data does
not support a sketch-first framing.

## Abstract

One paragraph. Must state the low-coverage result *and* the bias caveat; an
abstract built on correlation alone would misrepresent Fig 2.

## Introduction

- PTR as a culture-independent growth proxy (Korem 2015; CoPTR; Pilea).
- The alignment-free turn: sketching makes GTDB-scale profiling tractable.
- Gap: a random sketch discards **position**. Sorted-rank regression is the
  workaround, and it is fragile — it rides on extreme order statistics.
- 2bRAD/Type IIB anchors are a sketch whose loci are motif-defined, hence known
  a priori and identical across samples. That makes coordinates usable.

## Results

### 1. sk2bGrow overview — **Fig 1** ✅ `fig1_overview`
(a) workflow schematic; (b) window rate against genome coordinate with the
fitted V; (c) the same windows under rank regression, where a single dropout
window moves the estimate from 1.86 to 2.51; (d) forest plot of the 16
per-enzyme estimates with the fused value and *I*². Panels b–d are real output
from one 2× sample, so the figure shows the method working, not a cartoon.
Full prose in `manuscript/methods.md` §1.

### 2. Accuracy on bacterial isolates — **Fig 2** ✅ + 🟡 `fig2_accuracy_vs_coverage`, `fig3_magnitude`
Zheng et al. 2020, E. coli K-12, 16 media (λ 0.40–1.72 h⁻¹), PRJNA615952.
- ✅ r = **0.981 at 1×**, 0.982 at 2×, 0.979 at 5× (n = 16). Accuracy is already
  saturated at 1×, matching Pilea's *published full-depth* r of 0.976.
- ✅ RMSE 0.157 at 1×, slope 0.78; 0.304 / 0.62 at 0.5×.
- ✅ Pilea at defaults returns **no estimate below 10×**.
- 🟡 Full-depth run to match Pilea's own Fig 2 directly — needs ~55 GB, deferred.
- ❌ Assembly-quality sensitivity (45,529 *Escherichia* assemblies × ANI × N50).

**[NEW] 2b. Negative control — `fig4_negative_control`.** Pilea *excluded* the
run-out samples. We use one: a replication run-out must give log₂PTR ≈ 0.
sk2bGrow 0.046–0.260; sorted regression returns **2.21** at 0.5×. ✅

**[NEW] 2d. Which gate suppresses Pilea, and is it right to? ✅**
Re-applying each default threshold to the gates-off output (exact, no extra
runs) shows `--min-cove 5` is the sole cause of **68/68** suppressed *E. coli*
estimates and **167/172** in the simulation. Because 150 bp gives 120 31-mers,
k-mer coverage 5 ≈ read coverage 6.6×. **The gate is partly right**: at 0.5×
Pilea's own estimator is degenerate (coverage exactly 1.000, log₂PTR exactly
0.000, all 17 samples). But at 1–5× it is informative (r 0.89/0.95/0.95) and
still suppressed. The honest claim is *"the gate is calibrated to where the
FracMinHash estimator breaks, and that point is 6.6× — sk2bGrow's is below 1×"*,
not *"Pilea is over-conservative"*.

**[NEW] 2c. Coverage titration.** Pilea ran full depth only. ✅

### 3. Multi-strain communities — **Fig 6** ✅ (laptop scale) `fig6_simulation`
Pilea's Fig-3 design at reduced scale: 16 reference genomes (incl. *E. coli*
K-12 / O157:H7 / *Shigella* as a shared-anchor stress test), 4/8/16 strains ×
1/2/4/8×, V-shaped profiles, log₂PTR ~ U[0,2].
- ✅ Recall 1.000 for sk2bGrow at every cell; Pilea-defaults 0.224.
- ✅ sk2bGrow RMSE **0.134** vs Pilea-gates-off 0.265 aggregate.
- ✅ Opposite bias, and sk2bGrow's is now negligible: **−0.013** vs Pilea **+0.168**.
- ❌ Scale to Pilea's full grid (32 strains, 32×, 400 samples) — needs a cluster.

### 4. Attribution — sketch or estimator? — **Fig 5** ✅ `fig5_attribution`
The full 2 × 2 (both sketches × both estimators) on the same subsampled reads.
Pearson r against measured growth rate:

| coverage | anchors + V-fit | FracMinHash + V-fit | anchors + rank | FracMinHash + rank (Pilea) |
|---|---|---|---|---|
| 0.5× | **0.913** | 0.724 | 0.164 | — (degenerate) |
| 1× | **0.981** | 0.940 | 0.683 | 0.889 |
| 2× | 0.982 | **0.984** | 0.756 | 0.947 |
| 5× | **0.979** | 0.977 | 0.914 | 0.954 |
| 10× | 0.968 | 0.942 | 0.913 | **0.971** |

- ✅ **The two factors interact; neither carries the result alone.** At 1× the
  coordinate fit is worth +0.30 r on anchors but only +0.05 on a FracMinHash
  sketch (interaction +0.25), and the *sketch* effect changes sign with the
  estimator: anchors are +0.04 ahead under the V-fit and −0.21 behind under rank
  regression. The earlier three-arm reading — "the gain is the estimator, not the
  sketch" — was an artefact of the missing cell.
- ✅ **On magnitude the estimator is the dominant factor.** RMSE at 1×: 0.157
  (A) and 0.213 (E) under the V-fit, against 1.027 (B) and 0.397 (C) under rank
  regression. Rank regression is biased upward on anchors by +0.97 log₂ units.
- ✅ **At 0.5× only the combination survives** (r = 0.913): the V-fit on a
  FracMinHash sketch drops to 0.724 and Pilea's own arm is fully degenerate.
  Deterministic anchors are what keeps windows populated at that depth; the
  coordinate fit is what turns them into an unbiased slope.
- Arm E is built by rewriting Pilea's sketch into sk2bGrow's count-table format
  and running the *unmodified* estimator, so the only difference from arm A is
  which loci are counted (Methods §4).

### 5. How many enzymes are needed? — **Fig 7 / Table 7** ✅
Rank all 16 by standalone accuracy (each `per_enzyme.tsv` is already an
independent V-shape fit, so the ranking is free — **Table 6**), then sweep
k ∈ {2,4,8,12,16}. Ranking is dominated by anchor yield: the top six are the six
with > 1,900 anchors; the bottom four have < 450 and are near-zero or negative at
0.5×.

**This analysis is what found the counting bug.** The first pass said k = 2 beat
the full panel by a wide margin at low depth. That was not a property of the
panel: adding Bsp24I (at k = 8) was corrupting CjePI through the double-count
described in Methods §1.2. The sweep was then repeated with the corrected
ranking, so Table 6 and Table 7 agree.

Final picture, averaged over ≤2×:

| enzymes | anchors | r | RMSE | run-out bias | s/sample |
|---|---|---|---|---|---|
| 2 | 17,055 | 0.951 | 0.188 | 0.073 | 3.2 |
| **4** | 24,753 | 0.959 | **0.165** | 0.074 | 4.4 |
| **8** | 37,232 | **0.969** | 0.171 | 0.123 | 6.5 |
| 12 | 41,662 | 0.969 | 0.183 | 0.132 | 7.8 |
| 16 | 43,735 | 0.959 | 0.196 | 0.145 | 8.2 |

**4–8 enzymes, not 16.** r peaks at 8; RMSE at 4; the four sparsest enzymes buy
nothing and double the run-out bias. Cost is linear in k. Note that a 2-enzyme
panel is 17,055 anchors against Pilea's 18,261 sketch k-mers — the like-for-like
sketch-size comparison, and it still gives r = 0.960 at 1× where Pilea's own
gates report nothing.

⚠️ The mechanism is **not established**. Two candidate explanations were tested
and both refuted: sparse enzymes are not more biased than dense ones
(mean bias at 0.5×: −0.319 vs −0.312), and forcing fixed-effect fusion weights
does not recover the loss (RMSE 0.581 vs 0.506 at 0.5×). Reported as an
empirical result, not a mechanism.

### 6. Computational efficiency — **Table 7** ✅
Measured the same way for both tools (`/usr/bin/time -l`, 8 threads, k-mer
counting inside the timed region). On the 85 *E. coli* cells: sk2bGrow 8.6 s /
192 MB at k = 16, 7.0 s at k = 8, 3.2 s at k = 2; Pilea **11.6 s / 157 MB** with
gates off. **The earlier "Pilea is ~7× faster" claim does not survive this
measurement** — it came from the simulation, where Pilea's shipped gates skip
fitting entirely below 5× and it reports nothing. Pilea also has a striking
non-monotonic cost profile: 5.6 s at 0.5×, **21.1 s at 2×**, 5.6 s at 10×,
because its ZTP-mixture EM (`--max-iter` defaults to infinity) is slowest where
the mixture is least identifiable — exactly the depth band this paper is about.

### 7. Application ❌
Pilea used a rotating biological contactor. We have no application dataset yet.

## Discussion

- Where the gain comes from — and where it does not (see below).
- Limits: single-strain-per-species; relic DNA; multi-fork; plasmids.
- Bsp24I ⊂ CjePI: the panel is ~15 independent strata, not 16 — and the fusion
  does **not** yet act on this (`enzyme::CONTAINMENTS` is declared but unused).
- The panel should probably be ~8 enzymes, not 16: the four sparsest contribute
  no accuracy and double the negative-control bias.
- **Cochran's Q earned its keep.** It was rejecting in 56–69 % of samples; that
  was not biology, it was our double-count. After the fix, 6–24 %. A design
  whose QC can detect its own implementation defects is worth the complexity —
  and this is the honest way to present that, not as a clean-room result.

## Methods ✅ → `manuscript/methods.md`

Written in full, mirroring Pilea's headings: algorithm (§1), enzyme panel (§2),
datasets (§3), comparison configuration and the gate ablation (§4), metric
definitions (§5), environment (§6), availability (§7).

---

## What we can and cannot claim

**Can:**
- Usable estimates at 1–2× where Pilea's defaults report nothing.
- Recall 1.00 in multi-strain communities down to 1×.
- A stationary-phase control that sorted regression fails badly.
- Enzyme definitions verified against the original Perl on three genomes.

**Cannot — and must not imply:**
- *That deterministic anchors are the advance.* Holding the estimator fixed,
  anchors are **behind** FracMinHash at 1× (r 0.61 vs 0.89). The contribution is
  the coordinate-aware fit the anchors make possible.
- *That sk2bGrow is unbiased.* It still compresses the range at low coverage
  (slope 0.62 at 0.5×, 0.78 at 1×) and underestimates; Pilea overestimates.
- *That sk2bGrow is uniformly faster.* Per sample on this dataset it is
  (8.6 s vs 11.6 s at k = 16, 3.2 s at k = 2), but Pilea **at its shipped
  defaults** is far cheaper below its gate because it does no fitting at all
  there — and returns nothing.
- *Anything about metagenomes yet.* Every real-data result is one strain against
  a complete reference.
