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
- ✅ r = 0.954 at 1×, 0.975 at 2×, 0.979 at 5× (n = 16).
- ✅ Pilea at defaults returns **no estimate below 10×**.
- 🟡 Full-depth run to match Pilea's own Fig 2 directly — needs ~55 GB, deferred.
- ❌ Assembly-quality sensitivity (45,529 *Escherichia* assemblies × ANI × N50).

**[NEW] 2b. Negative control — `fig4_negative_control`.** Pilea *excluded* the
run-out samples. We use one: a replication run-out must give log₂PTR ≈ 0.
sk2bGrow 0.045–0.113; sorted regression returns **2.17** at 0.5×. ✅

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
- ✅ Recall 1.00 for sk2bGrow at every cell; Pilea-defaults 0.185.
- ✅ sk2bGrow RMSE 0.181 vs Pilea-relaxed 0.293 aggregate.
- ✅ Opposite bias: sk2bGrow −0.082, Pilea +0.200.
- ❌ Scale to Pilea's full grid (32 strains, 32×, 400 samples) — needs a cluster.

### 4. Attribution — sketch or estimator? — **Fig 5** ✅ `fig5_attribution`
Holding the estimator fixed, anchors are *behind* FracMinHash at 1×. The gain is
the coordinate-aware fit, not the deterministic sketch.

### 5. How many enzymes are needed? — **Fig 7 / Table 5** 🟡 *(running)*
Rank all 16 by standalone accuracy (each `per_enzyme.tsv` is already an
independent V-shape fit, so the ranking is free), then sweep k ∈ {2,4,8,12,16}.
Ranking is dominated by anchor yield: the top 6 are exactly the 6 with > 2,500
anchors; the bottom 4 have < 450 and are anti-correlated at 0.5×.
Early partial result (k = 2 vs the full panel, 0.5×): r 0.917 vs 0.907,
**RMSE 0.279 vs 0.506**, slope 0.733 vs 0.519, at 2.3 s vs 14.8 s.
If this holds, the sparse enzymes are adding variance, not strata — a result
that argues against the paper's own 16-enzyme design and must be reported as
such.

### 6. Computational efficiency — **Table 2** ✅
Index, profile wall-clock, peak RSS, database size. **Pilea is faster** at equal
panel size; the panel sweep is what closes the gap.

### 7. Application ❌
Pilea used a rotating biological contactor. We have no application dataset yet.

## Discussion

- Where the gain comes from — and where it does not (see below).
- Limits: single-strain-per-species; relic DNA; multi-fork; plasmids.
- Bsp24I ⊂ CjePI: the panel is ~15 independent strata, not 16 — and the fusion
  does **not** yet act on this (`enzyme::CONTAINMENTS` is declared but unused).
- If §5 holds, the panel should be *smaller*, not larger.

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
- *That sk2bGrow is unbiased.* It compresses the PTR range at low coverage
  (slope 0.52 at 0.5×) and underestimates; Pilea overestimates.
- *That sk2bGrow is faster.* It is not. Pilea at defaults is ~7× faster.
- *Anything about metagenomes yet.* Every real-data result is one strain against
  a complete reference.
