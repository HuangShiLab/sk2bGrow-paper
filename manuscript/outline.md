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

### 1. sk2bGrow overview — **Fig 1** ❌
Workflow: digest → anchor DB (uniqueness masking) → count (≤2 mismatch) → EM
reassignment → per-enzyme adaptive windows → ZTP/NB window rates → per-enzyme
GC correction → shared-origin V-shape fit → inverse-variance fusion + χ² QC.
*Needs drawing.*

### 2. Accuracy on bacterial isolates — **Fig 2** ✅ + 🟡
Zheng et al. 2020, E. coli K-12, 16 media (λ 0.40–1.72 h⁻¹), PRJNA615952.
- ✅ r = 0.954 at 1×, 0.975 at 2×, 0.979 at 5× (n = 16).
- ✅ Pilea at defaults returns **no estimate below 10×**.
- 🟡 Full-depth run to match Pilea's own Fig 2 directly — needs ~55 GB, deferred.
- ❌ Assembly-quality sensitivity (45,529 *Escherichia* assemblies × ANI × N50).

**[NEW] 2b. Negative control.** Pilea *excluded* the three run-out samples. We
use them: a stationary culture must give log₂PTR ≈ 0. sk2bGrow 0.045–0.113;
sorted regression returns **2.17** at 0.5×. ✅

**[NEW] 2c. Coverage titration.** Pilea ran full depth only. ✅

### 3. Multi-strain communities — **Fig 3** ✅ (laptop scale)
Pilea's Fig-3 design at reduced scale: 16 reference genomes (incl. *E. coli*
K-12 / O157:H7 / *Shigella* as a shared-anchor stress test), 4/8/16 strains ×
1/2/4/8×, V-shaped profiles, log₂PTR ~ U[0,2].
- ✅ Recall 1.00 for sk2bGrow at every cell; Pilea-defaults 0.185.
- ✅ sk2bGrow RMSE 0.181 vs Pilea-relaxed 0.293 aggregate.
- ✅ Opposite bias: sk2bGrow −0.082, Pilea +0.200.
- ❌ Scale to Pilea's full grid (32 strains, 32×, 400 samples) — needs a cluster.

### 4. Computational efficiency — **Table 2** ✅
Index, profile wall-clock, peak RSS, database size. **Pilea is faster.**

### 5. Application ❌
Pilea used a rotating biological contactor. We have no application dataset yet.

## Discussion

- Where the gain comes from — and where it does not (see below).
- Limits: single-strain-per-species; relic DNA; multi-fork; plasmids.
- Bsp24I ⊂ CjePI: the panel is ~15 independent strata, not 16.

## Methods

Mirror Pilea's Methods headings. Enzyme panel transcribed from
`2bRADExtraction.pl`, verified against three genomes (46/48 cells within 3%).

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
