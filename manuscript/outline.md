> **SUPERSEDED 部分（2026-09-13）**：Table 2 已于 HPC C1 双端实例再生成（arm A 0.5–10×：0.923/0.912/0.958/0.971/0.970 (2026-09-17 signed fixed-origin)），形式交互检验仅 2× 显著为负。本文档中的旧数字与 "interaction" 主张以 manuscript.md 为准。详见 data/repro_check/ 与 data/repro_check/multiseed/INTERACTION_REPORT.md。

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

### 2. Accuracy on bacterial isolates — **Fig 2a–c**, **Table 2** ✅ + 🟡 `fig2_zheng_benchmark`
Zheng et al. 2020, E. coli K-12, 16 media (λ 0.40–1.72 h⁻¹), PRJNA615952.
- ✅ r = **0.981 at 1×**, 0.982 at 2×, 0.979 at 5× (n = 16). Accuracy is already
  saturated at 1×, matching Pilea's *published full-depth* r of 0.976.
- ✅ RMSE 0.157 at 1×, slope 0.78; 0.304 / 0.62 at 0.5×.
- ✅ Pilea at defaults returns **no estimate below 10×**.
- 🟡 Full-depth run to match Pilea's own Fig 2 directly — needs ~55 GB, deferred.
- ❌ Assembly-quality sensitivity (45,529 *Escherichia* assemblies × ANI × N50).

**[NEW] 2b. Negative control — **Fig 2d** (in `fig2_zheng_benchmark`).** Pilea *excluded* the
run-out samples. We use one: a replication run-out must give log₂PTR ≈ 0.
sk2bGrow 0.046–0.260; sorted regression returns **2.21** at 0.5×. ✅

**[NEW] 2d. Which gate suppresses Pilea, and is it right to? ✅**
Re-applying each default threshold to the gates-off output (exact, no extra
runs) shows `--min-cove 5` is the sole cause of **68/68** suppressed *E. coli*
estimates and **167/172** in the simulation. Because 150 bp gives 120 31-mers,
k-mer coverage 5 ≈ read coverage 6.25× (5 × 150/120). **The gate is partly right**: at 0.5×
Pilea's own estimator is degenerate (coverage exactly 1.000, log₂PTR exactly
0.000, all 17 samples). But at 1–5× it is informative (r 0.89/0.95/0.95) and
still suppressed. The honest claim is *"the gate is calibrated to where the
FracMinHash estimator breaks, and that point is 6.25× — sk2bGrow's is below 1×"*,
not *"Pilea is over-conservative"*.

**[NEW] 2c. Coverage titration.** Pilea ran full depth only. ✅

### 3. Multi-strain communities — **Fig 3** ✅ (laptop scale) `fig3_simulation`
Pilea's Fig-3 design at reduced scale: 16 reference genomes (incl. *E. coli*
K-12 / O157:H7 / *Shigella* as a shared-anchor stress test), 4/8/16 strains ×
1/2/4/8×, V-shaped profiles, log₂PTR ~ U[0,2].
- ✅ Recall 1.000 for sk2bGrow at every cell; Pilea-defaults 0.224.
- ✅ sk2bGrow RMSE **0.134** vs Pilea-gates-off 0.265 aggregate.
- ✅ Opposite bias, and sk2bGrow's is now negligible: **−0.013** vs Pilea **+0.168**.
- ❌ Scale to Pilea's full grid (32 strains, 32×, 400 samples) — needs a cluster.

### 4. What do the landmarks contribute? — **Fig 4** 🟡 rewritten post-F1
Three layers, every comparison density-stated (HPC, `data/f1_sketch/`):

- **(a) 2×2 at Pilea's operating point** (arm E at scale 250 = 2.4× lower
  density than the panel — the shipped comparison). The interaction stands as
  measured *at that operating point*:

| coverage | anchors + V-fit | FracMinHash + V-fit | anchors + rank | FracMinHash + rank (Pilea) |
|---|---|---|---|---|
| 0.5× | **0.913** | 0.724 | 0.164 | — (degenerate) |
| 1× | **0.981** | 0.940 | 0.683 | 0.889 |
| 2× | 0.982 | **0.984** | 0.756 | 0.947 |
| 5× | **0.979** | 0.977 | 0.914 | 0.954 |
| 10× | 0.968 | 0.942 | 0.913 | **0.971** |

- **(b) Density-matched (F1).** At matched landmarks/Mb (scale 104 vs k = 16:
  9,645 vs 9,422), the two landmark sources are **indistinguishable at ≥ 1×**
  (r within noise across 1–10×). At 0.5× the sketch is *ahead on survivors*
  (r 0.82–0.84 vs 0.50–0.57 in the same harness) but loses 6/17 cells at the
  gate. The 0.5× arm-E deficit in (a) was a density artefact, as
  HPC_TASKS_SKETCH_MODE.md §1.2 predicted.
- **(c) Mechanism.** Per-window observed landmarks and detected fractions are
  near-identical between modes at matched density (17.3–17.5 vs 14.1–21.3
  landmarks/window; detected fraction 0.175 vs 0.185). **Window population is
  a function of density, not of landmark determinism** — the sentence
  "deterministic anchors are what keeps windows populated at 0.5×" is retired
  everywhere. The panel's remaining low-depth robustness is **multi-strata
  fusion redundancy**: a single-stratum sketch that finds no downhill origin
  at 0.5× (6/17 cells) has no fusion to rescue it; sixteen strata do.
- Scope caveat for (a): the A-baseline in the F1 harness (r 0.568 at 0.5×)
  is harness-dependent (C8 review); cross-harness comparisons are invalid,
  within-harness A-vs-E comparisons are valid.
- Arm E is built by rewriting Pilea's sketch into sk2bGrow's count-table
  format and running the *unmodified* estimator, so the only difference from
  arm A is which loci are counted (Methods §4).

### 5. How many enzymes are needed, and over what GC range? — **Fig 5** ✅ (panel-size sweep a–d + GC sweep e, f in `fig5_panel_design`)
Rank all 16 by standalone accuracy (each `per_enzyme.tsv` is already an
independent V-shape fit, so the ranking is free — **Table 6**), then sweep
k ∈ {2,4,8,12,16}. Ranking is dominated by anchor yield, with exceptions: the top six by accuracy
all carry > 1,900 anchors (BslFI, 2,501 anchors, ranks ninth), and the bottom four by
accuracy (FalI/AloI/BplI/PsrI, 687/441/379/419 anchors) are near-zero or negative at
0.5× — but the sparsest enzyme, PpiI (336 anchors), holds r = 0.61 at 0.5×.

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

**[NEW] GC robustness (F2, HPC, `data/f2_gc_sweep/`).** The A6 recommendation
was derived at one GC value (*E. coli*, 50.8%). F2 measured density and
planted-gradient PTR accuracy on 18 complete genomes spanning GC 25.4–72.0%,
both landmark modes, at 0.5/1/2/5×:

- ✅ **A6 holds outside *E. coli*.** k8 ≈ k16 accuracy across the whole GC
  range at ≥ 1× (difference < 0.07 r; at 0.5× the gap grows to 0.12–0.21 r near the low-GC boundary); even at 72% GC, where k8 density
  falls to 0.64× k16, accuracy is unharmed. The 8-enzyme recommendation
  stands.
- ✅ **k16 density rises monotonically with GC** (~4.5k/Mb at 26% → 9.4k at
  50.8% → 12.4k/Mb at 72%); only the low-GC end is depressed. The
  "both-ends-fall" expectation from the three-genome table was wrong;
  Syn2b's non-monotonic peak was not observed.
- ✅ **FracMinHash is flat across GC** (matched within 4%) and ties the panel
  at matched density — the F1 verdict generalises out of *E. coli*.
- ⚠️ **Instrument boundary measured, not assumed:** at GC ≤ 30% *and* 0.5×
  depth, both modes fail (r −0.5–0.7; Buchnera, 0.65 Mb — genome-size
  confound flagged, not resolved). Stated range: GC ≳ 30%, depth ≳ 1×.
- Within the stated range (GC ≥ 30%, depth ≥ 1×) all cells r ≥ 0.85 at ≥ 1×; 2×/5× ≥ 0.95/0.98.

### 6. Fragmented references and MAGs — **Fig 6 / Table 8** ✅ `fig6_fragmentation`
Pilea's Fig 3 protocol on the Zheng data: complete chromosome vs 100 shuffled
lognormal contigs vs those contigs re-ordered by `sk2bgrow scaffold`, n = 16
media per cell, identical reads throughout (43,707 of 43,735 anchors survive the
cut, so only the coordinate changes).

| coverage | | complete | 100 contigs | scaffolded vs O157:H7 | Pilea on 100 contigs |
|---|---|---:|---:|---:|---:|
| 1× | r | 0.981 | 0.550 | 0.977 | 0.827 |
| | slope | 0.779 | 0.103 | 0.755 | 0.628 |
| 10× | r | 0.968 | 0.859 | 0.967 | 0.960 |
| | RMSE | 0.063 | 0.862 | 0.086 | 0.077 |
| | slope | 0.951 | 0.210 | 0.984 | 0.820 |
| | QC pass | 75% | **100%** | 88% | — |

- ✅ **Fragmentation removes the gradient rather than adding noise.** Every
  estimate lands at about a fifth of truth. r = 0.86 at 10× makes this look
  survivable; the slope of 0.21 says it is not. This is the clearest case in the
  paper for reporting slope beside correlation.
- ✅ **Pilea wins outright on unscaffolded contigs** (r 0.827 vs 0.550 at 1×);
  fragmentation costs its rank regression almost nothing, 0.889 to 0.827. This
  is a genuine architectural advantage of the sorted estimator and the paper
  should say so rather than only reporting the scaffolded column.
- ✅ **`scaffold` restores the complete-reference result exactly, across
  strains.** Contig order returns with Spearman 1.0000 against *E. coli*
  O157:H7, orientation 99/99 correct. The 712 kb raw placement error is a rigid
  rotation, invisible to a fit that searches for the origin instead of assuming
  it. The method is therefore not restricted to closed genomes.
- ⚠️ **The QC is anti-correlated with the failure**: 100% of fragmented
  estimates pass at 5–10× against 75% of correct ones. Cochran's Q asks whether
  the enzymes agree, and a destroyed coordinate makes all sixteen agree there is
  no gradient. This belongs in the Discussion as a stated limit of the QC, with
  a contig-count guard as the fix.
- ✅ **There is no safe contig count, and r cannot find one** (panel c, 10×,
  n = 16): across 1, 2, 5, 10, 20, 50 and 100 contigs the correlation stays
  0.86–0.97 while the slope falls 0.95 → 0.21 and the bias grows to −0.81. At 50
  contigs (N50 156 kb, a good draft) r reads 0.96 and every estimate is 44% of
  truth. Degradation is smooth and monotone, so the rule is "scaffold anything
  not closed", not "stay under N contigs". This is the paper's strongest single
  argument for reporting slope beside correlation.
- ✅ **The coordinate is not necessary in principle, only for our estimator.**
  Because log₂ coverage is uniform across the genome with width log₂(PTR), the
  spread alone determines PTR. Fitting that distribution with the per-window
  standard error in the model (panel d) takes RMSE on contigs from 0.87 to 0.14
  at 5× with no scaffolding, and returns 0.000 rather than Pilea's 1.153 on the
  stationary control at 1×. Below 5× it carries too little information (bias
  −0.98 at 1×) and it does not match Pilea at depth (RMSE 0.227 vs 0.077 at 10×).
  Reported as a prototype; the design it argues for is **two estimators chosen by
  whether a coordinate exists**, which also closes the QC blind spot.
- ❌ Genuinely incomplete MAGs (missing sequence, contamination) and more distant
  scaffolding references — needs the cluster.
- ❌ An overdispersion term for the spread estimator, which should fix both its
  low-depth shrinkage and its ~0.4 floor on the stationary control at 5–10×.

**[NEW post-F4] Landmark-agnosticism and the search-artifact caveat
(`data/f4_fragment_sketch/`).** Re-running the fragmentation protocol with
FracMinHash landmarks (F4):
- ✅ **The collapse is landmark-agnostic.** Enzyme and sketch arms give
  *identical* slope collapse (≈1.0 → 0.26–0.38 at 5×) and identical RMSE
  blow-up (0.02 → 0.84–0.87). The mechanism is coordinate loss, full stop.
- ⚠️ **The residual high r on fragmented references is partly a search
  artifact.** A permutation test (200 shuffled-coordinate refits) shows the
  ori grid search fits a fake V to a flat profile: shuffled profiles yield
  the same estimates, scaled with the planted value. The r column in the
  table above was itself inflated by this. Slope and RMSE against truth are
  the honest metrics on fragmented references.
- ✅ **`scaffold` is landmark-source-agnostic.** With FMH landmarks: self-arm
  ties the enzyme arm exactly (Spearman 1.0000, 5× RMSE 0.017–0.033);
  0.1%-divergence relative: 100/100 placed, full recovery; divergent
  relative: sketch slightly *better* (directional, n = 6). "Scaffolding
  needs the enzyme panel" is dropped as a route-A argument — it remains a
  route-B capability (enzyme-only data). Caveat: FMH placement used a 1:1
  Python port of `scaffold.rs` (validated by the perfect self-arm); the
  shipped CLI digests enzyme tags only.

### 7. Real metagenomes 🟡 **[NEW — the biological section, `data/sun_three_arm/`, `data/c5_review/`]**

**7a. Cross-method concordance in a fecal cohort (Sun, PRJNA689204) — **Fig 7 / Table 10**.**
The first real-community result: sk2bGrow (WGS) vs Pilea at its shipped
defaults on the same three fecal samples, common-denominator species × sample
pairs (n = 58/68/78 per sample).

- ✅ **The two methods agree.** Bias stable at +0.22–+0.25 log₂ (Pilea
  higher by ~17%), LoA ±0.4–0.9; r 0.34–0.51, CCC 0.25–0.37. The low CCC is
  deflation from a narrow dynamic range (observed span ~1.4), not broken
  agreement — every coefficient is printed beside its observed range.
- ✅ **Pilea's gate does real work (control arm).** With gates off, Pilea vs
  sk2bGrow r ≈ 0 — the same C5 lesson symmetrically confirmed on an
  independent dataset.
- Framed as concordance validation (n = 3 samples), not a biological
  discovery.

**7b. Real 2bRAD libraries, exact dedup, and per-anchor efficiency (Sun).**
Same cohort, in-silico WGS vs real 2bRAD library — same estimator, two
preparations (BcgI only; single-stratum caveat stated).

- ✅ **Exact deduplication before counting destroys the PTR signal.**
  raw vs B: S06 r 0.797 / CCC 0.728 / bias +0.60; deduped: r 0.30–0.60 /
  bias −0.6~−1.2 across samples. Cap-dedup sensitivity (cap ∈ {1,2,3,5,∞})
  shows no robust intermediate (cap 2 best in one sample, worst in another).
  **Process recommendation: do not exactly dedup 2bRAD libraries before
  counting.**
- ✅ **σ_eff ≈ 1.53 measured** (per-anchor capture-efficiency dispersion —
  no published measurement exists; Pilea cannot produce it, it needs paired
  2bRAD+WGS). ICC decomposition: ~0.75 correctable within batch; residual
  ~0.7 diluted by √99 window averaging to ~0.07–0.10 against gradients of
  1.0–1.7 — measured, correctable, negligible after correction. Wording
  limit: "within-batch reproducible"; Hou (independent batch) pending.

**7c. MAG-scale application and QC validity (RBC reactor, PRJNA974210;
522 MAGs × 9 samples) — **Fig 8a, b**.**

- ✅ **Recall under a common protocol:** sk2bGrow QC-pass 3.8–13.8%
  (20–72/522) vs Pilea 5.0–11.1% (26–58/522) — overlapping. The raw
  522/522 = 1.00 is a denominator artefact (this path has no output gate)
  and is never printed unqualified.
- ✅ **QC works on real data.** Pass rate anti-correlates with MAG
  fragmentation: ρ(qc_rate, n_contigs) = −0.41, −0.34 controlling coverage
  (p < 1e-4). This is R2's positive result on real data, and survives the
  coverage confound.
- ✅ **Cost at scale flips, and the flip is fixable.** 89.5–240.8× slower per
  sample than Pilea at MAG-panel scale (measured, `c5_cost_per_sample.tsv`); mm = 1 narrows it to ~10–25×
  (63× lookup speedup, 5.3% anchor loss, no genome lost); the M4 containment
  screen clears the absent-genome false-positive floor (4,048/4,700
  synthetic genomes hit in an unrestricted 8.3 Gb index; 0 with screen) —
  measured on the false-positive side only; at ~90% presence a direct A/B shows no
  speedup (0.99× wall) and the current default drops 65/522 real MAGs, so its value
  is a projection for database-scale presence rates (< 1%).

### 8. Computational efficiency and scale — **Fig 8c, d** ✅ + 🟡 post-F3/F5
Measured the same way for both tools (`/usr/bin/time -l`, 8 threads, k-mer
counting inside the timed region). Averaged over all five depths on the 85
*E. coli* cells: sk2bGrow 8.2 s / 195 MB at k = 16, 6.5 s at k = 8, 3.2 s at
k = 2; Pilea **11.5 s / 157 MB** with gates off. The average understates the
gap where it matters — at 1× and 2× Pilea takes 16.8 s and 21.0 s against our
6.8 s and 8.3 s, i.e. it is ~2.5× *slower* in the band this paper is about, and
cheaper only at 0.5× and 10×. **The earlier "Pilea is ~7× faster" claim does not survive this
measurement** — it came from the simulation, where Pilea's shipped gates skip
fitting entirely below 5× and it reports nothing. Pilea also has a striking
non-monotonic cost profile: 5.6 s at 0.5×, **21.0 s at 2×**, 5.6 s at 10×,
because its ZTP-mixture EM (`--max-iter` defaults to infinity) is slowest where
the mixture is least identifiable — exactly the depth band this paper is about.

**[NEW post-F3/F5] Scale outlook (C6).**
- **Mismatch tolerance is a budget decision, not an accuracy one (F3).**
  Enzyme mode: misassignment ≤1.4e-4 pooled at mm ≤ 2 (exact-first
  suppression + motif gate make tolerance nearly free) — default mm = 2
  stands; mm = 1 is safe if memory forces it (sensitivity loss 0.44%→0.42%
  unmatched). Sketch mode: **must lock mm = 0** — at scale ~100 the mm ≥ 1
  rescue channel credits 1–2% of off-sketch k-mers to a neighbour at the
  wrong coordinate, independent of sequencing error (the near-neighbour
  census cannot see this channel).
- **Index cost is a wash between modes at matched density (F5):** ~39
  B/anchor either way, equal build speed. A1's 752 GB was peak RSS; disk
  projection for GTDB species reps: enzyme k16 ~232 GB, k8 ~184 GB, FMH
  half-density (scale 200) ~124 GB.
- **Half density is the practical C6 answer:** scale 200 carries 0.53×
  landmarks at ≤0.007 r cost ≥1×; the only cost is at 0.5×, growing with
  GC (−0.03 r at 72%). Engineering gap stated: a merged `--mode
  fracminhash` index does not exist yet (the per-genome sketch pipeline is
  slower end-to-end despite equal per-landmark storage).

## Discussion

- Where the gain comes from — and where it does not (see below).
- Limits: single-strain-per-species; relic DNA; multi-fork; plasmids.
- Bsp24I ⊂ CjePI: measured on three genomes it is **100% containment**
  (1636/1636, 891/891, 2910/2910 — `enzymes.md:186`), with an additional
  partial Bsp24I p0 ⊂ CjeI p1 relation (48.4/47.4/50.9%). So "~15 independent
  strata, not 16" is an *optimistic* reading; and the fusion does **not** yet
  act on this (`enzyme::CONTAINMENTS` is declared but unused).
- The panel should probably be ~8 enzymes, not 16: the four sparsest contribute
  no accuracy and double the negative-control bias.
- **The panel is what works at low depth, not any single fit.** At 0.5× the mean
  per-enzyme V-fit has r² = −0.009 — individually worthless — yet fusing sixteen
  of them correlates 0.913 with measured growth rate. This is the clearest
  argument for a multi-enzyme panel over a single deterministic sketch, and it is
  measured rather than asserted.
- **Bsp24I ⊂ CjePI is a real dependency with no practical effect.** Re-fusing
  every sample with it dropped moves the estimate by 0.007–0.011 on average and
  0.050 at worst, and does not improve the Q rejection rate. Report the
  dependency; do not claim the correction matters.
- **An annotated origin is worth ~10% of RMSE below 2× and nothing above it**
  (0.304 → 0.271 at 0.5×, 0.157 → 0.143 at 1×, identical at 5–10×). It does not
  explain the low-coverage slope compression: with the origin exactly right the
  0.5× slope is still 0.74.
- **Cochran's Q earned its keep, and has a blind spot.** It was rejecting in
  56–69 % of samples; that was not biology, it was our double-count. After the
  fix, 6–24 %. A design whose QC can detect its own implementation defects is
  worth the complexity — and this is the honest way to present that, not as a
  clean-room result. But it tests *agreement between strata*, so it cannot see a
  failure that is identical across them: on a fragmented reference every enzyme
  agrees there is no gradient and 100 % of the wrong answers pass. A
  contig-count guard, not a better Q, is the fix.

### Landmarks: what the enzyme panel is, and what it is not **[NEW, post F1/C5 review]**

The thesis stays on the **interaction** (§4): sorted rank regression requires
uniform positional sampling, the coordinate-aware fit does not — and that
asymmetry, not landmark quality, is what makes a motif-constrained but
wet-lab-realizable landmark set usable. F1 cannot overturn this: at matched
density the two landmark sources are indistinguishable at ≥ 1×.

**Failure is loud, not silent — and loud is not a compliment.**
sk2bGrow has no coverage gate: on C5 it emitted estimates for 522/522 MAGs
(recall 1.00) including references the sample may not even contain; Pilea's
gate refuses silently. High yield therefore *includes wrong answers shipped*.
This belongs in Discussion as a first-class caveat (and drives the M4
containment pre-screen), not buried in a footnote.

**Two distinct SNP-sensitivity geometries.** One SNP in a 150 bp read destroys
31 of its 120 31-mers (~26 %) — FMH containment decays smoothly with
divergence. A SNP in a 6 bp recognition motif removes the whole site but
leaves neighbours untouched — blocky loss, buffered by the other enzymes.
The naive reading ("anchors tolerate divergence better") is **not** what the
measurement says: at 0.1 % substitutions, four-enzyme retention is 89.5 %
against FMH 94.6 %. Report both mechanism and number.

**Strata heterogeneity is the real QC asset.** A hash-prefix split of one FMH
sketch could also run a Cochran-style test — but hash strata share every
systematic bias of the sketch. The enzyme strata differ in motif, GC and
(methylatable) base composition, so cross-enzyme agreement tests bias
*heterogeneity*, not just sampling noise. This is also why the measured
Bsp24I ⊂ CjePI containment shrinks the effective independence.

**Coordinates make failure modes enumerable.** Wrong reference → the estimate
is emitted and wrong (C5). Destroyed coordinate → degradation is smooth,
monotone, detectable by scaffold length, and repairable (Fig 6). FMH's
analogous failure (dilution by an unreferenced strain) is silent and has no
coordinate to inspect. Do not overclaim: this is a difference in the *shape*
of failure, not in accuracy.

**Why landmarks must be genome-wide.** Single-copy markers (rpoB etc.) cover
< 1 % of a genome and cannot support a PTR fit at all; PTR forces a
genome-wide spread of loci, which in turn forces contact with within-species
(accessory-genome) variation. No landmark scheme resolves this tension; it is
inherent to the quantity being estimated.

**Dropped claims (checked and refuted):** "ori-adjacent gradient is steepest"
(contradicts the piecewise-linear model in fit.py — the V has constant slope
per arm); "the anchor set is constant across samples, unlike FMH" (an FMH
reference sketch is deterministic under the same hash-prefix rule; read
sampling is random for both); "enzyme sites are enumerable/filterable, FMH is
not" (both rules are deterministic and enumerable).

**Proposed Table (landmark-source comparison, one row per property, every row
carries its evidence; unmeasured rows marked 未测):**

| property | FracMinHash (Pilea) | enzyme anchors (ours) | evidence |
|---|---|---|---|
| landmark rule | hash-prefix subsample of all 31-mers | Type IIB recognition motifs, enumerated per reference | Table 1, manifest.json |
| coordinate enters estimator? | no — sorted ranks | yes — windowed V-fit | Fig 4 (2×2) |
| reference sketch determinism | deterministic | deterministic | — (both; claim dropped) |
| internal-consistency strata | hash-prefix split possible, shares all biases | ~15 motif strata, heterogeneous biases (motif/GC/methylation) | Cochran's Q; Bsp24I⊂CjePI containment |
| divergence geometry | smooth (SNP kills ~26 % of a read's 31-mers) | blocky (SNP kills one site, neighbours intact) | retention at 0.1 % subst: 94.6 % vs 89.5 % (4 enzymes) |
| wrong/absent reference | silent (gate refusal) | loud (estimate emitted, wrong) | C5 recall 1.00; gate ablation §2d |
| fragmented reference | rank regression barely affected (0.889→0.827) | gradient destroyed identically (landmark-agnostic), residual r partly search artifact, QC-blind, scaffold-repairable with either landmark set | Fig 6; F4 |
| mismatch tolerance | mm ≥ 1 opens a rescue channel crediting 1–2% of off-sketch k-mers to wrong coordinates (lock mm = 0) | mm = 2 nearly free (misassignment ≤1.4e-4; mm = 1 safe if budget forces) | F3 |
| index cost at matched density | ~39 B/landmark (per-genome pipeline; merged index not yet implemented) | ~39 B/anchor, merged index shipped | F5 |
| density control | scale parameter (content-agnostic; half density ≈ 124 GB at GTDB scale, ≤0.007 r cost ≥1×) | panel composition (biology-aware, per-clade tunable) | F5; Table 6; per-clade survey 未测 |
| wet-lab realization | computational only | physically produced by the 2bRAD protocol | route B |
| accuracy at matched density, ≥ 1× | identical | identical | F1 (HPC, this study) |
| 0.5× yield | 10–11/17 cells (single stratum; failures = no downhill origin) | 15–16/17 cells — rescued by multi-strata fusion redundancy, incl. wrong outputs (no gate) | F1 mechanism; C5 caveat |
| methylation-bias effect on anchor yield | n/a | 未测 | — |

## Methods ✅ → `manuscript/methods.md`

Written in full, mirroring Pilea's headings: algorithm (§1), enzyme panel (§2),
datasets (§3), comparison configuration and the gate ablation (§4), metric
definitions (§5), environment (§7), availability (§8).

---

## What we can and cannot claim

**Can:**
- Usable estimates at 1–2× where Pilea's defaults report nothing.
- Recall 1.00 in multi-strain communities down to 1×.
- A stationary-phase control that sorted regression fails badly.
- Enzyme definitions verified against the original Perl on three genomes.
- Cross-method PTR concordance with Pilea on real fecal metagenomes
  (bias +0.22–0.25 log2, n = 58–78 pairs per sample).
- A QC that demonstrably filters fragmented MAGs on real data
  (ρ = −0.41 with coverage controlled).
- A measured instrument boundary: GC ≳ 30% and depth ≳ 1×, both landmark modes.
- A process-level result others can cite: exact dedup before counting destroys
  the PTR signal on real 2bRAD libraries.

**Cannot — and must not imply:**
- *That deterministic anchors are the advance, or that the estimator is.*
  Neither factor carries the result: the coordinate fit is worth +0.30 r on
  anchors and +0.05 on a FracMinHash sketch, and the sketch effect changes sign
  with the estimator. The claim the data supports is about the **combination**,
  which is the only cell that works at 0.5×.
- *That sk2bGrow is unbiased.* It compresses the range at low coverage (slope
  0.62 at 0.5×, 0.78 at 1×) and mildly over-estimates at depth (bias +0.03 log₂ at 10×); Pilea over-estimates throughout. **And we cannot yet say why.** Three
  candidate mechanisms were tested and excluded — outlier trimming (no windows
  are trimmed), inverse-variance fusion (it improves the slope, 0.681 vs 0.497
  unweighted), and origin error (worth only 0.06 of the 0.38 shortfall). Do not
  repeat the errors-in-variables explanation: the predictor in that regression is
  effectively noise-free, so it does not apply.
- *That sk2bGrow is uniformly faster.* Per sample on this dataset it is on
  average (8.2 s vs 11.5 s at k = 16, 3.2 s at k = 2) and by 2.5× at 1–2×, but
  Pilea is cheaper at 0.5× and 10×, and **at its shipped defaults** it is far
  cheaper below its gate because it does no fitting at all there — and returns
  nothing.
- *That sk2bGrow handles draft assemblies out of the box.* On unscaffolded
  contigs Pilea is better, by a wide margin at 1×. The claim is that scaffolding
  fixes it, and scaffolding is a step the user has to take.
- *Anything beyond metagenome validation.* The real-data results are
  concordance and QC-validity on 3 fecal samples and one RBC reactor — not
  biological findings; the real 2bRAD arm additionally carries the single-
  enzyme (BcgI) and within-batch (σ_eff) qualifiers.
- *That the QC catches what goes wrong.* It catches enzyme disagreement, which
  is why it caught our own double-count; it is blind to a destroyed coordinate,
  which it passes at 100%. Say both.
- *That high yield is a clean win.* sk2bGrow reports an estimate for every
  reference it is given (C5: 522/522, recall 1.00) because this path has no
  gate; Pilea refuses silently. Yield therefore includes wrong answers shipped
  — pair every yield number with the C5-style caveat or the QC-pass rate.
