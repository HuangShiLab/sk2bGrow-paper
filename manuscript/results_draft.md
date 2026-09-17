> **SUPERSEDED 部分（2026-09-13）**：Table 2 已于 HPC C1 双端实例再生成（arm A 0.5–10×：0.941/0.919/0.959/0.971/0.970），形式交互检验仅 2× 显著为负。本文档中的旧数字与 "interaction" 主张以 manuscript.md 为准。详见 data/repro_check/ 与 data/repro_check/multiseed/INTERACTION_REPORT.md。

# Results (new sections, draft)

*Draft covering §4–§8. Sections §1–§3 (overview, isolate accuracy, multi-strain
communities) already exist as prose elsewhere. All numbers follow the latest
REVISION of the corresponding `data/*/REVIEW.md` and tsv files.*

## 4. What do the landmarks contribute?

The coordinate-aware estimator, not the deterministic landmark set, is what
carries the low-coverage result; at matched landmark density the two landmark
sources are indistinguishable at ≥ 1×, and the panel's residual advantage below
1× comes from multi-strata fusion redundancy, not from landmark determinism. We
tested this in three layers on the Zheng *E. coli* grid (HPC; `data/f1_sketch/`).

First, at Pilea's own operating point (FracMinHash arm E at scale 250, a 2.4×
lower landmark density than the panel — the shipped comparison), the
2 × 2 attribution (landmark source × estimator) stands as measured at that
operating point:

| coverage | anchors + V-fit | FracMinHash + V-fit | anchors + rank | FracMinHash + rank (Pilea) |
|---|---|---|---|---|
| 0.5× | **0.913** | 0.724 | 0.164 | — (degenerate) |
| 1× | **0.981** | 0.940 | 0.683 | 0.889 |
| 2× | 0.982 | **0.984** | 0.756 | 0.947 |
| 5× | **0.979** | 0.977 | 0.914 | 0.954 |
| 10× | 0.968 | 0.942 | 0.913 | **0.971** |

Arm E was constructed by rewriting Pilea's sketch into sk2bGrow's count-table
format and running the *unmodified* sk2bGrow estimator, so the only difference
from arm A is which loci are counted (Methods §4).

Second, matching density removes the sketch's apparent 0.5× deficit: at
9,645 (sketch, scale 104) versus 9,422 (panel, k = 16) landmarks/Mb — within
2% — the two landmark sources are indistinguishable at 1–10×: a paired media
bootstrap (10,000 resamples of the media) bounds the 95% CI of the accuracy
difference at ±0.16 r at every depth ≥ 1× (−0.00 [−0.16, +0.12] at 1×;
data/ci_bootstrap/f1_ci.tsv). The source × estimator interaction is formally
significant at 1×, 2× and 10× on the Fisher-z scale (not at 5×, where all four
cells are already close; data/ci_bootstrap/interaction.tsv). At 0.5× the single-stratum sketch is
actually *ahead of the panel on surviving cells* (r 0.82–0.84, n = 10–11, vs
panel r 0.50–0.57, n = 15–16, in the same harness), but it loses 6/17 cells at
the "no downhill origin" gate. The sketch's 0.5× r is computed on those
survivors,
so the comparison carries a survivor-selection effect; because the dropped
cells fail for absence of a fittable gradient — the hard cells, not easy ones —
the selection acts conservatively. Two harness caveats apply: the A and E arms
were run under their respective parity flags (A with `--windows`, E without),
so the harness is not fully symmetric (a symmetric rerun remains a check, not
a result); and cross-harness comparisons are invalid — the A-baseline measured
in the F1 harness (r 0.568 at 0.5×) differs from the committed arm-A number in
the C8 harness (0.913), so only within-harness A-versus-E contrasts are used
here.

Third, the mechanism is density, not determinism: per-window observed landmark
counts (17.3–17.5 sketch vs 14.1–21.3 panel) and detected fractions (0.175 vs
0.185) are near-identical at matched density. The sentence "deterministic
anchors are what keeps windows populated at 0.5×" is therefore retired. What
actually keeps the panel working at 0.5× is fusion redundancy: a single-stratum
sketch that finds no downhill origin has no second stratum to rescue it
(6/17 cells), whereas sixteen strata need only a subset of per-enzyme fits to
succeed. A property-by-property comparison of the two landmark sources,
with each row carrying its evidence, is given in Table 9.

*(Fig. 4)*

## 5. How many enzymes are needed, and over what GC range?

Four to eight enzymes, not sixteen, are sufficient; the recommendation holds
across the measured GC range of 25.4–72.0%, with a measured instrument
boundary at GC ≲ 30% combined with 0.5× depth.

Ranking the sixteen enzymes by standalone accuracy (each per-enzyme output is
already an independent V-shape fit, so the ranking is free), we find it
dominated by anchor yield without being dictated by it: the top six by
accuracy all carry > 1,900 anchors (BslFI, 2,501 anchors, is the yield exception
at rank nine), and the bottom four by accuracy (FalI, AloI, BplI, PsrI;
687/441/379/419 anchors) are near-zero or negative at 0.5×; the sparsest enzyme
of all, PpiI (336 anchors), is the accuracy exception at r = 0.61 at 0.5×. Sweeping panel size k ∈ {2, 4, 8, 12, 16} (accuracy averaged
over ≤ 2× on the Zheng grid) gives:

| enzymes | anchors | r | RMSE | run-out bias | s/sample |
|---|---|---|---|---|---|
| 2 | 17,055 | 0.951 | 0.188 | 0.073 | 3.2 |
| 4 | 24,753 | 0.959 | **0.165** | 0.074 | 4.4 |
| 8 | 37,232 | **0.969** | 0.171 | 0.123 | 6.5 |
| 12 | 41,662 | 0.969 | 0.183 | 0.132 | 7.8 |
| 16 | 43,735 | 0.959 | 0.196 | 0.145 | 8.2 |

accuracy is statistically indistinguishable across 4–12 enzymes (bootstrap
95% CIs overlap), with the point estimate peaking at k = 8 and RMSE at k = 4;
the four sparsest enzymes add no accuracy and double the run-out bias, while
cost is linear in k. A two-enzyme panel still
carries 17,055 anchors — comparable to Pilea's 18,261 sketch k-mers, the
like-for-like sketch-size comparison — and delivers r = 0.960 at 1×, where
Pilea's shipped gates report nothing. We report this as an empirical result,
not a mechanism: sparse enzymes are not more biased than dense ones (mean bias
at 0.5×, −0.319 vs −0.312), and forcing fixed-effect fusion weights does not
recover the loss (RMSE 0.581 vs 0.506 at 0.5×). This analysis is also what
exposed a double-counting bug (Bsp24I corrupting CjePI at k ≥ 8; Methods
§1.2); the sweep was repeated with the corrected ranking.

We then asked whether the 8-enzyme recommendation, derived at *E. coli* GC
(50.8%), holds across genomes. On 18 near-complete genomes spanning GC
25.4–72.0% (GTDB R232 metadata standing in for the planned R226 release; reads
were simulated with planted replication gradients, n = 8 per cell), three
results emerged. (i) Panel density rises *monotonically* with GC — ~4.5k
landmarks/Mb at 26% GC, ~9.4k at 50.8%, ~12.4k/Mb at 72% — with depression
only at the low-GC end; the both-ends-fall expectation and Syn2b's
non-monotonic peak were not observed. (ii) The recommendation holds: k8 and
k16 accuracy differ by < 0.07 r across the whole GC range at ≥ 1× (at 0.5× the gap grows to 0.12–0.21 r near the low-GC boundary), and even at 72% GC —
where k8 density falls to 0.64× k16 — accuracy is unharmed (r 0.98–0.99 for
both at ≥ 1×). (iii) At matched density (≤ 4% difference, recomputed per
genome) FracMinHash is flat across GC and ties the panel everywhere, so the F1
verdict generalises out of *E. coli*. The instrument boundary is measured, not
assumed: at GC ≤ 30% *and* 0.5× depth both modes fail (r −0.5–0.7, RMSE
0.6–0.9); the two failing genomes are 0.65 Mb *Buchnera*, so the boundary is
confounded with genome size and is flagged, not resolved. The stated operating
range is GC ≳ 30% at depth ≳ 1× (all cells at GC ≥ 30% and depth ≥ 1× give
r ≥ 0.85; at 2×/5× r ≥ 0.95/0.98).

*(Fig. 5)*

## 6. Fragmented references and MAGs

### 6.1 Fragmentation removes the coordinate, and the collapse is landmark-agnostic

Shuffling a closed chromosome into 100 lognormal contigs does not add noise to
the estimate; it removes the coordinate signal, and it does so identically for
enzyme and FracMinHash landmarks. On the Zheng protocol (complete chromosome vs
100 shuffled contigs vs scaffolded contigs; n = 16 media per cell, identical
reads throughout, 43,707 of 43,735 anchors surviving the cut so that only the
coordinate changes), fragmented references give r = 0.550 and slope 0.103 at
1× (complete: r 0.981, slope 0.779): every estimate lands at about a fifth of
truth, and the r = 0.859 at 10× (slope 0.210) makes the failure look far more
survivable than it is — the clearest case in the paper for reporting slope
beside correlation. Pilea, by contrast, wins outright on unscaffolded contigs
(r 0.827 vs 0.550 at 1×; its rank regression loses almost nothing,
0.889 → 0.827): a genuine architectural advantage of the sorted estimator that
we report as such. Re-running the protocol with matched-density FracMinHash
landmarks (`data/f4_fragment_sketch/`) gave an *identical* collapse — slope
≈ 1.0 → 0.26–0.38 and RMSE 0.02 → 0.84–0.87 at 5× in both arms — so the
mechanism is coordinate loss, full stop. Moreover, a permutation test (200
shuffled-coordinate refits) shows the residual high r on fragmented references
is partly a search artifact: the origin grid search fits a fake V to a flat
profile, and shuffled profiles return the same estimates, scaled with the
planted value. Slope and RMSE against truth are the honest metrics on
fragmented references; the r column of the original experiment was itself
inflated by this artifact. There is no safe contig count, and r cannot find
one: across 1–100 contigs at 10×, r stays 0.86–0.97 while the slope falls
0.95 → 0.21 and bias grows to −0.81; at 50 contigs (N50 156 kb, a good draft)
r reads 0.96 while every estimate is 44% of truth. Degradation is smooth and
monotone, so the rule is "scaffold anything not closed".

### 6.2 `scaffold` restores the complete-reference result, with either landmark source

`sk2bGrow scaffold` recovers the complete-reference accuracy essentially
exactly: contig order returns with Spearman 1.0000 against the *E. coli*
O157:H7 reference, orientation is 99/99 correct, and the 712 kb raw placement
error is a rigid rotation, invisible to a fit that searches for the origin
rather than assuming it (scaffolded vs complete: r 0.977 vs 0.981, slope 0.755
vs 0.779 at 1×; RMSE 0.086 vs 0.063 at 10×). The method is therefore not
restricted to closed genomes. The F4 experiment shows scaffolding is also
landmark-source-agnostic: with FracMinHash landmarks, self-scaffolding ties the
enzyme arm item-for-item (98–100/100 contigs placed, Spearman 1.0000, 5× RMSE
0.017–0.033 vs 0.017–0.030), and scaffolding against a 0.1%-divergence relative
places 100/100 with full PTR recovery in both arms. Against genuinely divergent
relatives, placement partially succeeds and the sketch arm is slightly *better*
than the enzyme arm (5× RMSE 0.223 vs 0.644 on *Bacillus*; directional only,
n = 6) — so "scaffolding needs the enzyme panel" is dropped as an argument for
the panel in the all-computational route; it remains a capability unique to the
enzyme route when only 2bRAD data exist. One caveat: the shipped `scaffold`
CLI digests enzyme tags only; the FracMinHash placement used a 1:1 Python port
of the Rust placement algorithm, validated by the near-perfect self-arm, so
sketch scaffolding is demonstrated, not shipped.

### 6.3 The coordinate is sufficient for our estimator, not necessary in principle — and the QC cannot see this failure

A spread-fit prototype recovers PTR from fragmented references without any
coordinate: because log₂ coverage is uniform across the genome with width
log₂(PTR), the spread alone determines the value, and fitting that distribution
with the per-window standard error in the model takes RMSE on contigs from
0.87 to 0.14 at 5×, returning 0.000 (rather than Pilea's 1.153) on the
stationary control at 1×. Below 5× it carries too little information
(bias −0.98 at 1×) and it does not match Pilea at depth (RMSE 0.227 vs 0.077
at 10×); we report it as a prototype whose main argument is architectural —
two estimators chosen by whether a coordinate exists. The fragmentation failure
also exposes a blind spot in the enzyme-consistency QC: 100% of fragmented
estimates pass at 5–10×, against 75% of correct complete-reference ones.
Cochran's Q asks whether the enzymes agree, and a destroyed coordinate makes
all sixteen agree that there is no gradient. A contig-count guard, not a
better Q, is the fix.

*(Fig. 6; Table 8)*

## 7. Real metagenomes

### 7a. Cross-method concordance in a fecal cohort

Two independently implemented methods agree on species-level growth rates in a
real community. On three fecal samples from Sun et al. (PRJNA689204; 125–135 Gb
per sample), comparing sk2bGrow (WGS) against Pilea at its shipped defaults on
common-denominator species × sample pairs (n = 58/68/78 per sample), the bias
is stable at +0.22–+0.25 log₂ (Pilea higher by ~17%) with limits of agreement
±0.4–0.9, r 0.34–0.51, and CCC 0.25–0.37. The low CCC is deflation from a
narrow observed dynamic range (span ~1.4 log₂ in both methods), not broken
agreement — every coefficient is printed beside its observed range. A gates-off
control arm shows Pilea's gate does real work: with its thresholds removed,
Pilea's outputs against sk2bGrow decorrelate entirely (r 0.00–0.10, LoA ±6),
symmetrically confirming on an independent dataset the recall caveat of §7c.
We frame this as concordance validation (n = 3 samples), not a biological
discovery.

### 7b. Real 2bRAD libraries: exact deduplication destroys the PTR signal

On the same cohort, comparing in-silico WGS against the real 2bRAD libraries
(same estimator, two preparations; BcgI only, so the single-stratum data
cannot exercise cross-enzyme fusion or the Q-based QC), exact deduplication
before counting destroys the PTR signal. Raw counts against the WGS arm reach
r = 0.797 (CCC 0.728, bias +0.60) in the best sample, while awk-exact dedup —
which collapses 95.7% of reads to unique sequences — drops agreement to
r 0.30–0.60 with bias −0.6 to −1.2 and flattens the dynamic range. A cap-dedup
sensitivity scan (cap ∈ {1, 2, 3, 5, ∞}) finds no robust intermediate: cap 2 is
the best deduplicating condition in one sample (r 0.831) and the worst in
another (r 0.205). The process recommendation is therefore to *not* exactly
deduplicate 2bRAD libraries before counting. The same paired data yield the
first measurement of per-anchor capture-efficiency dispersion in 2bRAD:
σ_eff ≈ 1.53 (Pilea cannot produce this quantity; it requires paired 2bRAD and
WGS). An ICC decomposition attributes a median 74.8% of the dispersion to
correctable within-batch effects; the residual (~0.7) is diluted by √99 window
averaging to ~0.07–0.10 against gradients of 1.0–1.7 — measured, correctable,
and negligible after correction, but so far *within-batch reproducible* only:
cross-laboratory transferability awaits an independent batch (Hou, in
preparation).

### 7c. MAG-scale application: recall under a common protocol, and a QC that works on real data

On a time-series RBC reactor (PRJNA974210; 522 MAGs × 9 samples), sk2bGrow's
QC-pass recall and Pilea's reported recall are overlapping under a common
denominator: sk2bGrow QC-pass 3.8–13.8% (20–72 of 522 MAGs) versus Pilea
5.0–11.1% (26–58 of 522), with every Pilea-passed MAG falling inside
sk2bGrow's output set. The headline raw figure of 522/522 (recall 1.00) is a
denominator artefact — this path has no output gate — and is never reported
unqualified; high yield includes wrong answers shipped. The enzyme-consistency
QC demonstrably filters fragmented references on real data: QC pass rate
anti-correlates with MAG contig count (ρ = −0.41; ρ = −0.34 controlling mean
coverage, p < 1e-4), with the relationship *strengthening* inside the
highest-coverage quartile (ρ = −0.432), so the association is not a power
artefact of low coverage — coverage is the strongest single predictor
(ρ = +0.516) and the two effects are partially entangled, which is the real
structure of this MAG population rather than a scoring artefact. After QC
filtering, cross-sample consistency of log₂PTR tightens ~4× (median std 0.030
vs 0.129). Cost, however, flips at this scale: sk2bGrow is 89.5–240.8× slower
per sample than Pilea at default on these deep samples (measured per sample,
`c5_cost_per_sample.tsv`), because the anchor-matching phase is 89.6% of count
time and scales as reads × reference anchors — both factors two orders of
magnitude larger than the isolate benchmark. Setting `--max-mismatch 1` narrows
the gap to ~10–25× (a 63× lookup speedup from avoiding the collapsed
middle-seed slot; 5.3% anchor loss on this low-divergence data, no genome
lost), and the M4 containment pre-screen clears the absent-genome
false-positive floor of the unrestricted index (4,048 of 4,700 guaranteed-absent
synthetic genomes receive nonzero counts unrestricted; zero with the screen),
with no measured speedup at this dataset's ~90% presence (0.99× wall A/B)
and a loss of 65/522 real MAGs at the current default threshold; its value is a
projection for database-scale presence rates (< 1%, ~21× anchor-set shrink). The
flip is therefore real and partially fixed: mismatch-1 is measured (subset), the
screen is measured on the false-positive side only.

*(Fig. 8)*

## 8. Computational efficiency and scale

Per sample on the isolate benchmark, measured identically for both tools
(`/usr/bin/time -l`, 8 threads, k-mer counting inside the timed region) and
averaged over the five depths on the 85 *E. coli* cells, sk2bGrow takes 8.2 s
and 195 MB at k = 16 (6.5 s at k = 8; 3.2 s at k = 2) against Pilea's 11.5 s
and 157 MB with gates off. The average understates the gap where it matters:
at 1× and 2× Pilea takes 16.8 s and 21.0 s against our 6.8 s and 8.3 s, i.e.
it is ~2.5× *slower* in the depth band this paper is about, and cheaper only at
0.5× and 10×. Pilea's non-monotonic profile (5.6 s at 0.5×, 21.0 s at 2×,
5.6 s at 10×) reflects its ZTP-mixture EM being slowest exactly where the
mixture is least identifiable. The earlier "Pilea is ~7× faster" claim did not
survive this measurement — it came from the simulation, where Pilea's shipped
gates skip fitting below 5× and it reports nothing.

At scale, mismatch tolerance is a budget decision, not an accuracy one. In
enzyme mode, empirical misassignment stays ≤ 1.4e-4 pooled at mm ≤ 2 (exact-
first suppression plus the motif gate make tolerance nearly free), so the
mm = 2 default stands and mm = 1 is safe if memory forces it (unmatched-window
loss 0.44% → 0.42% at 0.1% error). In sketch mode the answer is the opposite:
mm must lock at 0 — at scale ~100, the mm ≥ 1 rescue channel credits 1–2% of
off-sketch k-mers to a neighbouring sketch key at the wrong coordinate,
independent of sequencing error; a landmark–landmark census cannot see this
channel (the FracMinHash empirical arm here is a rule-based simulation, as the
sketch counting mode is not yet in the CLI; the enzyme arm is the measured Rust
counter, validated anchor-by-anchor against an independent replica). Index cost
is a wash between modes at matched density (~39 B per landmark either way,
equal build speed); the earlier 752 GB figure was peak RSS, with true disk
about a quarter of it. Projected to GTDB species representatives, a k16 enzyme
index is ~232 GB, k8 ~184 GB, and a half-density FracMinHash index (scale 200,
0.53× landmarks) ~124 GB — the practical budget answer, since half density
costs ≤ 0.007 r at ≥ 1×, with the only cost at 0.5× growing with GC (−0.03 r
at 72%). One engineering gap is stated honestly: a merged
`--mode fracminhash` index does not exist yet, so the half-density route is
projected, not shipped.

*(Fig. 8c, d)*
