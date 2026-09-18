> **SUPERSEDED 部分（2026-09-13）**：Table 2 已于 HPC C1 双端实例再生成（arm A 0.5–10×：0.923/0.912/0.958/0.971/0.970 (2026-09-17 signed fixed-origin)），形式交互检验仅 2× 显著为负。本文档中的旧数字与 "interaction" 主张以 manuscript.md 为准。详见 data/repro_check/ 与 data/repro_check/multiseed/INTERACTION_REPORT.md。

# Methods

Section order follows Pilea (Microbiome 2026) so the two papers can be read
side by side. Every number quoted here is reproduced by the scripts in
[sk2bGrow/benches](https://github.com/HuangShiLab/sk2bGrow) from the tables in
`data/`. Sections 3.5–3.6 and 6 report experiments that have no Pilea
counterpart — the landmark-source, robustness and HPC-scale questions are new
with this paper.

---

## 1. The sk2bGrow algorithm

sk2bGrow estimates the peak-to-trough ratio (PTR) of a genome from the
replication gradient in read coverage. For an exponentially growing population
with replication period *C* and doubling time τ,

$$\log_2 \mathrm{PTR} \;=\; C/\tau \;=\; \lambda C/\ln 2 ,$$

so log₂PTR is proportional to growth rate λ at fixed *C*. The pipeline
(**Fig. 1a**) has six stages.

### 1.1 Anchor database construction (offline, once)

Reference genomes are digested *in silico* with a panel of 16 Type IIB
restriction enzymes (§2). Each enzyme recognises a short degenerate motif and
cleaves on **both** sides of it, excising a fixed-length tag of 25–33 bp. A
**anchor** is one such tag together with its genomic coordinate.

The digest slides a window of the enzyme's tag length along both strands and
tests it against that enzyme's IUPAC pattern set. Pattern sets are
reverse-complement closed by construction (enforced by a unit test), so a locus
is found from either strand and is reported once. Windows containing a non-ACGT
base are rejected rather than resolved; soft-masked (lowercase) reference
regions are **kept**, which is a deliberate divergence from Fast2bRAD-M, whose
exact-byte matcher silently skips them.

Each anchor carries:

* its genome, contig and coordinate;
* the local GC fraction over ±250 bp, quantised to 0.5 % steps;
* uniqueness flags — unique within its genome, unique across the database,
  multi-copy, shared, non-chromosomal.

Multi-copy and shared anchors are flagged, not deleted; the flags are what the
EM step (§1.2) and the QC layer act on. For *E. coli* K-12 MG1655
(GCF_000005845.2) the 16-enzyme panel yields **43,735 anchors**, of which
**42,234 are usable** after masking: mean spacing 110 bp, median 59 bp, 99th
percentile 647 bp, largest gap 6,007 bp (3 gaps > 5 kb). That is **227 usable
anchors per 25 kb** against **98 k-mers per 25 kb** in a Pilea FracMinHash
sketch of the same genome (18,261 k-mers; k = 31, s = 250) — 2.3× the density,
before any difference in estimator.

The essential property is that anchor coordinates are **motif-defined and
therefore known before any read is seen**, and are identical in every sample.
A FracMinHash sketch is a random subset of k-mers whose positions carry no
usable structure, which is why sketch-based PTR estimators regress on window
*rank* rather than window *position*.

### 1.2 Read counting and anchor attribution

Reads are scanned exactly as reference genomes are: slide a window of each
enzyme's tag length and test the enzyme's patterns. A 2bRAD tag is not an
arbitrary k-mer but a fixed-length window under a motif constraint, so this is
both the correct model and cheaper than hashing every k-mer. It also makes
shotgun metagenome reads (route A) and bench-generated 2bRAD reads (route B) one
code path — in route B the read *is* the tag.

Matching is strand-canonical and tolerates up to **2 mismatches**. With a budget
of *m* mismatches a tag is split into *m*+1 contiguous seeds; by the pigeonhole
principle at least one seed survives intact, so probing all *m*+1 seed slots
cannot miss a true match. Candidates are verified by full Hamming distance
against the stored tag.

A tag is counted only when it lies **wholly inside** the read. For 150 bp reads
and a 32 bp tag this retains (150 − 32 + 1)/150 ≈ 0.79 of the local depth — the
same factor that applies to a k = 31 sketch, so the two methods are compared at
matched *effective* depth rather than matched nominal coverage.

A single physical read window can satisfy **two enzymes'** patterns, because the
panel is not disjoint: every Bsp24I site is also a CjePI site (§2). Such a window
is one observation and is credited **once to each stratum it belongs to**. The
read is therefore scanned once per distinct tag length, not once per enzyme —
scanning per enzyme visits the shared window twice and counts the locus twice for
both enzymes. That is not hypothetical: it was a defect in this implementation,
and it inflated CjePI's total count on *E. coli* K-12 by 21.5 %, concentrated at
1,352 loci hit exactly 2.00× of which 94.7 % were Bsp24I sites. Because the
inflation is *local* rather than uniform it distorts the coverage profile instead
of cancelling in the fit's intercept.

Anchors shared between reference genomes are split by expectation-maximisation,
following the shared-k-mer reassignment Pilea borrows from sylph. Each shared
count is apportioned in proportion to current genome abundances; abundances are
re-estimated from **unique anchors only** on every iteration, which keeps the
fixed point identifiable (a genome with no unique anchors cannot inflate itself
without bound). Convergence is a relative abundance change below 10⁻⁶ or 100
iterations.

### 1.3 Windows and window rates

Windows are cut **inside each (genome, enzyme) series**, in coordinate order,
holding a fixed *number of anchors* rather than a fixed number of base pairs.
Equal-anchor windows equalise statistical power along the genome; equal-bp
windows do not, because anchor density varies. Windows never straddle a contig
boundary — two anchors on different contigs have no defined genomic distance.

Anchor density varies ~20-fold across the panel (CjeI 1,962/Mb, PpiI 74/Mb on
*E. coli*), so a single window size cannot serve all sixteen enzymes: at a flat
100 anchors/window the three sparsest enzymes would get 3–4 windows each and
drop out. The size is therefore chosen per enzyme as

```
anchors_per_window = clamp(n_anchors / 25, 25, 100)
```

targeting ≈25 windows per enzyme. Sparse enzymes get smaller, noisier windows
and consequently *earn less weight* in the fusion (§1.6), which is the correct
way for them to count less.

Within a window, counts are modelled as a **zero-truncated Poisson (ZTP)
mixture** — parameters by EM, component count by BIC — and the highest-weight
component's rate is the window's expected coverage. Truncation matters because
systematic sequence divergence makes some reference anchors genuinely absent
from a sample; counting those as true zeros drags an untruncated fit downward.
A **zero-truncated negative binomial (ZTNB)** branch is available for residual
overdispersion and is selected by BIC rather than by a coverage threshold; at
low coverage NB is unidentifiable and BIC falls back to ZTP automatically.

Every window rate is returned with a standard error. That is what makes the
inverse-variance fusion of §1.6 possible and is the piece a bootstrap has to
approximate.

*Numerical note.* The NB log-likelihood needs the rising factorial
lgamma(k+r) − lgamma(r) with r = 1/α. As α → 0 both terms diverge and their
difference loses all precision in double arithmetic; the implementation instead
accumulates Σᵢ log(r+i) directly. Without this, a degenerate α ≈ 0 fit gains
~12 spurious log-likelihood units and BIC selects NB everywhere.

### 1.4 GC bias correction

Every Type IIB recognition site has its own base composition, so the 16 enzymes
sample sixteen different, individually narrow GC neighbourhoods. A single global
correction curve averages them into a shape that fits none of them. sk2bGrow
therefore fits a loess curve of log₂ efficiency against GC **within each
enzyme**, and applies it **at anchor resolution** (each anchor's own ±250 bp GC),
averaging into the window only afterwards — correcting on a window's mean GC
would discard exactly the within-window variation the correction exists to
remove.

The correction is an **offset in log₂ space**, never a rescaling of counts: the
window models need integer counts, and a multiplicative fudge would silently
break the zero-truncation.

A constant per-enzyme efficiency factor cannot bias PTR, because each enzyme is
fitted separately and a constant offset is absorbed by that fit's intercept.
Only GC *slope within an enzyme* matters. Per-enzyme efficiency factors are
still reported, for QC.

We use this historical single-pass correction as the primary default. A residual
second pass was tested as a diagnostic but not adopted. Holding reads and
estimator fixed, it changed 0.5×/1× anchor-panel accuracy from r = 0.923/0.912
to 0.848/0.879 while leaving the density-matched sketch arm essentially
unchanged (0.883/0.912), the pattern expected when position-linked GC structure
absorbs part of the ori-ter gradient. The residual mode remains opt-in
(`--gc-residual-pass`) and is not used in primary benchmarks
(`data/m1_signed/gc_correction_diagnostic.tsv`).

### 1.5 Origin placement and V-shape fitting

The origin is a property of the chromosome, not of an enzyme. It is estimated
**once per genome** from all enzymes pooled, after median-centring each enzyme's
rates (which removes its efficiency offset while leaving the shared gradient
intact). Fitting an origin separately per enzyme both wastes power and injects
between-enzyme variance that has nothing to do with biology — two enzymes
landing on origins 200 kb apart would make the consistency test of §1.6 report
disagreement where the enzymes in fact agree. An externally supplied origin
(DoriC, Ori-Finder, *dnaA*) is used when given.

With coordinates in hand, log₂ coverage is regressed directly on circular
distance *d* from the origin:

$$\log_2\mu(x) \;=\; a - b_1\min(d,k) - b_2\max(0,\,d-k), \qquad d = \mathrm{dist}(x,\ \mathrm{ori})$$

With *b*₁ = *b*₂ this is the plain V that CoPTR-Ref shows to be the maximum
likelihood model, and log₂PTR is the fitted drop from origin to terminus. The
two-slope form exists for multi-fork replication, where overlapping rounds put a
genuine kink in the profile at PTR > 2. Which of the two is used is decided by
**BIC** on the weighted residual χ² (matching the weighted fit objective); the
segmented form is only offered when ≥30 windows are available and the preliminary
log₂PTR exceeds 2. Sign selection occurs only while searching for the origin,
where rejecting an uphill solution distinguishes origin from terminus. Once that
shared origin is fixed, each enzyme's slopes are allowed to be negative: a
stationary control can yield a small negative gradient, and truncating those
slopes at zero manufactures a positive control bias. The window standard errors
are re-expressed as a smooth function of the first-pass fitted profile in a
second weighted pass, so a window's weight is not correlated with its own
upward count fluctuation.

Standard errors are scaled by the reduced χ² of the fit. The window standard
errors of §1.3 describe counting noise only; anchor efficiency noise, residual
GC structure and profile misspecification all add scatter on top, and taking the
residuals at face value is what keeps the error bars honest.

For comparison and for the attribution analysis, the **sorted-rank regression**
of iRep/Pilea is implemented on the same windows: sort window log₂ rates, regress
value on rank, multiply the slope by the window count. It needs no coordinates,
which is why it works on fragmented references — and its estimate is set by the
largest and smallest window, so it rides on two extreme order statistics
(**Fig. 1c**).

### 1.6 Cross-enzyme fusion and quality control

Each of the 16 enzymes measures the same biological quantity through a different
set of loci, with its own digestion efficiency and GC neighbourhood. That is a
stratified design with real replication built in.

**Fused estimate.** Inverse-variance weighting is the minimum-variance linear
combination of unbiased estimates, so enzymes with more anchors and cleaner fits
count for more automatically. Enzymes with fewer than 30 usable anchors, no
finite estimate, or no usable standard error are excluded and named in the
output.

**Consistency test.** Under the null that every enzyme measures one common
value, Cochran's *Q* is χ² on *k* − 1 degrees of freedom. A significant *Q*
(α = 0.05) means the enzymes disagree — too few anchors, a methylation-blocked
site class, a mis-assembled region — and is a QC signal available at zero extra
sequencing cost. *I*² is reported alongside.

**Escalation.** When *Q* rejects, the fixed-effect standard error is known to be
too small, because it assumes the only scatter is sampling noise. The estimator
then switches to DerSimonian–Laird random-effects weights 1/(sᵢ² + τ̂²), with τ̂²
the moment estimator of the between-enzyme variance. Fig. 1d is a fixed-effect example (*I*² = 0.00) chosen to show the strata; random-effects escalation is exercised in the QC audit and reported with fusion model, *Q*, *I*² and τ² for every estimate.

A sample passes QC on estimated coverage, dispersion, the fraction of reference
windows covered, and EM containment; failures are reported with a reason rather
than dropped silently.

**What the consistency test caught.** Cross-enzyme heterogeneity was the signal
that exposed the double-counting defect described in §1.2: before the fix,
Cochran's *Q* rejected in 56–69 % of samples with mean *I*² = 0.34–0.42, and the
random-effects escalation was firing routinely. After the fix the same samples
give mean *I*² = 0.07–0.28 and *Q* rejects in 6–24 %. The enzymes were correctly
reporting that something was wrong with one of them.

**Known gap.** The fusion currently treats all 16 enzymes as independent strata.
They are not: Bsp24I's site set is *totally contained* in CjePI's (1,636/1,636,
891/891 and 2,910/2,910 tags on three genomes — the two patterns, written in
opposite strand orientations, differ at exactly one position), and about half of
Bsp24I's tags are also CjeI sites. The panel therefore offers at most ~15
independent strata, and *Q* is mildly conservative (positive dependence deflates *Q*, so cross-strata heterogeneity is under-detected rather than over-called). The containment
relations are recorded in the code (`enzyme::CONTAINMENTS`) but are not yet
consumed by the weighting.

### 1.7 Implementation

Digestion, counting, EM and windowing are Rust (rayon-parallel); the statistics
layer is Python (numpy/scipy/statsmodels). The two communicate through TSV count
tables, so either can be run alone — the A/B benchmarks in §4 exploit this to
hold one layer fixed and vary the other.

---

## 2. Enzyme panel

The panel is the 16 Type IIB enzymes of the 2bRAD-M lineage: AlfI, AloI, BaeI,
BcgI, BplI, BsaXI, BslFI, Bsp24I, CjeI, CjePI, CspCI, FalI, HaeIV, Hin4I, PpiI,
PsrI.

Recognition patterns and flank offsets were **transcribed from
`2bRADExtraction.pl`**, the reference implementation, rather than from secondary
tables; five of sixteen definitions we had initially taken from secondary
sources were wrong. Correctness was then checked by comparing measured anchor
density against the published per-genome densities on three genomes
(*E. coli* K-12, *B. subtilis* 168, *P. putida*): **46 of 48 cells agree within
3 %** (Table 1). HaeIV differs by exactly 2.00×, a locus-versus-window counting
convention; Hin4I is unreconciled — no IUPAC variant, spacer length or
deduplication convention we tried reproduces the published 1,650/2,010/1,057.

BslFI is a Type IIS enzyme (GGGAC, 10-11/14-15) that cuts on one side only;
it is retained because a fixed-length tag is still excised, but it is not a
Type IIB enzyme and is labelled as such.

---

## 3. Benchmark datasets

### 3.1 *E. coli* growth-rate panel (Zheng et al. 2020)

**Source.** BioProject PRJNA615952 — *Escherichia coli* K-12 MG1655 grown to
steady state in a panel of defined media spanning a ~4.3-fold range of growth
rate, each with an independently measured λ and an independently measured
C + D period. This is the same dataset Pilea uses for isolate accuracy.

**Samples used.** One run per medium, first replicate: 16 growth media
(λ = 0.399–1.721 h⁻¹, doubling time 0.40–1.74 h) plus one **run-out** sample.
Accessions are in `benches/zheng2020/picks.tsv`.

**Ground truth.** Two references are used and kept distinct:

* λ, the directly measured growth rate — used for correlation, because it is
  measured, not modelled;
* the *predicted* log₂PTR = λ·C/ln 2, using the C period reported by the same
  study — used for RMSE and slope, because a correlation cannot detect a
  compressed dynamic range.

  This second reference is **not fully independent of sequencing**: C was
  determined in part by marker-frequency analysis, which is itself an ori:ter
  ratio. RMSE and slope against it are therefore a *consistency check on
  magnitude*, not an independent validation. The correlation against λ, an
  optical growth measurement, is independent, and is the number to read as
  validation.

**Run-out sample.** Pilea excluded the run-out samples. We keep one
(`E.coli_gDNA_RUN_OUT`, sample alias CON-C) as a **negative control**. In a
replication run-out, ongoing rounds are allowed to finish while new initiation
is blocked, so every chromosome in the population ends fully replicated and
log₂PTR must be ≈ 0. This tests a failure mode that correlation against a growth
panel cannot see, and it costs nothing: a method that reports a large PTR here is
reading noise as a gradient.

**Coverage titration.** Both mates are 150 bp and the committed Zheng grid is paired-end for every arm. For each run the first 600,000 read pairs were retained (≈19× per end), then truncated to nominal **0.5×, 1×, 2×, 5× and 10×** at the pair level. SRA order is not random with respect to quality; a 2× check gave file-order r = 0.974 versus fixed-seed random r = 0.957 (`data/m1_rerun/fact_checks.txt`). All 16 media, plus the control, are present at every depth. The committed A/B/E estimates are the 2026-09-18 paired-end statistics rerun with signed fixed-origin fitting and single-pass GC correction; the regenerated summary is `data/m1_signed/hpc_singlepass_grid_summary.tsv`.

This titration is a **deviation from Pilea, which ran full depth only**. It is
the axis the paper is about: PTR at metagenomic per-strain depth, not at isolate
depth.

**Reference.** A single complete genome, GCF_000005845.2. Every real-data result
in this paper is therefore one strain against a complete reference; nothing here
establishes behaviour on a metagenome or a fragmented MAG.

### 3.2 Multi-strain simulation

Follows Pilea's Fig. 3 design, at a scale a laptop can run.

**Genomes.** 16 complete RefSeq bacterial chromosomes
(`benches/accs.txt`), deliberately including *E. coli* K-12, *E. coli* O157:H7
and *Shigella dysenteriae* as a shared-anchor stress test — those three share a
large fraction of their tag sets, so the EM attribution of §1.2 is exercised
rather than bypassed.

**Generative model** (Pilea's, verbatim): replication initiates at position zero
and terminates at the assembly midpoint; coverage decreases log₂-linearly from
origin to terminus; log₂PTR ~ Uniform[0, 2] independently per strain. Reads are
150 bp single-end, drawn multinomially over 1 kb position bins with weights
w(x) ∝ 2^(−log₂PTR · d(x)), d ∈ [0,1] the normalised distance from the origin,
and reverse-complemented with probability 0.5. No sequencing error is added —
which favours the exact-match k-mer method, not ours.

**Grid.** {4, 8, 16} strains × {1, 2, 4, 8}× per-strain coverage × 2
replicates = 24 cells and 224 genome-level truth values (Pilea's gates-off arm
returned 223 of them). Strains are sampled without replacement per cell.

**Deviations from Pilea's grid**, all in the direction of a smaller experiment:
Pilea used 120 genomes, up to 32 strains, up to 32× and 400 samples. The
laptop-scale grid does not establish behaviour at 32 strains or 32×. Three unfavourable variants on one 16-genome community were also run: 1% substitution errors (bias +0.001 versus control +0.001), a real GC gradient (bias −0.195), and near-neighbour O157:H7 at mismatch 0/1/2 (mean bias −0.068/−0.089/−0.080).

### 3.3 Reference fragmentation

The coordinate fit of §1.4 needs a genomic coordinate, which a draft assembly
does not have. Following Pilea's Fig 3, the *E. coli* reference is cut into 100
contigs with lognormal lengths (μ = 0, σ = 1, floor 500 bp — the indexer's own
minimum, so no contig is discarded), the order is shuffled and each contig is
independently reverse-complemented with probability 0.5.

Sequence content is untouched: the fragmented database holds 43,707 anchors
against the complete genome's 43,735, the 28 lost being tags that straddled a
cut. Reads are identical. **The coordinate is the only variable.**

Four reference conditions:

| condition | reference |
|---|---|
| complete | the finished chromosome |
| frag | 100 contigs |
| scafSelf | `frag` re-ordered by `sk2bgrow scaffold` against the same genome |
| scafRel | `frag` re-ordered against *E. coli* O157:H7 |

`scafSelf` is circular by construction and is reported only to separate
"scaffolding does not work" from "scaffolding does not work across strains".
`scafRel` is the case a MAG user faces.

Scaffolding places contigs by shared tags and is scored against the known
layout. Against O157:H7 it places 99 of 100 contigs, orients 99 of 99 correctly,
and recovers the contig **order** exactly (Spearman 1.0000 against the truth).
Its median start error of 712 kb is almost entirely a rigid rotation — O157:H7's
origin sits elsewhere — and §1.5 searches for the origin rather than assuming
it, so a rotation does not reach the estimate; after removing it the median
residual is 73 kb, 1.6% of the chromosome, from O157:H7's strain-specific
insertions.

The scaffolded contigs are re-emitted as a single pseudo-contig with each contig
written at its inferred start and the gaps filled with N, which the digester
rejects as ambiguous and so contributes no anchors. Where two contigs are placed
overlapping — possible when the coordinates come from another strain — the later
write wins; against O157:H7 this loses 89,597 bp, 1.93% of the draft. The
fragmentation results are therefore a slight *under*-estimate of what
scaffolding delivers.

### 3.4 An order-free estimator, as a control on §3.3

Fragmentation is only fatal if the coordinate is necessary. It is not, in
principle: since log₂(coverage) is a tent function of position with equal |slope|
on both replichores, over a uniformly tiled genome the coverage *values* are
uniform on [log₂ c_ter, log₂ c_ori] and the width of that uniform **is**
log₂(PTR). Position never enters, which is exactly why the sorted-rank estimator
of §4 is indifferent to fragmentation.

The reason not to adopt that estimator wholesale is that it reads the *observed*
spread directly. Sorted values of any sample rise, so sampling noise is credited
as growth; on the stationary control it returns log₂PTR = 1.15 at 1×.

As a control we therefore fit the same distributional model with the noise
included. Window *w* contributes *y_w* = *μ_w* + *e_w* with *μ_w* ~ U(*a*, *a*+*W*)
and *e_w* ~ N(0, *s_w*²), where *s_w* is the standard error §1.3 already reports.
Marginalising *μ_w*,

    p(y_w | a, W) = [Φ((a + W − y_w)/s_w) − Φ((a − y_w)/s_w)] / W

is maximised over (*a*, *W*) per enzyme and the per-enzyme widths are fused by
the same inverse-variance scheme as §1.6. The estimator is order-free by
construction and cannot manufacture a gradient from its own sampling error.

It is reported as a **prototype**, not as part of the method: it recovers most of
what fragmentation destroys above 5× but carries too little information below,
and it does not match the sorted-rank estimator at depth. Its role here is to
separate "the coordinate is necessary" from "our estimator needs the coordinate".

A contig-count sweep (2, 5, 10, 20, 50, 100 contigs at 10×, same 16 media)
locates the threshold rather than assuming one. It finds that there is not one:
degradation is smooth and monotone in bias from two contigs upward, and Pearson
r stays between 0.86 and 0.97 across the whole range while the fitted slope falls
from 0.88 to 0.21. Correlation is therefore not a usable diagnostic for
fragmentation, which is the concrete case behind §5's warning that r is
invariant to affine transformation.

### 3.5 Faecal microbiome (Sun et al.)

**Source.** BioProject PRJNA689204 — deep shotgun sequencing of three human
stool samples (SRR13371683, SRR13371682, SRR13371681; 2 × 150 bp,
125–135 Gb per sample, ≈394 Gb in total) with bench-generated 2bRAD libraries
(BcgI) prepared from the same samples. This is the paper's only real route-B
data, and it bounds what the dataset can show: with a single enzyme, the
cross-enzyme fusion, Cochran's *Q* and the panel-level QC of §1.6 are not
testable on it.

**Observation unit.** Species × sample; reads are attributed to species-level
references with the EM step of §1.2. For the Pilea arm the abundance floor is
set to **p ≥ 0.015 %** so that its coverage gate lands at the same effective
depth as in shallower studies (the floor scales with the measured 125–135 Gb
per-sample depth).

**Three arms.** (Arm labels here are local to this experiment and distinct
from the attribution arms A/B/E of §4.)

| arm | input | method |
|---|---|---|
| A | WGS reads | Pilea (§4), default gates and gates off |
| B | WGS reads | sk2bGrow |
| C | real 2bRAD library | sk2bGrow |

Arm A vs B is the cross-method comparison on identical input; B vs C is the
in-silico vs real-library comparison for one method.

**Exact dedup.** In the published 2bRAD data, reads were collapsed by exact
sequence identity (an awk exact-match dedup; 95.7 % of reads fold onto a
unique sequence). Because dedup status drives the arm-C result, arm C is run
on the deduplicated counts, on the raw pre-dedup counts, and on a
capped-dedup sensitivity series C′ = min(C, cap) with cap ∈ {1, 2, 3, 5}.

**Per-anchor efficiency.** For each species with ≥30 shared anchors observed
in all three samples, the per-anchor relative counts define an empirical
capture-efficiency profile *eᵢ*; its dispersion σ_eff and its between/within
decomposition are defined in §5.

### 3.6 Rotating biological contactor metagenome (PRJNA974210, "C5")

**Source.** BioProject PRJNA974210 — 9 biofilm samples from the flowpath of a rotating biological contactor (Hong Kong), the dataset of Pilea's MAG-scale benchmark
(67–79 M read pairs each). References: **522 MAGs**. The project's NCBI
esummary lists 525 assemblies; 523 were collected as `.fna`, and 3 archaeal
MAGs were removed by the domain filter, leaving 522 bacterial MAGs (24.1 M
anchors at k16). MAG quality (completeness, contamination) is scored by
CheckM2 1.1.0; contig counts and N50 are computed from the FASTA directly.

**Reporting conventions.** sk2bGrow has no output gate — every MAG gets a row
— so recall is reported three ways: *reported_fraction* (rows returned / 522,
always 1.00), *estimate_fraction* (rows with a finite PTR estimate / 522) and
*qc_recall* (QC passes / 522, the denominator-matched counterpart of Pilea's
default-gate output). Where Pilea's gates-off arm crashes (a `min_samples`
failure on 6 of 9 samples, stably reproducible) it is reported as a property
of the arm, not censored.

**Cost instrumentation.** Wall-clock and peak RSS per stage from
`/usr/bin/time -v`; the count stage is further decomposed by phase timing on a
1 M read-pair subset (database load, index build, match, window write-out).

**Mismatch-1 arm.** `--max-mismatch` is a runtime parameter (§1.2), so the
mm = 1 arm reruns counting against the same index. A posting-list diagnostic
(`seed_hist`: key count and mean/hit posting length per seed slot) is used to
explain where lookup time goes.

**Containment screen and dilution experiment.** To test whether a containment
prefilter makes unrestricted-database counting feasible at GTDB scale, a
diluted database was built: the 522 real MAGs plus **4,700 synthetic random
genomes** (3.5 Mb each, fixed seed, absent by construction), ≈242 M anchors
in total, indexed with `--screen-scale 2000`. On an identical 1 M read-pair
subset, counting is run unrestricted and with `--screen`. The synthetic
genomes are a false-positive census (truth = zero counts); the real MAGs are
an equivalence check (screened counts must be a subset of the unrestricted
counts, with per-anchor ratios ≈ 1).

### 3.7 GC-spanning genome set

Eighteen complete NCBI assemblies selected from GTDB R232 metadata (a local
GTDB R226 snapshot was substituted where R232 metadata was unavailable on the
analysis machine; all measured quantities are FASTA-derived and unaffected)
with: ≤5 contigs, CheckM2 completeness ≥98 %, contamination ≤1 %, 1–2
representatives per genus, and measured GC spanning **25.4–72.0 %** in
roughly even bands (25.4, 26.1, 30.6, 33.0–33.5, 43.5, 50.7–50.8, 57.0,
61.5, 65.3, 66.0–66.1, 72.0). The set feeds the F2–F4 experiments (§6). At
its low-GC end the set contains two *Buchnera* genomes of only 0.65 Mb, so
low GC and small genome size are confounded there and are read as a joint
boundary, not separated.

---

## 4. Comparison methods

**Pilea** (v1.3.8, installed from source) was run on identical
input files, single-end, 8 threads, with the database built from the same
reference genomes. Sketch parameters are Pilea's defaults: k = 31, s = 250,
w = 25,000 bp.

Pilea is run in **two configurations**, and both are reported everywhere:

* **defaults** — `--min-cove 5 --min-frac 0.75 --min-cont 0.25`;
* **gates off** — `-x 0 -z 0 -c 0`.

The second configuration exists because at its shipped defaults Pilea returns
*no estimate at all* below 10× nominal coverage, and a benchmark of
"estimate vs no estimate" is not a comparison of accuracy. Reporting only the
gates-off arm would be equally misleading in the other direction, since the
gates are part of the shipped tool.

**Which gate binds.** The gates-off run still reports the three gated statistics
per genome, so the ablation is exact and needs no extra runs: re-applying each
default threshold to the gates-off output shows that `--min-cove` is responsible
for **68 of 68** suppressed *E. coli* estimates and **167 of 172** in the
simulation (sole cause in 145). `--min-frac` and `--min-cont` are essentially
never binding at these depths. Because a 150 bp read yields 120 31-mers, Pilea's
threshold of 5 on *k-mer* coverage corresponds to ≈6.25× *read* coverage (5 × 150/120); at 5×
nominal the measured median is 3.93 and the gate fires.

The gate is not arbitrary conservatism. At 0.5× nominal, Pilea's gates-off
output is fully degenerate — measured coverage exactly 1.000 and log₂PTR exactly
0.000 for all 17 samples — so the gate is correctly refusing to report. Between
1× and 5×, however, the gates-off estimator is informative (r = 0.89, 0.95, 0.95
against λ) and is nevertheless suppressed.

**Attribution arms.** sk2bGrow differs from Pilea in two places at once — what
is counted (deterministic 2bRAD anchors vs a FracMinHash sketch) and how the
gradient is fitted (a coordinate V-fit vs sorted-rank regression). Comparing the
two end-to-end pipelines cannot say which change is responsible, so the two
factors are crossed into a full 2 × 2 and all four cells are run on the same
subsampled reads:

| | coordinate V-fit | sorted-rank regression |
|---|---|---|
| **2bRAD anchors** | sk2bGrow (arm A) | arm B |
| **FracMinHash** | arm E | Pilea, gates off (arm C) |

Arm B takes sk2bGrow's anchor counts and estimates from them exactly as Pilea
does (25 kb fixed-bp windows, sorted-rank regression with RANSAC). Arm E is the
transpose: Pilea's own sketch output is rewritten into sk2bGrow's count-table
format — one row per hashed locus, with its genome coordinate — and then passed
through the unmodified sk2bGrow estimator, so it receives the same zero-truncated
mixture, the same GC and outlier handling, the same per-enzyme inverse-variance
fusion, and the same standard errors. Building arm E by feeding raw window rates
into the V-fit would *not* be a controlled comparison: an earlier attempt at
that returned log₂PTR = 2.65 against a measured 1.73, because the surrounding
machinery, not the fitting geometry, was what had been removed.

The 2 × 2 is what licenses any statement of the form "the gain comes from X".

**Faecal three-arm comparison.** The Sun dataset (§3.5) adds a second,
independent arm structure: A = Pilea on WGS, B = sk2bGrow on WGS, C =
sk2bGrow on the real 2bRAD library. These labels are local to that experiment
and unrelated to the attribution arms above. Pilea is run at both its default
gates and gates off (thresholds as above); sk2bGrow is run on the
deduplicated, raw and capped-dedup count tables (§3.5), because the effect of
exact dedup on the 2bRAD route is itself one of the questions.

---

## 5. Evaluation metrics

Let *ŷᵢ* be an estimated log₂PTR and *yᵢ* the corresponding reference value, over
the *n* genome-estimates a method actually returned.

| metric | definition | what it detects | what it misses |
|---|---|---|---|
| **Pearson r** | corr(*ŷ*, λ) across media | whether the ranking of growth rates is recovered | a compressed or inflated scale — r is invariant to affine transformation |
| **Spearman ρ** | rank correlation | monotone recovery, robust to outliers | magnitude, as above |
| **RMSE** | √(mean(*ŷᵢ* − *yᵢ*)²) | magnitude error, in log₂PTR units | the direction of the error |
| **bias** | mean(*ŷᵢ* − *yᵢ*) | systematic over- or under-estimation | scatter (a method wrong by ±1 in equal measure has bias 0) |
| **slope** | OLS slope of *ŷ* on λ | dynamic-range compression | offset |
| **recall** | *n* reported / *n* truly present | silence — a method that reports nothing scores no error | whether what *was* reported is right |
| **L2** | ‖*ŷ* − *y*‖₂ over reported genomes | Pilea's own accuracy metric | it *rewards silence*: fewer terms, smaller norm |
| **spurious** | genomes reported that are not present | false positives | — |

Three points of method:

1. **Recall and error must be read together.** L2 (Pilea's metric) is computed
   only over the genomes a tool chose to report, so a tool that reports 2 of 16
   strains accurately scores better than one that reports all 16 with small
   errors. Every simulation table here reports recall beside RMSE for that
   reason, and RMSE (per-genome) is preferred to L2 (per-sample) because L2
   grows with the number of genomes reported.

2. **RMSE and bias are computed against the *predicted* log₂PTR** (λ·C/ln 2,
   from the source paper's own measured C + D), not against λ, since λ and
   log₂PTR differ by a unit-bearing constant. Correlation uses λ directly.

3. **A constant output has no correlation.** At 0.5× Pilea's gates-off arm
   returns log₂PTR = 0 for every sample; Pearson r is undefined there and is
   reported as absent rather than as zero, which would read as "scored badly"
   instead of "produced no signal".

**Agreement between methods (faecal three-arm).** Where two methods estimate
the same genome in the same sample, agreement is quantified by:

* **Bland–Altman** — bias = mean(*ŷ* − *y*) and the 95 % limits of agreement
  bias ± 1.96·SD of the pairwise differences. Detects systematic offset; a
  narrow observed dynamic range depresses correlation-based agreement without
  implying mis-ranking.
* **CCC** — Lin's concordance correlation coefficient,
  ρ_c = 2*r*σ_ŷσ_y / (σ_ŷ² + σ_y² + (μ_ŷ − μ_y)²): the Pearson correlation
  penalised for deviation from the 45° line. Low CCC together with tight
  limits of agreement indicates dynamic-range compression rather than
  disagreement.

**Per-anchor efficiency dispersion σ_eff and ICC decomposition.** On the Sun
data, per-anchor relative counts *eᵢ* (observed counts normalised by the
species mean, per sample) carry a capture-efficiency component beyond Poisson
noise; its magnitude is σ_eff, the CV of *eᵢ* after Poisson correction,
estimated per species (§3.5). A variance-components decomposition of *eᵢ*
across the three samples splits the dispersion into a between-anchor component
(a stable locus property, correctable by an empirical *e_prior*) and a
within-anchor component (sample-specific residual); the intraclass correlation
ICC = σ²_between / (σ²_between + σ²_within) is the correctable fraction.
Because *eᵢ* is spatially uncorrelated along the genome, window averaging
divides the residual by √(anchors per window): σ_eff bounds the window-level
excess coefficient of variation, not the fitted PTR directly.

**Computational cost.** Wall-clock and peak resident set size are taken from
`/usr/bin/time -l`. On macOS this propagates a grandchild's `ru_maxrss`
(verified explicitly), so the figure includes the Python statistics layer that
the Rust binary spawns. Both tools were given 8 threads. Cost is reported per
sample, with the one-off index/database build reported separately.

**Uncertainty quantification and analysis provenance.** Unless stated
otherwise, correlations and fitted slopes are reported with 95% confidence
intervals from a media bootstrap: the media are resampled with replacement
10,000 times, the statistic is recomputed on each resample, and the interval is
the 2.5–97.5% percentile range. Comparisons between arms evaluated on the same
cells (e.g. the two landmark sources in F1) use a paired bootstrap resampling
the same media for both arms, and are reported as a difference with its
bootstrap interval rather than as a hypothesis test. The source × estimator
interaction of the Results is tested on the Fisher-z scale as
(z_A − z_B) − (z_E − z_C) under paired resampling. The primary inference uses
three independent subsampling instances (48 media × instance units): the
interaction is significantly negative only at 2×, null at 1×/5×/10×, and
untestable at 0.5× (`data/repro_check/multiseed/INTERACTION_REPORT.md`). The Zheng
titration, multi-strain simulation, fragmentation protocol and Sun three-arm
analysis were pre-specified; the F1–F5 robustness experiments are post-hoc
analyses undertaken after internal review, and are labelled as such at first
mention.

---

## 6. Landmark-source and robustness experiments (F1–F5)

All five experiments isolate the **landmark source** — deterministic 2bRAD
anchors vs a density-matched FracMinHash sketch — from the estimator around
it. F1 and F2 ask whether accuracy differences survive density matching; F3
quantifies the misassignment risk of mismatch-tolerant counting; F4 asks
whether the fragmentation results of §3.3 are landmark-agnostic; F5 costs the
two sources at GTDB-like scale. The historical flag convention is retained for provenance but is functionally a
no-op (the current statistics layer does not consume `--windows`; symmetric
reruns are byte-identical). A matched pair is only ever compared inside the same
harness, and cross-harness numbers are never quoted against each other.

### 6.1 Density matching on the Zheng grid (F1)

**Purpose.** Test whether the anchor-vs-sketch accuracy gap survives once
landmark density is matched — i.e. whether the panel's advantage is density or
structure.

**Inputs.** The Zheng grid (§3.1) at 0.5–10× nominal coverage.

**Arms.** Panel arms use nested enzyme subsets of increasing density
(k2/k4/k8/k16); sketch arms use Pilea FracMinHash with the scale chosen so
that landmarks per Mb agree with the matched panel arm within ≈2 %
(s268/s184/s123/s104 against k2/k4/k8/k16), plus two over-density references
(s30/s60, 3.5×/1.8× panel density) as controls. FMH counts are rewritten into
armE count tables and pass through the sk2bGrow estimator (§4).

**Outputs.** Pearson r per (arm, depth) on the four matched pairs;
per-window landmark counts and detection fraction; per-cell failure cause.

### 6.2 GC sweep (F2)

**Purpose.** Measure panel anchor density and PTR accuracy as functions of
genome GC content.

**Inputs.** The 18 genomes of §3.7. Reads come from a purpose-built
planted-gradient simulator: the reference simulator is count-level Monte
Carlo with no read output, so a read-level simulator was written to the same
model — 150 bp single-end, no sequencing error, GC-neutral, 50 %
reverse-complement; log₂ copy(*p*) = −log₂PTR·dist(*p*, ori)/(*L*/2) with the
origin planted at *L*/2; depth *d* is the mean per-base depth over the
terminus half; start positions are drawn copy-weighted, an inhomogeneous
Poisson process with shot noise included.

**Grid.** planted log₂PTR ∈ {0.5, 1.0, 1.5, 2.0} × 2 seeds (n = 8 truth
replicates per cell) × depths {0.5, 1, 2, 5}× × arms {panel k8, panel k16,
FMH matched to k8, FMH matched to k16}. The matching scale is recomputed per
genome from its measured anchor density, scale = round(*L*/n_anchors), and
the matched FMH density agrees with the panel within ≤4 %.

**Outputs.** r/RMSE/slope per GC band and depth; density per arm per genome
(landmarks/Mb).

### 6.3 Mismatch tolerance and misassignment (F3)

**Purpose.** Quantify the risk that mismatch-tolerant counting credits an
observation to the wrong landmark.

**Census.** For each of the 18 F2 genomes, every pair of indexed anchors in
the same tag-length group at Hamming distance 1 or 2 is enumerated
(distance-0 pairs are the multi-copy families themselves, §1.1). Multi-copy,
shared and non-chromosomal anchors are included, because the counting layer
can retrieve them; each pair is classified uu / um / mm by uniqueness flags
(both unique / exactly one unique / neither).

**Empirical, enzyme panel.** Simulated reads — the F2 simulator with the
gradient removed (misassignment needs no gradient), 0.5×, i.i.d. substitutions
at 0, 0.1 and 1 %, no indels, no quality bias — are counted by the real Rust
counting layer at `--max-mismatch 0/1/2`. A Python replica of the full
counting rules (exact-first suppression, seed ranges, best-distance tie
handling, multi-mapper retention) attributes each recorded observation — read
names encode the true origin — as correct / wrong_tie / wrong_only. The
replica is validated anchor-exact against the Rust output on every cell
before use.

**Sketch arm, stated plainly.** sk2bGrow has no sketch counting mode in the
CLI, so the FMH arm is a **rule-level simulation**: the same counting rules
applied to the k16-matched sketch keys over a full *k* = 31 sliding window,
with no motif gate — mirroring Pilea's per-read k-mer counting. It answers
what *would* happen were sketch counting to accept *m* mismatches, not what
any shipped tool does.

**Outputs.** Census pairs per 1,000 landmarks (uu/um/mm × d ≤ 1 / d ≤ 2);
misassigned fraction of recorded observations per (mode, error rate, *m*).

### 6.4 Fragmentation × sketch (F4)

**Purpose.** Two questions: does the coordinate-fit collapse of §3.3 recur
identically on sketch landmarks, and can contigs be scaffolded on FMH
landmarks as well as on enzyme tags?

**Inputs.** Three genomes from the F2 set (*E. coli* 50.8 % GC, *B. subtilis*
43.5 %, *P. aeruginosa* 65.3 %). The F2 genomes are 3–4 contigs and are first
concatenated into a single pseudo-chromosome (the index's coordinate system),
then fragmented with the §3.3 protocol (100 lognormal contigs, shuffled,
independently reverse-complemented, seed 0). *P. aeruginosa* has no close
sister in the set, so its cross-strain scaffold reference is a **synthetic
0.1 %-divergence SNP mutant** (7,534 SNPs, fixed seed, labelled synthetic);
the other two use real sister species with large-scale rearrangements. Reads
are the F2 simulated FASTQ reused: planted log₂PTR ∈ {0.5, 1.0, 2.0} × 2
seeds × {1, 5}×, n = 6 per (genome, condition, depth).

**Arms and conditions.** Enzyme arm: sk2bGrow counting at C1 parity; sketch
arm: Pilea-interpreted FMH counts at the k16-matched scale at armE parity.
Reference conditions as in §3.3: complete, frag, scafSelf, scafRel.

**Sketch scaffolding, stated plainly.** The shipped `sk2bgrow scaffold`
digests enzyme tags only; the FMH arm therefore runs a **1:1 Python port of
the placement algorithm** (`scaffold.rs`) onto Pilea FMH keys — same
shared-tag mapping, Kendall-vote orientation, median start placement,
min_tags = 3, min_concordance = 0.8, same output schema — so downstream
scoring and pseudo-contig rebuild are unchanged. Sketch scaffolding here is a
port, not the shipped binary; its correctness is supported by the near-perfect
self-arm placements.

**Permutation null.** Whether an estimate on fragmented references exceeds
chance is judged by shuffling window positions and refitting the identical
V-fit **200 times** per cell.

**Outputs.** PTR recovery (r/RMSE/slope) per condition and arm; placement
rate, orientation accuracy and order Spearman; permutation null distributions.

### 6.5 Index cost at GTDB-like scale (F5)

**Purpose.** Cost the two landmark sources at matched density, feeding the
GTDB-scale budget decision of §3.6.

**Genomes.** 300 genomes, 1–20 Mb, GC 24.3–73.2 %: 282 drawn at random from
the local GTDB pool (2,200 random draws thinned to a GC-uniform 282) plus the
18 genomes of §3.7; merged 300-genome indexes are built in a single
invocation.

**Arms.** sk2bGrow `index` at k16 (full panel) and k8 (top-8); an FMH arm of
Pilea-style reference sketches (`sketch_loci`, *k* = 31, armE hash semantics)
at the F2-matched scale, plus a **half-density reference arm** (s200, ≈0.53×
the k16-matched landmark count). The binary index does not yet accept FMH
anchors, so the FMH arm's costs are those of the current Python
implementation — a gap that is itself reported.

**Measurement.** Builds wrapped in `/usr/bin/time -v` (wall-clock, peak RSS);
disk by `du -sb`; profile cost on the first 2 M single-end reads of two real
WGS samples (PRJNA1280254), 8 threads, `--no-stats`, matching the F2
profile-cost convention.

**Half-density accuracy.** 3 genomes (50.8/65.3/72.0 % GC) × 4 arms (panel
k16/k8, FMH matched, FMH s200) × 4 depths (0.5–5×) × 8 truth replicates.

**Outputs.** Bytes per landmark; build wall-clock/RSS and on-disk size per
index; per-sample profile wall-clock/RSS; half-density accuracy against
matched density.

---

## 7. Computational environment

Apple M3 Max, 16 cores, 48 GB RAM, macOS 14.7; Rust 1.92.0, Python 3.12.4.
The initial benchmark grid was run on a laptop; the final C1 paired-end A/B/E statistics grid and F1–F5/C5-scale experiments used the internal HPC (SLURM; Intel and AMD partitions). The final primary grid used 8 CPU threads on the AMD partition under sk2bGrow review-final commit `929f4c2` (job 4076237). Runs deferred for scale —
full-depth *E. coli*, the 45,529-assembly quality sweep, Pilea's full
32-strain/32× grid, and the marine metagenome application — are marked as such
where they appear. HPC runs used a SLURM cluster with Lustre storage.

---

## 8. Data and code availability

sk2bGrow: <https://github.com/HuangShiLab/sk2bGrow>.
Manuscript, figures and figure code: <https://github.com/HuangShiLab/sk2bGrow-paper>.
Figures are a pure function of the tables in `data/`; `python3 figures/make_figures.py`
regenerates all of them with no network access and no recomputation from reads.
Review-response provenance is in `data/m1_signed/`: the final primary A/B/E grid
(`hpc_singlepass_grid_results.tsv`), its summary
(`hpc_singlepass_grid_summary.tsv`), and the rejected residual-GC diagnostic
(`gc_correction_diagnostic.tsv`).
Sequencing data: PRJNA615952 (Zheng *E. coli* panel), PRJNA689204 (Sun faecal
study), PRJNA974210 (rotating biological contactor metagenome).
