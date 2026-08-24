# Methods

Section order follows Pilea (Microbiome 2026) so the two papers can be read
side by side. Every number quoted here is reproduced by the scripts in
[sk2bGrow/benches](https://github.com/HuangShiLab/sk2bGrow) from the tables in
`data/`.

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
**BIC**, and the segmented form is only offered when ≥30 windows are available.

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
the moment estimator of the between-enzyme variance. **Fig. 1d** shows this
firing on a real 2× sample (*I*² = 0.71).

A sample passes QC on estimated coverage, dispersion, the fraction of reference
windows covered, and EM containment; failures are reported with a reason rather
than dropped silently.

**Known gap.** The fusion currently treats all 16 enzymes as independent strata.
They are not: Bsp24I's site set is *totally contained* in CjePI's (1,636/1,636,
891/891 and 2,910/2,910 tags on three genomes — the two patterns, written in
opposite strand orientations, differ at exactly one position), and about half of
Bsp24I's tags are also CjeI sites. The panel therefore offers at most ~15
independent strata, and *Q* is mildly anti-conservative. The containment
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
* the *predicted* log₂PTR = λ·C/ln 2 from the paper's own measured C — used for
  RMSE and slope, because a correlation cannot detect a compressed dynamic
  range.

**Run-out sample.** Pilea excluded the run-out samples. We keep one
(`E.coli_gDNA_RUN_OUT`, sample alias CON-C) as a **negative control**. In a
replication run-out, ongoing rounds are allowed to finish while new initiation
is blocked, so every chromosome in the population ends fully replicated and
log₂PTR must be ≈ 0. This tests a failure mode that correlation against a growth
panel cannot see, and it costs nothing: a method that reports a large PTR here is
reading noise as a gradient.

**Coverage titration.** The first 600,000 R1 reads per run were downloaded
(≈19× of the 4.64 Mb genome); SRA preserves flowcell order, which is random with
respect to genome position. Each sample was then truncated to nominal
**0.5×, 1×, 2×, 5× and 10×** (n = ⌈cov·L/150⌉ reads). One medium (M13) is
missing at 10× — its download finished after the subsampling pass had already
measured the file — so n = 15 at 10× and 16 at every other depth.

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
laptop-scale grid does not establish behaviour at 32 strains or 32×.

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
threshold of 5 on *k-mer* coverage corresponds to ≈6.6× *read* coverage; at 5×
nominal the measured median is 3.93 and the gate fires.

The gate is not arbitrary conservatism. At 0.5× nominal, Pilea's gates-off
output is fully degenerate — measured coverage exactly 1.000 and log₂PTR exactly
0.000 for all 17 samples — so the gate is correctly refusing to report. Between
1× and 5×, however, the gates-off estimator is informative (r = 0.89, 0.95, 0.95
against λ) and is nevertheless suppressed.

**Attribution arm.** To separate "deterministic anchors" from "coordinate-aware
estimator", a third arm runs sk2bGrow's own anchors through Pilea's windowing
and estimator (25 kb fixed-bp windows, sorted-rank regression with RANSAC). Any
difference between this arm and full sk2bGrow is attributable to the estimator;
any difference between this arm and Pilea is attributable to the sketch.

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

**Computational cost.** Wall-clock and peak resident set size are taken from
`/usr/bin/time -l`. On macOS this propagates a grandchild's `ru_maxrss`
(verified explicitly), so the figure includes the Python statistics layer that
the Rust binary spawns. Both tools were given 8 threads. Cost is reported per
sample, with the one-off index/database build reported separately.

---

## 6. Computational environment

Apple M3 Max, 16 cores, 48 GB RAM, macOS 14.7; Rust 1.92.0, Python 3.12.4.
All benchmarks in this paper run on a single laptop. Runs deferred for scale —
full-depth *E. coli*, the 45,529-assembly quality sweep, Pilea's full
32-strain/32× grid, and the marine metagenome application — are marked as such
where they appear.

---

## 7. Data and code availability

sk2bGrow: <https://github.com/HuangShiLab/sk2bGrow>.
Manuscript, figures and figure code: <https://github.com/HuangShiLab/sk2bGrow-paper>.
Figures are a pure function of the tables in `data/`; `python3 figures/make_figures.py`
regenerates all of them with no network access and no recomputation from reads.
