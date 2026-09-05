# Framing: what we can claim, against what evidence

Working note for the Introduction and Discussion. Every row carries its source.
Rows marked **未测** are hypotheses, not results — they must not be written as
claims until an experiment exists.

Target: *Microbiome*-level. §6 says plainly what is still missing for that.

---

## 1. The thesis, in one sentence

**Order-free rank regression requires landmarks that sample genomic position
uniformly; a coordinate-aware estimator does not — and that is precisely what
makes a motif-constrained, wet-lab-realisable landmark set usable at all.**

The whole 2×2 falls out of it. Pearson r against measured growth rate, Zheng
*E. coli*, 1× subsampling, n = 16:

| | coordinate V-fit | sorted-rank regression |
|---|---|---|
| **2bRAD anchors** | **0.981** | 0.683 |
| **FracMinHash** | 0.940 | 0.889 (Pilea) |

Read the columns, not the rows. Rank regression works on FracMinHash (0.889) and
collapses on anchors (0.683), because its validity rests on the tent function's
*values* being uniformly distributed — which follows from uniform *position*
sampling. A hash threshold samples positions uniformly by construction; Type IIB
sites cluster on their recognition motifs. Supplying the coordinate removes that
requirement, and only then are enzyme anchors competitive.

This framing survives whichever way the pending density-matched experiment (F1)
goes, because it is a statement about estimators, not about landmark quality.

**Do not** frame the paper as "our landmarks beat FracMinHash." Arm E ran at
scale 250 = 3,934 landmarks/Mb against our panel's 9,422/Mb — a 2.4× density
disadvantage — so any landmark-quality claim can be overturned by one run.

---

## 2. The comparison table

Half a figure's worth. Source in every row; **未测** where we have no measurement.

| property | FracMinHash | Type IIB anchors | source |
|---|---|---|---|
| selection rule | `h(canonical kmer) < MAX/scale` | sequence matches a recognition site | both deterministic |
| context-free | yes | yes | Syn2b §1 |
| density control | continuous (`scale`) | 4 steps: 1/2/4/16 enzymes | Syn2b §3.2 |
| density, *E. coli* | 3,934 /Mb at s=250 | 9,422 /Mb at 16 enzymes | our index |
| GC dependence, density max/min | **1.04×** | 1.16× (16 enz), 2.9× (4), 9.7× (BcgI) | Syn2b §3.1; `data/density_3genomes.tsv` |
| GC compensation vs panel size | n/a | 1.51× (k=2) → 1.16× (k=16) | `data/density_3genomes.tsv`, 3 genomes GC 43.5–61.5% |
| position uniformity | uniform by construction | clustered on motifs | Syn2b §3.4: max breakpoint error 1,031 vs **3,671** bp at matched density |
| retention at 0.1% divergence | **94.6%** | 89.5% | Syn2b §3.5 |
| near-duplicate risk (1 sub from a multi-copy family) | **0.00%** | 0.34% (4 enz) | Syn2b §3.3 |
| genuine multi-copy families carried | 17 | 13 (BcgI), 38 (4 enz) | Syn2b §3.3 — a genome property, not a selection property |
| SNP failure geometry | smooth: one SNP kills ~26% of a 150 bp read's 31-mers | blocky: one SNP kills one locus, neighbours intact | arithmetic (120 31-mers per read, 31 destroyed) |
| independent strata | one pool; splittable by hash, but every split shares the same systematic bias | ~15 enzymes with **different** motif, GC and methylation biases | `fusion.py`; Bsp24I ⊂ CjePI is 100% containment (`docs/enzymes.md`), so 16 strata are not 16 |
| **wet-lab realisable** | **no** — needs WGS plus an assembly | **yes** — it is a sequencing protocol | categorical |
| failure when the reference is wrong | **silent** (containment falls, no output) | **loud** (an estimate is produced at a wrong coordinate) | our C5: recall 522/522 with no coverage gate anywhere |
| PTR at 1×, coordinate V-fit | 0.940 | **0.981** | `benches/zheng2020/RESULTS.txt` — at a 2.4× density disadvantage |
| PTR at 1×, rank regression | **0.889** | 0.683 | same |
| PTR at 0.5×, coordinate V-fit | 0.724 | **0.913** | same; density caveat applies with more force |
| behaviour on unscaffolded contigs | **0.827** | 0.550 | `data/fragmentation.tsv`, 1× |

Two rows decide the paper: **wet-lab realisable**, and the two PTR rows read as a
pair.

---

## 3. What we can claim

**C1 — A coordinate-aware estimator makes non-uniform landmark sets usable.**
The 2×2 above. The interaction is the result, not either main effect.

**C2 — Cross-enzyme fusion is stratified replication with heterogeneous bias.**
Sixteen enzymes carry *different* motif, GC and methylation biases, so agreement
across them tests sampling noise **and** systematic-bias heterogeneity. Splitting
one FracMinHash sketch into k strata tests only sampling noise, because every
split shares the same selection bias. State the limit honestly in the same
breath: above 1× per-species depth, BcgI alone is within **0.011–0.029 r** of the
full panel (`data/per_enzyme_zheng.tsv`), so the panel's contribution is
concentrated at 0.5× and in the QC.

**C3 — 2bRAD is a data-generation strategy, not a computational filter.**
FracMinHash needs an assembly and WGS reads. 2bRAD delivers landmarks from a
mixed sample directly, at a fraction of WGS library cost. **PTR is a dynamics
measurement**, and for dynamics more timepoints at modest depth beat fewer at
high depth. Framing 2bRAD-for-PTR as a *temporal-resolution* argument is, as far
as we can tell, new.

**C4 — Loud failure versus silent failure, and why loud is not simply worse.**
When the reference is wrong, containment-gated methods go quiet; anchor-based
methods emit an estimate at a wrong coordinate. Our own C5 is the demonstration:
recall 522/522 with no coverage gate anywhere in `fit` or `fuse`. High yield is
not a virtue on its own — it includes the wrong answers. This must be written by
us, not discovered by a reviewer.

**C5 — The QC has an enumerable blind spot.** Cochran's Q asks whether the
enzymes agree, and a destroyed coordinate makes all sixteen agree there is no
gradient: 100% of fragmented estimates pass against 75% of correct ones
(`data/fragmentation.tsv`). Because the coordinate is explicit, the failure is
*nameable, detectable and fixable* (`scaffold`). A containment-dilution failure
has no coordinate to inspect.

---

## 4. What we must not claim

- **Not** that our landmarks are more species-specific. Specificity comes from
  the reference set and the uniqueness filter, not the selection rule, and our
  near-duplicate exposure is *higher* (0.34% vs 0.00%) because the motif
  constraint compresses landmarks into a small subspace — and we run at
  `--max-mismatch 2`, above the distance that measurement used.
- **Not** that our landmarks span the genome better. At matched density they
  span it *worse*: max breakpoint localisation error 3,671 bp vs 1,031 bp,
  because adding enzymes cannot fill a gap that has no sites.
- **Not** that the anchor set is uniquely enumerable a priori. A hash threshold
  is equally deterministic and equally enumerable.
- **Not** that the same loci recur across samples in a way FracMinHash cannot
  match. A FracMinHash sketch of a fixed reference is deterministic and constant.
- **Not** that a single sketch cannot support a within-sample consistency test.
  It can be split by hash prefix. The real distinction is bias heterogeneity (C2).
- **Not** that the replication gradient is steepest near *ori*. Under the plain V
  the slope is constant on each replichore; only multi-fork replication (PTR > 2)
  introduces a kink, and BIC decides whether it is used.
- **Not** the low-coverage advantage, until F1 reports. Arm E's 0.5× deficit was
  measured at 2.4× lower landmark density.

---

## 5. The estimator problem — status and solution

Two estimators, neither dominant:

| | complete reference | unscaffolded contigs | stationary control |
|---|---|---|---|
| coordinate V-fit | 0.981 at 1× | slope 0.21, bias −0.81 | 0.077 at 1× |
| order-free spread MLE | — | RMSE 0.87 → 0.14 at 5× | **0.000** |
| Pilea rank regression | 0.889 at 1× | **0.827** at 1× | 1.153 |

The order-free mode is the fragmentation-proof half and it has two defects. They
have different causes and different fixes.

### Defect A — no random-effects escalation. Fixed by porting 15 lines.

`benches/fragmentation/spread_estimator.py::fuse` computes Cochran's Q and then
never uses it. It is pure inverse-variance weighting, so the reported standard
error assumes the only scatter is sampling noise:

| depth | mean Q (dof 15) | reported SE |
|---:|---:|---:|
| 0.5× | 1,352 | 0.006 |
| 1× | 6,384 | 0.005 |
| 2× | 2,632 | 0.019 |
| 10× | 177 | 0.013 |

Q = 37 is already p < 0.002. These SEs are absurdly tight around values the
enzymes violently disagree about — confidently wrong, the worst failure mode a
QC can have. **The production `fusion.py` already implements the fix**
(DerSimonian–Laird escalation when Q rejects, `s_eff² = s_w² + τ²`). This is not
research; it is porting. Do it before any further spread-mode benchmarking.

### Defect B — collapse to exactly zero at low coverage. Probably upstream.

At 0.5× the spread MLE returns exactly 0.0000 in **15 of 17** conditions, and 14
of 17 at 1×. Random effects will not fix this: DL inflates the interval, it does
not move a point estimate off a boundary.

The V-fit compresses at the same depths — slope 0.615 at 0.5×, 0.779 at 1×,
rising to 0.951 at 10× — but degrades *gracefully* rather than collapsing. Both
estimators lose signal at the same depths and in the same direction, which points
at a cause upstream of both: the window rates themselves, not the fitting.

That is exactly open question **R4/A4** — whether low-count windows are pulled
toward the mean before any fit sees them — and its decisive experiment is already
specified: simulate with known λ per window and compare recovered rates to truth
at 0.5–1×.

**So A4 and R5 are one defect seen through two estimators.** One cheap simulation
settles both. If window rates are being shrunk, no change to either fitting
routine can help, and the fix belongs in `ztp.py`.

And the asymmetry is itself evidence for the thesis in §1: the V-fit survives
upstream compression because coordinates tell it which windows *should* be high;
the order-free estimator sees only the value distribution, so when that
distribution is compressed there is nothing left to recover. **Coordinates are
redundancy, and redundancy is what buys robustness.**

---

## 6. What is still missing for *Microbiome*

Stated plainly, because the gap is real.

1. **No biological result yet.** Everything above is method validation. A
   *Microbiome* paper needs a finding about a community, not only a benchmark
   against a competitor. The paired 2bRAD/WGS datasets (Sun 2022, Hou 2025) are
   the most likely route: growth rates across a real gut community, from a data
   type nobody has used for PTR.
2. **σ_eff is unmeasured.** Roadmap R1 calls it the largest open risk, and every
   route-B claim rests on an assumed 0.3–0.6. It is measurable from data already
   in hand — see `benches/HPC_TASKS_PAIRED_2BRAD.md` task 1.
3. **C5 has no ground truth**, so recall is its only metric, and recall alone
   cannot be falsified — returning noise scores 1.00. See `benches/C5_REVIEW.txt`.
4. **The cost reversal is unexplained.** 16–27 h per sample against Pilea's
   6–15 min, having been *faster* on isolates (8.2 s vs 11.5 s). The M4 two-tier
   screen exists to address this; it needs measuring at community scale.
5. **The multi-enzyme architecture is unvalidated on real community data.** Both
   paired datasets are BcgI-only, so cross-enzyme fusion and the Q-based QC —
   C2 above — cannot be tested on them. That still needs the roadmap's M3
   experiment.
