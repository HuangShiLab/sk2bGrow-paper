> **SUPERSEDED 部分（2026-09-13）**：Table 2 已于 HPC C1 双端实例再生成（arm A 0.5–10×：0.941/0.919/0.959/0.971/0.970），形式交互检验仅 2× 显著为负。本文档中的旧数字与 "interaction" 主张以 manuscript.md 为准。详见 data/repro_check/ 与 data/repro_check/multiseed/INTERACTION_REPORT.md。

# Story plan — reorganization for Microbiome submission (2026-09-08)

Scope: how to integrate the HPC results (C5, Sun/σ_eff, C8, F1, F2; F3–F5
pending) into the manuscript, what story it tells, and what must change in
existing figures/tables. Builds on `outline.md` (which stays the section-level
contract); this file is the delta plan.

Target: *Microbiome* (Pilea's venue) or above. The bar that matters:
**a real biological result**, not only benchmarks. Before this round the
outline's §8 Application was ❌ ("no application dataset yet"). That gap is now
fillable — this plan's center of gravity moves there.

---

## 1. What the story is now

**One-sentence claim.** A coordinate-aware estimator turns sparse, motif-defined
2bRAD anchors into PTR estimates at 1–2× read depth — an order of magnitude
below the coverage gate of the current sketch-based standard — and the two
independent methods agree on real fecal metagenomes.

**Why this is a story and not a benchmark suite.** The method claim (low-depth
accuracy) was already in hand. What was missing, and what the HPC round adds:

1. **Biological validation.** Sun fecal cohort: sk2bGrow-WGS vs Pilea on the
   same samples agree (bias +0.22–0.25 log2, n = 58–78 species×sample pairs).
   This is the first real-community result in the project and the strongest
   candidate for the paper's biological anchor.
2. **A stated instrument boundary.** F2: usable range GC ≳ 30%, depth ≳ 1×,
   in both landmark modes. Reviewers reward a measured boundary over an
   implied universal one.
3. **An honest landmark verdict.** F1: at matched density the two landmark
   sources are indistinguishable ≥1×; the panel's low-depth edge is fusion
   redundancy across ~15 heterogeneous strata, not landmark quality. The
   paper becomes *narrower and more defensible* — the contribution is the
   estimator + wet-lab realizability + strata redundancy, not "deterministic
   anchors are better".
4. **Process-level findings that others will cite.** Exact dedup before
   counting destroys the PTR signal (Sun §6–7); QC pass rate anti-correlates
   with fragmentation on real MAGs (C5); cost-at-scale flip and its mm=1 /
   containment-screen fixes (C5 §4e–4g).
5. **A novel measurement.** σ_eff ≈ 1.53 per-anchor capture-efficiency
   dispersion, mostly correctable within batch — no published measurement
   exists; it is route-B-specific and Pilea cannot produce it.

**Title/abstract consequence.** Working title (coordinate-aware estimation)
survives F1 — the coordinate fit remains the contribution. Abstract must now
also carry the concordance result and the GC boundary; correlation-only
language is doubly insufficient.

---

## 2. Proposed Results structure (section → figures/tables → data status)

Order rationale: method and isolate benchmark first (unchanged), landmark
question resolved early because everything downstream (panel size,
fragmentation, real data) references it, then the real-data payoff, then cost.

### §1 Overview — Fig 1 ✅ unchanged
(exemplar_M1_2x; note `growth_rates.tsv` pred is derived from the same data —
flag in Methods, not new)

### §2 Isolate accuracy — Fig 2–4 ✅ unchanged
Zheng grid, 0.5–10×. Numbers stand (Part 9 voided by C8 closeout; committed
arm A numbers are the valid ones). **Check before submission:** sim_results
aggregate RMSE reads 0.134 in `data/sim_results.tsv` vs 0.168 in README —
reconcile.

### §3 What do the landmarks contribute? — **Fig 6 rewritten** 🟡
Was: 2×2 attribution at Pilea's operating point (scale 250 ≠ panel density).
Now three layers, all density-stated:
- (a) 2×2 at unmatched density (arm E as shipped) — the interaction stands as
  measured **at Pilea's operating point**; claim explicitly scoped.
- (b) Density-matched (F1): FMH ≈ panel at ≥1× (r within noise), sketch ahead
  on survivors at 0.5× but 6/17 cells lost at the no-gradient gate.
- (c) Mechanism (F1_mechanism): per-window landmark counts and detection
  fractions near-identical → window population is density, not determinism;
  panel's low-depth robustness = multi-strata fusion redundancy (0/1-stratum
  sketch cannot be rescued; 16-strata panel can).
- **Retire the sentence "deterministic anchors are what keeps windows
  populated at 0.5×"** everywhere (Fig 6 caption, old §4 text).
- New headline sentence for §3: *at matched density the estimator is
  landmark-agnostic ≥1×; below 1× the multi-enzyme fusion, not the landmark
  set, is what keeps the instrument working.*

### §4 Panel size and GC robustness — Fig 7 ✅ + **new Fig 8 GC** 🟡
- Panel sweep (existing Table 7 / Fig 7): 4–8 enzymes.
- **[NEW] GC sweep (F2):** ~18 genomes, GC 25–72%, density + planted-gradient
  accuracy, both modes. Messages: (i) k16 density rises monotonically with GC;
  low-GC depression only; (ii) k8 ≈ k16 accuracy across GC → **A6 holds
  outside E. coli**; (iii) instrument boundary GC ≲ 30% + 0.5× (both modes;
  Buchnera confounded with genome size — state it); (iv) FMH flat and tied.
- This section also absorbs the density table: **Table 1 extended to the
  F2 genome set** (density_3genomes.tsv → F2_density.tsv; keep 3-genome
  enzyme-vs-report audit as supplement or first block).

### §5 Fragmented references — Fig 9 ✅ caption fixed 🟡
Substantive numbers unchanged. Caption edits per F1: the collapse is about the
coordinate, not the landmark set; QC blind spot stays. **F4 landed
(2026-09-08, data/f4_fragment_sketch/):**
- V-fit collapse is landmark-agnostic — enzyme and sketch arms give identical
  slope collapse (≈1.0 → 0.26–0.38 at 5×) and identical RMSE blow-up (0.02 →
  0.84–0.87).
- **Winner's-curse refinement to the Fig 9 story (important).** F4's
  permutation test (200 shuffled-coordinate refits) shows the residual
  r ≈ 1.0 / r = 0.86-on-contigs is partly a *search artifact*: the ori grid
  search fits a fake V to a flat profile, and shuffled profiles produce the
  same proportional estimates. The honest metrics on fragmented references
  are slope and RMSE against truth; r overstates survival. This strengthens —
  and sharpens — the existing "report slope beside correlation" argument: the
  r column in the old Fig 9 was itself inflated by the same artifact.
- **scaffold is landmark-source-agnostic.** Self-arm: sketch ties enzyme
  exactly (E. coli 98/100, Bacillus/Pseudomonas 100/100 placed, Spearman
  1.0000, 5× RMSE 0.017–0.033). Rel-arm: 0.1% divergence relative both
  100/100 with full recovery; divergent-relative sketch is slightly *better*
  than enzyme (RMSE 0.223 vs 0.644, Bacillus, n=6 directional only).
  Consequence: "scaffolding needs the enzyme panel" is dropped as a route-A
  argument; scaffolding stays a route-B (enzyme-only data) capability.
- Mechanism sentence for §5: fragmentation = coordinate-fit death + search
  artifact residue, not noise addition.
- Caveat carried from F4 REVIEW: shipped `scaffold` digests enzyme tags only;
  the FMH placement was a 1:1 Python port of scaffold.rs (validated by the
  perfect self-arm), so sketch scaffolding is demonstrated, not shipped.

### §6 Real metagenomes **[NEW — the biological section]** ❌→🟡
Two datasets, two distinct questions.

- **§6a Fecal cohort cross-method concordance (Sun).** WGS only: sk2bGrow vs
  Pilea-defaults, common denominator, BA + CCC with observed ranges printed
  beside every coefficient (dynamic range is narrow; CCC deflation expected
  and explained). Bias +0.22–0.25 log2 (Pilea higher), LoA ±0.4–0.9.
  Gates-off arm reported as the control showing Pilea's gate does real work
  (r ≈ 0). **This is the paper's biological anchor** — frame as concordance
  validation in a real community, n = 3 samples, species×sample pairs
  58/68/78, not as a biological discovery.
- **§6b Real 2bRAD arm and the dedup result.** Same cohort, in-silico vs
  real 2bRAD library (same method, two preparations). Headline: exact
  deduplication before counting flattens the PTR dynamic range (raw r 0.73–0.80
  vs deduped 0.30–0.60); cap-dedup sensitivity shows no robust intermediate
  (cap 2 best in one sample, worst in another) → **process recommendation:
  do not exactly dedup 2bRAD libraries before counting**. Include σ_eff here
  as the route-B characterization: σ_eff ≈ 1.53, ICC-correctable share
  ~0.75 (within-batch), window averaging dilutes the residual to
  ~0.07–0.10 on the gradient scale (ICC-based lower bound ~0.07,
  structure-test upper ~0.10) — i.e. measured, correctable, negligible
  after correction. Single-enzyme (BcgI) caveat stated up front.
- **§6c MAG-scale application and QC validity (C5, RBC reactor).** Three
  messages: (i) recall under a common protocol — QC-pass 3.8–13.8% vs Pilea
  5.0–11.1%, overlapping (the 1.00 headline is a denominator artefact and is
  never printed unqualified); (ii) **QC works on real data** — pass rate
  anti-correlates with MAG fragmentation (ρ = −0.41; −0.34 controlling
  coverage); (iii) cost at scale flips (89.5–240.8× → **~10–25× after mm=1**,
  cite the revised number only) and the M4 containment screen clears the
  absent-genome false-positive floor (4,048/4,700 synthetic genomes hit in
  the unrestricted index, 0 with screen; at ~90% presence no speedup measured (0.99× wall) and 65/522 MAGs dropped; value projected for < 1% presence).
- Outline's old ❌ "Application — we have no dataset" is deleted; §6 replaces
  it. Pilea's RBC parallel gives us a reviewer-ready comparison frame.

### §7 Cost and scale — Table (existing §7) + C6 outlook 🟡
Existing per-sample costs stand. **F3 landed (2026-09-09, data/f3_mismatch/):**
- Default **mm=2 stands for the enzyme mode** — misassignment is a non-issue
  there (≤1.4e-4 pooled, worst single genome 6.1e-4; exact-first suppression +
  motif gate make tolerance nearly free). If GTDB memory forces the cut,
  **mm=1 is safe**: misassignment identical to mm2 at realistic error rates,
  sensitivity loss negligible (0.44%→0.42% unmatched). The mm decision is
  purely budget-vs-sensitivity; the misassignment risk is removed from it.
- **Sketch mode must lock mm=0.** The naive "FMH has ~66× fewer near
  neighbors ⇒ it can afford higher tolerance" is **refuted**: at scale ~100,
  ~99% of read k-mers are off-sketch, and mm≥1 opens a rescue channel that
  credits ~1% (mm1) / ~2% (mm2) of those to a neighboring sketch key at the
  wrong coordinate — independent of sequencing error. Landmark–landmark
  census cannot see this channel; only the empirical read-level test can.
  (Caveat carried: FMH empirical arm is a rule-based simulation — the sketch
  counting mode is not in the CLI — declared in F3 REVIEW §0.)
- These feed the Discussion landmark table (mismatch-tolerance row) and the
  C6/mm=1 sentence in §6c.
- **F5 landed (2026-09-09, data/f5_index_cost/) — the C6 answer:**
  (i) **Index cost is a wash between modes at matched density** — ~39 B/anchor
  either way, equal build speed; mode choice is accuracy/wet-lab (F1/F2:
  tied), not cost. (ii) **A1's 752 GB was peak RSS; the disk figure is ~¼ of
  it.** GTDB species-rep projection: enzyme k16 ~232 GB, k8 ~184 GB,
  **FMH half-density (scale 200) ~124 GB**. (iii) **Half density is the
  practical C6 answer**: s200 carries 0.53× the landmarks and loses ≤0.007 r
  at ≥1× (RMSE sometimes better); the only cost appears at 0.5× and grows
  with GC (−0.03 r at 72% GC). (iv) Engineering gaps stated honestly:
  `sketch index --mode fracminhash` (merged FMH index) does not exist — F5's
  FMH arm is the per-genome Python sketch pipeline, which rescans reads per
  genome and is therefore *slower end-to-end* than the merged enzyme index
  (8 s/genome vs one 530 s pass) despite equal per-landmark storage.
- §7 C6 outlook paragraph: enzyme k16 with the M4 containment screen, or
  sketch half-density at ~124 GB — both viable; screen and mm=1 from §6c/§7
  apply to the enzyme route.

### Discussion — existing content + four updates
- Landmark positioning table (already drafted in outline Discussion): fill F1
  (matched-density row), F2 (GC row), F4 (scaffold row, pending), F5
  (density-control row: sketch has a knob; panel composition is
  biology-aware — the one remaining qualitative panel advantage).
- Loud vs silent failure paragraph: keep, now strengthened by §6c numbers.
- **New paragraph: dedup recommendation** (process finding, cite-able).
- **New paragraph: instrument boundary** (GC ≳ 30%, ≥1×; per-clade panel
  tuning as future work).
- Keep QC blind spot + Bsp24I⊂CjePI + annotated-origin paragraphs unchanged.

---

## 3. Figure/table delta list

**Must fix (conflicts with new data):**
1. `fig5_attribution` — caption/text rewrite per §3 above; add density-matched
   panels b–c from `f1_sketch/F1_results.tsv` + `F1_mechanism.tsv`.
2. `fig8_fragmentation` — caption sentence "anchors keep windows populated"
   removed; coordinate-loss framing only.
3. README figure table — numbering is stale (fig5_simulation etc.), refresh
   after renumbering.

**New figures:**
4. `fig6_gc_sweep` (F2): density-vs-GC curves (panel k∈{2,8,16}, FMH matched)
   + accuracy-vs-GC heatmap at 1× (and 0.5× showing the dead zone). Data:
   `f2_gc_sweep/F2_density.tsv`, `F2_accuracy.tsv`.
5. `fig9_metagenome` (Sun): BA scatter per sample, sk2bGrow vs Pilea, with
   ranges annotated; inset or second panel: deduped-vs-raw C-arm compression.
   Data: `sun_three_arm/p2_review_agreement.tsv`, REVIEW §6.
6. `fig10_mag_qc` (C5): QC pass rate vs contig count (binned, coverage-controlled)
   + recall three-ways bars. Data: `c5_review/c5_qc_x_ncontigs_binned.tsv`,
   `c5_recall_three_ways.tsv`.
7. (Supplement) σ_eff: per-anchor e distribution + ICC decomposition. Data:
   `sun_three_arm/sigma_icc.tsv`, `sigma_eff_e_prior_sun.tsv`.
8. (Supplement) cost-at-scale waterfall: 89.5–240.8× → mm=1 → screen. Data:
   `c5_review/c5_cost_per_sample.tsv`, `c5_count_breakdown.tsv`.

**New/updated tables:**
9. T1 → extended density table (F2 18 genomes; keep 3-genome audit block).
10. New T-landmark (Discussion): the property×mode comparison table from
    outline.md, with F1/F2/F4/F5 rows filled and 未测 rows marked.
11. New T-metagenome: Sun concordance per sample (n, r, CCC, bias, SD, LoA,
    ranges) + C5 recall three ways + cost revision.

**Before-submission data hygiene (from audit):**
- Reconcile sim RMSE 0.134 (tsv) vs 0.168 (README).
- `growth_rates.tsv` pred is derived, not independent — Methods must say so.
- Dedup `sun_three_arm/REVIEW.md` (§6/§7 pasted 3×); same for any
  double-appended sections in `c5_review/REVIEW.md`.
- Never quote: C5 suspicious=0 (structurally vacuous), Sun [D] bias–coverage
  correlation (outlier-driven, refuted by winsorization), §4d's 65–277× (revised to 89.5–240.8× per tsv)
  (superseded by mm=1 revision), F1 cross-harness comparisons.

---

## 4. Sequencing and ownership

- Now (data in hand): §3 rewrite inputs, §4 GC figure, §6a/§6c figures,
  table updates, caption fixes — all local, no HPC dependency.
- When agents return: nothing left — F3/F4/F5 all landed (see §7, §5).
  F7 design note: F3's mm=0 lock for sketch strata means a
  hybrid landmark set runs enzyme strata at mm=2 and the sketch stratum at
  mm=0 — technically straightforward, one more reason F7 is viable.
- Hou data (external batch σ_eff validation) remains the one dataset-level
  gap; if it lands before submission it strengthens §6b's within-batch
  qualifier to a cross-batch statement.
