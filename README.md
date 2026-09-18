# sk2bGrow — paper

Figures and figure-producing code for the sk2bGrow manuscript. The tool itself
lives in the separate [`sk2bGrow`](../sk2bGrow) repository.

## Layout

```
*.docx        the two source design documents (algorithm report, architecture)
manuscript/   outline and `manuscript.md`, the authoritative manuscript
manuscript/submission/  formatted review draft for journal submission
data/         committed benchmark outputs (small TSVs) — the only input to figures & tables
figures/      make_figures.py, make_tables.py, style.py
figures/out/  generated PNG + PDF (regenerable; safe to delete)
tables/       generated markdown + TSV (regenerable)
```

`manuscript/manuscript.md` is authoritative. It follows the BMC/*Microbiome*
Research-article order (Background, Methods, Results, Discussion, Conclusions,
Declarations) and embeds the submission figures and tables. The current review
draft is
`manuscript/submission/sk2bGrow_Microbiome_submission_draft.docx`; the
cover-letter and submission states are recorded in `COVER_LETTER.md` and
`SUBMISSION_CHECKLIST.md`.

The two Chinese-language design documents are kept here rather than in the code
repository: they are the provenance of the method and belong with the
manuscript, not with the installable tool.

## Regenerating

```bash
python3 figures/make_figures.py
```

Reads only `data/*.tsv`. No network, no recomputation — the figures are a pure
function of the committed data, so a reviewer can reproduce every panel without
rerunning the benchmark or installing the tool. To regenerate the *data*, see
`benches/zheng2020/` in the code repo.

## Provenance

Every figure traces to one benchmark: **Zheng et al. 2020**
([BioProject PRJNA615952](https://www.ebi.ac.uk/ena/browser/view/PRJNA615952)),
E. coli K-12 MG1655 across 16 growth media (λ = 0.40–1.72 h⁻¹) plus a
stationary-phase control, against **Pilea v1.3.8** (bioconda) on byte-identical
inputs at 0.5/1/2/5/10×.

Record with any submitted version: the sk2bGrow commit hash, the Pilea version,
and the subsampling read counts. `data/results_raw.tsv` carries one row per
(arm, medium, coverage).

## Figures

| figure | claim |
|---|---|

| `fig1_overview` (Fig 1) | the method, on real output from one 2× sample |
| `fig2_zheng_benchmark` (Fig 2) | Zheng isolate benchmark, four panels: (a) sk2bGrow is at or above Pilea at every coverage — the margin is in the 1–2× band; (b, c) **both tools are biased at 1×, in opposite directions** — see below; (d) the sorted-regression estimator invents a gradient on a non-growing culture |
| `fig3_simulation` (Fig 3) | multi-strain communities: recall, accuracy and cost — **read the three panels together** |
| `fig4_attribution` (Fig 4) | the estimator × landmark-source interaction, at the Pilea operating point and density-matched |
| `fig5_panel_design` (Fig 5) | 4–8 enzymes suffice; the four sparsest add nothing (a–d); panel density rises with GC, instrument boundary at GC ≲ 30% × 0.5× (e, f) |
| `fig6_fragmentation` (Fig 6) | fragmentation destroys the coordinate; `scaffold` restores the complete-reference result |
| `fig7_metagenome` (Fig 7) | Sun cohort: cross-method concordance + the exact-dedup ablation |
| `fig8_mag_qc_cost` (Fig 8) | C5: QC pass rate falls with MAG fragmentation and recall under a common protocol (a, b); cost at scale, measured: 89.5–240.8× → mm=1 → containment screen (c, d) |


## Multi-strain simulation (Pilea's Fig-3 design, laptop scale)

16 reference genomes (including *E. coli* K-12, *E. coli* O157:H7 and *Shigella
dysenteriae* as a deliberate shared-anchor stress test), 4/8/16 strains per
sample at 1/2/4/8x, V-shaped profiles, log2PTR ~ U[0,2], 2 replicates.

| method | recall | RMSE | bias | sec | peak RSS |
|---|---:|---:|---:|---:|---:|

| sk2bGrow | **1.000** | **0.134** | **−0.013** | 16.4 | **187 MB** |

| Pilea (defaults) | 0.224 | 0.083 | −0.045 | **2.8** | 223 MB |
| Pilea (gates off) | 0.997 | 0.265 | +0.168 | 17.7 | 222 MB |

Neither tool reported a genome that was not in the sample (spurious = 0).

**Recall and RMSE have to be read together.** Pilea at its shipped defaults has
the best RMSE in this table, but it earns it by answering only 22% of cases —
the high-coverage ones. At matched recall (gates off) its RMSE is 0.265 against
sk2bGrow's 0.134. Quoting either RMSE alone misrepresents the comparison.

**Cost crosses over at ~4×.** sk2bGrow's anchor scan is linear in read count
(5.0 s at 1×, 9.4 s at 2×, 17.3 s at 4×, 33.7 s at 8×), while Pilea with gates
off is *non-monotonic* and most expensive exactly where this paper is aimed
(22.6 s at 1×, 24.7 s at 2×, 17.0 s at 4×, 6.6 s at 8×) because its ZTP-mixture
EM is slowest where the mixture is least identifiable. So sk2bGrow is 4.5×
cheaper at 1× and 2.6× at 2×, level at 4×, and 5.1× more expensive at 8×;
averaged over the grid, 16.4 s against 17.7 s. At its shipped defaults Pilea is
far cheaper than either (2.8 s) because below its gate it does no fitting at all
— and returns nothing. Memory is modestly lower throughout (187 vs 222 MB).

## Two things not to overclaim

**The headline correlation hides a bias.** Fig 2 is the honest panel. At 1×
sk2bGrow ties the density-matched sketch on ranking (r 0.912 vs 0.912) but sits
below the identity line: the slope against predicted log₂(PTR) is 0.58, against
Pilea gates-off 0.63. Better or equal ranking, not better magnitude. Neither
tool is unbiased at 1×; only by 10× does sk2bGrow track y = x closely
(slope 1.07 vs Pilea 0.92). Any abstract sentence built on r alone is
misleading.

**The landmark source is not the winner.** Fig 4 compares landmark source and
estimator. At matched density the two landmark sources are statistically
indistinguishable; the coordinate-aware V-fit is what separates the method from
sorted-rank regression. The multi-enzyme panel contributes wet-lab
realizability, heterogeneous QC strata and fusion redundancy, not a magical
deterministic sketch. The old “interaction” narrative has been withdrawn after
the multi-instance test and the signed-origin fix.
