# sk2bGrow — paper

Figures and figure-producing code for the sk2bGrow manuscript. The tool itself
lives in the separate [`sk2bGrow`](../sk2bGrow) repository.

## Layout

```
*.docx        the two source design documents (algorithm report, architecture)
manuscript/   outline; the manuscript itself lives here
data/         committed benchmark outputs (small TSVs) — the only input to figures & tables
figures/      make_figures.py, make_tables.py, style.py
figures/out/  generated PNG + PDF (regenerable; safe to delete)
tables/       generated markdown + TSV (regenerable)
```

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
| `fig1_overview` | the method working on one real 2× sample: fitted V, rank-regression fragility, per-enzyme forest |
| `fig2_accuracy_vs_coverage` | sk2bGrow is at or above Pilea at every coverage; the margin is in the 1–2× band |
| `fig3_magnitude` | **both tools are biased at 1×, in opposite directions** — see below |
| `fig4_negative_control` | the sorted-regression estimator invents a gradient on a non-growing culture |
| `fig5_attribution` | sketch and estimator interact; neither carries the result alone — see below |
| `fig6_simulation` | multi-strain communities: recall, accuracy and cost — **read the three panels together** |
| `fig7_panel_size` | 4–8 enzymes, not 16 |
| `fig8_fragmentation` | fragmentation removes the gradient; `scaffold` restores it |

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

**The headline correlation hides a bias.** Fig 3 is the honest panel. At 1×
sk2bGrow correlates better than Pilea (0.981 vs 0.889) but sits systematically
*below* the identity line — mean −0.13 log₂ units, growing to −0.23 at the top
of the range because the slope is 0.86, not 1 — while Pilea sits systematically
*above* it (mean +0.35, slope 0.63). Better ranking, not better magnitude. Neither tool is
unbiased at 1×; by 10× sk2bGrow slightly overshoots (slope 1.07 vs Pilea 0.94).
Any abstract sentence built on r alone is misleading, and a reviewer will find
this panel.

**Neither the sketch nor the estimator carries the result alone.** Fig 5 is the
full 2 × 2 — both sketches × both estimators on the same subsampled reads. At 1×
the coordinate fit is worth +0.30 r on anchors but only +0.05 on a FracMinHash
sketch, and the *sketch* effect changes sign with the estimator: anchors are
+0.04 ahead under the V-fit (0.981 vs 0.940) and 0.21 behind under rank
regression (0.683 vs 0.889). The two factors interact, and at 0.5× only the
combination survives at all (0.913, against 0.724 for the V-fit on a
FracMinHash sketch and a fully degenerate Pilea arm). Deterministic anchors keep
the windows populated at that depth; the coordinate fit turns them into an
unbiased slope. The paper must claim the **combination**, not either factor.

> An earlier three-arm reading of this figure concluded "the estimator, not the
> sketch, carries the gain". That was an artefact of the missing fourth cell
> (FracMinHash + V-fit) and has been withdrawn.

## Ground truth

`data/growth_rates.tsv` from the paper's supplementary source data. Two targets,
**not equally independent**:

- `growth_rate` (λ, h⁻¹) — measured by OD/microscopy. Independent. This is the
  real test and what Pilea's published r = 0.9764 refers to.
- `pred_log2ptr` = λC/ln2 — theory-predicted log₂(PTR). **Not independent**: the
  authors derived C from this same sequencing by marker-frequency analysis.
  Useful for checking magnitude; circular if quoted as accuracy.

## Scope

One organism, one strain, a complete single-contig reference, no community.
This is the easiest possible case and says nothing yet about the metagenomic
setting, which is where PTR estimation is actually hard.

## Style

`figures/style.py` fixes the categorical hues per *arm*, assigned in fixed order
and never cycled, so a figure that drops an arm does not repaint the survivors.
The palette is validated for colour-vision deficiency (worst adjacent pair
ΔE 9.1 protan, 22.9 normal vision).
# sk2bGrow-paper
