# sk2bGrow — paper

Figures and figure-producing code for the sk2bGrow manuscript. The tool itself
lives in the separate [`sk2bGrow`](../sk2bGrow) repository.

## Layout

```
data/       committed benchmark outputs (small TSVs) — the figures' only input
figures/    make_figures.py, style.py
figures/out/ generated PNG + PDF (regenerable; safe to delete)
```

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
| `fig1_accuracy_vs_coverage` | sk2bGrow is at or above Pilea at every coverage; the margin is in the 1–2× band |
| `fig2_magnitude` | **both tools are biased at 1×, in opposite directions** — see below |
| `fig3_negative_control` | the sorted-regression estimator invents a gradient on a non-growing culture |
| `fig4_attribution` | the estimator, not the deterministic sketch, carries the gain |

## Two things not to overclaim

**The headline correlation hides a bias.** Fig 2 is the honest panel. At 1×
sk2bGrow correlates better than Pilea (0.954 vs 0.889) but sits systematically
*below* the identity line — it underestimates log₂(PTR) by roughly 0.3 — while
Pilea sits systematically *above* it. Better ranking, not better magnitude.
Neither tool is unbiased at 1×; only by 10× does sk2bGrow track y = x
(slope 1.06 vs Pilea 0.94). Any abstract sentence built on r alone is
misleading, and a reviewer will find this panel.

**The deterministic sketch is not what wins.** Fig 4: holding the sketch fixed
and swapping the estimator moves accuracy far more than the reverse, and under a
sorted-regression estimator the anchors are *behind* FracMinHash at 1×
(0.61 vs 0.89). The contribution is the coordinate-aware V-shape fit that the
anchors make possible — not the anchors as a sketch. The paper should claim
that, because it is what the data shows.

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
