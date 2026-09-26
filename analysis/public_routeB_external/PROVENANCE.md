# Run provenance

## Data

| Label | Dataset | Runs profiled | Raw FASTQ source | Labels |
|---|---|---:|---|---|
| Bladder | PRJNA946904 | 22 | ENA | NMIBC (7) versus MIBC (15) |
| Ovary | PRJNA1005427 | 20 | ENA | benign B-1–B-10 versus cancer C-1–C-10 |
| Maternal/meconium | PRJCA030517 / CRA019236 | 63 | NGDC GSA | maternal faeces MST (33) versus infant meconium MEC (30) |

The local raw-data mirror is `/Volumes/MoneyCat/Data/2bRAD_public`. HPC raw FASTQ intermediates were removed after profiling; per-sample sk2bGrow summaries remain on HPC and are copied here.

## Profiling

- Binary: sk2bGrow 0.1.0, HPC `src_m4` release build.
- Database: GTDB species-representative BcgI database (970 representatives).
- Route: route B, real 2bRAD; BcgI; mismatch 0.
- Statistics: no residual two-pass GC correction; GLM window-rate fitting.
- Batch jobs: maternal NGDC set job `4175494`; ENA set job `4175496`.
- Status: 63/63 maternal and 42/42 ENA samples completed; 105 output tables total.

## Analysis

- Outputs were aggregated with `scripts/public_routeB_summarize.py`.
- Robust exploratory PTR filter: finite log2 PTR 0–3, at least 5 windows and 50 anchors, origin confidence at least 0.5, PTR SE at most 1.
- Classification: leave-one-out CV; feature filtering/imputation inside each training fold; robust PTR missingness represented explicitly and median imputation used.
- Support arm: coverage/support fields inferred by the same route-B count layer (no independently published abundance table was available for these runs).
- Statistical uncertainty: 5,000-sample bootstrap confidence intervals for AUROC and 5,000-permutation paired score tests.
- Multiple testing: BH q values within each group contrast.
- Figure: `figures/public_routeB_story.py`; collision audit passed with 0 failures and 0 warnings.

## Caveats

Only one strict-QC estimate passed; all disease/body-site findings are exploratory. The bladder and ovarian cohorts are small, and no independently published abundance tables were available for these raw public runs.
