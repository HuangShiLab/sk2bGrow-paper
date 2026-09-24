# Lim ORPI denture route-B benchmark

## Dataset

- 97 denture samples with clean BcgI 2bRAD-M FASTQs in
  `/Volumes/MoneyCat/Data/Lim_ORPI/clean_data`.
- Labels: 55 `Clean_Denture` and 42 `Unclean_Denture`.
- Vendor species abundance table: `external/Abundance_Stat.all.xls`
  (1,601 species rows).
- Two `Control*.BcgI.fq.gz` files were present but were excluded because they
  are absent from `meta.txt`.

## sk2bGrow route B

All 97 samples were run on HKU HPC using the GTDB species-representative BcgI
database used for the route-B benchmark:

```bash
sk2bgrow profile <sample>.BcgI.fq.gz \
  --db db_gtdb_bcgI_sp --mode 2brad --max-mismatch 0 \
  --threads 8 --no-stats
python -m sk2bgrow.cli profile ... --no-gc-correct --method glm
```

The large count/window intermediates were deleted after statistics. Per-sample
finite PTR yields are in
`analysis/lim_orpi_denture_routeB/per_sample_feature_yield.tsv`.

The route-A-derived QC gate passed for 0/970 genome-sample estimates. Finite
PTRs therefore enter the ML benchmark as exploratory features rather than
high-confidence biological PTRs. Missingness is represented explicitly.

## ML benchmark

`scripts/denture_routeB_ml_benchmark.py` uses repeated stratified 10-fold CV
(100 repeats) with 5,000 Random Forest trees per fold. The positive class is
`Unclean_Denture`; feature filtering, transformations, imputation, and scaling
are fit inside each training fold. The three arms are:

* A: vendor species abundance.
* B: sk2bGrow log2 PTR, anchor support, coverage, and missing indicators.
* C: A plus sk2bGrow log2 PTR.

Results:

```text
analysis/lim_orpi_denture_routeB/cv_metrics_summary.tsv
analysis/lim_orpi_denture_routeB/paired_auroc_deltas.tsv
```
