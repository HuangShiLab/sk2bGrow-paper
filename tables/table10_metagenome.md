**Table 10. Metagenome results: Sun fecal-cohort concordance, the real-2bRAD dedup finding, and the C5 RBC application**

**(a) Sun fecal cohort (PRJNA689204): cross-method agreement on the WGS arms**

| comparison      | sample   |   n |   pearson_r |   ccc |   bias |   loa_lo |   loa_hi | sk2bgrow_wgs_range   | comparator_range   |
|:----------------|:---------|----:|------------:|------:|-------:|---------:|---------:|:---------------------|:-------------------|
| A_default_vs_B  | S01      |  58 |       0.431 | 0.313 |  0.247 |   -0.418 |    0.913 | 0.04..1.45           | 0.32..1.74         |
| A_default_vs_B  | S06      |  68 |       0.344 | 0.253 |  0.222 |   -0.445 |    0.888 | 0.01..1.45           | 0.28..1.23         |
| A_default_vs_B  | S07      |  78 |       0.505 | 0.368 |  0.225 |   -0.433 |    0.882 | 0.03..1.49           | 0.30..1.36         |
| A_gatesoff_vs_B | S01      | 381 |       0.095 | 0.070 | -0.932 |   -6.363 |    4.498 | 0.00..15.86          | 0.00..13.70        |
| A_gatesoff_vs_B | S06      | 459 |       0.005 | 0.003 | -1.153 |   -7.812 |    5.505 | 0.01..35.46          | 0.00..8.00         |
| A_gatesoff_vs_B | S07      | 403 |       0.046 | 0.026 | -1.238 |   -7.561 |    5.087 | 0.00..27.21          | 0.00..10.85        |

**(b) Same cohort, in-silico vs real 2bRAD library (same method, two preparations)**

| arm                 | sample   |   n |   pearson_r |   ccc |   bias |   loa_lo |   loa_hi | sk2bgrow_wgs_range   | comparator_range   |
|:--------------------|:---------|----:|------------:|------:|-------:|---------:|---------:|:---------------------|:-------------------|
| real 2bRAD, deduped | S01      |  18 |       0.297 | 0.186 | -0.730 |   -4.446 |    2.986 | 0.08..7.73           | 0.04..2.55         |
| real 2bRAD, deduped | S06      |  18 |       0.435 | 0.192 | -1.210 |   -6.452 |    4.030 | 0.14..10.84          | 0.05..2.33         |
| real 2bRAD, deduped | S07      |  20 |       0.604 | 0.308 | -0.631 |   -5.568 |    4.306 | 0.03..10.45          | 0.17..3.16         |
| real 2bRAD, deduped | ALL      |  56 |       0.455 | 0.238 | -0.849 |   -5.478 |    3.780 | 0.03..10.84          | 0.04..3.16         |
| real 2bRAD, raw     | S01      |   8 |       0.449 | 0.440 | -0.230 |  nan     |  nan     | nan                  | nan                |
| real 2bRAD, raw     | S06      |  10 |       0.797 | 0.728 |  0.600 |  nan     |  nan     | nan                  | nan                |
| real 2bRAD, raw     | S07      |  16 |       0.734 | 0.566 |  0.620 |  nan     |  nan     | nan                  | nan                |

**(b, cont.) Per-anchor capture efficiency (route-B characterization)**

| metric                                   | value               | note                                                                                                  |
|:-----------------------------------------|:--------------------|:------------------------------------------------------------------------------------------------------|
| sigma_eff, raw per-anchor dispersion     | 1.53                | first per-anchor capture-efficiency measurement for 2bRAD                                             |
| ICC-correctable fraction, median (IQR)   | 0.748 (0.679-0.816) | between/(between+within), Poisson-corrected; 21 species with >=30 shared anchors across the 3 samples |
| residual sigma_eff after e_prior, median | 0.699               | within-batch; window averaging (/sqrt(99)) dilutes it to ~0.07-0.10 on the gradient scale             |

**(c) C5 rotating biological contactor (PRJNA974210, 9 samples, 522 MAGs)**

| metric                       | value                                                                                               | note                                                                                                                |
|:-----------------------------|:----------------------------------------------------------------------------------------------------|:--------------------------------------------------------------------------------------------------------------------|
| sk2bgrow reported_fraction   | 1.00 (522/522)                                                                                      | no output gate; a denominator artifact, never a performance claim                                                   |
| sk2bgrow estimate_fraction   | 0.989-1.000                                                                                         | 19 rows carry no PTR estimate (coverage NaN, 0 QC pass)                                                             |
| sk2bgrow qc_recall           | 3.8-13.8% (20-72 of 522)                                                                            | the same-denominator comparator to Pilea default                                                                    |
| Pilea default reported       | 5.0-11.1% (26-58 of 522)                                                                            | overlaps sk2bgrow qc_recall almost exactly                                                                          |
| common denominator           | Pilea-default MAGs are a strict subset of sk2bgrow outputs; sk2bgrow QC passes 15-36 of 26-58 on it | c5_common_denominator.tsv                                                                                           |
| QC pass vs MAG fragmentation | Spearman rho -0.41 vs n_contigs; -0.34 controlling mean coverage                                    | QC works on real data; coverage is the strongest single predictor but the fragmentation effect survives the control |
| cost vs Pilea, per sample    | 89.5-240.8x (measured, c5_cost_per_sample.tsv) -> ~10-25x after --max-mismatch 1                    | mm=1: lookup 63x faster (collapsed middle seed), 5.3% anchors lost, 0 genomes lost; C5 REVIEW §3/§4f                |

Agreement coefficients are over species x sample units. In (a) the comparator is Pilea at its shipped gates; the observed ranges span only ~1.4 log2, so the low CCC (0.25-0.37) is dynamic-range deflation, not broken concordance; the gates-off rows (r ~0) are the control showing Pilea's gate does real work. In (b) exact read-level deduplication before counting flattens the PTR dynamic range (deduped r 0.30-0.60 vs raw 0.45-0.80); the cap-dedup sensitivity sweep found no robust intermediate (cap=2 is best in S07, worst in S06), so the process recommendation is to not exactly dedup 2bRAD libraries before counting; a GLM on deduped counts does not recover the signal (r -0.03-0.41). sigma_eff is within-batch (3 libraries, one study/one centre; cross-lab transfer awaits the Hou cohort) and measured on the single enzyme (BcgI) present in the Sun libraries. In (c) the 1.00 recall headline is reported_fraction (no gate), not accuracy.
