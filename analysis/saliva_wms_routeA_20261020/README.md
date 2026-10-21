
## Completion status, 2026-10-21

All 32 matched saliva WMS samples completed route-A profiling. One sample is
labelled `skipped` in the status file because its pilot result already existed
and was reused; all 32 final `output.tsv` files are present.

Protocol:

- route A / shotgun WMS;
- same GTDB species-representative BcgI anchor database as route B;
- `--mode wms --max-mismatch 2 --window-anchors 100`;
- `--no-gc-correct`;
- both `--method auto` (`noGC_auto`) and `--method glm` (`noGC_glm`) were run
  from the same count table.

Final paired analyses are in `results/paired_analysis/`.

Key yield:

| Assay / method | species-sample rows | finite log2 PTR | exploratory robust | strict QC |
|---|---:|---:|---:|---:|
| WMS auto | 31,040 | 1,690 | 1,132 | 76 |
| WMS glm | 31,040 | 13,123 | 9,112 | 361 |
| matched 2bRAD glm | 31,040 | 14,040 | 8,785 | 26 |

Across 8 hosts and 4 timepoints, all 32 samples had at least one exploratory
robust estimate for WMS-auto, WMS-glm and 2bRAD-glm.

Preliminary paired log2-PTR correlations:

- WMS-glm vs 2bRAD-glm, all shared finite: `n=11,343`,
  Pearson `r=0.460`, Spearman `rho=0.436`.
- Both exploratory robust: `n=6,434`,
  Pearson `r=0.418`, Spearman `rho=0.403`.
- WMS-auto correlations were weak.

Sample mapping correction: the leading code in `S11_1` is the timepoint
(`11AM`); the trailing code is the host (`1`). The final paired tables use 8
hosts and 4 timepoints.
