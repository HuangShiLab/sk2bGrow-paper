# Public real-2bRAD route-B story check

## Run setting

- Public raw FASTQs were profiled with sk2bGrow route B against the GTDB species-representative BcgI database.
- Real-2bRAD settings: `--mode 2brad`, `--max-mismatch 0`, no residual two-pass GC correction, `glm` statistics.
- Three label-complete datasets were analyzed: bladder cancer (`PRJNA946904`, 22), ovarian cancer versus benign tissue (`PRJNA1005427`, 20), and maternal faeces versus infant meconium (`PRJCA030517/CRA019236`, 63).
- The IC/urinary dataset was not used as a benchmark because only four of the paper's 22 samples were exposed through ENA when checked.

## Yield

| Dataset | Samples | Median finite PTRs | Median robust PTRs |
|---|---:|---:|---:|
| Bladder | 22 | 12 | 10 |
| Ovary | 20 | 12 | 9 |
| Maternal faeces/meconium | 63 | 47 | 24 |

Only one of 61,110 genome-sample estimates passed the original strict route-A-style QC. The “robust” exploratory estimates therefore use a route-B-aware support filter (`n_windows >= 5`, `n_anchors >= 50`, origin confidence `>= 0.5`, PTR SE `<= 1`, PTR 1–8). They should not be equated with strict-QC PTRs.

## Classification result

Feature arms were processed inside leave-one-out cross-validation. Robust PTR missingness was represented by explicit missing indicators; median—not zero—imputation was used.

| Dataset | Feature arm | AUROC | AUROC bootstrap 95% CI | AUPRC |
|---|---|---:|---:|---:|
| Ovarian cancer vs benign | Taxonomic support | 0.780 | 0.515–1.000 | 0.730 |
| Ovarian cancer vs benign | Growth PTR | 0.930 | 0.786–1.000 | 0.929 |
| Ovarian cancer vs benign | Combined | 0.810 | 0.560–1.000 | 0.739 |
| Bladder NMIBC vs MIBC | Taxonomic support | 0.990 | 0.946–1.000 | 0.982 |
| Bladder NMIBC vs MIBC | Growth PTR | 0.957 | 0.857–1.000 | 0.913 |
| Bladder NMIBC vs MIBC | Combined | 0.981 | 0.914–1.000 | 0.968 |
| Maternal faeces vs meconium | Taxonomic support | 0.617 | 0.476–0.750 | 0.665 |
| Maternal faeces vs meconium | Growth PTR | 0.538 | 0.393–0.683 | 0.581 |
| Maternal faeces vs meconium | Combined | 0.594 | 0.454–0.734 | 0.649 |

Paired LOOCV score comparisons did not establish a significant improvement of growth PTR over taxonomic support: ovarian `ΔAUC(B−A)=0.150`, permutation `p=0.099`; bladder `ΔAUC(B−A)=−0.033`, `p=0.554`; maternal `ΔAUC(B−A)=−0.079`, `p=0.012`.

## Exploratory group effects

No contrast survived BH FDR `q < 0.1`, but ovarian cancer showed the clearest direction:

- `Ralstonia sp007997035`: median log2 PTR effect `+0.420` in cancer, raw `p=0.014`, `q=0.126`.
- `Sphingomonas paucimobilis`: `−0.552`, raw `p=0.038`, `q=0.147`.
- `Ralstonia mannitolilytica`: `+0.403`, raw `p=0.049`, `q=0.147`.

For meconium versus maternal faeces, `Parabacteroides distasonis` had the largest effect (`+0.576`; raw `p=0.012`, `q=0.235`).

## Story decision

There is a defensible exploratory story, but it should not claim universal superiority. The strongest framing is:

> In real 2bRAD-M public data, sk2bGrow can recover a usable growth dimension even from very sparse single-enzyme data. In the ovarian-tissue cohort, robust PTR features alone gave the best point estimate for cancer versus benign classification, whereas in the bladder cohort taxonomic composition was already near-perfect and growth features were complementary rather than necessary.

This is consistent with the earlier positioning: sk2bGrow adds growth dynamics in shallow 2bRAD data, but does not universally outperform FracMinHash/Pilea or abundance-only analysis.

## Required caveats

- Cohorts are small; ovarian and bladder confidence intervals are wide and LOOCV can be optimistic.
- Most estimates fail strict QC; group tests are exploratory.
- The BcgI-only species database has 970 GTDB representatives and does not represent the full GTDB/strain space.
- No published abundance table was available for these public runs, so “taxonomic support” is the route-B coverage/support layer rather than an independent 2bRAD-M abundance baseline.
- Breast milk and the incomplete IC dataset were not included in this first external story check.
