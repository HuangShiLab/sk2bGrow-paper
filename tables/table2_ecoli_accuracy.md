**Table 2. Accuracy on the Zheng et al. E. coli dataset, by coverage**

|   coverage | method                              |   n |   pearson_r |   rmse_vs_predicted |    bias |   slope |    r_lo |    r_hi |   slope_lo |   slope_hi |
|-----------:|:------------------------------------|----:|------------:|--------------------:|--------:|--------:|--------:|--------:|-----------:|-----------:|
|        0.5 | sk2bGrow: anchors + V-fit           |  16 |       0.941 |               0.407 |  -0.370 |   0.617 |   0.882 |   0.975 |      0.513 |      0.736 |
|        0.5 | FracMinHash + V-fit                 |  16 |       0.883 |               0.522 |  -0.473 |   0.483 |   0.766 |   0.959 |      0.351 |      0.637 |
|        0.5 | anchors + rank regression           |  16 |       0.517 |               0.933 |   0.871 |   0.162 |   0.008 |   0.788 |      0.002 |      0.292 |
|        0.5 | Pilea gates off: FracMinHash + rank |  16 |     nan     |             nan     | nan     | nan     | nan     | nan     |      0.000 |      0.000 |
|        0.5 | Pilea (defaults)                    |   0 |     nan     |             nan     | nan     | nan     | nan     | nan     |    nan     |    nan     |
|        1.0 | sk2bGrow: anchors + V-fit           |  16 |       0.919 |               0.346 |  -0.312 |   0.795 |   0.819 |   0.972 |      0.578 |      0.952 |
|        1.0 | FracMinHash + V-fit                 |  16 |       0.912 |               0.286 |  -0.232 |   0.675 |   0.797 |   0.975 |      0.458 |      0.811 |
|        1.0 | anchors + rank regression           |  16 |       0.775 |               0.395 |   0.303 |   0.361 |   0.419 |   0.936 |      0.174 |      0.485 |
|        1.0 | Pilea gates off: FracMinHash + rank |  16 |       0.778 |               0.555 |   0.509 |   0.487 |   0.478 |   0.940 |      0.279 |      0.662 |
|        1.0 | Pilea (defaults)                    |   0 |     nan     |             nan     | nan     | nan     | nan     | nan     |    nan     |    nan     |
|        2.0 | sk2bGrow: anchors + V-fit           |  16 |       0.959 |               0.148 |  -0.109 |   0.788 |   0.906 |   0.986 |      0.639 |      0.902 |
|        2.0 | FracMinHash + V-fit                 |  16 |       0.958 |               0.182 |  -0.147 |   0.783 |   0.891 |   0.987 |      0.657 |      0.892 |
|        2.0 | anchors + rank regression           |  16 |       0.925 |               0.162 |   0.074 |   0.652 |   0.816 |   0.978 |      0.474 |      0.772 |
|        2.0 | Pilea gates off: FracMinHash + rank |  16 |       0.848 |               0.243 |   0.173 |   0.696 |   0.646 |   0.942 |      0.412 |      0.896 |
|        2.0 | Pilea (defaults)                    |   0 |     nan     |             nan     | nan     | nan     | nan     | nan     |    nan     |    nan     |
|        5.0 | sk2bGrow: anchors + V-fit           |  16 |       0.971 |               0.060 |  -0.006 |   0.924 |   0.945 |   0.990 |      0.809 |      1.050 |
|        5.0 | FracMinHash + V-fit                 |  16 |       0.990 |               0.056 |  -0.013 |   0.930 |   0.976 |   0.997 |      0.872 |      0.995 |
|        5.0 | anchors + rank regression           |  16 |       0.947 |               0.089 |   0.049 |   0.868 |   0.880 |   0.985 |      0.714 |      1.017 |
|        5.0 | Pilea gates off: FracMinHash + rank |  16 |       0.988 |               0.113 |   0.091 |   0.800 |   0.970 |   0.997 |      0.743 |      0.858 |
|        5.0 | Pilea (defaults)                    |   0 |     nan     |             nan     | nan     | nan     | nan     | nan     |    nan     |    nan     |
|       10.0 | sk2bGrow: anchors + V-fit           |  16 |       0.970 |               0.049 |   0.022 |   0.951 |   0.933 |   0.990 |      0.833 |      1.067 |
|       10.0 | FracMinHash + V-fit                 |  16 |       0.987 |               0.040 |   0.013 |   0.943 |   0.969 |   0.995 |      0.870 |      1.014 |
|       10.0 | anchors + rank regression           |  16 |       0.938 |               0.105 |   0.062 |   0.900 |   0.828 |   0.985 |      0.725 |      1.047 |
|       10.0 | Pilea gates off: FracMinHash + rank |  16 |       0.973 |               0.117 |  -0.104 |   0.822 |   0.941 |   0.992 |      0.707 |      0.912 |
|       10.0 | Pilea (defaults)                    |  16 |       0.973 |               0.117 |  -0.104 |   0.822 | nan     | nan     |    nan     |    nan     |

Pearson r is against independently measured growth rate. RMSE is against the lambda*C-derived prediction, which is NOT independent (derived from the same reads by marker-frequency analysis). Blank rows: the method returned no estimate, or a constant. The first four methods are the 2x2 of sketch (2bRAD anchors vs FracMinHash) by estimator (coordinate V-fit vs sorted-rank regression), all run on the same subsampled reads. pearson_r and slope carry 95% CIs from a media bootstrap (10,000 resamples of the 16 media; data/ci_bootstrap/table2_ci.tsv).
