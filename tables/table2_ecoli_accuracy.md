**Table 2. Accuracy on the Zheng et al. E. coli dataset, by coverage**

|   coverage | method                            |   n |   pearson_r |   rmse_vs_predicted |   slope |
|-----------:|:----------------------------------|----:|------------:|--------------------:|--------:|
|      0.500 | sk2bGrow                          |  16 |       0.907 |               0.506 |   0.519 |
|      0.500 | sk2bGrow (Pilea-parity estimator) |  16 |       0.445 |               1.129 |   0.115 |
|      0.500 | Pilea (defaults)                  |   0 |     nan     |             nan     | nan     |
|      0.500 | Pilea (gates off)                 |  16 |     nan     |             nan     | nan     |
|      1.000 | sk2bGrow                          |  16 |       0.954 |               0.339 |   0.720 |
|      1.000 | sk2bGrow (Pilea-parity estimator) |  16 |       0.605 |               0.972 |   0.130 |
|      1.000 | Pilea (defaults)                  |   0 |     nan     |             nan     | nan     |
|      1.000 | Pilea (gates off)                 |  16 |       0.888 |               0.397 |   0.550 |
|      2.000 | sk2bGrow                          |  16 |       0.975 |               0.185 |   0.755 |
|      2.000 | sk2bGrow (Pilea-parity estimator) |  16 |       0.818 |               0.685 |   0.311 |
|      2.000 | Pilea (defaults)                  |   0 |     nan     |             nan     | nan     |
|      2.000 | Pilea (gates off)                 |  16 |       0.947 |               0.259 |   0.755 |
|      5.000 | sk2bGrow                          |  16 |       0.979 |               0.034 |   0.918 |
|      5.000 | sk2bGrow (Pilea-parity estimator) |  16 |       0.912 |               0.371 |   0.584 |
|      5.000 | Pilea (defaults)                  |   0 |     nan     |             nan     | nan     |
|      5.000 | Pilea (gates off)                 |  16 |       0.953 |               0.103 |   0.909 |
|     10.000 | sk2bGrow                          |  15 |       0.974 |               0.056 |   0.945 |
|     10.000 | sk2bGrow (Pilea-parity estimator) |  15 |       0.898 |               0.257 |   0.715 |
|     10.000 | Pilea (defaults)                  |  15 |       0.972 |               0.079 |   0.842 |
|     10.000 | Pilea (gates off)                 |  15 |       0.972 |               0.079 |   0.842 |

Pearson r is against independently measured growth rate. RMSE is against the lambda*C-derived prediction, which is NOT independent (derived from the same reads by marker-frequency analysis). Blank rows: the method returned no estimate, or a constant.
