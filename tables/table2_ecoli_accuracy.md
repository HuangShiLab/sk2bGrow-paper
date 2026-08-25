**Table 2. Accuracy on the Zheng et al. E. coli dataset, by coverage**

|   coverage | method                            |   n |   pearson_r |   rmse_vs_predicted |   slope |
|-----------:|:----------------------------------|----:|------------:|--------------------:|--------:|
|      0.500 | sk2bGrow                          |  16 |       0.913 |               0.304 |   0.615 |
|      0.500 | sk2bGrow (Pilea-parity estimator) |  16 |       0.164 |               1.188 |   0.058 |
|      0.500 | Pilea (defaults)                  |   0 |     nan     |             nan     | nan     |
|      0.500 | Pilea (gates off)                 |  16 |     nan     |             nan     | nan     |
|      1.000 | sk2bGrow                          |  16 |       0.981 |               0.157 |   0.779 |
|      1.000 | sk2bGrow (Pilea-parity estimator) |  16 |       0.683 |               1.027 |   0.123 |
|      1.000 | Pilea (defaults)                  |   0 |     nan     |             nan     | nan     |
|      1.000 | Pilea (gates off)                 |  16 |       0.888 |               0.397 |   0.550 |
|      2.000 | sk2bGrow                          |  16 |       0.982 |               0.128 |   0.840 |
|      2.000 | sk2bGrow (Pilea-parity estimator) |  16 |       0.756 |               0.667 |   0.291 |
|      2.000 | Pilea (defaults)                  |   0 |     nan     |             nan     | nan     |
|      2.000 | Pilea (gates off)                 |  16 |       0.947 |               0.259 |   0.755 |
|      5.000 | sk2bGrow                          |  16 |       0.979 |               0.039 |   0.920 |
|      5.000 | sk2bGrow (Pilea-parity estimator) |  16 |       0.914 |               0.346 |   0.610 |
|      5.000 | Pilea (defaults)                  |   0 |     nan     |             nan     | nan     |
|      5.000 | Pilea (gates off)                 |  16 |       0.953 |               0.103 |   0.909 |
|     10.000 | sk2bGrow                          |  16 |       0.968 |               0.063 |   0.951 |
|     10.000 | sk2bGrow (Pilea-parity estimator) |  16 |       0.913 |               0.230 |   0.760 |
|     10.000 | Pilea (defaults)                  |  16 |       0.971 |               0.077 |   0.841 |
|     10.000 | Pilea (gates off)                 |  16 |       0.971 |               0.077 |   0.841 |

Pearson r is against independently measured growth rate. RMSE is against the lambda*C-derived prediction, which is NOT independent (derived from the same reads by marker-frequency analysis). Blank rows: the method returned no estimate, or a constant.
