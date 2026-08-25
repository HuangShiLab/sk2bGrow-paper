**Table 2. Accuracy on the Zheng et al. E. coli dataset, by coverage**

|   coverage | method                              |   n |   pearson_r |   rmse_vs_predicted |    bias |   slope |
|-----------:|:------------------------------------|----:|------------:|--------------------:|--------:|--------:|
|        0.5 | sk2bGrow: anchors + V-fit           |  16 |       0.913 |               0.304 |  -0.252 |   0.615 |
|        0.5 | FracMinHash + V-fit                 |  16 |       0.724 |               0.382 |  -0.275 |   0.564 |
|        0.5 | anchors + rank regression           |  16 |       0.164 |               1.188 |   1.119 |   0.058 |
|        0.5 | Pilea gates off: FracMinHash + rank |  16 |     nan     |             nan     | nan     | nan     |
|        0.5 | Pilea (defaults)                    |   0 |     nan     |             nan     | nan     | nan     |
|        1.0 | sk2bGrow: anchors + V-fit           |  16 |       0.981 |               0.157 |  -0.130 |   0.779 |
|        1.0 | FracMinHash + V-fit                 |  16 |       0.940 |               0.213 |  -0.153 |   0.930 |
|        1.0 | anchors + rank regression           |  16 |       0.683 |               1.027 |   0.968 |   0.123 |
|        1.0 | Pilea gates off: FracMinHash + rank |  16 |       0.888 |               0.397 |   0.354 |   0.550 |
|        1.0 | Pilea (defaults)                    |   0 |     nan     |             nan     | nan     | nan     |
|        2.0 | sk2bGrow: anchors + V-fit           |  16 |       0.982 |               0.128 |  -0.078 |   0.840 |
|        2.0 | FracMinHash + V-fit                 |  16 |       0.984 |               0.093 |  -0.039 |   0.875 |
|        2.0 | anchors + rank regression           |  16 |       0.756 |               0.667 |   0.608 |   0.291 |
|        2.0 | Pilea gates off: FracMinHash + rank |  16 |       0.947 |               0.259 |   0.241 |   0.755 |
|        2.0 | Pilea (defaults)                    |   0 |     nan     |             nan     | nan     | nan     |
|        5.0 | sk2bGrow: anchors + V-fit           |  16 |       0.979 |               0.039 |   0.001 |   0.920 |
|        5.0 | FracMinHash + V-fit                 |  16 |       0.977 |               0.051 |   0.031 |   0.950 |
|        5.0 | anchors + rank regression           |  16 |       0.914 |               0.346 |   0.315 |   0.610 |
|        5.0 | Pilea gates off: FracMinHash + rank |  16 |       0.953 |               0.103 |   0.024 |   0.909 |
|        5.0 | Pilea (defaults)                    |   0 |     nan     |             nan     | nan     | nan     |
|       10.0 | sk2bGrow: anchors + V-fit           |  16 |       0.968 |               0.063 |   0.031 |   0.951 |
|       10.0 | FracMinHash + V-fit                 |  16 |       0.942 |               0.149 |   0.096 |   1.018 |
|       10.0 | anchors + rank regression           |  16 |       0.913 |               0.230 |   0.203 |   0.760 |
|       10.0 | Pilea gates off: FracMinHash + rank |  16 |       0.971 |               0.077 |  -0.048 |   0.841 |
|       10.0 | Pilea (defaults)                    |  16 |       0.971 |               0.077 |  -0.048 |   0.841 |

Pearson r is against independently measured growth rate. RMSE is against the lambda*C-derived prediction, which is NOT independent (derived from the same reads by marker-frequency analysis). Blank rows: the method returned no estimate, or a constant. The first four methods are the 2x2 of sketch (2bRAD anchors vs FracMinHash) by estimator (coordinate V-fit vs sorted-rank regression), all run on the same subsampled reads.
