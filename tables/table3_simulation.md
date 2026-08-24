**Table 3. Multi-strain simulation: accuracy and computational cost**

| arm               |   n_strains |   coverage |   recall |    rmse |    bias |   seconds |   peak_rss_mb |
|:------------------|------------:|-----------:|---------:|--------:|--------:|----------:|--------------:|
| Pilea (defaults)  |           4 |      1.000 |    0.000 | nan     | nan     |     1.210 |       181.674 |
| Pilea (defaults)  |           4 |      2.000 |    0.000 | nan     | nan     |     1.160 |       192.471 |
| Pilea (defaults)  |           4 |      4.000 |    0.000 | nan     | nan     |     1.300 |       195.445 |
| Pilea (defaults)  |           4 |      8.000 |    1.000 |   0.075 |  -0.038 |     4.865 |       199.590 |
| Pilea (defaults)  |           8 |      1.000 |    0.000 | nan     | nan     |     1.185 |       201.015 |
| Pilea (defaults)  |           8 |      2.000 |    0.000 | nan     | nan     |     1.300 |       205.234 |
| Pilea (defaults)  |           8 |      4.000 |    0.000 | nan     | nan     |     1.850 |       212.140 |
| Pilea (defaults)  |           8 |      8.000 |    0.938 |   0.080 |  -0.038 |     6.740 |       229.229 |
| Pilea (defaults)  |          16 |      1.000 |    0.000 | nan     | nan     |     1.550 |       229.736 |
| Pilea (defaults)  |          16 |      2.000 |    0.000 | nan     | nan     |     1.970 |       246.948 |
| Pilea (defaults)  |          16 |      4.000 |    0.000 | nan     | nan     |     2.545 |       288.670 |
| Pilea (defaults)  |          16 |      8.000 |    0.750 |   0.095 |  -0.059 |     8.340 |       296.116 |
| Pilea (gates off) |           4 |      1.000 |    1.000 |   0.386 |   0.226 |    21.125 |       181.166 |
| Pilea (gates off) |           4 |      2.000 |    1.000 |   0.297 |   0.271 |    21.060 |       191.734 |
| Pilea (gates off) |           4 |      4.000 |    1.000 |   0.119 |   0.030 |    14.515 |       195.994 |
| Pilea (gates off) |           4 |      8.000 |    1.000 |   0.075 |  -0.038 |     4.660 |       193.749 |
| Pilea (gates off) |           8 |      1.000 |    1.000 |   0.535 |   0.479 |    22.040 |       201.990 |
| Pilea (gates off) |           8 |      2.000 |    1.000 |   0.340 |   0.291 |    21.610 |       206.340 |
| Pilea (gates off) |           8 |      4.000 |    1.000 |   0.170 |   0.100 |    11.690 |       211.780 |
| Pilea (gates off) |           8 |      8.000 |    1.000 |   0.076 |  -0.035 |     6.230 |       229.097 |
| Pilea (gates off) |          16 |      1.000 |    0.969 |   0.636 |   0.439 |    24.720 |       228.418 |
| Pilea (gates off) |          16 |      2.000 |    1.000 |   0.293 |   0.236 |    31.420 |       249.422 |
| Pilea (gates off) |          16 |      4.000 |    1.000 |   0.159 |   0.079 |    24.755 |       285.090 |
| Pilea (gates off) |          16 |      8.000 |    1.000 |   0.090 |  -0.057 |     8.955 |       287.900 |
| sk2bGrow          |           4 |      1.000 |    1.000 |   0.389 |  -0.298 |     3.160 |       188.809 |
| sk2bGrow          |           4 |      2.000 |    1.000 |   0.125 |  -0.057 |     6.255 |       187.146 |
| sk2bGrow          |           4 |      4.000 |    1.000 |   0.101 |   0.014 |    10.380 |       190.169 |
| sk2bGrow          |           4 |      8.000 |    1.000 |   0.049 |   0.011 |    20.320 |       190.882 |
| sk2bGrow          |           8 |      1.000 |    1.000 |   0.244 |  -0.166 |     6.185 |       187.343 |
| sk2bGrow          |           8 |      2.000 |    1.000 |   0.138 |  -0.069 |    11.010 |       190.882 |
| sk2bGrow          |           8 |      4.000 |    1.000 |   0.127 |   0.049 |    19.455 |       189.587 |
| sk2bGrow          |           8 |      8.000 |    1.000 |   0.159 |   0.090 |    38.130 |       182.764 |
| sk2bGrow          |          16 |      1.000 |    1.000 |   0.371 |  -0.309 |    10.925 |       190.816 |
| sk2bGrow          |          16 |      2.000 |    1.000 |   0.162 |  -0.116 |    21.085 |       188.826 |
| sk2bGrow          |          16 |      4.000 |    1.000 |   0.060 |  -0.020 |    40.915 |       188.621 |
| sk2bGrow          |          16 |      8.000 |    1.000 |   0.091 |   0.027 |    80.385 |       185.418 |

**recall** is the fraction of truly-present strains for which the method returned any estimate, and must be read alongside RMSE — a method that reports only the easy cases earns a flattering RMSE. Pilea at its shipped defaults returns nothing below 8x.
