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
| sk2bGrow          |           4 |      1.000 |    1.000 |   0.269 |  -0.170 |     2.470 |       185.672 |
| sk2bGrow          |           4 |      2.000 |    1.000 |   0.080 |   0.016 |     4.700 |       190.611 |
| sk2bGrow          |           4 |      4.000 |    1.000 |   0.095 |   0.028 |     7.640 |       189.342 |
| sk2bGrow          |           4 |      8.000 |    1.000 |   0.050 |   0.018 |    14.960 |       190.235 |
| sk2bGrow          |           8 |      1.000 |    1.000 |   0.177 |  -0.080 |     4.670 |       186.425 |
| sk2bGrow          |           8 |      2.000 |    1.000 |   0.105 |   0.003 |     8.130 |       188.596 |
| sk2bGrow          |           8 |      4.000 |    1.000 |   0.142 |   0.073 |    14.395 |       182.567 |
| sk2bGrow          |           8 |      8.000 |    1.000 |   0.164 |   0.101 |    27.795 |       182.018 |
| sk2bGrow          |          16 |      1.000 |    1.000 |   0.251 |  -0.122 |     7.980 |       185.377 |
| sk2bGrow          |          16 |      2.000 |    1.000 |   0.133 |  -0.059 |    15.330 |       191.111 |
| sk2bGrow          |          16 |      4.000 |    1.000 |   0.051 |  -0.001 |    29.780 |       191.267 |
| sk2bGrow          |          16 |      8.000 |    1.000 |   0.094 |   0.037 |    58.395 |       183.624 |

**recall** is the fraction of truly-present strains for which the method returned any estimate, and must be read alongside RMSE — a method that reports only the easy cases earns a flattering RMSE. Pilea at its shipped defaults returns nothing below 8x.
