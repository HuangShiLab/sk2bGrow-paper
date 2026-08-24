**Table 4. Multi-strain simulation, aggregate over the whole grid**

| arm               |   recall |   rmse |   bias |   spurious |   seconds |   peak_rss_mb |
|:------------------|---------:|-------:|-------:|-----------:|----------:|--------------:|
| Pilea (defaults)  |    0.224 |  0.083 | -0.045 |      0.000 |     2.835 |       223.189 |
| Pilea (gates off) |    0.997 |  0.265 |  0.168 |      0.000 |    17.732 |       221.890 |
| sk2bGrow          |    1.000 |  0.168 | -0.070 |      0.000 |    22.350 |       188.439 |

**spurious** counts genomes reported that were not in the sample (false positives); both methods scored zero.
