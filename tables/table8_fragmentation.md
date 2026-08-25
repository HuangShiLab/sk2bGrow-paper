**Table 8. Reference fragmentation: the same reads against a complete, a fragmented, and a scaffolded reference**

|   coverage | reference                             |   n |   pearson_r |   rmse |   bias |   slope |   qc_pass |   stationary_control |
|-----------:|:--------------------------------------|----:|------------:|-------:|-------:|--------:|----------:|---------------------:|
|        0.5 | complete chromosome                   |  16 |       0.913 |  0.304 | -0.252 |   0.615 |     0.000 |                0.260 |
|        0.5 | 100 contigs                           |  16 |       0.382 |  0.833 | -0.750 |   0.082 |     0.000 |                0.150 |
|        0.5 | 100 contigs, scaffolded vs itself     |  16 |       0.922 |  0.304 | -0.255 |   0.632 |     0.000 |                0.257 |
|        0.5 | 100 contigs, scaffolded vs a relative |  16 |       0.921 |  0.266 | -0.227 |   0.685 |     0.000 |                0.299 |
|        0.5 | 100 contigs, Pilea                    |  16 |     nan     |  1.117 | -1.045 | nan     |   nan     |                0.000 |
|        0.5 | 100 contigs, order-free spread MLE    |  16 |       0.308 |  1.067 | -0.994 |   0.143 |   nan     |                0.000 |
|        1.0 | complete chromosome                   |  16 |       0.981 |  0.157 | -0.130 |   0.779 |     0.000 |                0.077 |
|        1.0 | 100 contigs                           |  16 |       0.550 |  0.890 | -0.817 |   0.103 |     0.000 |                0.192 |
|        1.0 | 100 contigs, scaffolded vs itself     |  16 |       0.983 |  0.155 | -0.128 |   0.772 |     0.000 |                0.078 |
|        1.0 | 100 contigs, scaffolded vs a relative |  16 |       0.977 |  0.152 | -0.117 |   0.755 |     0.000 |                0.113 |
|        1.0 | 100 contigs, Pilea                    |  16 |       0.827 |  0.545 |  0.508 |   0.628 |   nan     |                1.153 |
|        1.0 | 100 contigs, order-free spread MLE    |  16 |      -0.273 |  1.108 | -0.977 |  -0.167 |   nan     |                0.000 |
|        2.0 | complete chromosome                   |  16 |       0.982 |  0.128 | -0.078 |   0.840 |    50.000 |                0.097 |
|        2.0 | 100 contigs                           |  16 |       0.774 |  0.886 | -0.829 |   0.208 |    75.000 |                0.135 |
|        2.0 | 100 contigs, scaffolded vs itself     |  16 |       0.986 |  0.123 | -0.074 |   0.837 |    50.000 |                0.104 |
|        2.0 | 100 contigs, scaffolded vs a relative |  16 |       0.982 |  0.099 | -0.046 |   0.856 |    62.500 |                0.123 |
|        2.0 | 100 contigs, Pilea                    |  16 |       0.944 |  0.347 |  0.334 |   0.790 |   nan     |                0.621 |
|        2.0 | 100 contigs, order-free spread MLE    |  16 |       0.384 |  0.914 | -0.712 |   0.537 |   nan     |                0.000 |
|        5.0 | complete chromosome                   |  16 |       0.979 |  0.039 |  0.001 |   0.920 |    75.000 |                0.060 |
|        5.0 | 100 contigs                           |  16 |       0.871 |  0.871 | -0.813 |   0.187 |   100.000 |                0.050 |
|        5.0 | 100 contigs, scaffolded vs itself     |  16 |       0.981 |  0.042 | -0.002 |   0.910 |    62.500 |                0.055 |
|        5.0 | 100 contigs, scaffolded vs a relative |  16 |       0.978 |  0.053 |  0.023 |   0.946 |    81.250 |                0.064 |
|        5.0 | 100 contigs, Pilea                    |  16 |       0.948 |  0.116 |  0.053 |   0.888 |   nan     |                0.255 |
|        5.0 | 100 contigs, order-free spread MLE    |  16 |       0.917 |  0.144 |  0.080 |   0.922 |   nan     |                0.383 |
|       10.0 | complete chromosome                   |  16 |       0.968 |  0.063 |  0.031 |   0.951 |    75.000 |                0.046 |
|       10.0 | 100 contigs                           |  16 |       0.859 |  0.862 | -0.808 |   0.210 |   100.000 |                0.059 |
|       10.0 | 100 contigs, scaffolded vs itself     |  16 |       0.973 |  0.056 |  0.026 |   0.948 |    68.750 |                0.059 |
|       10.0 | 100 contigs, scaffolded vs a relative |  16 |       0.967 |  0.086 |  0.052 |   0.984 |    87.500 |                0.052 |
|       10.0 | 100 contigs, Pilea                    |  16 |       0.960 |  0.077 | -0.031 |   0.820 |   nan     |                0.203 |
|       10.0 | 100 contigs, order-free spread MLE    |  16 |       0.871 |  0.227 |  0.149 |   0.919 |   nan     |                0.458 |

The 100-contig reference holds 43,707 of the complete genome's 43,735 anchors, so the genomic coordinate is the only variable. qc_pass is the share of estimates the fusion QC accepts -- note that it is HIGHEST where the estimates are worst. stationary_control is the estimate for the RUN_OUT sample, whose true log2(PTR) is ~0. Blank r or slope: the method returned nothing, or a constant.
