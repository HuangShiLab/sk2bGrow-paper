**Table 6. Single-enzyme PTR accuracy on the E. coli panel**

|   rank | enzyme   |   anchors |   windows |   fit_rate |   r@0.5x |   r@1x |   r@2x |   r@5x |   r@10x |   r_low |   rmse_low |   ctl_bias |   score |
|-------:|:---------|----------:|----------:|-----------:|---------:|-------:|-------:|-------:|--------:|--------:|-----------:|-----------:|--------:|
|      1 | CjeI     |      8797 |        88 |      0.988 |    0.873 |  0.959 |  0.979 |  0.975 |   0.962 |   0.937 |      0.197 |      0.042 |   0.877 |
|      2 | CjePI    |      7466 |        74 |      0.940 |    0.789 |  0.945 |  0.972 |  0.958 |   0.900 |   0.902 |      0.496 |      0.015 |   0.774 |
|      3 | AlfI     |      1989 |        25 |      0.917 |    0.808 |  0.897 |  0.881 |  0.936 |   0.952 |   0.862 |      0.289 |      0.094 |   0.767 |
|      4 | HaeIV    |      6435 |        64 |      0.976 |    0.622 |  0.912 |  0.965 |  0.967 |   0.941 |   0.833 |      0.272 |      0.099 |   0.740 |
|      5 | BcgI     |      2872 |        29 |      0.952 |    0.682 |  0.952 |  0.952 |  0.962 |   0.957 |   0.862 |      0.267 |      0.302 |   0.719 |
|      6 | Hin4I    |      5230 |        52 |      0.964 |    0.525 |  0.883 |  0.967 |  0.965 |   0.971 |   0.792 |      0.342 |      0.063 |   0.691 |
|      7 | Bsp24I   |      1587 |        25 |      0.976 |    0.760 |  0.807 |  0.866 |  0.885 |   0.958 |   0.811 |      0.700 |      0.058 |   0.621 |
|      8 | BslFI    |      2501 |        25 |      0.929 |    0.380 |  0.705 |  0.900 |  0.925 |   0.978 |   0.661 |      0.309 |      0.020 |   0.579 |
|      9 | CspCI    |       563 |        23 |      0.869 |    0.308 |  0.708 |  0.935 |  0.946 |   0.948 |   0.650 |      0.460 |      0.193 |   0.487 |
|     10 | BsaXI    |       932 |        25 |      0.905 |    0.565 |  0.877 |  0.757 |  0.879 |   0.934 |   0.733 |      0.698 |      0.318 |   0.479 |
|     11 | BaeI     |       729 |        24 |      0.869 |    0.559 |  0.453 |  0.805 |  0.912 |   0.953 |   0.605 |      0.514 |      0.336 |   0.393 |
|     12 | FalI     |       687 |        25 |      0.881 |    0.014 |  0.459 |  0.869 |  0.960 |   0.923 |   0.447 |      0.523 |      0.297 |   0.242 |
|     13 | BplI     |       366 |        14 |      0.905 |   -0.196 |  0.294 |  0.745 |  0.852 |   0.893 |   0.281 |      0.580 |      0.287 |   0.065 |
|     14 | PsrI     |       419 |        17 |      0.917 |   -0.150 |  0.290 |  0.611 |  0.913 |   0.969 |   0.250 |      0.602 |      0.143 |   0.064 |
|     15 | PpiI     |       311 |        12 |      0.762 |   -0.185 | -0.069 |  0.488 |  0.906 |   0.850 |   0.078 |      0.891 |      0.181 |  -0.190 |
|     16 | AloI     |       416 |        17 |      0.881 |   -0.367 | -0.110 |  0.440 |  0.881 |   0.888 |  -0.013 |      0.878 |      0.346 |  -0.319 |

Each enzyme fitted alone, before fusion; r is against measured growth rate across 16 media. r_low and rmse_low average the 0.5/1/2x depths; ctl_bias is mean |log2PTR| on the run-out control, which should be 0. The ranking tracks anchor yield almost exactly: the top six are the six enzymes with more than 2,500 anchors, the bottom four have fewer than 450 and are anti-correlated at 0.5x. These fits share an origin estimated from all 16 enzymes pooled, so the table measures how informative an enzyme is given a good origin, not how a panel of that enzyme alone would behave -- Table 7 tests that directly.
