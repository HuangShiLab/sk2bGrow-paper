**Table 6. Single-enzyme PTR accuracy on the E. coli panel**

|   rank | enzyme   |   anchors |   windows |   fit_rate |   r@0.5x |   r@1x |   r@2x |   r@5x |   r@10x |   r_low |   rmse_low |   ctl_bias |   score |
|-------:|:---------|----------:|----------:|-----------:|---------:|-------:|-------:|-------:|--------:|--------:|-----------:|-----------:|--------:|
|      1 | CjePI    |      7616 |        76 |      0.953 |    0.901 |  0.934 |  0.970 |  0.963 |   0.939 |   0.935 |      0.193 |      0.029 |   0.879 |
|      2 | CjeI     |      8797 |        88 |      1.000 |    0.863 |  0.958 |  0.975 |  0.975 |   0.962 |   0.932 |      0.186 |      0.058 |   0.871 |
|      3 | AlfI     |      1989 |        25 |      0.953 |    0.821 |  0.914 |  0.884 |  0.936 |   0.937 |   0.873 |      0.268 |      0.058 |   0.792 |
|      4 | Hin4I    |      5330 |        53 |      1.000 |    0.697 |  0.949 |  0.971 |  0.961 |   0.952 |   0.872 |      0.227 |      0.109 |   0.788 |
|      5 | HaeIV    |      6435 |        64 |      1.000 |    0.711 |  0.920 |  0.964 |  0.968 |   0.939 |   0.865 |      0.259 |      0.101 |   0.775 |
|      6 | BcgI     |      2872 |        29 |      0.953 |    0.651 |  0.956 |  0.952 |  0.961 |   0.957 |   0.853 |      0.258 |      0.232 |   0.731 |
|      7 | Bsp24I   |      1524 |        24 |      0.965 |    0.675 |  0.717 |  0.833 |  0.894 |   0.923 |   0.742 |      0.362 |      0.055 |   0.637 |
|      8 | BsaXI    |       970 |        26 |      0.953 |    0.327 |  0.917 |  0.883 |  0.976 |   0.965 |   0.709 |      0.370 |      0.074 |   0.598 |
|      9 | BslFI    |      2501 |        25 |      0.941 |    0.342 |  0.684 |  0.904 |  0.924 |   0.977 |   0.644 |      0.316 |      0.047 |   0.553 |
|     10 | CspCI    |       563 |        23 |      0.882 |    0.358 |  0.720 |  0.932 |  0.947 |   0.944 |   0.670 |      0.450 |      0.205 |   0.506 |
|     11 | PpiI     |       336 |        13 |      0.871 |    0.609 |  0.314 |  0.814 |  0.863 |   0.869 |   0.579 |      0.566 |      0.073 |   0.419 |
|     12 | BaeI     |       729 |        24 |      0.871 |    0.557 |  0.497 |  0.815 |  0.910 |   0.949 |   0.623 |      0.507 |      0.400 |   0.396 |
|     13 | FalI     |       687 |        25 |      0.906 |   -0.056 |  0.528 |  0.867 |  0.960 |   0.890 |   0.446 |      0.537 |      0.246 |   0.250 |
|     14 | AloI     |       441 |        18 |      0.906 |   -0.029 |  0.591 |  0.916 |  0.934 |   0.933 |   0.493 |      0.488 |      0.510 |   0.243 |
|     15 | BplI     |       379 |        15 |      0.906 |    0.089 |  0.178 |  0.754 |  0.851 |   0.884 |   0.340 |      0.568 |      0.240 |   0.138 |
|     16 | PsrI     |       419 |        17 |      0.918 |   -0.069 |  0.326 |  0.620 |  0.913 |   0.967 |   0.292 |      0.582 |      0.153 |   0.109 |

Each enzyme fitted alone, before fusion; r is against measured growth rate across 16 media. r_low and rmse_low average the 0.5/1/2x depths; ctl_bias is mean |log2PTR| on the run-out control, which should be 0. The ranking tracks anchor yield almost exactly: the top six are the six enzymes with more than 2,500 anchors, the bottom four have fewer than 450 and are anti-correlated at 0.5x. These fits share an origin estimated from all 16 enzymes pooled, so the table measures how informative an enzyme is given a good origin, not how a panel of that enzyme alone would behave -- Table 7 tests that directly.
