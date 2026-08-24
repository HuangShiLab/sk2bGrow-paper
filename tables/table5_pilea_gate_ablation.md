**Table 5. Which of Pilea's quality gates suppresses the estimate**

| dataset    |   depth |   n |   median k-mer coverage |   min-cove (-x 5) |   min-frac (-z 0.75) |   min-cont (-c 0.25) |   reported by defaults |
|:-----------|--------:|----:|------------------------:|------------------:|---------------------:|---------------------:|-----------------------:|
| ecoli      |   0.500 |  17 |                   1.000 |             0.000 |              100.000 |              100.000 |                  0.000 |
| ecoli      |   1.000 |  17 |                   1.454 |             0.000 |              100.000 |              100.000 |                  0.000 |
| ecoli      |   2.000 |  17 |                   1.995 |             0.000 |              100.000 |              100.000 |                  0.000 |
| ecoli      |   5.000 |  17 |                   3.933 |             0.000 |              100.000 |              100.000 |                  0.000 |
| ecoli      |  10.000 |  16 |                   7.573 |           100.000 |              100.000 |              100.000 |                100.000 |
| simulation |   1.000 |  55 |                   1.439 |             0.000 |               96.364 |               89.091 |                  0.000 |
| simulation |   2.000 |  56 |                   1.969 |             0.000 |               92.857 |               94.643 |                  0.000 |
| simulation |   4.000 |  56 |                   3.292 |             0.000 |               92.857 |               94.643 |                  0.000 |
| simulation |   8.000 |  56 |                   6.416 |           100.000 |               94.643 |               96.429 |                 91.071 |

Percentages are the share of genome-estimates passing that gate alone; "reported by defaults" is the share passing all three. Of the 240 estimates the defaults discard, the sole cause is min-cove (-x 5) for 213, min-frac (-z 0.75) for 3, min-cont (-c 0.25) for 2. A 150 bp read yields 120 31-mers, so the k-mer-coverage threshold of 5 corresponds to about 6.6x read coverage.
