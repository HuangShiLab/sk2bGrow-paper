| WMS fitting | Comparison | Species-time pairs | Pearson r | Spearman rho | Median abs. delta log2 PTR |
|---|---|---:|---:|---:|---:|
| auto | Shared finite estimates | 1,549 | -0.021 | 0.055 | 0.639 |
| auto | Both exploratory robust | 820 | 0.135 | 0.152 | 0.363 |
| glm | Shared finite estimates | 11,343 | 0.460 | 0.436 | 0.556 |
| glm | Both exploratory robust | 6,434 | 0.418 | 0.403 | 0.361 |
| glm | Both strict QC | 2 | — | — | — |

WMS-auto denotes the shipped coordinate/sorted fallback policy; WMS-glm and
2bRAD use the coordinate-aware Poisson GLM fit on the same BcgI anchor
database. "Both exploratory robust" is a pre-specified yield filter, not a
claim of absolute PTR accuracy.
