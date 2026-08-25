**Table 7. Enzyme-panel size against accuracy and cost (E. coli, 16 media)**

|   enzymes |   anchors |   r (<=2x) |   RMSE (<=2x) |   slope (<=2x) |   run-out bias (<=2x) |   index_s |   profile_s |   peak_RSS_MB |   speedup_vs_16 |
|----------:|----------:|-----------:|--------------:|---------------:|----------------------:|----------:|------------:|--------------:|----------------:|
|         2 |     17055 |      0.951 |         0.188 |          0.808 |                 0.073 |     0.260 |       3.243 |           162 |           2.533 |
|         4 |     24753 |      0.959 |         0.165 |          0.825 |                 0.074 |     0.430 |       4.431 |           173 |           1.854 |
|         8 |     37232 |      0.969 |         0.171 |          0.809 |                 0.123 |     0.820 |       6.493 |           186 |           1.265 |
|        12 |     41662 |      0.969 |         0.183 |          0.784 |                 0.132 |     1.220 |       7.762 |           192 |           1.058 |
|        16 |     43735 |      0.959 |         0.196 |          0.745 |                 0.145 |     1.530 |       8.215 |           195 |           1.000 |

Subsets are the top k of Table 6. Accuracy columns average the 0.5/1/2x depths, the regime the panel exists for; cost columns average all five. Index cost is one-off per reference. Accuracy peaks at 4-8 enzymes and decays to 16, while cost is linear in k -- the four sparsest enzymes (< 450 anchors) buy nothing and double the bias on the replication run-out control. We did not establish the mechanism: two candidate explanations were tested and both refuted (sparse enzymes are not more biased than dense ones, and forcing fixed-effect weights does not recover the loss). Same cells, Pilea, gates off 11.5 s / 157 MB; Pilea, defaults 1.9 s / 157 MB.
