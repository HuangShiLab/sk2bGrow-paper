**Table 12. C5 deployment benchmark: current 16-enzyme/mismatch-2 policy versus 8-enzyme/mismatch-1 fast mode.**

| arm | enzymes | mismatch | screen | threads | wall_s | wall_h | peak_rss_gb | n_finite_ptr | n_qc_pass | median_log2ptr_qc | median_abs_log2ptr_common_qc |
|:---|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|
| k16_mm2_current | 16 | 2 | no | 8 | 76218.63 | 21.1718 | 15.3719 | 8 | 4 | 0.126238 | NA |
| k8_mm1_fast | 8 | 1 | no | 8 | 4551.00 | 1.2642 | 10.3765 | 8 | 3 | 0.148049 | 0.019429 |

Current policy is the C5 refusion arm under the conservative fragmented-reference rule; its runtime combines the original k16/mismatch-2 count stage with the retained-window refit. The fast arm is end-to-end and uses the top-8 ranked enzymes and mismatch 1 without a containment screen. `median_abs_log2ptr_common_qc` is the median absolute log2 PTR difference for the three genomes that pass QC in both arms. Across all eight finite calls, Pearson r was 0.291 and median absolute difference was 0.0924. SLURM job 4095386; sample SRR28338156; 67,423,986 reads.
