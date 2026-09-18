# Signed-origin rerun and analysis decision

## Provenance

* Source grid: HPC C1 paired-end Zheng instance, 16 media plus `RUN_OUT`.
* Per-enzyme fits: `../m1_rerun/per_enzyme_zheng.tsv`.
* Rule: the shared origin is selected once, with the usual uphill-solution
  rejection. At that fixed origin, slopes are **not** censored at zero.
* Fusion: production `fuse_table` (inverse-variance, Cochran-Q/DL escalation);
  all 16 enzyme strata survive at every cell.
* The negative-r2 sensitivity below is not the paper default. Reviewer code
  suggestion M1 asked for signed fixed-origin slopes; the separate positive-r2
  gate was evaluated as a sensitivity because many low-depth true gradients have
  adjusted r2 < 0 under over-optimistic window SEs.

## Main decision (used in Table 2 and Fig 2)

Signed fixed-origin estimates, no positive-r2 fusion gate, are the primary arm.
This restores the correct null semantics for the run-out control:

| depth | RUN_OUT log2PTR |
|---:|---:|
| 0.5x | 0.2375 |
| 1x | 0.0178 |
| 2x | 0.0381 |
| 5x | 0.0408 |
| 10x | 0.0260 |

Accuracy remains strong: r vs measured growth rate is
0.9229/0.9116/0.9577/0.9705/0.9701 at 0.5/1/2/5/10x (n = 16). Thus the
positive-only version's apparent 0.94 at 0.5x was partly truncation, but the
method's usable 0.5x performance is not an artefact.

## Positive-r2 sensitivity

`r2gate_sensitivity.tsv` shows the reviewer-suggested gate. It drops one 0.5x
sample and reduces r from 0.9229 to 0.6266; other depths are 0.902/0.957/0.971/
0.970. It also makes the stationary 0.5x estimate worse (0.590). We therefore do
not adopt positive r2 as a fusion gate. A better low-depth uncertainty model is
future work.

## GC-correction default decision (2026-09-18)

We do **not** adopt the residual two-pass GC correction as the default for the
primary 2bRAD benchmark. A controlled HPC diagnostic
(`gc_correction_diagnostic.tsv`) held the reads and estimator fixed and varied
only the correction:

| mode | arm | 0.5x r | 1x r |
|---|---|---:|---:|
| historical single pass | anchors + V-fit | 0.9229 | 0.9116 |
| historical single pass | FracMinHash + V-fit | 0.8831 | 0.9121 |
| no GC correction | anchors + V-fit | 0.9177 | 0.9080 |
| no GC correction | FracMinHash + V-fit | 0.8635 | 0.9027 |
| residual two pass | anchors + V-fit | 0.8485 | 0.8785 |
| residual two pass | FracMinHash + V-fit | 0.9008 | 0.9139 |

The residual pass absorbs part of the position-associated replication gradient
in the motif-structured anchor panel, especially at 0.5–1x. It is retained as an
opt-in diagnostic (`--gc-residual-pass`), not a default. The primary benchmark
uses the historical single-pass GC correction.

After making the single-pass default explicit in code (HPC `review-final`
commit `929f4c2`), a full A/B/E statistics-only rerun confirmed the primary
values: anchors + V-fit gives r =
0.9229/0.9116/0.9577/0.9705/0.9701 at 0.5/1/2/5/10x
(`hpc_singlepass_grid_summary.tsv`; SLURM job `4076237`).

## Macstudio versus HPC

* Completed locally on Mac Studio: refusion of committed per-enzyme windows,
  bootstrap CIs, tables/figures, LOMO panel check, GC-test aggregation, tests and
  code changes.
* Completed on HPC: paired-end count/window production, arm-B merged-25kb and
  adaptive variants, unfriendly simulation, Sun shallow reads, C5 scale.
* Future full-grid change: statistics-only reruns can be done locally from
  `windows.rates.tsv` if archived. Raw FASTQ recounting, sketch construction and
  C5/community-scale sweeps should remain on HPC.
