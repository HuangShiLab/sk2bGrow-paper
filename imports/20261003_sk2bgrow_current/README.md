# sk2bGrow-paper

Independent manuscript and benchmark repository for sk2bGrow.

This repository was split from the sk2bGrow implementation repository at code
commit `72b5bbe1d62fb347232e600f56a0de525842de88`. The implementation lives in `sk2bGrow`; this repository keeps
paper-facing results, figures, manuscript text, validation summaries and the
Microbiome submission package.

## Layout

```
benches/                                  benchmark scripts, summaries and result archives
docs/paper/                               manuscript, supplement, figures, submission package
docs/PAPER_RESULTS.md                     current results and provenance
docs/PAPER_CLAIMS.md                      claim-calibration rules
docs/PEER_REVIEW_RESPONSE.md              internal adversarial review response
docs/BENCHMARK_EVALUATION.md              benchmark significance assessment
```

The two repositories are expected as siblings by some legacy HPC scripts:

```
parent/
  sk2bGrow/
  sk2bGrow-paper/
```

## Primary checks

Run from this repository:

```bash
../sk2bGrow/.venv/bin/python benches/refresh_20260930/validate_manuscript_numbers.py
../sk2bGrow/.venv/bin/python benches/refresh_20260930/make_paper_figures.py
../sk2bGrow/.venv/bin/python benches/refresh_20260930/make_supplementary_tables.py
```

Rebuild the Microbiome package:

```bash
../sk2bGrow/.venv/bin/python docs/paper/build_microbiome_submission.py
```
