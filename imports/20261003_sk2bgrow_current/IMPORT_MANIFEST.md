# Imported sk2bGrow current-refresh snapshot

Import date: 2026-10-03

## Purpose

This directory is a safe, non-overwriting import of the current paper-facing
results and manuscript assets that were temporarily separated from the
`sk2bGrow` implementation repository. It is a snapshot: the authoritative
manuscript remains `manuscript/manuscript.md` in the parent repository until a
separate integration review is performed.

## Source and provenance

- Temporary source worktree: `/Users/macstudio/Downloads/sk2bGrow-paper`
- Source paper snapshot commit: `f4e775a472a8337578ad38b8da1e3d320b6b098a`
- `sk2bGrow` implementation commit after split: `e28bbb438f756987247ff3ce13a8d1eb7fb52d15`
- `sk2bGrow` implementation commit before split: `56d0e368659686043625dc8a2e071466763db6c9`
- Destination repository base commit: `ff2ab6752087362e2e0e03f877645c592f2eb062`
- Imported files: 282

## Scope

- `benches/`: benchmark scripts, summaries, logs and small result archives
- `docs/paper/`: refreshed manuscript, supplement, figures, tables and
  Microbiome/BMC submission package
- `docs/PAPER_*.md`, `docs/BENCHMARK_EVALUATION.md`, and
  `docs/PEER_REVIEW_RESPONSE.md`: current claims, provenance and review notes

Large ignored intermediates, genomes, FASTQs, database directories and
`__pycache__` directories were not imported.

## Relationship to the existing manuscript

This snapshot is deliberately isolated under `imports/20261003_sk2bgrow_current`. It does not
replace `manuscript/manuscript.md`, `manuscript/submission/`, `figures/`,
`tables/`, `data/`, or any other existing formal asset. A later editorial pass
should compare the two manuscript versions and selectively integrate the
current-refresh results.
