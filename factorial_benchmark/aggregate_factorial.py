#!/usr/bin/env python3
"""Merge factorial benchmark parts and emit a wide estimator-contrast table."""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", type=Path, required=True)
    ap.add_argument("--long-output", type=Path, required=True)
    ap.add_argument("--contrast-output", type=Path, required=True)
    args = ap.parse_args()

    parts = sorted(args.results_dir.glob("part-*.tsv"))
    if not parts:
        raise SystemExit("no part-*.tsv files found")
    long = pd.concat([pd.read_csv(path, sep="\t") for path in parts], ignore_index=True)
    long.to_csv(args.long_output, sep="\t", index=False, float_format="%.8g")

    keys = [
        "n_strains", "depth", "shared_fraction", "similarity", "abundance",
        "landmark_source",
    ]
    wide = long.pivot(index=keys, columns="estimator")
    wide.columns = [f"{metric}_{est}" for metric, est in wide.columns]
    wide = wide.reset_index()
    wide["bias_sorted_minus_coord"] = wide["bias_sorted_rank"] - wide["bias_coordinate_vfit"]
    wide["rmse_sorted_minus_coord"] = wide["rmse_sorted_rank"] - wide["rmse_coordinate_vfit"]
    wide.to_csv(args.contrast_output, sep="\t", index=False, float_format="%.8g")
    print(f"merged {len(parts)} parts: {len(long)} rows; contrasts: {len(wide)} rows")


if __name__ == "__main__":
    main()
