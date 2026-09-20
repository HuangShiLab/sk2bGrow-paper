#!/usr/bin/env python3
"""Count-level factorial simulator for coordinate and sorted PTR estimators.

The simulator deliberately does not write FASTQ.  It generates expected
landmark counts under a V-shaped replication profile, models shared anchors as
ambiguous observations, aggregates assigned counts to fixed windows, and then
applies two estimators to the same count table.  The output is aggregate only.
"""
from __future__ import annotations

import argparse
import itertools
import math
import os
import time
from pathlib import Path

import numpy as np
import pandas as pd

from sk2bgrow.fit import fit_sorted_ransac, fit_v_shape


N_LANDMARKS = 5_000
GENOME_LEN = 3_000_000
WINDOW_BP = 5_000
READ_LEN = 150
N_WINDOWS = GENOME_LEN // WINDOW_BP

N_STRAINS = (4, 8, 16, 32)
DEPTHS = (0.5, 1.0, 2.0, 4.0, 8.0)
SHARED = ((0.00, "low"), (0.05, "medium"), (0.15, "high"))
ABUNDANCE = ("even", "10:1", "100:1")
LANDMARKS = ("coordinate_grid", "random_sketch")


def abundance_weights(label: str, m: int, rng: np.random.Generator) -> np.ndarray:
    """Return ratio weights with geometric mean one, then randomize rank."""
    if label == "even":
        raw = np.ones(m)
    elif label == "10:1":
        raw = np.linspace(1.0, 10.0, m)
    elif label == "100:1":
        raw = np.linspace(1.0, 100.0, m)
    else:
        raise ValueError(label)
    geo_mean = float(np.exp(np.mean(np.log(raw))))
    return rng.permutation(raw / geo_mean)


def landmark_coordinates(source: str, n_shared: int, rng: np.random.Generator):
    """Return common coordinate scaffold split into shared/private identities.

    Coordinates are held constant across strains. Identity, rather than
    coordinate, determines whether an observation is assignable: private
    landmarks at the same coordinate are distinct across strains, while shared
    landmarks represent the same sequence/hash identity in every strain.
    """
    n_private = N_LANDMARKS - n_shared
    if source == "coordinate_grid":
        coords = (np.arange(N_LANDMARKS) + 0.5) * GENOME_LEN / N_LANDMARKS
        shared_idx = np.linspace(0, N_LANDMARKS - 1, n_shared, dtype=int)
        shared = coords[shared_idx]
        private = np.delete(coords, shared_idx)
    elif source == "random_sketch":
        coords = np.sort(rng.choice(GENOME_LEN, size=N_LANDMARKS, replace=False)).astype(float)
        shared = coords[:n_shared]
        private = coords[n_shared:]
    else:
        raise ValueError(source)
    return shared.astype(float), private.astype(float)


def v_profile(position: np.ndarray, log2_ptr: float) -> np.ndarray:
    """Pilea-style linear-in-log2 decrease from ori=0 to ter=L/2."""
    d = np.minimum(position, GENOME_LEN - position)
    relative = d / (GENOME_LEN / 2.0)
    weights = np.exp2(-float(log2_ptr) * relative)
    return weights / weights.sum()


def assigned_window_counts(
    m: int,
    depth: float,
    shared_fraction: float,
    source: str,
    log2_ptr: np.ndarray,
    weights: np.ndarray,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    """Simulate one community and return assigned counts plus landmark exposure."""
    n_shared = int(round(shared_fraction * N_LANDMARKS))
    counts = np.zeros((m, N_WINDOWS), dtype=np.int32)
    strain_depth = depth * weights

    shared_coords, private_coords = landmark_coordinates(source, n_shared, rng)
    all_coords = np.concatenate([shared_coords, private_coords])
    landmark_exposure = np.bincount(
        np.floor(all_coords / WINDOW_BP).astype(int) % N_WINDOWS,
        minlength=N_WINDOWS,
    ).astype(float)
    shared_contrib: list[np.ndarray] = []

    # Private identities are source-specific and therefore assigned directly.
    for i in range(m):
        p_private = v_profile(private_coords, log2_ptr[i])
        p_shared = v_profile(shared_coords, log2_ptr[i])
        private_lambda = strain_depth[i] * GENOME_LEN / READ_LEN * (1.0 - shared_fraction) * p_private
        shared_lambda = strain_depth[i] * GENOME_LEN / READ_LEN * shared_fraction * p_shared
        obs = rng.poisson(private_lambda)
        win = np.floor(private_coords / WINDOW_BP).astype(int) % N_WINDOWS
        np.add.at(counts[i], win, obs)
        shared_contrib.append(rng.poisson(shared_lambda))

    # Shared identities are ambiguous. The total observation is assigned to
    # one strain with probability proportional to each strain's latent expected
    # contribution. This is a conservative model of k-mer/anchor confusion.
    if n_shared:
        contrib = np.vstack(shared_contrib)
        totals = contrib.sum(axis=0)
        win = np.floor(shared_coords / WINDOW_BP).astype(int) % N_WINDOWS
        for j in np.nonzero(totals)[0]:
            probs = contrib[:, j].astype(float)
            source_i = rng.choice(m, p=probs / probs.sum())
            counts[source_i, win[j]] += int(totals[j])

    return counts, landmark_exposure


def window_log2_rates(
    counts: np.ndarray, landmark_exposure: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Normalize each window by its number of sampled landmarks.

    Without this exposure term, a random sketch would appear worse only because
    sparse windows have larger Poisson noise. That is a density artifact, not a
    landmark-source effect.
    """
    per_landmark = max(float(counts.sum()) / N_LANDMARKS, 1e-6)
    expected = landmark_exposure * per_landmark
    rate = np.log2((counts + 0.5) / (expected + 0.5))
    # Delta-method standard error for log2 of a Poisson count.  Coordinate
    # V-fit needs these weights when windows have unequal landmark exposure;
    # otherwise sparse random-sketch windows create a Jensen-bias artifact.
    se = 1.0 / (math.log(2.0) * np.sqrt(np.maximum(expected, 1.0)))
    return rate, se


def paired_delta(metric_by_est: dict[str, list[float]], estimator: str, metric: str) -> float:
    """Cell-level estimator-minus-coordinate difference."""
    coord = np.concatenate([rep[metric] for rep in metric_by_est["coordinate_vfit"]]).astype(float)
    other = np.concatenate([rep[metric] for rep in metric_by_est[estimator]]).astype(float)
    return float(np.mean(other - coord))


def score_cell(rows: list[dict], factor_row: dict, reps: list[dict]) -> None:
    by_est: dict[str, list[dict]] = {"coordinate_vfit": [], "sorted_rank": []}
    for rep in reps:
        for est, est_rep in rep["estimators"].items():
            by_est[est].append(est_rep)

    for estimator, vals in by_est.items():
        err = np.concatenate([x["error"] for x in vals]).astype(float)
        ok = np.isfinite(err)
        err = err[ok]
        truth = np.concatenate([x["true"] for x in vals]).astype(float)[ok]
        estimate = np.concatenate([x["estimate"] for x in vals]).astype(float)[ok]
        corr = np.corrcoef(truth, estimate)[0, 1] if len(err) > 2 and np.std(truth) > 0 else np.nan
        row = dict(factor_row)
        row.update(
            estimator=estimator,
            n=int(len(np.concatenate([x["error"] for x in vals]))),
            n_ok=int(ok.sum()),
            reported_fraction=float(ok.mean()),
            bias=float(np.mean(err)) if len(err) else np.nan,
            rmse=float(np.sqrt(np.mean(err**2))) if len(err) else np.nan,
            mae=float(np.mean(np.abs(err))) if len(err) else np.nan,
            correlation=float(corr),
            mean_true_log2ptr=float(np.mean(truth)) if len(truth) else np.nan,
            mean_est_log2ptr=float(np.mean(estimate)) if len(estimate) else np.nan,
            delta_bias_vs_coord=paired_delta(by_est, estimator, "error"),
            delta_rmse_vs_coord=paired_delta(by_est, estimator, "abs_error"),
        )
        rows.append(row)


def run_task(task_id: int, reps: int, seed_base: int) -> pd.DataFrame:
    partition = list(itertools.product(N_STRAINS, ABUNDANCE, LANDMARKS))
    n_strains, abundance, landmark = partition[task_id]
    rows: list[dict] = []
    all_reps: list[dict] = []
    started = time.time()

    for depth, (shared_fraction, shared_label), rep in itertools.product(
        DEPTHS, SHARED, range(reps)
    ):
        if rep == 0:
            all_reps = []
        seed = seed_base + task_id * 1_000_003 + int(depth * 10) * 10_009 + int(shared_fraction * 1000) * 101 + rep
        rng = np.random.default_rng(seed)
        m = n_strains
        weights = abundance_weights(abundance, m, rng)
        truth = rng.uniform(0.0, 2.0, size=m)
        counts, landmark_exposure = assigned_window_counts(
            m, depth, shared_fraction, landmark, truth, weights, rng
        )
        rate_pairs = [
            window_log2_rates(counts[i], landmark_exposure) for i in range(m)
        ]
        rates = np.vstack([pair[0] for pair in rate_pairs])
        window_ses = np.vstack([pair[1] for pair in rate_pairs])
        positions = (np.arange(N_WINDOWS) + 0.5) * WINDOW_BP
        rep_out: dict[str, dict] = {}

        for estimator in ("coordinate_vfit", "sorted_rank"):
            estimates = np.full(m, np.nan)
            for i in range(m):
                y = rates[i]
                if np.isfinite(y).sum() < 5:
                    continue
                if estimator == "coordinate_vfit":
                    fit = fit_v_shape(
                        positions,
                        y,
                        GENOME_LEN,
                        se=window_ses[i],
                        ori=0.0,
                        allow_segmented=False,
                        refine=False,
                    )
                else:
                    fit = fit_sorted_ransac(y, seed=int(seed % 2**31) + i)
                estimates[i] = fit.log2_ptr

            rep_out[estimator] = {
                "true": truth,
                "estimate": estimates,
                "error": estimates - truth,
                "abs_error": np.abs(estimates - truth),
            }

        all_reps.append({"estimators": rep_out})

        if rep == reps - 1:
            score_cell(
                rows,
                {
                    "n_strains": m,
                    "depth": depth,
                    "shared_fraction": shared_fraction,
                    "similarity": shared_label,
                    "abundance": abundance,
                    "landmark_source": landmark,
                },
                all_reps,
            )
    runtime = time.time() - started
    for row in rows:
        row["runtime_sec"] = runtime
    agg = pd.DataFrame(rows)
    return agg


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--task-id", type=int, required=True)
    ap.add_argument("--n-tasks", type=int, default=24)
    ap.add_argument("--reps", type=int, default=10)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--seed-base", type=int, default=20260920)
    args = ap.parse_args()
    if not 0 <= args.task_id < args.n_tasks:
        raise ValueError("task id outside partition")
    if args.n_tasks != len(N_STRAINS) * len(ABUNDANCE) * len(LANDMARKS):
        raise ValueError("n-tasks must be 24")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    frame = run_task(args.task_id, args.reps, args.seed_base)
    out = args.out_dir / f"part-{args.task_id:03d}.tsv"
    frame.to_csv(out, sep="\t", index=False, float_format="%.8g")
    print(f"wrote {out} with {len(frame)} aggregate rows")


if __name__ == "__main__":
    main()
