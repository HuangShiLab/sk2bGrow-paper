#!/usr/bin/env python3
"""Nested-free, leakage-safe CV benchmark for ECC diagnosis.

The published baseline is the BcgI 2bRAD-M species-abundance table used for
Fig. 7c-f of Ruan et al. (2025).  Feature construction is deliberately simple
and is applied inside every training fold:

1. keep features present (non-zero) in at least `--min-prevalence` of the
   training samples;
2. log-transform proportions with a half-minimum positive pseudocount;
3. standardize continuous features.

For sk2bGrow PTR features, missing values are represented by an explicit
`PTR_missing` indicator, never by zero imputation.  The same folds are used for
every feature arm.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    f1_score,
    roc_auc_score,
)
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


def read_published(path: Path) -> tuple[pd.DataFrame, pd.Series]:
    x = pd.read_csv(path, sep="\t", index_col=0)
    meta_path = path.with_name("published_2b_metadata.tsv")
    meta = pd.read_csv(meta_path, sep="\t", index_col=0)
    shared = x.index.intersection(meta.index)
    x = x.loc[shared]
    y = meta.loc[shared, "status"].map({"H": "healthy", "Health": "healthy", "ECC": "ecc"}).astype("category")
    return x, y


def read_ptr_summary(path: Path) -> pd.DataFrame:
    d = pd.read_csv(path, sep="\t", low_memory=False)
    # Route-A QC is intentionally not used as a feature filter.  The pilot showed
    # that a gut-biased UHGG reference can yield finite but QC-failed PTRs.
    # Feature prevalence and internal CV guard against unstable features.
    d["genome"] = d["genome"].astype(str)
    for col in ["PTR", "log2(PTR)", "coverage", "n_anchors", "n_windows", "se", "ori_confidence"]:
        if col not in d:
            d[col] = np.nan
    d = d[pd.to_numeric(d["log2(PTR)"], errors="coerce").notna()].copy()
    return d


def make_ptr_features(summaries: dict[str, Path], samples: list[str], max_species_per_sample: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return PTR/log2PTR/support features and matching missing indicators.

    A genome must be QC-passing in at least `max_species_per_sample` samples to
    enter the feature matrix.  This is a population-level prevalence filter, not
    a class-aware selection.  Sample-specific values remain missing when absent.
    """
    frames = []
    for sample, path in summaries.items():
        d = read_ptr_summary(path)
        d["sample"] = sample
        frames.append(d)
    long = pd.concat(frames, ignore_index=True)
    if "taxonomy" in long.columns:
        extracted = long["taxonomy"].astype(str).str.extract(r"(s__[^;]+)")
        long["species"] = extracted[0].fillna(long["genome"])
    else:
        long["species"] = long["genome"]

    long = (
        long.groupby(["sample", "species"], as_index=False)
        .agg({
            "log2(PTR)": "median",
            "n_anchors": "sum",
            "coverage": "median",
        })
        .rename(columns={"species": "genome"})
    )
    counts = long.groupby("genome")["sample"].nunique()
    keep = counts[counts >= max_species_per_sample].index
    long = long[long["genome"].isin(keep)].copy()

    values = long.pivot(index="sample", columns="genome", values="log2(PTR)")
    support = long.pivot(index="sample", columns="genome", values="n_anchors")
    coverage = long.pivot(index="sample", columns="genome", values="coverage")
    support.columns = [f"PTR_support__{c}" for c in support.columns]
    coverage.columns = [f"PTR_coverage__{c}" for c in coverage.columns]
    values.index = values.index.astype(str)
    support.index = support.index.astype(str)
    coverage.index = coverage.index.astype(str)
    support = support.reindex(samples).fillna(0)
    coverage = coverage.reindex(samples).fillna(0)
    values = values.reindex(samples)
    values = values.loc[:, ~values.columns.duplicated()]
    support = support.loc[:, ~support.columns.duplicated()]
    coverage = coverage.loc[:, ~coverage.columns.duplicated()]
    missing = values.isna().astype(int)
    missing = missing.loc[:, ~missing.columns.duplicated()]
    return pd.concat([values, support, coverage], axis=1), missing


def feature_preprocess(x_train: pd.DataFrame, x_test: pd.DataFrame, min_prev: float, rng_seed: int) -> tuple[np.ndarray, np.ndarray, list[str]]:
    prevalence = (x_train != 0).mean(axis=0)
    keep = prevalence >= min_prev
    if isinstance(x_train.iloc[0, 0], (bool, int, float)):
        keep |= x_train.columns.str.startswith(("PTR_support__", "PTR_coverage__", "missing__"))
    x_train = x_train.loc[:, keep]
    x_test = x_test.loc[:, keep]

    # PTR values and support/coverage are continuous, but have different scales;
    # all are standardized after median imputation.
    positive = x_train.to_numpy(dtype=float)
    finite = positive[np.isfinite(positive) & (positive != 0)]
    pseudo = float(np.nanmin(finite)) / 2 if finite.size else 1e-8
    train_values = np.where(np.isfinite(positive), positive, np.nan)
    test_values = np.where(np.isfinite(x_test.to_numpy(dtype=float)), x_test.to_numpy(dtype=float), np.nan)
    # For abundance features, missing never occurs. For PTR features, explicit
    # indicator columns were added upstream; median imputation only handles
    # sporadic missingness in retained support features.
    imp = SimpleImputer(strategy="median")
    cols = list(x_train.columns)
    train_imp = imp.fit_transform(train_values)
    test_imp = imp.transform(test_values)
    if train_imp.shape[1] != len(cols):
        # SimpleImputer drops all-missing columns; the matching missing
        # indicators preserve their information downstream.
        keep_positions = [j for j, name in enumerate(cols) if not name.startswith("missing__")]
        cols = [cols[j] for j in keep_positions]
        train_imp = train_imp[:, keep_positions]
        test_imp = test_imp[:, keep_positions]
    log_train = np.sign(train_imp) * np.log1p(np.abs(train_imp))
    log_test = np.sign(test_imp) * np.log1p(np.abs(test_imp))
    sc = StandardScaler()
    return sc.fit_transform(log_train), sc.transform(log_test), cols


def permutation_delta(y: pd.Series, prob_a: pd.Series, prob_b: pd.Series, n_perm: int, seed: int) -> tuple[float, float]:
    delta = roc_auc_score(y, prob_b) - roc_auc_score(y, prob_a)
    rng = np.random.default_rng(seed)
    obs = prob_b.to_numpy() - prob_a.to_numpy()
    null = []
    yy = (y == "ecc").astype(int).to_numpy()
    for _ in range(n_perm):
        signs = rng.choice([-1.0, 1.0], size=len(obs))
        null.append(abs(roc_auc_score(yy, prob_b.to_numpy() + signs * obs) - roc_auc_score(yy, prob_a.to_numpy() + signs * obs)))
    return float(delta), float((np.sum(np.abs(null) >= abs(delta)) + 1) / (n_perm + 1))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--published", type=Path, default=Path("data/ecc_routeB/external/published_2b_abundance.tsv"))
    ap.add_argument("--outdir", type=Path, default=Path("analysis/ecc_routeB"))
    ap.add_argument("--summaries-glob", default=None, help="optional glob of sk2bGrow output.tsv files")
    ap.add_argument("--min-prevalence", type=float, default=0.1)
    ap.add_argument("--ptr-min-prevalence", type=int, default=3)
    ap.add_argument("--n-trees", type=int, default=5000)
    ap.add_argument("--repeats", type=int, default=20)
    ap.add_argument("--permutations", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=20260923)
    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    x0, y = read_published(args.published)
    if args.summaries_glob:
        import glob
        summaries = {}
        for p in map(Path, glob.glob(args.summaries_glob)):
            sample = p.name.removesuffix(".output.tsv")
            summaries[sample] = p
        ptr, missing = make_ptr_features(summaries, list(x0.index), args.ptr_min_prevalence)
        arms = {
            "A_published": x0,
            "B_sk2bgrow": pd.concat([ptr, missing.add_prefix("missing__")], axis=1),
            "C_published_PTR": pd.concat([x0, ptr, missing.add_prefix("missing__")], axis=1),
        }
    else:
        arms = {"A_published": x0}

    cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=args.seed)
    all_metrics = []
    predictions = []
    deltas = []
    for rep in range(args.repeats):
        seed = args.seed + rep
        cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=seed)
        fold_array = np.empty(len(y), dtype=int)
        for fold, (_, test) in enumerate(cv.split(np.zeros(len(y)), y)):
            fold_array[test] = fold
        rep_probs = {}
        for arm, raw in arms.items():
            # Matrix placeholder; preprocessing is fold-specific below.
            pred, imp, metric = _cv_one(raw, y, fold_array, args.min_prevalence, args.n_trees, seed)
            metric.update({"arm": arm, "repeat": rep, "n_features": raw.shape[1]})
            all_metrics.append(metric)
            pred["repeat"] = rep
            pred["arm"] = arm
            predictions.append(pred)
            rep_probs[arm] = pred.set_index("sample")["prob_ecc"]
        if "B_sk2bgrow" in rep_probs:
            d1, p1 = permutation_delta(y, rep_probs["A_published"], rep_probs["B_sk2bgrow"], args.permutations, seed)
            d2, p2 = permutation_delta(y, rep_probs["A_published"], rep_probs["C_published_PTR"], args.permutations, seed)
            deltas.extend([{"repeat": rep, "contrast": "B-A", "auroc_delta": d1, "p_permutation": p1}, {"repeat": rep, "contrast": "C-A", "auroc_delta": d2, "p_permutation": p2}])
    pd.DataFrame(all_metrics).to_csv(args.outdir / "cv_metrics_by_repeat.tsv", sep="\t", index=False)
    pd.concat(predictions).to_csv(args.outdir / "cv_predictions.tsv", sep="\t", index=False)
    metrics_df = pd.DataFrame(all_metrics)
    summary = metrics_df.groupby("arm").agg(
        n=("auroc", "size"), auroc_mean=("auroc", "mean"), auroc_sd=("auroc", "std"),
        auprc_mean=("auprc", "mean"), auprc_sd=("auprc", "std"),
        balacc_mean=("balanced_accuracy", "mean"), balacc_sd=("balanced_accuracy", "std"),
        macrof1_mean=("macro_f1", "mean"), macrof1_sd=("macro_f1", "std"),
    )
    summary.to_csv(args.outdir / "cv_metrics_summary.tsv", sep="\t")
    if deltas:
        pd.DataFrame(deltas).to_csv(args.outdir / "paired_auroc_deltas.tsv", sep="\t", index=False)


def _cv_one(raw: pd.DataFrame, y: pd.Series, folds: np.ndarray, min_prev: float, n_trees: int, seed: int):
    raw = raw.loc[:, ~raw.columns.duplicated()]
    classes = list(y.cat.categories)
    pos = classes.index("ecc")
    pred = pd.Series(index=y.index, dtype=object)
    prob = pd.Series(index=y.index, dtype=float)
    importances = []
    for fold in sorted(np.unique(folds)):
        test = folds == fold
        train = ~test
        x_train_num, x_test_num, cols = feature_preprocess(raw.loc[train], raw.loc[test], min_prev, seed)
        x_train = pd.DataFrame(x_train_num, index=raw.index[train], columns=cols)
        x_test = pd.DataFrame(x_test_num, index=raw.index[test], columns=cols)
        model = RandomForestClassifier(n_estimators=n_trees, class_weight="balanced_subsample", random_state=seed + fold, n_jobs=-1)
        model.fit(x_train, y.loc[train])
        pred.loc[test] = model.predict(x_test)
        prob.loc[test] = model.predict_proba(x_test)[:, pos]
        importances.append(pd.Series(model.feature_importances_, index=cols))
    yy = (y == "ecc").astype(int)
    metrics = {
        "auroc": roc_auc_score(yy, prob),
        "auprc": average_precision_score(yy, prob),
        "balanced_accuracy": balanced_accuracy_score(y, pred),
        "macro_f1": f1_score(y, pred, average="macro"),
    }
    imp_series = pd.concat(importances, axis=1)
    imp = imp_series.groupby(level=0).mean().iloc[:, 0].sort_values(ascending=False)
    out = pd.DataFrame({"sample": y.index, "y_true": y.astype(str), "y_pred": pred.astype(str), "prob_ecc": prob})
    return out, imp, metrics


if __name__ == "__main__":
    main()
