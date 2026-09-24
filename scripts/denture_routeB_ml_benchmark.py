#!/usr/bin/env python3
"""Repeated CV benchmark for Lim ORPI clean/unclean denture diagnosis.

Class coding: Unclean_Denture is the positive class, so AUROC quantifies
detection of unclean denture.  Abundance and sk2bGrow feature construction are
fit inside every training fold.  PTR missingness is represented by an explicit
indicator and never imputed as zero.
"""
from __future__ import annotations

import argparse
import glob
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.metrics import average_precision_score, balanced_accuracy_score, f1_score, roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler


def read_abundance(path: Path, meta_path: Path) -> tuple[pd.DataFrame, pd.Series]:
    a = pd.read_csv(path, sep="\t")
    a = a.rename(columns={"#Kingdom": "Kingdom"})
    # Keep one row per species; this vendor table has unique Species labels.
    a = a.drop_duplicates("Species").set_index("Species")
    x = a[[c for c in a.columns if c.startswith("SF")]].T.astype(float)
    x.index.name = "sample-id"
    meta = pd.read_csv(meta_path, sep="\t", index_col=0)
    shared = x.index.intersection(meta.index)
    y = meta.loc[shared, "Group"].map({"Clean_Denture": "clean", "Unclean_Denture": "unclean"}).astype("category")
    return x.loc[shared], y


def read_ptr(path: Path) -> pd.DataFrame:
    d = pd.read_csv(path, sep="\t", low_memory=False)
    d["log2_ptr"] = pd.to_numeric(d["log2(PTR)"], errors="coerce")
    d = d[d["log2_ptr"].notna()].copy()
    tax = d["taxonomy"].astype(str) if "taxonomy" in d else pd.Series("", index=d.index)
    extracted = tax.str.extract(r"(s__[^;]+)")
    d["species"] = extracted[0].str.removeprefix("s__").fillna(d["genome"].astype(str))
    return d


def ptr_matrix(summaries: dict[str, Path], samples: list[str], min_prevalence: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    long = pd.concat([read_ptr(p).assign(sample=s) for s, p in summaries.items()], ignore_index=True)
    agg = (
        long.groupby(["sample", "species"], as_index=False)
        .agg(log2_ptr=("log2_ptr", "median"), n_anchors=("n_anchors", "sum"), coverage=("coverage", "median"))
    )
    counts = agg.groupby("species")["sample"].nunique()
    keep = set(counts[counts >= min_prevalence].index)
    agg = agg[agg["species"].isin(keep)]
    values = agg.pivot(index="sample", columns="species", values="log2_ptr")
    support = agg.pivot(index="sample", columns="species", values="n_anchors")
    coverage = agg.pivot(index="sample", columns="species", values="coverage")
    for d in [values, support, coverage]:
        d.index = d.index.astype(str)
    values = values.reindex(samples)
    support = support.reindex(samples).fillna(0)
    coverage = coverage.reindex(samples).fillna(0)
    values.columns = [f"PTR_value__{c}" for c in values.columns]
    support.columns = [f"PTR_support__{c}" for c in support.columns]
    coverage.columns = [f"PTR_coverage__{c}" for c in coverage.columns]
    missing = values.isna().astype(int)
    missing.columns = [c.replace("PTR_value__", "PTR_missing__") for c in values.columns]
    return pd.concat([values, support, coverage, missing], axis=1), values.isna().astype(int)


def preprocess(train: pd.DataFrame, test: pd.DataFrame, min_prev: float) -> tuple[np.ndarray, np.ndarray, list[str]]:
    train = train.loc[:, ~train.columns.duplicated()]
    test = test.loc[:, ~test.columns.duplicated()]
    prevalence = (train != 0).mean(axis=0)
    keep = prevalence >= min_prev
    keep |= train.columns.str.startswith("PTR_missing__")
    train = train.loc[:, keep]
    test = test.loc[:, keep]

    transformed = []
    transformed_test = []
    names = []
    for col in train.columns:
        a, b = train[col].astype(float), test[col].astype(float)
        if col.startswith("abundance__"):
            finite = a[np.isfinite(a) & (a != 0)]
            pseudo = (float(np.nanmin(finite)) / 2) if finite.size else 1e-10
            a, b = np.log1p(np.maximum(a, 0) + pseudo), np.log1p(np.maximum(b, 0) + pseudo)
        elif col.startswith(("PTR_support__", "PTR_coverage__")):
            a, b = np.log1p(np.maximum(a, 0)), np.log1p(np.maximum(b, 0))
        # PTR_value__ is already log2-scaled; PTR_missing__ is binary.
        transformed.append(a.to_numpy())
        transformed_test.append(b.to_numpy())
        names.append(col)
    a = np.column_stack(transformed)
    b = np.column_stack(transformed_test)
    imp = SimpleImputer(strategy="median")
    a = imp.fit_transform(a)
    b = imp.transform(b)
    sc = StandardScaler()
    return sc.fit_transform(a), sc.transform(b), names


def one_cv(raw: pd.DataFrame, y: pd.Series, folds: np.ndarray, min_prev: float, n_trees: int, seed: int):
    positive = "unclean"
    pred = pd.Series(index=y.index, dtype=object)
    prob = pd.Series(index=y.index, dtype=float)
    imps = []
    for fold in sorted(np.unique(folds)):
        test = folds == fold
        xtr, xte, names = preprocess(raw.loc[~test], raw.loc[test], min_prev)
        model = RandomForestClassifier(
            n_estimators=n_trees, class_weight="balanced_subsample", random_state=seed + fold, n_jobs=-1
        )
        model.fit(pd.DataFrame(xtr, index=raw.index[~test], columns=names), y.loc[~test])
        xx = pd.DataFrame(xte, index=raw.index[test], columns=names)
        pred.loc[test] = model.predict(xx)
        prob.loc[test] = model.predict_proba(xx)[:, 1]
        imps.append(pd.Series(model.feature_importances_, index=names))
    yy = (y == positive).astype(int)
    metrics = {
        "auroc": roc_auc_score(yy, prob),
        "auprc": average_precision_score(yy, prob),
        "balanced_accuracy": balanced_accuracy_score(y, pred),
        "macro_f1": f1_score(y, pred, average="macro"),
    }
    importance = pd.concat(imps, axis=1).groupby(level=0).mean().iloc[:, 0].sort_values(ascending=False)
    out = pd.DataFrame({"sample": y.index, "y_true": y.astype(str), "y_pred": pred.astype(str), "prob_unclean": prob})
    return out, importance, metrics


def paired_auc_delta(y: pd.Series, pa: pd.Series, pb: pd.Series, n_perm: int, seed: int) -> tuple[float, float]:
    yy = (y == "unclean").astype(int).to_numpy()
    delta = roc_auc_score(yy, pb) - roc_auc_score(yy, pa)
    d = pb.to_numpy() - pa.to_numpy()
    rng = np.random.default_rng(seed)
    null = []
    for _ in range(n_perm):
        signs = rng.choice([-1.0, 1.0], size=len(d))
        null.append(abs(roc_auc_score(yy, pb + signs * d) - roc_auc_score(yy, pa + signs * d)))
    return float(delta), float((np.sum(np.asarray(null) >= abs(delta)) + 1) / (n_perm + 1))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--abundance", type=Path, required=True)
    ap.add_argument("--metadata", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--summaries-glob")
    ap.add_argument("--min-prevalence", type=float, default=0.1)
    ap.add_argument("--ptr-min-samples", type=int, default=10)
    ap.add_argument("--n-trees", type=int, default=5000)
    ap.add_argument("--repeats", type=int, default=20)
    ap.add_argument("--permutations", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=20260924)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    x0, y = read_abundance(args.abundance, args.metadata)
    arms = {"A_species_abundance": x0.add_prefix("abundance__")}
    if args.summaries_glob:
        summaries = {p.name.removesuffix(".output.tsv"): p for p in map(Path, glob.glob(args.summaries_glob))}
        ptr, missing = ptr_matrix(summaries, list(x0.index), args.ptr_min_samples)
        ptr = ptr.loc[:, ~ptr.columns.duplicated()]
        missing = missing.loc[:, ~missing.columns.duplicated()]
        arms["B_sk2bgrow"] = pd.concat([ptr, missing.add_prefix("dup_missing__")], axis=1)
        arms["C_abundance_plus_sk2bgrow"] = pd.concat([arms["A_species_abundance"], ptr], axis=1)

    all_metrics, predictions, deltas, importances = [], [], [], []
    for rep in range(args.repeats):
        seed = args.seed + rep
        folds = np.empty(len(y), dtype=int)
        cv = StratifiedKFold(10, shuffle=True, random_state=seed)
        for fold, (_, test) in enumerate(cv.split(np.zeros(len(y)), y)):
            folds[test] = fold
        probs = {}
        for arm, raw in arms.items():
            out, imp, metric = one_cv(raw, y, folds, args.min_prevalence, args.n_trees, seed)
            metric.update({"arm": arm, "repeat": rep, "n_features": raw.shape[1]})
            all_metrics.append(metric)
            out["repeat"] = rep
            out["arm"] = arm
            predictions.append(out)
            importances.append(imp.rename(arm).to_frame("importance").assign(repeat=rep, arm=arm))
            probs[arm] = out.set_index("sample")["prob_unclean"]
        if len(arms) > 1:
            for contrast, pb_name in [("B-A", "B_sk2bgrow"), ("C-A", "C_abundance_plus_sk2bgrow")]:
                if pb_name in probs:
                    d, p = paired_auc_delta(y, probs["A_species_abundance"], probs[pb_name], args.permutations, seed)
                    deltas.append({"repeat": rep, "contrast": contrast, "auroc_delta": d, "p_permutation": p})

    pd.DataFrame(all_metrics).to_csv(args.outdir / "cv_metrics_by_repeat.tsv", sep="\t", index=False)
    pd.concat(predictions).to_csv(args.outdir / "cv_predictions.tsv", sep="\t", index=False)
    pd.concat(importances).to_csv(args.outdir / "cv_feature_importances.tsv.gz", sep="\t", index=True, compression="gzip")
    md = pd.DataFrame(all_metrics)
    summary = md.groupby("arm").agg(
        n=("auroc", "size"), auroc_mean=("auroc", "mean"), auroc_sd=("auroc", "std"),
        **{
            "auroc_p2.5": ("auroc", lambda z: np.quantile(z, 0.025)),
            "auroc_p97.5": ("auroc", lambda z: np.quantile(z, 0.975)),
        },
        auprc_mean=("auprc", "mean"), auprc_sd=("auprc", "std"),
        balacc_mean=("balanced_accuracy", "mean"), balacc_sd=("balanced_accuracy", "std"),
        macrof1_mean=("macro_f1", "mean"), macrof1_sd=("macro_f1", "std"),
    )
    summary.to_csv(args.outdir / "cv_metrics_summary.tsv", sep="\t")
    if deltas:
        pd.DataFrame(deltas).to_csv(args.outdir / "paired_auroc_deltas.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
