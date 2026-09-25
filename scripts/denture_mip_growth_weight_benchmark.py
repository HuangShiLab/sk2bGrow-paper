#!/usr/bin/env python3
"""Evaluate a growth-weighted MIP candidate for Lim ORPI denture data.

Candidate formula (before log transformation):

    MIP_growth = sum_p abundance_p * 2^clip(log2PTR_p, 0, 3)

The clip keeps PTR in [1, 8].  A pathogen without a robust PTR receives
growth-neutral weight 1 in `all_missing1`; the observed-only version omits its
abundance contribution.  The pathogen set is an explicit, reproducible oral/
opportunistic list used only for this sensitivity analysis; it is not a claim
that every member is a confirmed pathogen.
"""
from __future__ import annotations

import argparse
import glob
import importlib.util
import os
from pathlib import Path

import numpy as np
import pandas as pd


PATHOGEN_PATTERNS = [
    "Porphyromonas_gingivalis", "Tannerella_forsythia", "Treponema_.*_denticola",
    "Filifactor_alocis", "Fusobacterium_nucleatum", "Prevotella_intermedia",
    "Prevotella_nigrescens", "Eikenella_corrodens",
    "Aggregatibacter_actinomycetemcomitans", "Selenomonas_sputigena",
    "Streptococcus_mutans", "Streptococcus_sobrinus", "Streptococcus_anginosus",
    "Bifidobacterium_dentium", "Scardovia_wiggsiae", "Parascardovia_denticolens",
    "Actinomyces_israelii", "Actinomyces_dentalis", "Granulicatella_adiacens",
    "Rothia_mucilaginosa", "Rothia_aeria", "Capnocytophaga",
    "Enterobacter_hormaechei", "Klebsiella", "Pseudomonas_aeruginosa",
    "Staphylococcus", "Enterococcus", "Peptidiphaga_gingivicola", "Kingella_oralis",
]


def load_benchmark_module(path: Path):
    spec = importlib.util.spec_from_file_location("denture_benchmark", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser()
    ap.add_argument("--abundance", default=root / "data/lim_orpi_denture/external/Abundance_Stat.all.xls")
    ap.add_argument("--metadata", default=root / "data/lim_orpi_denture/external/meta.tsv")
    ap.add_argument("--summaries-glob", default=root / "data/lim_orpi_denture/sk2bgrow_sp/summaries_glm/*.output.tsv")
    ap.add_argument("--benchmark-script", default=root / "scripts/denture_routeB_ml_benchmark.py")
    ap.add_argument("--outdir", default=root / "analysis/lim_orpi_denture_routeB")
    ap.add_argument("--repeats", type=int, default=20)
    ap.add_argument("--n-trees", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=20260925)
    args = ap.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    bench = load_benchmark_module(Path(args.benchmark_script))
    abundance, y = bench.read_abundance(args.abundance, args.metadata)
    a = pd.read_csv(args.abundance, sep="\t").drop_duplicates("Species").set_index("Species")
    pathogen_abundance = a[a.index.str.contains("|".join(PATHOGEN_PATTERNS), case=True, regex=True)]
    pathogen_abundance = pathogen_abundance[list(abundance.index)].T.fillna(0)

    frames = []
    for f in map(Path, glob.glob(str(args.summaries_glob))):
        d = pd.read_csv(f, sep="\t", low_memory=False)
        d["sample_id"] = os.path.basename(f).removesuffix(".output.tsv")
        frames.append(d)
    ptr = pd.concat(frames, ignore_index=True)
    for col in ["PTR", "log2(PTR)", "n_anchors", "n_windows", "ori_confidence", "se"]:
        ptr[col] = pd.to_numeric(ptr[col], errors="coerce")
    ptr = ptr[
        ptr["log2(PTR)"].between(0, 3)
        & (ptr["n_windows"] >= 5)
        & (ptr["n_anchors"] >= 50)
        & (ptr["ori_confidence"] >= 0.5)
        & (ptr["se"] <= 1)
    ].copy()
    ptr["species"] = ptr["taxonomy"].astype(str).str.extract(r"(s__[^;]+)")[0].str.replace("^s__", "", regex=True)
    ptr_matrix = ptr.pivot(index="sample_id", columns="species", values="PTR").reindex(
        index=abundance.index, columns=pathogen_abundance.columns
    )

    growth_weight = ptr_matrix.clip(1, 8)
    abundance_matrix = pathogen_abundance.reindex(columns=ptr_matrix.columns).fillna(0)
    candidates = pd.DataFrame(index=abundance.index)
    candidates["MIP_log_pathogen_abundance"] = np.log1p(abundance_matrix.sum(axis=1))
    candidates["MIP_log_growth_weighted_all_missing1"] = np.log1p(
        (abundance_matrix * growth_weight.fillna(1)).sum(axis=1)
    )
    candidates["MIP_log_growth_weighted_observed"] = np.log1p(
        (abundance_matrix * growth_weight).sum(axis=1)
    )
    candidates["MIP_robust_pathogen_count"] = ptr_matrix.notna().sum(axis=1)
    candidates.to_csv(outdir / "mip_per_sample_candidate.tsv", sep="\t")

    abundance_features = abundance.add_prefix("abundance__")
    arms = {
        "A_abundance": abundance_features,
        "A_plus_growth_MIP": pd.concat([abundance_features, candidates], axis=1),
        "growth_MIP_only": candidates,
    }
    all_metrics, all_predictions = [], []
    for rep in range(args.repeats):
        seed = args.seed + rep
        folds = np.empty(len(y), dtype=int)
        cv = bench.StratifiedKFold(10, shuffle=True, random_state=seed)
        for fold, (_, test) in enumerate(cv.split(np.zeros(len(y)), y)):
            folds[test] = fold
        for arm, raw in arms.items():
            pred, _, metric = bench.one_cv(raw, y, folds, 0.1, args.n_trees, seed)
            metric.update({"arm": arm, "repeat": rep, "n_features": raw.shape[1]})
            all_metrics.append(metric)
            pred["repeat"] = rep
            pred["arm"] = arm
            all_predictions.append(pred)

    metrics = pd.DataFrame(all_metrics)
    predictions = pd.concat(all_predictions)
    summary = metrics.groupby("arm").agg(
        n=("auroc", "size"), auroc_mean=("auroc", "mean"), auroc_sd=("auroc", "std"),
        auprc_mean=("auprc", "mean"), auprc_sd=("auprc", "std"),
        balacc_mean=("balanced_accuracy", "mean"), balacc_sd=("balanced_accuracy", "std"),
        macrof1_mean=("macro_f1", "mean"), macrof1_sd=("macro_f1", "std"),
    )
    metrics.to_csv(outdir / "mip_cv_metrics_by_repeat.tsv", sep="\t", index=False)
    predictions.to_csv(outdir / "mip_cv_predictions.tsv", sep="\t", index=False)
    summary.to_csv(outdir / "mip_cv_metrics_summary.tsv", sep="\t")


if __name__ == "__main__":
    main()
