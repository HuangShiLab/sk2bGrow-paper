#!/usr/bin/env python3
"""Summarize robust sk2bGrow growth activity in the Lim ORPI denture data.

The estimates are exploratory because none passed route-A-style containment QC.
This summary therefore applies strict support/plausibility filters before
ranking species:

* finite log2 PTR in the biologically plausible interval [0, 3] (PTR 1--8);
* >=5 modelled windows and >=50 anchors;
* ori confidence >=0.5;
* PTR SE <=1 on log2 scale.

A species must have robust estimates in at least `--min-samples` samples to be
ranked. Group differences are descriptive (Mann-Whitney U + BH FDR), not a
causal or diagnostic claim.
"""
from __future__ import annotations

import argparse
import glob
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu


PATHOGEN_PATTERNS = [
    "Porphyromonas_gingivalis", "Tannerella_forsythia", "Treponema_.*_denticola",
    "Filifactor_alocis", "Fusobacterium_nucleatum", "Prevotella_intermedia",
    "Prevotella_nigrescens", "Eikenella_corrodens", "Aggregatibacter_actinomycetemcomitans",
    "Selenomonas_sputigena", "Streptococcus_mutans", "Streptococcus_sobrinus",
    "Streptococcus_anginosus", "Bifidobacterium_dentium", "Scardovia_wiggsiae",
    "Parascardovia_denticolens", "Actinomyces_israelii", "Actinomyces_dentalis",
    "Granulicatella_adiacens", "Rothia_mucilaginosa", "Rothia_aeria",
    "Capnocytophaga", "Enterobacter_hormaechei", "Klebsiella", "Pseudomonas_aeruginosa",
    "Staphylococcus", "Enterococcus", "Peptidiphaga_gingivicola", "Kingella_oralis"
]


def read_one(path: Path) -> pd.DataFrame:
    d = pd.read_csv(path, sep="\t", low_memory=False)
    d["sample_id"] = os.path.basename(path).removesuffix(".output.tsv")
    for c in ["PTR", "log2(PTR)", "containment", "n_anchors", "n_windows", "ori_confidence", "se", "coverage"]:
        d[c + "__num"] = pd.to_numeric(d[c], errors="coerce")
    return d


def bh_adjust(p: pd.Series) -> pd.Series:
    p = p.astype(float)
    n = len(p)
    order = p.sort_values().index
    ranked = p.loc[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked.values[::-1])[::-1]
    adj = pd.Series(np.minimum(ranked, 1), index=order).sort_index()
    return adj


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--summaries-glob", default="data/lim_orpi_denture/sk2bgrow_sp/summaries_glm/*.output.tsv")
    ap.add_argument("--metadata", default="data/lim_orpi_denture/external/meta.tsv")
    ap.add_argument("--abundance", default="data/lim_orpi_denture/external/Abundance_Stat.all.xls")
    ap.add_argument("--outdir", default="analysis/lim_orpi_denture_routeB/growth_activity")
    ap.add_argument("--min-samples", type=int, default=5)
    args = ap.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    x = pd.concat([read_one(Path(f)) for f in glob.glob(args.summaries_glob)], ignore_index=True)
    x["species"] = x["taxonomy"].astype(str).str.extract(r"(s__[^;]+)")[0].str.replace("^s__", "", regex=True)
    qual = x[
        x["log2(PTR)__num"].between(0, 3, inclusive="both")
        & (x["n_windows__num"] >= 5)
        & (x["n_anchors__num"] >= 50)
        & (x["ori_confidence__num"] >= 0.5)
        & (x["se__num"] <= 1)
    ].copy()
    qual["clean"] = qual["sample_id"].map(
        pd.read_csv(args.metadata, sep="\t", index_col=0)["Group"].to_dict()
    )

    def q25(s):
        return s.quantile(0.25)

    def q75(s):
        return s.quantile(0.75)

    g = qual.groupby("species").agg(
        n_samples=("sample_id", "nunique"), n_observations=("sample_id", "size"),
        median_log2PTR=("log2(PTR)__num", "median"), q1_log2PTR=("log2(PTR)__num", q25),
        q3_log2PTR=("log2(PTR)__num", q75), median_PTR=("PTR__num", "median"),
        median_coverage=("coverage__num", "median"), median_containment=("containment__num", "median"),
        median_windows=("n_windows__num", "median"), median_anchors=("n_anchors__num", "median"),
        median_ori_confidence=("ori_confidence__num", "median"),
    )

    group_med = qual.groupby(["species", "clean"])["log2(PTR)__num"].median().unstack()
    for col in ["Clean_Denture", "Unclean_Denture"]:
        if col in group_med:
            g[f"{col.lower()}_median_log2PTR"] = group_med[col]
    g["unclean_minus_clean_log2PTR"] = g.get("unclean_denture_median_log2PTR", np.nan) - g.get("clean_denture_median_log2PTR", np.nan)

    pvals = []
    for sp, d in qual.groupby("species"):
        a = d.loc[d["clean"] == "Clean_Denture", "log2(PTR)__num"]
        b = d.loc[d["clean"] == "Unclean_Denture", "log2(PTR)__num"]
        if len(a) >= 3 and len(b) >= 3:
            pvals.append((sp, mannwhitneyu(b, a, alternative="two-sided").pvalue))
    if pvals:
        pm = pd.DataFrame(pvals, columns=["species", "mannwhitney_p"]).set_index("species")
        pm["mannwhitney_fdr_bh"] = bh_adjust(pm["mannwhitney_p"])
        g = g.join(pm)

    # Add abundance context from the vendor species table.
    a = pd.read_csv(args.abundance, sep="\t").drop_duplicates("Species").set_index("Species")
    sample_cols = [c for c in a.columns if c.startswith("SF")]
    abd_median_all = a[sample_cols].median(axis=1)
    abd_median_nonzero = a[sample_cols].replace(0, np.nan).median(axis=1)
    g["median_relative_abundance_all_samples"] = g.index.map(abd_median_all)
    g["median_nonzero_relative_abundance"] = g.index.map(abd_median_nonzero)
    g = g.sort_values(["median_log2PTR", "n_samples"], ascending=[False, False])
    g.to_csv(outdir / "species_growth_activity_quality_filtered.tsv", sep="\t")

    top = g[g["n_samples"] >= args.min_samples].copy()
    top.to_csv(outdir / f"species_growth_activity_n{args.min_samples}.tsv", sep="\t")
    top.head(30).to_csv(outdir / f"top30_growth_activity_n{args.min_samples}.tsv", sep="\t")

    path = top[top.index.str.contains("|".join(PATHOGEN_PATTERNS), case=False, regex=True)].copy()
    path.insert(0, "pathogen_category", "oral/opportunistic/pathogen-associated")
    path.to_csv(outdir / f"pathogen_associated_growth_activity_n{args.min_samples}.tsv", sep="\t")

    with open(outdir / "README.txt", "w") as w:
        w.write(
            "Quality filters: finite log2PTR in [0,3]; n_windows>=5; n_anchors>=50; "
            "ori_confidence>=0.5; SE<=1.\n"
            f"Species ranking requires robust estimates in at least {args.min_samples} samples.\n"
            "Unclean_Denture is the comparison numerator in unclean_minus_clean_log2PTR.\n"
            "These are exploratory route-B estimates; no genome-sample estimate passed route-A-style QC.\n"
        )


if __name__ == "__main__":
    main()
