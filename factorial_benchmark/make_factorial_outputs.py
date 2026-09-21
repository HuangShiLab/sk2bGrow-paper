#!/usr/bin/env python3
"""Summarize the count-level factorial and make its manuscript figure/table.

    python3 factorial_benchmark/make_factorial_outputs.py

Reads only factorial_benchmark/factorial_long.tsv; writes Table 11 and Fig. 5.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "figures"))
from style import apply, grid  # noqa: E402

LONG = ROOT / "factorial_benchmark" / "factorial_long.tsv"
OUT_TABLE = ROOT / "tables" / "table11_factorial_count_level"
FIG_STEM = ROOT / "figures" / "out" / "fig5_factorial_mechanism"

ESTIMATOR_LABEL = {
    "coordinate_vfit": "Coordinate V-fit",
    "sorted_rank": "Sorted-rank regression",
}
SOURCE_LABEL = {
    "coordinate_grid": "Regular coordinates (2bRAD-like)",
    "random_sketch": "Random coordinates (FracMinHash-like)",
}


def aggregate(long: pd.DataFrame) -> pd.DataFrame:
    keys = [
        "estimator", "landmark_source", "n_strains", "depth",
        "shared_fraction", "similarity", "abundance",
    ]
    return (
        long.groupby(keys, as_index=False)
        .agg(
            n=("n", "sum"),
            n_ok=("n_ok", "sum"),
            reported_fraction=("reported_fraction", "mean"),
            bias=("bias", "mean"),
            rmse=("rmse", "mean"),
            mae=("mae", "mean"),
            correlation=("correlation", "mean"),
        )
        .sort_values(keys)
        .reset_index(drop=True)
    )


def source_contrast(agg: pd.DataFrame) -> pd.DataFrame:
    keys = ["estimator", "n_strains", "depth", "shared_fraction", "similarity", "abundance"]
    grid = agg[agg["landmark_source"] == "coordinate_grid"][
        keys + ["bias", "rmse", "correlation"]
    ].rename(columns={"bias": "bias_grid", "rmse": "rmse_grid", "correlation": "r_grid"})
    random = agg[agg["landmark_source"] == "random_sketch"][
        keys + ["bias", "rmse", "correlation"]
    ].rename(columns={"bias": "bias_random", "rmse": "rmse_random", "correlation": "r_random"})
    paired = grid.merge(random, on=keys, validate="one_to_one")
    paired["rmse_grid_minus_random"] = paired["rmse_grid"] - paired["rmse_random"]
    return paired


def write_table(agg: pd.DataFrame) -> None:
    coord = agg[
        agg["depth"].isin([0.5, 1.0, 8.0])
        & agg["landmark_source"].eq("coordinate_grid")
    ]
    show = (
        coord.groupby(["estimator", "shared_fraction", "depth"], as_index=False)[
            ["bias", "rmse", "correlation"]
        ]
        .mean()
    )
    show["estimator"] = show["estimator"].map(ESTIMATOR_LABEL)
    show["shared_anchors"] = (show["shared_fraction"] * 100).round(0).astype(int).astype(str) + "%"
    show["depth"] = show["depth"].map({0.5: "0.5×", 1.0: "1×", 8.0: "8×"})
    # Avoid a signed-zero artifact when rounding near-zero means.
    show.loc[show["bias"].abs() < 5e-4, "bias"] = 0.0
    show = show[
        ["estimator", "shared_anchors", "depth", "bias", "rmse", "correlation"]
    ].rename(
        columns={
            "bias": "Mean bias (log2)",
            "rmse": "Mean RMSE (log2)",
            "correlation": "Mean Pearson r",
        }
    )
    tsv = show.to_csv(sep="\t", index=False, float_format="%.3f")
    md = show.to_markdown(index=False, floatfmt=".3f")
    OUT_TABLE.with_suffix(".tsv").write_text(tsv)
    OUT_TABLE.with_suffix(".md").write_text(
        "**Table 11. Count-level factorial decomposition of estimator, depth and "
        "shared-anchor ambiguity.**\n\n"
        "Values are means over 4/8/16/32-strain communities, even/10:1/100:1 "
        "abundance ratios and 10 replicates, shown for the regular coordinate "
        "arm. The random FracMinHash-like arm and paired source "
        "contrasts are in `factorial_benchmark/factorial_long.tsv` and "
        "`factorial_benchmark/factorial_source_contrasts.tsv`. Bias and RMSE "
        "are in log2(PTR) units.\n\n" + md + "\n"
    )


def make_figure(agg: pd.DataFrame, paired: pd.DataFrame) -> None:
    apply()
    fig = plt.figure(figsize=(12.4, 3.55))
    axes = fig.subplots(1, 3, gridspec_kw={"wspace": 0.38})

    # (a) estimator effect at the private-anchor control
    ax = axes[0]
    control = (
        agg[agg["shared_fraction"].eq(0.0)]
        .groupby(["estimator", "depth"], as_index=False)["rmse"]
        .mean()
    )
    for estimator, color in [
        ("coordinate_vfit", "#2a78d6"),
        ("sorted_rank", "#eb6834"),
    ]:
        x = control[control["estimator"].eq(estimator)]
        ax.plot(x["depth"], x["rmse"], marker="o", color=color,
                label=ESTIMATOR_LABEL[estimator])
    ax.set_xlabel("Per-strain depth (×; log2 scale)")
    ax.set_ylabel("Mean RMSE (log2PTR)")
    ax.set_title("a  Private-anchor control")
    ax.set_xticks([0.5, 1, 2, 4, 8])
    ax.set_xscale("log", base=2)
    ax.set_ylim(0, 0.70)
    ax.legend(loc="upper right")
    grid(ax)

    # (b) ambiguity increases coordinate-fit bias, and depth does not remove it
    ax = axes[1]
    coord = (
        agg[agg["estimator"].eq("coordinate_vfit")]
        .groupby(["shared_fraction", "depth"], as_index=False)["bias"]
        .mean()
    )
    palette = {"0.0": "#2a78d6", "0.05": "#eda100", "0.15": "#d0342c"}
    labels = {"0.0": "0% shared", "0.05": "5% shared", "0.15": "15% shared"}
    for shared, group in coord.groupby("shared_fraction"):
        ax.plot(group["depth"], group["bias"], marker="o",
                color=palette[str(shared)], label=labels[str(shared)])
    ax.axhline(0.0, color="#8a8a8a", linewidth=0.8)
    ax.set_xlabel("Per-strain depth (×; log2 scale)")
    ax.set_ylabel("Coordinate V-fit bias (log2PTR)")
    ax.set_title("b  Shared-anchor ambiguity")
    ax.set_xticks([0.5, 1, 2, 4, 8])
    ax.set_xscale("log", base=2)
    ax.set_ylim(-1.4, 0.1)
    ax.legend(loc="lower left")
    grid(ax)

    # (c) source × estimator × ambiguity interaction at 8x
    ax = axes[2]
    deep = (
        agg[agg["depth"].eq(8.0)]
        .groupby(["estimator", "landmark_source", "shared_fraction"], as_index=False)["rmse"]
        .mean()
    )
    source_color = {
        "coordinate_grid": "#2a78d6",
        "random_sketch": "#1baf7a",
    }
    for (estimator, source), group in deep.groupby(["estimator", "landmark_source"]):
        group = group.sort_values("shared_fraction")
        ax.plot(
            group["shared_fraction"], group["rmse"], marker="o",
            color=source_color[source],
            linestyle="-" if estimator == "coordinate_vfit" else "--",
            label=f"{SOURCE_LABEL[source]} · {ESTIMATOR_LABEL[estimator]}",
        )
    ax.set_xlabel("Shared-anchor fraction")
    ax.set_ylabel("Mean RMSE at 8× (log2PTR)")
    ax.set_title("c  Ambiguity, placement and estimator")
    ax.set_xticks([0, 0.05, 0.15])
    ax.set_ylim(0, 2.62)
    ax.legend(loc="upper left", fontsize=7)
    grid(ax)

    for ext in ("png", "pdf"):
        fig.savefig(FIG_STEM.with_suffix(f".{ext}"))
    plt.close(fig)


def main() -> None:
    long = pd.read_csv(LONG, sep="\t")
    agg = aggregate(long)
    paired = source_contrast(agg)
    OUT_TABLE.parent.mkdir(parents=True, exist_ok=True)
    write_table(agg)
    make_figure(agg, paired)
    # Persist the paired source contrasts used in the caption and Discussion.
    paired.to_csv(ROOT / "factorial_benchmark" / "factorial_source_contrasts.tsv",
                  sep="\t", index=False, float_format="%.8g")
    print(f"wrote {OUT_TABLE}.tsv/.md, {FIG_STEM}.png/.pdf, and source contrasts")


if __name__ == "__main__":
    main()
