#!/usr/bin/env python3
"""Per-sample full-sample anchor-retention analysis: mm=1 (this task) vs mm=2
(baseline), WITHOUT rerunning mm=2.

mm=2 per-anchor counts were deleted post-stats by 01_cell.sh, but the baseline
stats stage's windows.rates.tsv retains, for every (genome, enzyme, window):
n_anchors / n_positive computed from those mm=2 counts.  Windows partition the
usable anchors of each (genome, enzyme, contig), so

    D(g) = sum of n_positive over all windows of genome g

equals the number of usable enzyme-anchors of g with count >= 1.
The usable set is count-independent (derived from anchor flags/GC), and both
arms go through the identical stats invocation, so D_mm1(g) vs D_mm2(g) is a
symmetric full-sample retention comparison.

Outputs per sample:
  <out>/mm1_vs_mm2_per_genome.tsv   D2, D1, ratio per genome
  <out>/analysis.json                summary metrics for the cost table
"""
import json
import sys
import time

import pandas as pd

BASE = "/lustre1/g/aos_shihuang/sk2bgrow-hpc"
C5 = f"{BASE}/bench/C5"
E2E = f"{BASE}/bench/mm1_e2e"


def load_detected(rates_path):
    """Per-genome usable detected anchors D(g) from a windows.rates.tsv."""
    df = pd.read_csv(rates_path, sep="\t",
                     usecols=["genome_id", "genome", "n_positive"])
    g = df.groupby("genome_id").agg(
        genome=("genome", "first"), detected=("n_positive", "sum"))
    return g


def main():
    run = sys.argv[1]
    mm1_dir = f"{E2E}/res/A_{run}"
    mm2_dir = f"{C5}/res/A_{run}"
    res = {"run": run, "analyzed_at": time.strftime("%Y-%m-%d %H:%M:%S")}

    d1 = load_detected(f"{mm1_dir}/windows.rates.tsv")
    d2 = load_detected(f"{mm2_dir}/windows.rates.tsv")
    m = d2.join(d1, on="genome_id", how="outer",
                lsuffix="_mm2", rsuffix="_mm1").fillna(0)
    m["ratio"] = m["detected_mm1"] / m["detected_mm2"].where(m["detected_mm2"] > 0)
    m = m.rename(columns={"detected_mm2": "D_mm2", "detected_mm1": "D_mm1"})
    m.to_csv(f"{mm1_dir}/mm1_vs_mm2_per_genome.tsv", sep="\t")

    res["anchors_nonzero_mm2"] = int(m["D_mm2"].sum())
    res["anchors_nonzero_mm1"] = int(m["D_mm1"].sum())
    lost = m[(m["D_mm2"] > 0) & (m["D_mm1"] == 0)]
    res["genomes_lost"] = int(len(lost))
    res["genomes_total"] = int(len(m))
    res["genomes_detected_mm2"] = int((m["D_mm2"] > 0).sum())
    shared = m[m["D_mm2"] > 0]
    res["median_count_ratio"] = float(shared["ratio"].median())
    res["anchor_loss_frac"] = float(1 - m["D_mm1"].sum() / m["D_mm2"].sum())
    res["genomes_lost_list"] = [str(x) for x in lost.index.tolist()[:20]]

    # raw (no usable filter) nonzero anchors in the mm1 counts table, reference
    try:
        c = pd.read_csv(f"{mm1_dir}/{run}.counts.tsv", sep="\t",
                        usecols=["genome_id", "count"])
        res["raw_nonzero_mm1"] = int((c["count"] > 0).sum())
        res["raw_rows_mm1"] = int(len(c))
    except Exception as e:
        res["raw_nonzero_mm1"] = None
        res["counts_read_error"] = str(e)

    # mismatch_hist comparison (full-sample, from stats.json of both arms)
    def counting(path):
        try:
            return json.load(open(path)).get("counting", {})
        except Exception:
            return {}
    c2 = counting(f"{mm2_dir}/{run}.stats.json")
    c1 = counting(f"{mm1_dir}/{run}.stats.json")
    res["mm2_mismatch_hist"] = c2.get("mismatch_hist")
    res["mm1_mismatch_hist"] = c1.get("mismatch_hist")
    if c2.get("mismatch_hist") and c1.get("mismatch_hist"):
        res["match_retention_mm1_over_mm2"] = float(
            sum(c1["mismatch_hist"]) / sum(c2["mismatch_hist"]))
    res["mm2_reads_total"] = c2.get("reads_total")
    res["mm1_reads_total"] = c1.get("reads_total")
    res["mm2_reads_with_anchor"] = c2.get("reads_with_anchor")
    res["mm1_reads_with_anchor"] = c1.get("reads_with_anchor")

    with open(f"{mm1_dir}/analysis.json", "w") as fh:
        json.dump(res, fh, indent=2)
    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
