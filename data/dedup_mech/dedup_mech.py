#!/usr/bin/env python3
"""Dedup mechanism regression (P2 follow-up).

Hypothesis: PCR preferentially duplicates high-copy templates, so anchors near
the replication origin (multicopy in a replicating genome pool) carry a higher
duplicate fraction; exact dedup therefore removes the very gradient sk2bGrow
needs. Predicted pattern: per-anchor duplication fraction df = 1 - dedup/raw
decreases with circular distance from ori  ->  slope(df ~ dist) < 0,
equivalently slope(df ~ proximity=0.5-dist) > 0.

Data (all pre-existing, P2_sigma):
  raw counts   : counts_C/Sxx/Sxx_*.counts.tsv        (sk2bgrow profile, raw fq)
  dedup counts : counts_C2/Sxx/Sxx_*.dedup.fq.gz.counts.tsv (same db_bcgI index)
  ori calls    : counts_C/Sxx/stats/output.tsv (stats run on the RAW counts)
  genome length: db_bcgI/manifest.json (sum of contig lengths)

Unit of observation = species (genome) x sample. Regression over anchors with
raw > 0 (df undefined otherwise), m >= 20. Unweighted OLS of df on circular
distance fraction (0 at ori, 0.5 opposite). n=3 samples: sign consistency only,
no inferential testing.

Usage: dedup_mech.py out_prefix
"""
import json
import sys

import numpy as np
import pandas as pd

P2 = "/lustre1/g/aos_shihuang/sk2bgrow-hpc/bench/P2_sigma"
SAMPLES = ["S01", "S06", "S07"]
RAW = {s: f"{P2}/counts_C/{s}" for s in SAMPLES}
DEDUP = {s: f"{P2}/counts_C2/{s}" for s in SAMPLES}
MIN_M = 20

import glob

def counts_path(d, dedup):
    pat = "*.dedup.fq.gz.counts.tsv" if dedup else "*.counts.tsv"
    return glob.glob(f"{d}/{pat}")[0]

# genome length from manifest (sum of contig lengths)
with open(f"{P2}/db_bcgI/manifest.json") as fh:
    manifest = json.load(fh)
genome_len = {}
for g in manifest["genomes"]:
    genome_len[g["name"]] = sum(c["length"] for c in g["contigs"])

rows = []
for s in SAMPLES:
    raw = pd.read_csv(counts_path(RAW[s], False), sep="\t",
                      usecols=["genome_id", "genome", "contig_id", "position",
                               "global_position", "count"])
    ded = pd.read_csv(counts_path(DEDUP[s], True), sep="\t",
                      usecols=["genome_id", "contig_id", "position", "count"])
    m = raw.merge(ded, on=["genome_id", "contig_id", "position"],
                  suffixes=("_raw", "_ded"))
    viol = int((m["count_ded"] > m["count_raw"]).sum())
    print(f"[{s}] merged anchors: {len(m)}, dedup>raw violations: {viol}")

    st = pd.read_csv(f"{RAW[s]}/stats/output.tsv", sep="\t",
                     usecols=["genome", "ori"])
    ori = {r["genome"]: float(r["ori"]) for _, r in st.iterrows()
           if str(r["ori"]) not in ("NA", "", "nan")}

    m = m[m["count_raw"] > 0].copy()
    m["df"] = 1.0 - m["count_ded"] / m["count_raw"]
    m["df"] = m["df"].clip(lower=0.0, upper=1.0)

    for gname, g in m.groupby("genome"):
        if gname not in ori or gname not in genome_len:
            continue
        L = genome_len[gname]
        d = (g["global_position"] - ori[gname]).abs() % L
        d = np.minimum(d, L - d) / L  # circular distance fraction, 0..0.5
        n = len(g)
        if n < MIN_M:
            continue
        x = d.values.astype(float)
        y = g["df"].values.astype(float)
        xc, yc = x - x.mean(), y - y.mean()
        var_x = (xc * xc).sum()
        if var_x <= 0:
            continue
        slope = float((xc * yc).sum() / var_x)
        intercept = float(y.mean() - slope * x.mean())
        ss_res = float(((y - (intercept + slope * x)) ** 2).sum())
        ss_tot = float((yc * yc).sum())
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
        rows.append({
            "species": gname, "sample": s,
            "slope": round(slope, 6), "slope_prox": round(-slope, 6),
            "intercept": round(intercept, 6),
            "n_anchors": n, "r2": round(r2, 4),
            "median_df": round(float(np.median(y)), 6),
        })
        print(f"[{s}] {gname}: n={n} slope={slope:+.5f} r2={r2:.3f} "
              f"median_df={np.median(y):.3f}")

res = pd.DataFrame(rows).sort_values(["sample", "species"])
out = sys.argv[1]
res.to_csv(out + ".tsv", sep="\t", index=False)

print("\n=== sign consistency (mechanism predicts slope<0, slope_prox>0) ===")
for s in SAMPLES:
    r = res[res["sample"] == s]
    neg = int((r["slope"] < 0).sum())
    print(f"{s}: species={len(r)}  slope<0: {neg}  slope>0: {len(r)-neg}  "
          f"frac_neg={neg/len(r):.2f}" if len(r) else f"{s}: no species")
print("\nwrote", out + ".tsv")
