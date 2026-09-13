#!/usr/bin/env python3
"""Build mm1_full_sample_cost.tsv from per-sample time.txt + analysis.json,
plus the baseline c5_cost_per_sample.tsv (mm=2). Also summarizes the existing
1M-pair subset A/B (SRR28338159, res/_mm1_bench) as the per-anchor validation.
Run on the HPC after samples complete (any subset of the 9 may be present).
"""
import glob
import json
import os
import re

import pandas as pd

BASE = "/lustre1/g/aos_shihuang/sk2bgrow-hpc"
C5 = f"{BASE}/bench/C5"
E2E = f"{BASE}/bench/mm1_e2e"
OUT = f"{E2E}/mm1_full_sample_cost.tsv"

MG = ["SRR28338172", "SRR28338171", "SRR28338162", "SRR28338161", "SRR28338160",
      "SRR28338159", "SRR28338158", "SRR28338157", "SRR28338156"]


def parse_time(path):
    if not os.path.exists(path):
        return {}
    txt = open(path, errors="replace").read()
    out = {}
    m = re.search(r"Elapsed \(wall clock\) time.*?:\s*(?:(\d+)-)?(?:(\d+):)?(\d+):([\d.]+)", txt)
    if m:
        d = int(m.group(1) or 0); h = int(m.group(2) or 0)
        out["wall_s"] = d * 86400 + h * 3600 + int(m.group(3)) * 60 + float(m.group(4))
    r = re.search(r"Maximum resident set size \(kbytes\): (\d+)", txt)
    if r:
        out["rss_kb"] = int(r.group(1))
    c = re.search(r"Percent of CPU this job got: (\d+)%", txt)
    if c:
        out["cpu_pct"] = int(c.group(1))
    out["exit_ok"] = "Exit status: 0" in txt
    return out


base = pd.read_csv(f"{C5}/review/c5_cost_per_sample.tsv", sep="\t")
base = base.set_index("run")

rows = []
for run in MG:
    mm1_dir = f"{E2E}/res/A_{run}"
    if not os.path.isdir(mm1_dir):
        continue
    ct = parse_time(f"{mm1_dir}/time.txt")
    st = parse_time(f"{mm1_dir}/stats.time.txt")
    aj = {}
    if os.path.exists(f"{mm1_dir}/analysis.json"):
        aj = json.load(open(f"{mm1_dir}/analysis.json"))
    b = base.loc[run] if run in base.index else {}
    count_wall = ct.get("wall_s")
    stats_wall = st.get("wall_s")
    e2e = (count_wall or 0) + (stats_wall or 0) if count_wall and stats_wall else None
    base_count = float(b["count_wall_s"]) if hasattr(b, "get") and b.get("count_wall_s") else None
    base_stats = float(b["stats_wall_s"]) if hasattr(b, "get") and b.get("stats_wall_s") else None
    base_e2e = (base_count + base_stats) if base_count and base_stats else None
    rows.append({
        "sample": run,
        "wall_s": e2e,
        "peak_rss_kb": max([x for x in [ct.get("rss_kb"), st.get("rss_kb")] if x] or [None]),
        "count_wall_s": count_wall,
        "stats_wall_s": stats_wall,
        "count_rss_kb": ct.get("rss_kb"),
        "stats_rss_kb": st.get("rss_kb"),
        "count_cpu_pct": ct.get("cpu_pct"),
        "count_exit_ok": ct.get("exit_ok"),
        "stats_exit_ok": st.get("exit_ok"),
        "anchors_nonzero_mm2": aj.get("anchors_nonzero_mm2"),
        "anchors_nonzero_mm1": aj.get("anchors_nonzero_mm1"),
        "genomes_lost": aj.get("genomes_lost"),
        "median_count_ratio": aj.get("median_count_ratio"),
        "anchor_loss_frac": aj.get("anchor_loss_frac"),
        "match_retention_mm1_over_mm2": aj.get("match_retention_mm1_over_mm2"),
        "mm1_reads_total": aj.get("mm1_reads_total"),
        "mm1_reads_with_anchor": aj.get("mm1_reads_with_anchor"),
        "mm2_reads_with_anchor": aj.get("mm2_reads_with_anchor"),
        "baseline_count_wall_s": base_count,
        "baseline_stats_wall_s": base_stats,
        "baseline_e2e_wall_s": base_e2e,
        "speedup_count": (base_count / count_wall) if base_count and count_wall else None,
        "speedup_e2e": (base_e2e / e2e) if base_e2e and e2e else None,
    })

df = pd.DataFrame(rows)
df.to_csv(OUT, sep="\t", index=False)
print(df.to_string(index=False))
print(f"\nwrote {OUT} ({len(df)} samples)")

# ---- 1M-pair subset A/B summary (SRR28338159, existing _mm1_bench) ----
RES = f"{C5}/res/_mm1_bench"
try:
    a = pd.read_csv(glob.glob(f"{RES}/mm2ns/*.counts.tsv")[0], sep="\t")
    b1 = pd.read_csv(glob.glob(f"{RES}/mm1/*.counts.tsv")[0], sep="\t")
    key = ["genome_id", "contig_id", "position"]
    m = a.merge(b1, on=key, suffixes=("_2", "_1"))
    nz2 = m["count_2"] > 0; nz1 = m["count_1"] > 0
    mm = m[nz2 & nz1]
    ratio = mm["count_1"] / mm["count_2"]
    ga = a.groupby("genome_id")["count"].sum()
    gb = b1.groupby("genome_id")["count"].sum()
    g = pd.concat([ga, gb], axis=1, keys=["mm2", "mm1"]).fillna(0)
    lost_g = int(((g["mm2"] > 0) & (g["mm1"] == 0)).sum())
    t2 = parse_time(f"{RES}/mm2ns.time.txt"); t1 = parse_time(f"{RES}/mm1.time.txt")
    summ = pd.DataFrame([{
        "sample": "SRR28338159", "subset": "1M_pairs",
        "anchors_compared": len(m),
        "anchors_nonzero_mm2": int(nz2.sum()),
        "anchors_nonzero_mm1": int(nz1.sum()),
        "anchors_lost_by_mm1": int((nz2 & ~nz1).sum()),
        "anchor_loss_frac": float((nz2 & ~nz1).sum() / nz2.sum()),
        "median_count_ratio_shared_anchors": float(ratio.median()),
        "genomes_nonzero_mm2": int((g["mm2"] > 0).sum()),
        "genomes_lost_entirely": lost_g,
        "count_wall_mm2_s": t2.get("wall_s"), "count_wall_mm1_s": t1.get("wall_s"),
        "speedup": t2.get("wall_s") / t1.get("wall_s") if t2.get("wall_s") and t1.get("wall_s") else None,
    }])
    summ.to_csv(f"{E2E}/mm1_subset_1M_ab_summary.tsv", sep="\t", index=False)
    print("\n1M subset A/B:"); print(summ.to_string(index=False))
except Exception as e:
    print(f"subset summary skipped: {e}")
