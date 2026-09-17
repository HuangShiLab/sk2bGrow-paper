#!/usr/bin/env python3
"""Extract one multiseed instance into multiseed/regen_seed<S>.tsv.
cols: instance, medium, srr, cov, arm, log2ptr, est_cov, growth_rate"""
import csv, os, sys

BASE = "/lustre1/g/aos_shihuang/sk2bgrow-hpc"
RC = f"{BASE}/bench/repro_check"
S = sys.argv[1]
ROOT = f"{RC}/multiseed/seed{S}"
MEDIA = ["M6","M3","M2","RUN_OUT","M27","M25","M24","M22","M23","M19",
         "M18","M17","M13","M12","M1","M10","M4"]
COVS = ["0.5", "1", "2", "5", "10"]

srr_of = {}
for line in open(f"{BASE}/bench/C1/run_medium.tsv"):
    s, m = line.split()
    srr_of.setdefault(m, s)
gr_of = {}
for row in csv.DictReader(open(f"{BASE}/bench/C1/growth_rates.tsv"), delimiter="\t"):
    gr_of[row["medium"]] = row["growth_rate"]

def clean(v):
    return "" if v in (None, "", "NA", "nan", "NaN") else v

def read_est(path):
    """-> (log2ptr, est_cov) or ('','') if absent/gate-rejected."""
    if not os.path.exists(path):
        return "", ""
    with open(path) as fh:
        r = next(csv.DictReader(fh, delimiter="\t"), None)
    if r is None:
        return "", ""
    return clean(r.get("log2(PTR)")), clean(r.get("coverage"))

rows, missing = [], []
for cov in COVS:
    for m in MEDIA:
        srr = srr_of[m]
        for arm, p in [
            ("A",        f"{ROOT}/out_A/{srr}.{cov}x/output.tsv"),
            ("B",        f"{ROOT}/out_B/{srr}.{cov}x/output.tsv"),
            ("E",        f"{ROOT}/out_E/{srr}.{cov}x/output.tsv"),
            ("C_default",  f"{ROOT}/pilea/default/{srr}.{cov}x/output.tsv"),
            ("C_relaxed",  f"{ROOT}/pilea/relaxed/{srr}.{cov}x/output.tsv"),
        ]:
            lp, ec = read_est(p)
            if not os.path.exists(p):
                missing.append((arm, srr, cov))
            rows.append([S, m, srr, cov, arm, lp, ec, gr_of.get(m, "")])

out = f"{RC}/multiseed/regen_seed{S}.tsv"
with open(out, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["instance", "medium", "srr", "cov", "arm", "log2ptr",
                "est_cov", "growth_rate"])
    w.writerows(rows)
print(f"rows: {len(rows)}  missing outputs: {len(missing)}")
for arm in ["A", "B", "E", "C_default", "C_relaxed"]:
    n = sum(1 for r in rows if r[4] == arm and r[5] != "")
    print(f"  {arm}: {n}/85 with estimate")
if missing:
    print("missing:", missing[:10])
print("wrote", out)
