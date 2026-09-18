#!/usr/bin/env python3
"""Aggregate current-stats C5 refusion outputs and QC associations.

Run on HPC with the sk2bgrow environment. SciPy/statsmodels are deliberately
avoided because their imports deadlocked on the current Lustre-backed env.
"""
import glob
import os
from math import erfc, sqrt

import numpy as np
import pandas as pd

BASE = "/lustre1/g/aos_shihuang/sk2bgrow-hpc/bench/C5"
NEW = "/lustre1/g/aos_shihuang/sk2bgrow-hpc/bench/C5_refuse_current"
OUT = NEW + "/review"
os.makedirs(OUT, exist_ok=True)
rng = np.random.default_rng(20260918)


def corr_p(x, y, nperm=5000):
    d = pd.DataFrame({"x": x, "y": y}).dropna()
    rx = d.x.rank().to_numpy()
    ry = d.y.rank().to_numpy()
    rx = (rx - rx.mean()) / rx.std()
    ry = (ry - ry.mean()) / ry.std()
    obs = float(np.sum(rx * ry) / len(rx))
    ge = sum(
        abs(float(np.sum(rx * rng.permutation(ry)) / len(rx))) >= abs(obs)
        for _ in range(nperm)
    )
    return obs, (ge + 1) / (nperm + 1), len(d)


def partial_corr_p(x, y, z, nperm=5000):
    d = pd.DataFrame({"x": x, "y": y, "z": z}).dropna()
    x, y, z = (d[col].to_numpy(float) for col in ["x", "y", "z"])

    def resid(a, b):
        A = np.column_stack([np.ones(len(a)), b])
        return a - A @ np.linalg.lstsq(A, a, rcond=None)[0]

    rx, ry, rz = (pd.Series(v).rank().to_numpy() for v in [x, y, z])
    ex, ey = resid(rx, rz), resid(ry, rz)
    ex = (ex - ex.mean()) / ex.std()
    ey = (ey - ey.mean()) / ey.std()
    obs = float(np.sum(ex * ey) / len(ex))
    ge = 0
    for _ in range(nperm):
        ey2 = resid(pd.Series(rng.permutation(y)).rank().to_numpy(), rz)
        ey2 = (ey2 - ey2.mean()) / ey2.std()
        ge += abs(float(np.sum(ex * ey2) / len(ex))) >= abs(obs)
    return obs, (ge + 1) / (nperm + 1), len(d)


# Aggregate the nine newly refused samples.
frames = []
for path in sorted(glob.glob(NEW + "/res/A_*/output.tsv")):
    frame = pd.read_csv(path, sep="\t")
    frame["run"] = os.path.basename(os.path.dirname(path)).removeprefix("A_")
    frames.append(frame)
cur = pd.concat(frames, ignore_index=True)
cur["pass_qc"] = cur.pass_qc.astype(str).isin(["True", "true", "1"])

# The refusion-only job intentionally starts from retained windows.rates files;
# reuse the count-phase coverage from the existing C5 run, as the current code
# no longer emits that field from an absent stats object.
old_full = pd.read_csv(BASE + "/res/c5_results_raw.tsv", sep="\t")
coverage_map = old_full[old_full.arm == "sk2bgrow"][["run", "mag", "coverage"]].rename(
    columns={"mag": "genome"})
cur = cur.drop(columns=["coverage"]).merge(coverage_map, on=["run", "genome"],
                                            how="left", validate="many_to_one")
cur.to_csv(OUT + "/c5_current_results.tsv", sep="\t", index=False)

summ = cur.groupby("run").agg(
    n=("genome", "size"), n_qc=("pass_qc", "sum"), qc_rate=("pass_qc", "mean"),
    mean_enzymes=("n_enzymes", "mean"), mean_fit_rate=("enzyme_fit_rate", "mean")
).reset_index()
summ.to_csv(OUT + "/c5_current_sample_summary.tsv", sep="\t", index=False)

q = cur.groupby("genome").agg(
    n_qc_pass=("pass_qc", "sum"), n=("pass_qc", "size"),
    qc_rate_current=("pass_qc", "mean"), mean_cov=("coverage", "mean"),
    med_cov=("coverage", "median"), max_cov=("coverage", "max")
).reset_index()
qual = pd.read_csv(BASE + "/review/c5_mag_quality.tsv", sep="\t")
m = qual.merge(q, left_on="mag", right_on="genome", how="inner")
m.to_csv(OUT + "/c5_current_qc_x_quality.tsv", sep="\t", index=False)

rows = []
for feature in ["n_contigs", "n50", "Completeness", "Contamination"]:
    rho0, p0, n = corr_p(m.qc_rate_current, m[feature])
    rhoc, pc, _ = partial_corr_p(m.qc_rate_current, m[feature], m.mean_cov)
    rows.append([feature, rho0, p0, rhoc, pc, n])
assoc = pd.DataFrame(rows, columns=[
    "feature", "rho_raw", "p_raw", "rho_partial_cov", "p_partial", "n"])
assoc.to_csv(OUT + "/c5_current_coverage_control.tsv", sep="\t", index=False)

m["cov_q"] = pd.qcut(m.mean_cov, 4, labels=["Q1", "Q2", "Q3", "Q4"])
m["contig_bin"] = pd.cut(
    m.n_contigs, [0, 10, 25, 50, 100, 200, 1e9],
    labels=["<=10", "11-25", "26-50", "51-100", "101-200", ">200"])
with open(OUT + "/c5_current_coverage_control.txt", "w") as fh:
    fh.write(f"n merged = {len(m)}\n\n")
    fh.write("=== raw vs coverage-partial Spearman (current qc_rate) ===\n")
    fh.write(assoc.to_string(index=False, float_format=lambda x: f"{x:.4g}") + "\n\n")
    fh.write("=== bivariate Spearman ===\n")
    for feature in ["n_contigs", "n50", "Completeness", "Contamination", "qc_rate_current"]:
        rho, p, n = corr_p(m[feature], m.mean_cov)
        fh.write(f"rho({feature}, mean_cov)={rho:+.3f}, p={p:.3g}, n={n}\n")
    fh.write("\n=== within coverage quartile: Spearman(qc_rate_current,n_contigs) ===\n")
    for label, group in m.groupby("cov_q", observed=True):
        rho, p, _ = corr_p(group.qc_rate_current, group.n_contigs)
        fh.write(
            f"{label}: n={len(group)}, cov {group.mean_cov.min():.2f}-"
            f"{group.mean_cov.max():.2f}, rho={rho:+.3f}, p={p:.4g}\n")
    fh.write("\n=== mean current qc rate: coverage x contig bin ===\n")
    pivot = m.pivot_table(index="cov_q", columns="contig_bin", values="qc_rate_current",
                          aggfunc="mean", observed=True)
    fh.write(pivot.to_string(float_format=lambda x: f"{x:.3f}") + "\n")

cons = cur.dropna(subset=["log2(PTR)"]).groupby("genome")["log2(PTR)"].agg(
    std_log2ptr="std", n="size").reset_index()
mc = cons.merge(m[["mag", "mean_cov", "qc_rate_current"]],
                left_on="genome", right_on="mag", how="inner")
mc["qc_any"] = mc.qc_rate_current > 0
mc["cov_q"] = pd.qcut(mc.mean_cov, 4, labels=["Q1", "Q2", "Q3", "Q4"])
mc.to_csv(OUT + "/c5_current_cross_sample_consistency.tsv", sep="\t", index=False)
d = mc.dropna(subset=["mean_cov", "std_log2ptr"]).copy()
d["log_cov"] = np.log10(d.mean_cov)
beta = se = pvals = None
if d.qc_any.nunique() == 2:
    X = np.column_stack([np.ones(len(d)), d.qc_any.astype(float), d.log_cov])
    y = d.std_log2ptr.to_numpy()
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    resid = y - X @ beta
    dof = len(y) - X.shape[1]
    mse = np.dot(resid, resid) / dof
    covb = mse * np.linalg.inv(X.T @ X)
    se = np.sqrt(np.diag(covb))
    pvals = [erfc(abs(b / s) / sqrt(2)) for b, s in zip(beta, se)]

with open(OUT + "/c5_current_coverage_control.txt", "a") as fh:
    fh.write("\n=== consistency: median std log2PTR, QC-any vs none ===\n")
    for label, group in [("ALL", mc)] + list(mc.groupby("cov_q", observed=True)):
        a = group.loc[group.qc_any, "std_log2ptr"]
        b = group.loc[~group.qc_any, "std_log2ptr"]
        fh.write(
            f"{label}: n_qc={len(a)}, median={a.median():.3f}; "
            f"n_no_qc={len(b)}, median={b.median():.3f}\n")
    if beta is None:
        fh.write("\nOLS omitted: the current refusal policy leaves no "
                 "within-genome cross-sample variance contrast.\n")
    else:
        fh.write("\nOLS: std_log2ptr ~ qc_any + log10(cov)\n")
        for term, b, s, p in zip(["intercept", "qc_any", "log_cov"], beta, se, pvals):
            fh.write(f"  beta_{term}={b:.4f}, p={p:.3g}\n")

# Direct check against the retained prior fusion.
old = old_full[old_full.arm == "sk2bgrow"][
    ["run", "mag", "log2ptr", "pass_qc", "n_enzymes"]
].rename(columns={
    "mag": "genome", "log2ptr": "old_log2ptr", "pass_qc": "old_pass_qc",
    "n_enzymes": "old_n_enzymes"})
old["old_pass_qc"] = old.old_pass_qc.astype(str).isin(["True", "true", "1"])
cmp = cur[["run", "genome", "log2(PTR)", "pass_qc", "n_enzymes"]].merge(
    old, on=["run", "genome"], how="outer", indicator=True)
cmp["d_log2ptr"] = cmp["log2(PTR)"] - cmp.old_log2ptr
cmp["same_pass_qc"] = cmp.pass_qc == cmp.old_pass_qc
cmp.to_csv(OUT + "/c5_current_vs_previous.tsv", sep="\t", index=False)
both = cmp[cmp._merge == "both"]
with open(OUT + "/c5_current_vs_previous.txt", "w") as fh:
    fh.write(f"n both={len(both)}, current_only={(cmp._merge == 'left_only').sum()}, "
             f"old_only={(cmp._merge == 'right_only').sum()}\n")
    fh.write(f"max abs delta log2PTR={both.d_log2ptr.abs().max():.6g}, "
             f"median abs delta={both.d_log2ptr.abs().median():.6g}\n")
    fh.write(f"QC disagreements={(~both.same_pass_qc).sum()} / {len(both)}\n")
    fh.write(f"old total QC={old.old_pass_qc.sum()}, current total QC={cur.pass_qc.sum()}\n")

print(open(OUT + "/c5_current_coverage_control.txt").read())
print(open(OUT + "/c5_current_vs_previous.txt").read())
print(summ.to_string(index=False))
