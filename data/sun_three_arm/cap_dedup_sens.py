#!/usr/bin/env python3
"""Cap-dedup sensitivity: does an intermediate dedup cap recover the C arm?

Transform deduped counts C2: C'_i = min(C_i, cap) for cap in {2,3,5}
(cap=1 == current C2, cap=inf == raw C). Rerun stats on the transformed
tables (no re-profiling), compare each vs B with the same gates as the
diagnosis (fraction>=0.5, dispersion<5, both arms).

Also recomputes [D] bias-vs-coverage with B-arm outliers winsorized
(|ba_diff| <= 3) — the uncontrolled version from c_arm_diag is not citable.

Usage: cap_dedup_sens.py            (all samples, all caps)
"""
import os
import subprocess
import sys

import numpy as np
import pandas as pd

P2 = "/lustre1/g/aos_shihuang/sk2bgrow-hpc/bench/P2_sigma"
PY = "/lustre1/g/aos_shihuang/sk2bgrow-hpc/micromamba/envs/sk2bgrow/bin/python"
ENV = dict(os.environ,
           PYTHONPATH="/lustre1/g/aos_shihuang/sk2bgrow-hpc/src/python"
           + os.pathsep + os.environ.get("PYTHONPATH", ""))
PAIRS = [("S01", "SRR13371683"), ("S06", "SRR13371682"), ("S07", "SRR13371681")]
CAPS = [2, 3, 5]

def base_counts(s_c):
    import glob
    hits = [p for p in glob.glob(f"{P2}/counts_C2/{s_c}/*.counts.tsv")
            if ".cap" not in p and ".epooled" not in p and ".eprior" not in p]
    assert len(hits) == 1, f"{s_c}: {hits}"
    return hits[0]

def run_stats(counts_path, out_dir, windows):
    os.makedirs(out_dir, exist_ok=True)
    if os.path.exists(f"{out_dir}/output.tsv"):
        print(f"skip {out_dir} (exists)")
        return
    cmd = [PY, "-m", "sk2bgrow.cli", "profile", counts_path,
           "--db", f"{P2}/db_bcgI", "--output", out_dir,
           "--windows", windows, "--use-rust-windows", "--count-model", "ztp"]
    with open(f"{out_dir}/stats.log", "w") as lg:
        subprocess.run(cmd, check=True, stdout=lg, stderr=lg, env=ENV)

def cmp_arm(tc, tb, label, s_c):
    m = tc.merge(tb, on="genome", suffixes=("_C", "_B"))
    m = m[(m.fraction_C >= 0.5) & (m.fraction_B >= 0.5)]
    m = m[(m.dispersion_C < 5) & (m.dispersion_B < 5)]
    if len(m) < 3:
        print(f"{label:22s} {s_c}: n={len(m)} (too few)")
        return None
    x = np.log2(m.PTR_B.values); y = np.log2(m.PTR_C.values)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    r = np.corrcoef(x, y)[0, 1]
    ccc = 2 * r * x.std() * y.std() / (x.std()**2 + y.std()**2 + (x.mean() - y.mean())**2)
    d = y - x
    print(f"{label:22s} {s_c}: n={len(m)} r={r:.3f} CCC={ccc:.3f} "
          f"bias={d.mean():+.3f} SD={d.std(ddof=1):.3f} "
          f"C_range=[{y.min():.2f},{y.max():.2f}] P(|d|<0.263)={np.mean(np.abs(d) < 0.263):.1%}")
    return m

rows = []
for s_c, s_b in PAIRS:
    tb = pd.read_csv(f"{P2}/counts_B/{s_b}/stats/output.tsv", sep="\t")
    win = f"{P2}/counts_C2/{s_c}/windows.tsv"
    df0 = pd.read_csv(base_counts(s_c), sep="\t")

    # cap transform + stats
    for cap in CAPS:
        out_counts = f"{P2}/counts_C2/{s_c}/{s_c}.cap{cap}.counts.tsv"
        if not os.path.exists(out_counts):
            dfc = df0.copy()
            dfc["count"] = dfc["count"].clip(upper=cap)
            dfc.to_csv(out_counts, sep="\t", index=False)
        out_dir = f"{P2}/counts_C2/{s_c}/stats_cap{cap}"
        run_stats(out_counts, out_dir, win)
        tc = pd.read_csv(f"{out_dir}/output.tsv", sep="\t")
        m = cmp_arm(tc, tb, f"C_cap{cap}_vs_B", s_c)
        if m is not None:
            rows.append((s_c, cap, len(m)))

    # bookends from existing results
    for label, path in [("C_cap1(dedup)_vs_B", f"{P2}/counts_C2/{s_c}/stats/output.tsv"),
                        ("C_capinf(raw)_vs_B", f"{P2}/counts_C/{s_c}/stats/output.tsv")]:
        tc = pd.read_csv(path, sep="\t")
        cmp_arm(tc, tb, label, s_c)

# [D] redo with winsorization, pooled across samples, per cap
print("\n[D-winsorized] Spearman(bias, log10 coverage_C), |bias|<=3")
for label, pathfun in [
    ("cap1", lambda s: f"{P2}/counts_C2/{s}/stats/output.tsv"),
    ("capinf", lambda s: f"{P2}/counts_C/{s}/stats/output.tsv"),
] + [(f"cap{c}", (lambda c: lambda s: f"{P2}/counts_C2/{s}/stats_cap{c}/output.tsv")(c))
     for c in CAPS]:
    bs, cs = [], []
    for s_c, s_b in PAIRS:
        tb = pd.read_csv(f"{P2}/counts_B/{s_b}/stats/output.tsv", sep="\t")
        tc = pd.read_csv(pathfun(s_c), sep="\t")
        m = tc.merge(tb, on="genome", suffixes=("_C", "_B"))
        m = m[(m.fraction_C >= 0.5) & (m.fraction_B >= 0.5)]
        m = m[(m.dispersion_C < 5) & (m.dispersion_B < 5)]
        d = np.log2(m.PTR_C) - np.log2(m.PTR_B)
        keep = np.isfinite(d) & (np.abs(d) <= 3)
        bs += d[keep].tolist()
        cs += np.log10(m.coverage_C.values[keep]).tolist()
    if len(bs) >= 8:
        rho = pd.Series(bs).corr(pd.Series(cs), method="spearman")
        print(f"  {label:8s} rho={rho:+.3f} n={len(bs)}")
print("CAP_DEDUP_SENS_DONE")
