#!/usr/bin/env python3
"""Merge instance 1 (C1 regen tsvs) with multiseed instances into
multiseed/all_instances.tsv; per-(arm,cov,instance) Pearson r; paired
cluster bootstrap of the 2x2 interaction per depth.

interaction_z = (z_A - z_B) - (z_E - z_C_relaxed), Fisher z of Pearson r
over the 16 media x 3 instances (RUN_OUT excluded); bootstrap resamples
media with replacement (all instances of a drawn medium kept together),
B=10000, rng seed 20260912.
"""
import csv, math, random
from collections import defaultdict

BASE = "/lustre1/g/aos_shihuang/sk2bgrow-hpc"
RC = f"{BASE}/bench/repro_check"
MS = f"{RC}/multiseed"
COVS = ["0.5", "1", "2", "5", "10"]
INSTANCES = ["1", "20260912", "20260913"]
ARMS = ["A", "B", "E", "C_default", "C_relaxed"]
B_BOOT, SEED_BOOT = 10000, 20260912

def clean(v):
    return "" if v in (None, "", "NA", "nan", "NaN") else v

rows = []
# instance 1 from existing regen tsvs (cols medium,srr,cov,log2ptr,est_cov[,pass_qc],growth_rate)
for arm, path, split in [
    ("A", f"{RC}/regen_armA_cells.tsv", None),
    ("B", f"{RC}/regen_armB_cells.tsv", None),
    ("E", f"{RC}/regen_armE_PE_cells.tsv", None),
    ("C_default", f"{RC}/regen_armC_cells.tsv", "default"),
    ("C_relaxed", f"{RC}/regen_armC_cells.tsv", "relaxed"),
]:
    for r in csv.DictReader(open(path), delimiter="\t"):
        if split and r["mode"] != split:
            continue
        rows.append({"instance": "1", "medium": r["medium"], "srr": r["srr"],
                     "cov": r["cov"], "arm": arm, "log2ptr": clean(r["log2ptr"]),
                     "est_cov": clean(r["est_cov"]),
                     "growth_rate": clean(r["growth_rate"])})
# instances 2/3
for inst in INSTANCES[1:]:
    for r in csv.DictReader(open(f"{MS}/regen_seed{inst}.tsv"), delimiter="\t"):
        rows.append({"instance": inst, "medium": r["medium"], "srr": r["srr"],
                     "cov": r["cov"], "arm": r["arm"], "log2ptr": clean(r["log2ptr"]),
                     "est_cov": clean(r["est_cov"]),
                     "growth_rate": clean(r["growth_rate"])})

cols = ["instance", "medium", "srr", "cov", "arm", "log2ptr", "est_cov", "growth_rate"]
with open(f"{MS}/all_instances.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, delimiter="\t", fieldnames=cols)
    w.writeheader()
    w.writerows(rows)

def pearson(pts):
    n = len(pts)
    if n < 3:
        return None
    xs = [x for x, _ in pts]; ys = [y for _, y in pts]
    mx, my = sum(xs) / n, sum(ys) / n
    sx = math.sqrt(sum((x - mx) ** 2 for x in xs))
    sy = math.sqrt(sum((y - my) ** 2 for y in ys))
    if sx == 0 or sy == 0:
        return None
    return sum((x - mx) * (y - my) for x, y in pts) / (sx * sy)

def fisher(r):
    return math.atanh(max(min(r, 0.999999), -0.999999)) if r is not None else None

# index: {(arm,cov,instance): {medium: (log2ptr, growth_rate)}}
idx = defaultdict(dict)
for r in rows:
    if r["medium"] == "RUN_OUT" or not r["log2ptr"] or not r["growth_rate"]:
        continue
    idx[(r["arm"], r["cov"], r["instance"])][r["medium"]] = \
        (float(r["log2ptr"]), float(r["growth_rate"]))

MEDIA = sorted({r["medium"] for r in rows if r["medium"] != "RUN_OUT"})

print("=== Pearson r per (arm, cov, instance) ===")
print(f"{'arm':11s} {'cov':>4s} " + " ".join(f"inst{i:>10s}" for i in INSTANCES) + "   pooled48")
rstore = {}
for arm in ARMS:
    for cov in COVS:
        vals = []
        for inst in INSTANCES:
            pts = list(idx[(arm, cov, inst)].values())
            r = pearson(pts)
            rstore[(arm, cov, inst)] = r
            vals.append(f"{r:>10.4f}" if r is not None else "        NA")
        pts48 = [v for inst in INSTANCES for v in idx[(arm, cov, inst)].values()]
        r48 = pearson(pts48)
        rstore[(arm, cov, "pooled")] = r48
        print(f"{arm:11s} {cov:>4s} " + " ".join(vals) +
              (f"   {r48:.4f}" if r48 is not None else "   NA"))

def interaction(media_set, cov):
    zs = {}
    for arm in ["A", "B", "E", "C_relaxed"]:
        pts = [idx[(arm, cov, inst)][m] for inst in INSTANCES
               for m in media_set if m in idx[(arm, cov, inst)]]
        r = pearson(pts)
        if r is None:
            return None
        zs[arm] = fisher(r)
    return (zs["A"] - zs["B"]) - (zs["E"] - zs["C_relaxed"])

rng = random.Random(SEED_BOOT)
print("\n=== interaction bootstrap (cluster by medium, B=10000, seed %d) ===" % SEED_BOOT)
print(f"{'cov':>4s} {'z_obs':>8s} {'boot_mean':>9s} {'boot_sd':>8s} {'CI95':>18s} {'valid':>6s}")
inter_res = {}
for cov in COVS:
    z_obs = interaction(MEDIA, cov)
    reps = []
    for _ in range(B_BOOT):
        sample = [rng.choice(MEDIA) for _ in range(len(MEDIA))]
        z = interaction(sample, cov)
        if z is not None:
            reps.append(z)
    if z_obs is None or len(reps) < B_BOOT * 0.5:
        inter_res[cov] = (z_obs, None, None, None, len(reps))
        print(f"{cov:>4s} {'NA' if z_obs is None else f'{z_obs:8.3f}'} "
              f"{'':>9s} {'':>8s} {'NA (unstable)':>18s} {len(reps):>6d}")
        continue
    reps.sort()
    lo = reps[int(0.025 * len(reps))]
    hi = reps[int(0.975 * len(reps))]
    mean = sum(reps) / len(reps)
    sd = math.sqrt(sum((x - mean) ** 2 for x in reps) / (len(reps) - 1))
    znorm = z_obs / sd if sd > 0 else float("nan")
    inter_res[cov] = (z_obs, mean, sd, (lo, hi, znorm), len(reps))
    print(f"{cov:>4s} {z_obs:8.3f} {mean:9.3f} {sd:8.3f} "
          f"[{lo:+.3f},{hi:+.3f}] {len(reps):>6d}  z_norm={znorm:+.2f}")

print("\n=== per-instance interaction decomposition ===")
for inst in INSTANCES:
    out = []
    for cov in COVS:
        zs = {}
        ok = True
        for arm in ["A", "B", "E", "C_relaxed"]:
            pts = list(idx[(arm, cov, inst)].values())
            r = pearson(pts)
            if r is None:
                ok = False; break
            zs[arm] = fisher(r)
        z = (zs["A"] - zs["B"]) - (zs["E"] - zs["C_relaxed"]) if ok else None
        out.append(f"{cov}x:{'%+.3f' % z if z is not None else 'NA'}")
    print(f"instance {inst:>9s}: " + "  ".join(out))

print("\nwrote", f"{MS}/all_instances.tsv")
