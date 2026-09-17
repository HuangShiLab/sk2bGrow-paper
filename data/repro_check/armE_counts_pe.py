#!/usr/bin/env python3
"""Paired-end variant of armE_counts.py: pools both mates' k-mer counts into
ONE sk2bgrow count table (C1/Pilea PE convention), instead of one table per
fq. Reuses armE_counts' constants and Pilea's count64. Run with the Pilea
interpreter. Writes <outdir>/<name>.counts.tsv; skip-if-exists."""
import argparse, os, pickle, sys

sys.path.insert(0, "/lustre1/g/aos_shihuang/sk2bgrow-hpc/src/benches/zheng2020")
from armE_counts import FLAG_UNIQUE, FLAG_MULTICOPY, HEADER  # noqa: E402
from pilea.io import parse_fastx_file                      # noqa: E402
from pilea.kmc import count64                              # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("reads", nargs="+", help="paired mates (_1 then _2)")
ap.add_argument("-k", type=int, default=31)
ap.add_argument("-s", type=int, default=104)
ap.add_argument("--cache", required=True, help="reference sketch pickle")
ap.add_argument("-o", "--outdir", required=True)
ap.add_argument("-g", "--genome", default="ecoli")
ap.add_argument("--name", required=True, help="sample base for output/table")
a = ap.parse_args()

maxhash = ((1 << 64) - 1) // a.s
loci, lens = pickle.load(open(a.cache, "rb"))
os.makedirs(a.outdir, exist_ok=True)
out = os.path.join(a.outdir, f"{a.name}.counts.tsv")
if os.path.exists(out):
    print(f"skip {out}", file=sys.stderr)
    sys.exit(0)

offsets, run = [], 0
for length in lens:
    offsets.append(run)
    run += length

kmc = {}
for fq in a.reads:
    for record in parse_fastx_file(fq):
        count64(record, None, a.k, maxhash, kmc)

with open(out, "w") as fh:
    fh.write(HEADER + "\n")
    for key, places in loci.items():
        flags = FLAG_MULTICOPY if len(places) > 1 else FLAG_UNIQUE
        cnt = int(kmc.get(key, 0))
        for cid, pos, gc in places:
            fh.write(f"{a.name}\t0\t{a.genome}\t{cid}\t{pos}"
                     f"\t{offsets[cid] + pos}\tFracMinHash\t+\t{flags}"
                     f"\t{min(200, round(gc / 50))}\t4294967295\t{cnt}\n")
print(f"  {a.name}: {sum(1 for k in loci if kmc.get(k)):,} observed",
      file=sys.stderr)
