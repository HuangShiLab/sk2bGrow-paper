#!/usr/bin/env python3
"""Debug: reproduce main() flow for g08 and dump per-(err,mm) count-row diffs."""
import json
import sys

sys.path.insert(0, "/lustre1/g/aos_shihuang/sk2bgrow-hpc/bench/F3/scripts")
import f3_13_replica as R
from f3_common import F2, F3, read_fasta, load_manifest

gid = "g08_Escherichia_coli"
mode = R.EnzymeMode(gid)
man = json.load(open(f"{F2}/db/{gid}_k16/manifest.json"))
seqs, order = read_fasta(f"{F2}/genomes/{gid}.fna")
mode.load_ref(seqs, order, [c["name"] for c in man["genomes"][0]["contigs"]])

for err in ["e0", "e0001"]:
    fq = f"{F3}/fq/{gid}.{err}.fq.gz"
    per_read = []
    for h, seq in R.read_fq(fq):
        _g, s, rc, _i = h.split("|")
        per_read.append((int(s), int(rc), seq))
    for mm in [0, 1, 2]:
        acc = R.new_acc(mode)
        for s, rc, seq in per_read:
            acc["reads"] += 1
            mode.process(list(mode.scan(s, rc, seq)), mm, acc)
        name = f"{gid}.{err}"
        cdir = f"{F3}/cells/{gid}/{err}/mm{mm}"
        rust = {}
        with open(f"{cdir}/{name}.counts.tsv") as fh:
            next(fh)
            for line in fh:
                f = line.rstrip("\n").split("\t")
                rust[(int(f[3]), int(f[4]), f[6])] = int(f[11])
        mine = {}
        for i, c in enumerate(acc["counts"]):
            if c:
                cid, pos = mode.locus[i]
                mine[(cid, pos, mode.enz_name[mode.enz_idx[i]])] = c
        diffs = [(k, rust.get(k, 0), mine.get(k, 0))
                 for k in set(rust) | set(mine)
                 if rust.get(k, 0) != mine.get(k, 0)]
        print(f"{err} mm{mm}: diffs={len(diffs)} "
              f"rust_total={sum(rust.values())} mine_total={sum(mine.values())}",
              flush=True)
        for d in diffs[:5]:
            print("   ", d, flush=True)
print("DEBUG_DONE")
