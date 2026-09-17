#!/usr/bin/env python3
"""Reproduce the replica main() flow and dump count diffs in the same process
that also runs the FMH mode (suspect: cross-mode state)."""
import sys

sys.path.insert(0, "/lustre1/g/aos_shihuang/sk2bgrow-hpc/bench/F3/scripts")
import f3_13_replica as R
from f3_common import F2, F3, read_fasta, load_manifest
import json

gid = "g08_Escherichia_coli"
ERRS = ["e0", "e0001", "e01"]
MMS = [0, 1, 2]


def diffs_for(mode, err, mm, acc):
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
    return diffs, sum(rust.values()), sum(mine.values())


gc = 50.8
fasta_seqs, fasta_order = read_fasta(f"{F2}/genomes/{gid}.fna")
enzyme = R.EnzymeMode(gid)
man = json.load(open(f"{F2}/db/{gid}_k16/manifest.json"))
enzyme.load_ref(fasta_seqs, fasta_order,
                [c["name"] for c in man["genomes"][0]["contigs"]])
modes = [enzyme, R.FmhMode(gid, fasta_seqs, fasta_order)]
for mode in modes:
    for err in ERRS:
        fq = f"{F3}/fq/{gid}.{err}.fq.gz"
        per_read = []
        for h, seq in R.read_fq(fq):
            _g, s, rc, _i = h.split("|")
            per_read.append((int(s), int(rc), seq))
        for mm in MMS:
            acc = R.new_acc(mode)
            for s, rc, seq in per_read:
                acc["reads"] += 1
                scanned = list(mode.scan(s, rc, seq))
                if mode.name.startswith("enzyme"):
                    mode.process(scanned, mm, acc)
                else:
                    mode.process(scanned, s, rc, mm, acc)
            if mode.name.startswith("enzyme"):
                d, rt, mt = diffs_for(mode, err, mm, acc)
                print(f"MAINFLOW {err} mm{mm}: diffs={len(d)} "
                      f"rust={rt} mine={mt}", flush=True)
                for x in d[:5]:
                    print("   ", x, flush=True)
print("MAINFLOW_DONE")
