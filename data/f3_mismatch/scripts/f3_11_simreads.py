#!/usr/bin/env python3
"""F3 step 2: simulate 150 bp single-end reads from known coordinates.

Unlike the F2 simulator there is no planted PTR gradient (misassignment does
not need one): start positions are uniform on the concatenated sequence, reads
wrap at the end like F2, 50% reverse-complement, quality 'I'*150. Per-base
substitution errors at rate err are drawn i.i.d. (substituted base is a uniform
ACGT draw, possibly the same base). The true start coordinate and orientation
are encoded in the read name as @{gid}|{start}|{rc}|{i} — this is what the
replica uses to attribute credited coordinates.

Depth 0.5x at the read level: n = round(0.5 * L / 150).

Idempotent: existing outputs are skipped.
"""
import argparse
import gzip
import hashlib
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from f3_common import F2, F3, genomes, read_fasta  # noqa: E402

READLEN = 150
DEPTH = 0.5
ERRS = [("e0", 0.0), ("e0001", 0.001), ("e01", 0.01)]
BASES = np.frombuffer(b"ACGT", dtype=np.uint8)


def simulate(path, gid, outdir):
    seqs, order = read_fasta(path)
    seq = "".join(seqs[n] for n in order)
    L = len(seq)
    arr = np.frombuffer((seq + seq[:READLEN]).encode(), dtype=np.uint8)
    n = int(round(DEPTH * L / READLEN))
    comp = np.arange(256, dtype=np.uint8)
    for a, b in zip(b"ACGT", b"TGCA"):
        comp[a] = b
    for errtag, err in ERRS:
        out = os.path.join(outdir, f"{gid}.{errtag}.fq.gz")
        if os.path.exists(out):
            print(f"skip {gid} {errtag}", flush=True)
            continue
        seed = int(hashlib.sha256(f"F3|{gid}|{errtag}".encode()).hexdigest()[:16], 16)
        rng = np.random.default_rng(seed)
        starts = rng.integers(0, L, size=n)
        rc = rng.random(n) < 0.5
        tmp = out + ".tmp"
        with gzip.open(tmp, "wt") as fh:
            for i, s in enumerate(starts):
                read = arr[s:s + READLEN].copy()
                if err > 0:
                    m = rng.random(READLEN) < err
                    if m.any():
                        read[m] = BASES[rng.integers(0, 4, size=int(m.sum()))]
                rcs = bool(rc[i])
                if rcs:
                    read = comp[read][::-1]
                fh.write(f"@{gid}|{s}|{1 if rcs else 0}|{i}\n"
                         f"{read.tobytes().decode()}\n+\n{'I' * READLEN}\n")
        os.rename(tmp, out)
        print(f"ok {gid} {errtag} nreads={n}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default=None)
    a = ap.parse_args()
    outdir = f"{F3}/fq"
    os.makedirs(outdir, exist_ok=True)
    for gid, _gc in genomes():
        if a.only and gid not in a.only.split(","):
            continue
        path = f"{F2}/genomes/{gid}.fna"
        if not os.path.exists(path):
            print(f"MISSING {path}", file=sys.stderr)
            continue
        simulate(path, gid, outdir)
    print("SIMREADS_DONE", flush=True)


if __name__ == "__main__":
    main()
