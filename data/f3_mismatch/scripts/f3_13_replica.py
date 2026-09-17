#!/usr/bin/env python3
"""F3 step 4: coordinate attribution (replica) for one genome.

The sk2bgrow counting layer does not emit the coordinates a read window was
credited to, only per-anchor counts. This script replicates the full counting
rule of sk2bgrow-core count.rs — every read offset at every tag length is
tested against the enzyme patterns (IUPAC, ported from enzyme.rs), exact
canonical lookup first (an exact hit suppresses the near-neighbour search),
then seed-based search with the same seed_ranges / best-distance /
record-all-best-anchors semantics, keep_multimappers=true — on the
known-source simulated reads. Every recorded observation is attributed as:

  correct     all credited landmark tags are identical to the reference tag
              at the true locus (Hamming 0; includes identical-tag
              multi-copy / multi-enzyme ties, which are masked at index time
              and handled by the EM layer, not an F3 near-collision risk)
  wrong_tie   the true tag is among the credited tags, but a tag at Hamming
              >= 1 from it is credited at the same best distance
  wrong_only  the true tag is not credited at all (a neighbour is strictly
              closer to the read window)

Enzyme mode is validated against the real Rust runs (bench/F3/cells): read
and tag counters, mismatch histogram, and the full per-anchor count table
must match exactly (parts/validate_{gid}.tsv). FMH mode has no Rust counter —
the sketch mode does not exist in the CLI yet — so it is emulated with the
same counting rule over the FracMinHash keys at the k16-matched scale; that
is the rule a future `sk2bgrow index --mode fracminhash` would run. Caveat is
reported in REVIEW.md.

Idempotent: existing parts are skipped.
"""
import bisect
import gzip
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from f3_common import (  # noqa: E402
    F2, F3, IUPAC, build_seeds, canonical, genomes, load_anchors,
    load_manifest, load_fmh_loci, lookup_replica, parse_panel, read_fasta,
)

ERRS = ["e0", "e0001", "e01"]
MMS = [0, 1, 2]
READLEN = 150
KFMH = 31

FIELDS = [
    "genome", "gc", "mode", "err", "mm",
    "n_windows", "n_motif_passed", "n_recorded", "n_correct",
    "n_wrong_tie", "n_wrong_only", "frac_wrong_recorded",
]


def read_fq(path):
    with gzip.open(path, "rt") as fh:
        while True:
            h = fh.readline().rstrip("\n")
            if not h:
                return
            seq = fh.readline().rstrip("\n")
            fh.readline()
            fh.readline()
            yield h[1:], seq


def make_motif_regexes(panel, tl):
    """One full-length regex per pattern of every enzyme with tag length tl:
    '.' anywhere, [allowed] at motif positions. Mirrors match_window."""
    out = []
    for _idx, _nm, etl, pats in panel:
        if etl != tl:
            continue
        for pat in pats:
            chars = ["."] * tl
            for off, bases in pat:
                for k, b in enumerate(bases):
                    chars[off + k] = "[" + IUPAC[b] + "]"
            out.append(re.compile("".join(chars)))
    return out


def wrap_slice(ref, q, tl):
    L = len(ref)
    if q + tl <= L:
        return ref[q:q + tl]
    return ref[q:] + ref[:q + tl - L]


class EnzymeMode:
    def __init__(self, gid):
        self.name = "enzyme_k16"
        man = load_manifest(f"{F2}/db/{gid}_k16")
        contigs = man["genomes"][0]["contigs"]
        self.panel = parse_panel()
        self.tag_len = {idx: tl for idx, _n, tl, _p in self.panel}
        self.enz_name = {idx: nm for idx, nm, _tl, _p in self.panel}
        anchors, packed = load_anchors(f"{F2}/db/{gid}_k16")
        self.n_anchors = len(anchors)
        self.locus = [(a["contig_id"], a["position"]) for a in anchors]
        self.enz_idx = [a["enzyme_idx"] for a in anchors]
        self.tags_by_id = {i: packed[i][:self.tag_len[a["enzyme_idx"]]]
                           for i, a in enumerate(anchors)}
        self.offsets = [0]
        for c in contigs:
            self.offsets.append(self.offsets[-1] + c["length"])
        # motif regexes per tag length (count.rs by_len)
        self.regexes = {}
        for tl in sorted({tl for _i, _n, tl, _p in self.panel}):
            self.regexes[tl] = make_motif_regexes(self.panel, tl)
        self.indexes = {}

    def load_ref(self, fasta_seqs, fasta_order, contig_names):
        self.ref = "".join(fasta_seqs[n] for n in contig_names)

    def index(self, mm):
        if mm not in self.indexes:
            exact = {}
            for i, t in self.tags_by_id.items():
                exact.setdefault(canonical(t), []).append(i)
            lengths = {i: len(t) for i, t in self.tags_by_id.items()}
            seeds = build_seeds(self.tags_by_id, lengths, mm) if mm > 0 else None
            self.indexes[mm] = (exact, seeds)
        return self.indexes[mm]

    def scan(self, s, rc, seq):
        """Every motif-passing window, like count.rs: (tl, off, window, q)
        with q = concat coordinate of the window's reference position."""
        Lg = self.offsets[-1]
        for tl, regexes in self.regexes.items():
            for off in range(0, READLEN - tl + 1):
                w = seq[off:off + tl]
                for rx in regexes:
                    if rx.match(w):
                        break
                else:
                    continue
                if rc:
                    q = s + READLEN - tl - off
                else:
                    q = s + off
                if q >= Lg:
                    q -= Lg
                yield tl, off, w, q

    def process(self, scanned, mm, acc):
        exact, seeds = self.index(mm)
        for tl, off, w, q in scanned:
            acc["covered"] += 1
            acc["motif_passed"] += 1
            if any(b not in "ACGT" for b in w):
                acc["ambiguous"] += 1
                continue
            hits = lookup_replica(w, exact, seeds, self.tags_by_id, mm)
            if not hits:
                acc["unmatched"] += 1
                continue
            best = min(d for _i, d in hits)
            best_ids = [i for i, d in hits if d == best]
            loci = frozenset(self.locus[i] for i in best_ids)
            if len(loci) > 1:
                acc["multi_locus"] += 1
            if len(best_ids) > 1 and len(loci) == 1:
                acc["multi_enzyme"] += 1
            for i in best_ids:
                acc["counts"][i] += 1
            acc["recorded"] += 1
            acc["hist"][best] += 1
            true_canon = canonical(wrap_slice(self.ref, q, tl))
            acc["events"].append(
                (true_canon,
                 frozenset(canonical(self.tags_by_id[i]) for i in best_ids)))


class FmhMode:
    def __init__(self, gid, fasta_seqs, fasta_order):
        self.name = "fmh_k16matched"
        loci, lens, scale = load_fmh_loci(gid)
        self.scale = scale
        self.contig_seqs = [fasta_seqs[n] for n in fasta_order]
        self.ref = "".join(self.contig_seqs)
        self.offsets = [0]
        for L in lens:
            self.offsets.append(self.offsets[-1] + L)
        self.key_loci = {k: [(c, p) for c, p, _g in v] for k, v in loci.items()}
        self.tags_by_id = {}
        self.key_of = {}
        for j, (key, places) in enumerate(loci.items()):
            cid, pos, _g = places[0]
            if cid >= len(fasta_order):
                continue
            t = self.contig_seqs[cid][pos:pos + KFMH]
            if len(t) != KFMH or any(b not in "ACGT" for b in t):
                continue
            self.tags_by_id[j] = t
            self.key_of[j] = key
        self.indexes = {}

    def index(self, mm):
        if mm not in self.indexes:
            exact = {}
            for i, t in self.tags_by_id.items():
                exact.setdefault(canonical(t), []).append(i)
            lengths = {i: KFMH for i in self.tags_by_id}
            seeds = build_seeds(self.tags_by_id, lengths, mm) if mm > 0 else None
            self.indexes[mm] = (exact, seeds)
        return self.indexes[mm]

    def scan(self, s, rc, seq):
        Lg = self.offsets[-1]
        for off in range(0, READLEN - KFMH + 1):
            yield off, seq[off:off + KFMH], None

    def process(self, scanned, s, rc, mm, acc):
        exact, seeds = self.index(mm)
        Lg = self.offsets[-1]
        for off, w, _q in scanned:
            acc["windows"] += 1
            if any(b not in "ACGT" for b in w):
                continue
            hits = lookup_replica(w, exact, seeds, self.tags_by_id, mm)
            if not hits:
                acc["unmatched"] += 1
                continue
            best = min(d for _i, d in hits)
            best_ids = [i for i, d in hits if d == best]
            loci = set()
            for i in best_ids:
                loci |= set(self.key_loci[self.key_of[i]])
            if len(loci) > 1:
                acc["multi_locus"] += 1
            if rc:
                p = s + READLEN - off - KFMH
            else:
                p = s + off
            if p >= Lg:
                p -= Lg
            cid = bisect.bisect_right(self.offsets, p) - 1
            true_tag = wrap_slice(self.ref, p, KFMH)
            acc["recorded"] += 1
            acc["hist"][best] += 1
            acc["events"].append(
                (canonical(true_tag),
                 frozenset(canonical(self.tags_by_id[i]) for i in best_ids)))


def classify(events):
    """Wrong = a credited landmark carries a tag at Hamming >= 1 from the
    reference tag at the true locus. Identical-tag multi-copy ties are masked
    at index time and handled by EM — not an F3 near-collision."""
    n_correct = n_tie = n_only = 0
    for true_canon, credited in events:
        if credited == frozenset((true_canon,)):
            n_correct += 1
        elif true_canon in credited:
            n_tie += 1
        else:
            n_only += 1
    return n_correct, n_tie, n_only


def validate_rust(mode, gid, err, mm, acc):
    cdir = f"{F3}/cells/{gid}/{err}/mm{mm}"
    name = f"{gid}.{err}"
    stats = json.load(open(f"{cdir}/{name}.stats.json"))["counting"]
    rows = []

    def cmp(metric, rust_v, my_v):
        rows.append((gid, err, mm, metric, rust_v, my_v,
                     "OK" if rust_v == my_v else "MISMATCH"))

    cmp("reads_total", stats["reads_total"], acc["reads"])
    cmp("tag_matched", stats["tag_matched"], acc["recorded"])
    cmp("tag_multi_locus", stats["tag_multi_locus"], acc["multi_locus"])
    cmp("tag_multi_enzyme", stats["tag_multi_enzyme"], acc["multi_enzyme"])
    cmp("tag_unmatched", stats["tag_unmatched"], acc["unmatched"])
    cmp("tag_ambiguous", stats["tag_ambiguous"], acc["ambiguous"])
    for d in (0, 1, 2):
        cmp(f"mismatch_hist_{d}", stats["mismatch_hist"][d], acc["hist"][d])
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
    cmp("counts_total", sum(rust.values()), sum(mine.values()))
    diffkeys = [k for k in set(rust) | set(mine)
                if rust.get(k, 0) != mine.get(k, 0)]
    cmp("counts_diff_rows", len(diffkeys), 0)
    return rows


def new_acc(mode):
    return dict(reads=0, windows=0, covered=0, motif_passed=0, recorded=0,
                ambiguous=0, unmatched=0, multi_locus=0, multi_enzyme=0,
                hist={0: 0, 1: 0, 2: 0, 3: 0}, events=[],
                counts=([0] * mode.n_anchors
                        if mode.name.startswith("enzyme") else None))


def main():
    gid = sys.argv[1]
    gc = dict(genomes())[gid]
    fasta_seqs, fasta_order = read_fasta(f"{F2}/genomes/{gid}.fna")
    out = f"{F3}/parts/misassign_{gid}.tsv"
    if os.path.exists(out):
        print(f"skip {gid}")
        return
    enzyme = EnzymeMode(gid)
    man = load_manifest(f"{F2}/db/{gid}_k16")
    enzyme.load_ref(fasta_seqs, fasta_order,
                    [c["name"] for c in man["genomes"][0]["contigs"]])
    modes = [enzyme, FmhMode(gid, fasta_seqs, fasta_order)]
    rows, vrows = [], []
    for mode in modes:
        for err in ERRS:
            fq = f"{F3}/fq/{gid}.{err}.fq.gz"
            # scan once per (mode, err); lookups per mm reuse it
            per_read = []
            for h, seq in read_fq(fq):
                _g, s, rc, _i = h.split("|")
                per_read.append((int(s), int(rc), seq))
            for mm in MMS:
                acc = new_acc(mode)
                for s, rc, seq in per_read:
                    acc["reads"] += 1
                    scanned = list(mode.scan(s, rc, seq))
                    if mode.name.startswith("enzyme"):
                        mode.process(scanned, mm, acc)
                    else:
                        mode.process(scanned, s, rc, mm, acc)
                n_correct, n_tie, n_only = classify(acc["events"])
                n_rec = acc["recorded"]
                frac = ((n_tie + n_only) / n_rec) if n_rec else 0.0
                if mode.name.startswith("enzyme"):
                    vrows.extend(validate_rust(mode, gid, err, mm, acc))
                rows.append(dict(
                    genome=gid, gc=gc, mode=mode.name, err=err, mm=mm,
                    n_windows=acc["covered"] if mode.name.startswith("enzyme")
                    else acc["windows"],
                    n_motif_passed=acc["motif_passed"],
                    n_recorded=n_rec, n_correct=n_correct,
                    n_wrong_tie=n_tie, n_wrong_only=n_only,
                    frac_wrong_recorded=round(frac, 8),
                ))
                print(f"{gid} {mode.name} {err} mm{mm}: recorded={n_rec} "
                      f"wrong_tie={n_tie} wrong_only={n_only}", flush=True)
    with open(out, "w") as fh:
        fh.write("\t".join(FIELDS) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[f]) for f in FIELDS) + "\n")
    with open(f"{F3}/parts/validate_{gid}.tsv", "w") as fh:
        fh.write("genome\terr\tmm\tmetric\trust\treplica\tmatch\n")
        for v in vrows:
            fh.write("\t".join(map(str, v)) + "\n")
    print(f"REPLICA_PART_DONE {gid}", flush=True)


if __name__ == "__main__":
    main()
