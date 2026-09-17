#!/usr/bin/env python3
"""F3 step 1: near-collision census for one genome, both landmark modes.

Enzyme mode: anchors from the F2 k16 index db; tag sequences extracted from
the reference FASTA by anchor coordinate (contig[position:position+tag_len],
forward strand), cross-checked against the packed tags stored in anchors.bin.
Census per tag-length group (cross-length pairs are not a counting channel).

FracMinHash mode: landmark instances from the F2 sketch cache at the scale
matched to the k16 panel (f2_scales.tsv, match_panel=16); tag = reference
31-mer at each key locus; unique = key occurs at exactly one locus (armE
convention).

Output part: bench/F3/parts/nn_{gid}.tsv

Run with the micromamba python on a compute node.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from f3_common import (  # noqa: E402
    F2, F3, FLAG_MASKED_MULTICOPY, FLAG_UNIQUE_IN_GENOME, genomes, load_anchors,
    load_fmh_loci, load_manifest, near_census, parse_panel, read_fasta,
)

FIELDS = [
    "genome", "gc", "mode", "tag_len", "n_landmarks", "n_unique",
    "pairs_le1_uu", "pairs_le1_um", "pairs_le1_mm",
    "pairs_le2_uu", "pairs_le2_um", "pairs_le2_mm",
]


def census_groups(tags, is_unique):
    groups = {}
    for t, u in zip(tags, is_unique):
        groups.setdefault(len(t), []).append((t, u))
    rows = []
    for L in sorted(groups):
        gt = [x[0] for x in groups[L]]
        gu = [x[1] for x in groups[L]]
        c = near_census(gt, gu)
        row = dict(tag_len=L, n_landmarks=len(gt),
                   n_unique=sum(1 for x in gu if x))
        for cls in ("uu", "um", "mm"):
            d1 = c.get((1, cls), 0)
            row[f"pairs_le1_{cls}"] = d1
            row[f"pairs_le2_{cls}"] = d1 + c.get((2, cls), 0)
        rows.append(row)
    return rows


def enzyme_census(gid, fasta_seqs):
    dbdir = f"{F2}/db/{gid}_k16"
    man = load_manifest(dbdir)
    contigs = man["genomes"][0]["contigs"]
    assert len(man["genomes"]) == 1
    panel = parse_panel()
    tag_len = {idx: tl for idx, _n, tl, _p in panel}
    anchors, packed = load_anchors(dbdir)
    tags = []
    bad_packed = 0
    for a, pk in zip(anchors, packed):
        cname = contigs[a["contig_id"]]["name"]
        seq = fasta_seqs[cname]
        L = tag_len[a["enzyme_idx"]]
        t = seq[a["position"]:a["position"] + L]
        if t != pk[:L]:
            bad_packed += 1
        tags.append(t)
    if bad_packed:
        print(f"WARN {gid}: {bad_packed} tags disagree with packed tags",
              flush=True)
    is_unique = [bool(a["flags"] & FLAG_UNIQUE_IN_GENOME) for a in anchors]
    n_multi = sum(1 for a in anchors if a["flags"] & FLAG_MASKED_MULTICOPY)
    print(f"{gid} enzyme: {len(tags)} anchors ({n_multi} multi-copy)",
          flush=True)
    return census_groups(tags, is_unique)


def fmh_census(gid, fasta_seqs, fasta_order):
    loci, lens, scale = load_fmh_loci(gid)
    tags, is_unique, skipped = [], [], 0
    for key, places in loci.items():
        u = len(places) == 1
        for cid, pos, _gc in places:
            if cid >= len(fasta_order):
                skipped += 1
                continue
            t = fasta_seqs[fasta_order[cid]][pos:pos + 31]
            if len(t) != 31 or any(b not in "ACGT" for b in t):
                skipped += 1
                continue
            tags.append(t)
            is_unique.append(u)
    print(f"{gid} fmh(s={scale}): {len(tags)} landmark instances "
          f"({skipped} skipped)", flush=True)
    return census_groups(tags, is_unique)


def main():
    gid = sys.argv[1]
    gc = dict(genomes())[gid]
    fasta_seqs, fasta_order = read_fasta(f"{F2}/genomes/{gid}.fna")
    os.makedirs(f"{F3}/parts", exist_ok=True)
    out = f"{F3}/parts/nn_{gid}.tsv"
    rows = []
    for mode, rrows in (("enzyme_k16", enzyme_census(gid, fasta_seqs)),
                        ("fmh_k16matched", fmh_census(gid, fasta_seqs,
                                                      fasta_order))):
        for r in rrows:
            r.update(genome=gid, gc=gc, mode=mode)
            rows.append(r)
    with open(out, "w") as fh:
        fh.write("\t".join(FIELDS) + "\n")
        for r in rows:
            fh.write("\t".join(str(r[f]) for f in FIELDS) + "\n")
    print(f"NN_PART_DONE {gid}", flush=True)


if __name__ == "__main__":
    main()
