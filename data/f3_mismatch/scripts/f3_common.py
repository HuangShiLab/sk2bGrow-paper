#!/usr/bin/env python3
"""F3 shared helpers: enzyme panel parsing, anchors.bin I/O, tag census.

Conventions (all cross-checked against the Rust source):
  * anchor tag = contig[position : position + tag_len] forward strand
    (sk2bgrow-core digest.rs).
  * anchors.bin = bincode(Vec<Anchor>, Vec<[u8;12]>): u64 len + n*26 bytes
    per anchor (u64 seq_hash, u32 genome_id, u16 contig_id, u64 position,
    u8 enzyme_idx, u8 strand, u8 flags, u8 local_gc), then u64 len + n*12
    packed tag bytes (2-bit, A=0 C=1 G=2 T=3, MSB first — tgt.rs pack_bases).
  * flags: UNIQUE_IN_GENOME=1, UNIQUE_ACROSS_DB=2, MASKED_MULTICOPY=4,
    MASKED_SHARED=8, NON_CHROMOSOMAL=16, GC_UNDEFINED=32.
  * near-census counts pairs of landmark *instances* with Hamming distance
    d in [1, 2]; cross-tag-length pairs are not enumerated (a read tag always
    has exactly one length, so cross-length confusion is not a counting
    channel). split: uu both unique-in-genome, um exactly one, mm none.
"""
import bisect
import json
import pickle
import struct

BASE = "/lustre1/g/aos_shihuang/sk2bgrow-hpc"
F2 = f"{BASE}/bench/F2"
F3 = f"{BASE}/bench/F3"
ENZYME_RS = f"{BASE}/src/crates/sk2bgrow-core/src/enzyme.rs"

FLAG_UNIQUE_IN_GENOME = 1 << 0
FLAG_MASKED_MULTICOPY = 1 << 2

IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "Y": "CT", "R": "AG", "W": "AT", "S": "CG", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}
COMP = str.maketrans("ACGT", "TGCA")


def revcomp(s):
    return s.translate(COMP)[::-1]


def canonical(s):
    r = revcomp(s)
    return s if s <= r else r


def hamming(a, b):
    return sum(x != y for x, y in zip(a, b))


def hamming_within(tag, probe, budget):
    """count.rs hamming_within: distance if <= budget and same length."""
    if len(tag) != len(probe):
        return None
    d = 0
    for x, y in zip(tag, probe):
        if x != y:
            d += 1
            if d > budget:
                return None
    return d


def seed_ranges(length, n):
    """count.rs seed_ranges: split into n contiguous seeds, remainder to the
    trailing seeds."""
    base, rem = divmod(length, n)
    out, lo = [], 0
    for k in range(n):
        sz = base + (1 if k >= n - rem else 0)
        out.append((lo, lo + sz))
        lo += sz
    return out


def read_fasta(path):
    seqs = {}
    order = []
    with open(path) as fh:
        name, cur = None, []
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(cur)
                name = line[1:].split()[0]
                order.append(name)
                cur = []
            else:
                cur.append(line.strip().upper())
        if name is not None:
            seqs[name] = "".join(cur)
    return seqs, order


def parse_panel():
    """Parse enzyme.rs into [(idx, name, tag_len, patterns)] where patterns is
    a list of patterns, each a list of (offset, bases). Aborts loudly if the
    parse does not reproduce the expected 16-enzyme panel."""
    import re
    src = open(ENZYME_RS).read()
    pat_blocks = {}
    for m in re.finditer(
        r'static (\w+?)_P: \[Pattern; \d+\] = \[(.*?)\];', src, re.S
    ):
        name, body = m.group(1), m.group(2)
        pats = []
        for pm in re.finditer(r"motifs!\[(.*?)\]", body, re.S):
            motifs = [
                (int(o), b)
                for o, b in re.findall(r'\((\d+),\s*"([A-Z]+)"\)', pm.group(1))
            ]
            pats.append(motifs)
        pat_blocks[name] = pats
    panel = []
    for m in re.finditer(
        r'Enzyme \{\s*idx: (\d+),\s*name: "(\w+)",\s*tag_len: (\d+),\s*'
        r"patterns: &(\w+)_P,",
        src,
    ):
        idx, name, tl, ref = int(m.group(1)), m.group(2), int(m.group(3)), m.group(4)
        panel.append((idx, name, tl, pat_blocks[ref]))
    panel.sort()
    assert len(panel) == 16, f"panel parse failed: {len(panel)} enzymes"
    assert [p[0] for p in panel] == list(range(16))
    return panel


def window_matches_patterns(window, patterns):
    for pat in patterns:
        ok = True
        for off, bases in pat:
            sub = window[off:off + len(bases)]
            if len(sub) != len(bases):
                ok = False
                break
            for wb, mb in zip(sub, bases):
                if wb not in IUPAC[mb]:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            return True
    return False


def load_anchors(dbdir):
    """Parse anchors.bin -> (anchors list of dicts, tags list of str decoded
    from the packed representation)."""
    data = open(f"{dbdir}/anchors.bin", "rb").read()
    n = struct.unpack_from("<Q", data, 0)[0]
    off = 8
    anchors = []
    for _ in range(n):
        seq_hash, genome_id, contig_id, position, enzyme_idx, strand, flags, gc = (
            struct.unpack_from("<QIHQBBBB", data, off)
        )
        off += 26
        anchors.append(
            dict(contig_id=contig_id, position=position, enzyme_idx=enzyme_idx,
                 flags=flags)
        )
    m = struct.unpack_from("<Q", data, off)[0]
    off += 8
    assert m == n, "tag count mismatch"
    LUT = "ACGT"
    tags = []
    for i in range(n):
        packed = data[off + 12 * i: off + 12 * i + 12]
        tags.append("".join(
            LUT[(packed[j // 4] >> (6 - 2 * (j % 4))) & 0b11]
            for j in range(48)
        ))
    assert off + 12 * n == len(data), "trailing bytes in anchors.bin"
    return anchors, tags


def load_manifest(dbdir):
    return json.load(open(f"{dbdir}/manifest.json"))


def genomes():
    """-> [(gid, gc)] from the F2 genome table."""
    out = []
    with open(f"{F2}/genomes.tsv") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        gi, gci = header.index("genome_id"), header.index("gc")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            out.append((f[gi], float(f[gci])))
    return out


def fmh_scale_k16(gid):
    with open(f"{F2}/f2_scales.tsv") as fh:
        next(fh)
        for line in fh:
            g, panel, scale, _mb = line.rstrip("\n").split("\t")
            if g == gid and panel == "16":
                return int(scale)
    raise KeyError(gid)


def load_fmh_loci(gid):
    """-> (loci dict key->[(cid,pos,gc)...], lens list) from the F2 sketch
    cache at the k16-matched scale."""
    s = fmh_scale_k16(gid)
    with open(f"{F2}/sketch_cache/{gid}_s{s}.pkl", "rb") as fh:
        loci, lens = pickle.load(fh)
    return loci, lens, s


def near_census(tags, is_unique):
    """Enumerate all same-length pairs with 1 <= Hamming <= 2 via a 3-piece
    exact-share index (pigeonhole: two strings at distance <=2 over 3
    contiguous pieces must share one whole piece). Returns
    counts dict {(d, cls): n} with cls in uu/um/mm."""
    n = len(tags)
    counts = {}
    if n < 2:
        return counts
    L = len(tags[0])
    assert all(len(t) == L for t in tags), "near_census needs uniform length"
    pieces = seed_ranges(L, 3)
    index = {}
    for i, t in enumerate(tags):
        for lo, hi in pieces:
            index.setdefault(t[lo:hi], []).append(i)
    cands = set()
    for v in index.values():
        if len(v) > 1:
            v = sorted(v)
            for a in range(len(v) - 1):
                va = v[a]
                for b in range(a + 1, len(v)):
                    cands.add((va, v[b]))
    for i, j in cands:
        d = hamming(tags[i], tags[j])
        if 1 <= d <= 2:
            u = (1 if is_unique[i] else 0) + (1 if is_unique[j] else 0)
            cls = "uu" if u == 2 else ("um" if u == 1 else "mm")
            counts[(d, cls)] = counts.get((d, cls), 0) + 1
    return counts


def lookup_replica(query, exact, seeds, tags_by_id, budget):
    """Replicates count.rs AnchorIndex::lookup on canonical strings.

    exact: dict canon(tag) -> [id]. seeds: list (per slot) of dict
    substring -> [id] or None when budget == 0. tags_by_id: id -> tag str.
    Returns list of (id, dist)."""
    out = []
    c = canonical(query)
    for i in exact.get(c, ()):
        tag = tags_by_id[i]
        if tag == query or tag == revcomp(query):
            out.append((i, 0))
    if out or budget == 0:
        return out
    seen = set()
    rc = revcomp(query)
    for probe in (query, rc):
        for slot in seeds:
            for lo_hi, table in slot:
                sub = probe[lo_hi[0]:lo_hi[1]]
                for i in table.get(sub, ()):
                    if i in seen:
                        continue
                    seen.add(i)
                    d = hamming_within(tags_by_id[i], probe, budget)
                    if d is not None:
                        out.append((i, d))
    return out


def build_seeds(tags_by_id, lengths, budget):
    """Seeds for one mismatch budget: budget+1 slots, each a list of
    (range, dict substring->ids) mirroring count.rs (per-tag seed ranges)."""
    n_slots = budget + 1
    ranges = {L: seed_ranges(L, n_slots) for L in set(lengths.values())}
    out = []
    for k in range(n_slots):
        tables = {}
        per_len = {}
        for L, rs in ranges.items():
            per_len[L] = rs[k]
            tables[rs[k]] = {}
        for i, tag in tags_by_id.items():
            r = per_len[lengths[i]]
            tables[r].setdefault(tag[r[0]:r[1]], []).append(i)
        out.append([(r, tables[r]) for r in sorted(tables)])
    return out
