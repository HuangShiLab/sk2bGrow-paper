#!/usr/bin/env python3
"""F3 step 5: aggregate part tables into F3_nearneighbor.tsv,
F3_misassign.tsv, F3_validation.tsv. Prints a headline summary."""
import glob
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from f3_common import F3  # noqa: E402


def concat(pattern, out):
    files = sorted(glob.glob(pattern))
    header = None
    rows = []
    for f in files:
        with open(f) as fh:
            h = fh.readline().rstrip("\n")
            if header is None:
                header = h
            elif h != header:
                raise SystemExit(f"header mismatch: {f}")
            rows.extend(l.rstrip("\n") for l in fh if l.strip())
    with open(out, "w") as fh:
        fh.write(header + "\n")
        fh.write("\n".join(rows) + ("\n" if rows else ""))
    return len(rows)


def main():
    n_nn = concat(f"{F3}/parts/nn_*.tsv", f"{F3}/F3_nearneighbor.tsv")
    n_ma = concat(f"{F3}/parts/misassign_*.tsv", f"{F3}/F3_misassign.tsv")
    n_va = concat(f"{F3}/parts/validate_*.tsv", f"{F3}/F3_validation.tsv")
    print(f"nn rows={n_nn} misassign rows={n_ma} validate rows={n_va}")

    # validation summary
    bad = ok = 0
    with open(f"{F3}/F3_validation.tsv") as fh:
        next(fh)
        for line in fh:
            if line.rstrip("\n").split("\t")[-1] == "OK":
                ok += 1
            else:
                bad += 1
                print("MISMATCH:", line.rstrip("\n"))
    print(f"validation: {ok} OK, {bad} MISMATCH")

    # headline: near-neighbour totals per mode
    totals = {}
    with open(f"{F3}/F3_nearneighbor.tsv") as fh:
        header = next(fh).rstrip("\n").split("\t")
        mi = {h: i for i, h in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            t = totals.setdefault(f[mi["mode"]],
                                  dict(landmarks=0, uu1=0, um1=0, uu2=0,
                                       um2=0))
            t["landmarks"] += int(f[mi["n_landmarks"]])
            t["uu1"] += int(f[mi["pairs_le1_uu"]])
            t["um1"] += int(f[mi["pairs_le1_um"]])
            t["uu2"] += int(f[mi["pairs_le2_uu"]])
            t["um2"] += int(f[mi["pairs_le2_um"]])
    for mode, t in totals.items():
        print(f"{mode}: landmarks={t['landmarks']} "
              f"d1 uu={t['uu1']} um={t['um1']} | d2 uu={t['uu2']} um={t['um2']}")

    # headline: misassignment, pooled per (mode, err, mm)
    pool = {}
    with open(f"{F3}/F3_misassign.tsv") as fh:
        header = next(fh).rstrip("\n").split("\t")
        mi = {h: i for i, h in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            k = (f[mi["mode"]], f[mi["err"]], f[mi["mm"]])
            p = pool.setdefault(k, [0, 0])
            p[0] += int(f[mi["n_wrong_tie"]]) + int(f[mi["n_wrong_only"]])
            p[1] += int(f[mi["n_recorded"]])
    for (mode, err, mm), (w, r) in sorted(pool.items()):
        frac = w / r if r else 0.0
        print(f"{mode} {err} mm{mm}: wrong={w} recorded={r} "
              f"frac={frac:.3e}")
    print("AGGREGATE_DONE")


if __name__ == "__main__":
    main()
