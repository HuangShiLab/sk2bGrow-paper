# census of synthetic-genome counts table: unique keys, genomes, nonzero anchors
NR > 1 {
    key = $2 "\t" $3
    if (!(key in seen)) {
        seen[key] = 1
        gsum[$3] += $4
        if ($4 > 0) anch[$1 "\t" $2] = 1
    }
}
END {
    for (x in anch) na++
    c = 0
    for (x in gsum) if (gsum[x] > 0) c++
    print "unique_keys=" length(seen) " genomes=" length(gsum) " nonzero_anchors=" na " genomes_with_counts=" c
}
