# Cover letter draft

Dear Editor,

We submit “Coordinate-aware V-fitting enables shallow-depth PTR estimation from 2bRAD anchor panels” for consideration as a Research article in *Microbiome*.

Peak-to-trough ratio (PTR) methods infer bacterial replication from coverage gradients, but sparse per-species coverage often limits metagenomic applications. sk2bGrow uses coordinate-aware windowed V-fitting over motif-defined 2bRAD anchor panels. On the Zheng *E. coli* growth-rate panel, the method returned usable estimates at 1–2× nominal sequencing depth, with explicit uncertainty and bias caveats. The controlled attribution analysis shows that the principal gain comes from coordinate-aware fitting rather than landmark determinism: after density matching, 2bRAD anchors and FracMinHash landmarks were not significantly different in the paired bootstrap analysis. Relative to a computational sketch, however, the 2bRAD formulation provides a wet-lab-realizable reduced-representation protocol and cross-enzyme consistency strata.

A post-hoc count-level factorial further separates shallow-depth gradient extraction from multi-reference landmark assignment. It shows that shared-anchor ambiguity is a distinct failure mode that sequencing depth does not by itself remove, and motivates reporting shared-anchor fractions and per-reference unique coverage in future multi-strain studies. We present this analysis as mechanistic and hypothesis-generating rather than as a read-level validation of any production pipeline.

The paper also reports a signed fixed-origin negative control, a sensitivity analysis for GC correction, fragmentation and scaffolding experiments, real faecal-metagenome concordance, and a MAG-scale cost analysis. We deliberately distinguish what the benchmark establishes from what it does not: the method is not unbiased at low depth, the current implementation is slower than Pilea at large MAG-panel scale, and fragmented references require a conservative handling policy. These limitations are stated in the abstract, Results, and Discussion.

We believe the work is well suited to *Microbiome* because it addresses a broad methodological bottleneck—microbial growth inference from shallow metagenomes—while evaluating performance against real microbiome datasets and reproducible community-scale benchmarks. All figure-generating tables, provenance records, and analysis code are available in the repositories listed in the manuscript. The work has not been published elsewhere and is not under consideration by another journal. All authors have approved the submission.

Suggested reviewers are listed in `REVIEWER_SUGGESTIONS.md`. Institutional contact details will be confirmed immediately before submission.

Sincerely,  
[Corresponding author]  
[Affiliation and contact details]  
on behalf of all authors
