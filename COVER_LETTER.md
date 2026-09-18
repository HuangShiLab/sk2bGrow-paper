# Cover letter draft

Dear Editor,

We submit “Coordinate-aware V-fitting enables shallow-depth PTR estimation from 2bRAD anchor panels” for consideration as a methods article in *NAR Genomics and Bioinformatics*.

Peak-to-trough ratio (PTR) methods infer bacterial replication from coverage gradients, but coverage sparsity and per-species depth often limit metagenomic applications. sk2bGrow uses a coordinate-aware windowed V-fitting framework over motif-defined 2bRAD anchor panels. On the Zheng *E. coli* growth-rate panel, the method returned usable estimates at 1–2× nominal depth, while retaining explicit uncertainty and bias caveats. The controlled attribution analysis shows that the principal gain comes from coordinate-aware fitting rather than landmark determinism: after density matching, 2bRAD anchors and FracMinHash landmarks were not significantly different in the paired bootstrap analysis, although only the 2bRAD formulation directly corresponds to a wet-lab reduced-representation protocol and supplies cross-enzyme consistency strata.

The paper also reports a signed fixed-origin negative control, a sensitivity analysis for residual GC correction, fragmentation and scaffolding experiments, real faecal-metagenome concordance, and a MAG-scale cost analysis. We deliberately distinguish what the benchmark establishes from what it does not: the method is not unbiased at low depth, the current implementation is slower than Pilea at large MAG-panel scale, and fragmented references require a conservative handling policy. These limitations are stated in the abstract, Results, and Discussion.

All figure-generating tables, provenance records, and analysis code are available in the repositories listed in the manuscript. The work has not been published elsewhere and is not under consideration by another journal. All authors have approved the submission.

We thank you for your consideration.

Sincerely,  
[Corresponding author]  
[Affiliation and contact details]  
on behalf of all authors

