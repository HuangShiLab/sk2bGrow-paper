# ECC saliva route-B benchmark

## Goal

Test whether sk2bGrow-derived growth-state features improve the published BcgI
2bRAD-M classifier for early-childhood caries (ECC), beyond the published
species-abundance baseline.

## Inputs

* Published Fig. 7c-f normalized species table and labels:
  `external/published_2b_abundance.tsv`, `external/published_2b_metadata.tsv`.
* Vendor GTDB species count table used for cross-checking:
  `external/vendor_feature_table_species_counts.tsv`.
* Raw clean BcgI 2bRAD FASTQs:
  `/lustre1/g/aos_shihuang/Strain2b/data/saliva_data/ECC_saliva/2b/data/Clean_data`.

There are 38 samples (19 ECC, 19 Health). The published paper reports
AUC 0.92 from 10-fold cross-validation.

## Current design

1. Reproduce the published species-abundance baseline with 10-fold CV.
2. Run sk2bGrow route B on each raw FASTQ (`--mode 2brad`,
   `--max-mismatch 0`, one BcgI stratum).
3. Build feature arms:
   * A: published species abundance;
   * B: sk2bGrow PTR/support/coverage;
   * C: published abundance plus sk2bGrow PTR/support/coverage.
4. Use identical stratified folds in every repeat, with all feature filtering
   and scaling inside the training fold. Missing PTR gets an explicit indicator
   and is never imputed as zero.

## Implementation note and pilot result

The pilot used the available coordinate-bearing UHGG BcgI database
(`P2_sigma/db_bcgI`, 4,698 MGYG genomes), not a custom sk2bGrow index of the
original 2bRAD-M GTDB tag database. This can under-represent oral taxa and is
recorded as a major interpretation caveat. We are evaluating whether a filtered
GTDB tag-coordinate database can be imported without full genome FASTAs.

Pilot S_8123218 gave 177 finite single-enzyme PTR estimates, but none passed the
route-A-derived QC gates. Reasons include low containment against the gut-biased
reference, fragmented MAG references, and dispersion/coverage thresholds. These
finite-but-QC-failed estimates may still be useful as ML features, but must not
be described as high-confidence biological PTRs. For the full run we retain both
`glm` and `sorted` fits and report PTR feature yield.

## HPC jobs

Full-array job: `4101589` (38 samples, 12 concurrent tasks). Per sample:

```bash
sk2bgrow profile ... --mode 2brad --max-mismatch 0 --no-stats
python -m sk2bgrow.cli profile ... --no-gc-correct --method glm
python -m sk2bgrow.cli profile ... --no-gc-correct --method sorted
```

Large intermediate count/window files are deleted after both fits to keep the
Lustre footprint small.
