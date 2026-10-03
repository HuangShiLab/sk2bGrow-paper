# Raw-read ENA manifests, 2026-10-03

These files preserve the source-of-truth metadata needed to re-download all raw
paired-end FASTQs used by the sk2bGrow benchmarks.

- Source host path: `/lustre1/g/aos_shihuang/sk2bgrow-hpc/manifests/`
- BioProjects: PRJNA615952, PRJNA1280254, PRJNA551656, PRJNA974210
- The `*_counts.tsv` files contain read/base counts and library layout.
- The non-count `PRJNA*.tsv` files contain run accessions, ENA FTP URLs,
  expected byte counts, and MD5 checksums.

The public raw FASTQs are not stored in this repository. On 2026-10-03 the
complete final FASTQ set on HPC was 273.3 GB (254.56 GiB) across 206 files;
82 partial `.part` download fragments (4.67 GiB) were deleted. These manifests
are sufficient to regenerate the raw-read cache from ENA.
