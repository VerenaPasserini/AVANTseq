# AVANTseq v1.0.1 – bug fixes

## Fixed
- **Pipelines now start from the repository root**: `AVANTseq.smk` reads `config/config.yaml` + `config/samples_tumor.yaml`; `CreatePoN.smk` reads `config/config.yaml` + `config/samples_normal.yaml` (previously pointed to non-existent files).
- **CreatePoN.smk target paths**: `rule all` now requests `variants/mutect2/pon/...`, matching the outputs of `rules/pon.smk` (previously failed with MissingInputException).
- **Mutect2 now uses the merged PoN** (`merged_pon`, custom + public) in `AVANTseq.smk`, as documented. Previously it used the public PoN (`pon`) only.
- **Thread handling**: `trim`, `align`, `fastqc` and `multibamsummary` declare `threads:` and use `{threads}`, so Snakemake no longer oversubscribes CPUs.
- **Logging**: bwa messages are now captured in the alignment log; `qc_stats` and `vt_normalize_decompose` now write log files.
- **multiBamSummary** uses the declared BAM inputs instead of a `*.bam` wildcard, and writes the required `--outFileName` matrix (`CoverageSummary.npz`).
- **FastQC** now runs on the raw FASTQ files (R1 and R2) instead of aligned BAMs.
- **vt normalization**: temporary file and `.tbi` index are declared as outputs.

## Changed
- New config option `ref_version` (default `hg38`) used by Funcotator.
- Duplicate-marking metrics declared as an output of `mark_duplicates`.
- Removed unused imports; corrected header comments (no CNV step; tumor-only mode).

## Documentation
- README: workflow overview (rule graph), run order (CreatePoN first), updated requirements (deepTools, tabix; removed unused bcftools), Snakemake ≥ 8 badge.
- docs: corrected output paths and descriptions (FastQC per read, CalculateContamination segments, PoN paths), rule graph captions.
