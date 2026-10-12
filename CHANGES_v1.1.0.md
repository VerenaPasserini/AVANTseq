# AVANTseq v1.1.0 – GATK best-practice alignment

## New
- **Base quality score recalibration (BQSR)** in both pipelines: GATK BaseRecalibrator (restricted to the target regions + 50 bp) and ApplyBQSR. Requires the new `known_sites` list in `config/config.yaml` (e.g. dbSNP, Mills & 1000G indels).
- **Read orientation bias filtering**: Mutect2 now collects F1R2 counts, LearnReadOrientationModel learns orientation bias priors and FilterMutectCalls uses them (`--ob-priors`). Removes artefacts such as FFPE deamination and oxoG damage.
- **fastp replaces Atropos** for trimming: automatic paired-end adapter detection, 3' quality trimming and length filtering. Settings are configurable in `config/config.yaml` (`fastp_cut_tail_quality`, `fastp_min_length`, `fastp_extra`). fastp reports are included in MultiQC.
- **Coverage QC report**: CollectHsMetrics of all samples in one table and plot, with a PASS / FLAG / FAIL class per sample (% target bases > 30x; thresholds configurable). FAIL samples are left out of the variant report.
- **Final report**: the per-sample MAFs are filtered (covered regions from the BED file, depth, alt reads, VAF, gnomAD population frequency), merged and summarised with maftools (summary plot and oncoplot in `report/`). Filters are configurable in `config/config.yaml`.
- **Conda environment** with all required tools, including Snakemake: `bash install_env.sh` (see `envs/avantseq.yaml`).

## Fixed
- Adapters were not trimmed in previous versions (Atropos ran with `--no-default-adapters` and no adapter sequence).
- vt: multi-allelic records are now decomposed first (`vt decompose -s`, which also splits per-allele fields such as AD and AF) and then normalized, as recommended by vt.

## Changed
- CreateSomaticPanelOfNormals uses the gnomAD germline resource (`--germline-resource`), as in the GATK PoN workflow.
- GenomicsDBImport uses the same 50 bp interval padding as Mutect2 and `--merge-input-intervals true` (much faster with many panel intervals).
- Read groups now include a library tag (`LB`).
- MultiQC also collects MarkDuplicates and BQSR metrics; duplicate metrics moved to `alignment/qc/markdup/`.
- The `{sample}` wildcard is constrained to the names in the sample file.

## Documentation
- README and docs updated (fastp, BQSR, orientation model, new outputs and config keys), new "Notes and limitations" sections (tumor-only calling, merged PoN rationale, germline/PoN sites in the VCF, contamination on small panels).
- Workflow overview diagram and rule graphs updated.

## Upgrading from v1.0.x
- Add `known_sites` (and optionally the `fastp_*` settings) to your `config/config.yaml`.
- Install fastp; Atropos is no longer needed.
- Results differ from v1.0.x (BQSR, adapter trimming, orientation filtering), so re-run both pipelines, starting with CreatePoN.smk.
