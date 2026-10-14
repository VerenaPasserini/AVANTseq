# AVANTseq v1.1.0 – GATK best practices, QC and variant reports

## New
- **Base quality score recalibration (BQSR)** in both pipelines: GATK BaseRecalibrator (restricted to the target regions + 50 bp) and ApplyBQSR. Requires the new `known_sites` list in `config/config.yaml` (e.g. dbSNP, Mills & 1000G indels).
- **Read orientation bias filtering**: Mutect2 now collects F1R2 counts, LearnReadOrientationModel learns orientation bias priors and FilterMutectCalls uses them (`--ob-priors`). Removes artefacts such as FFPE deamination and oxoG damage.
- **fastp replaces Atropos** for trimming: automatic paired-end adapter detection, 3' quality trimming and length filtering. Settings are configurable in `config/config.yaml` (`fastp_cut_tail_quality`, `fastp_min_length`, `fastp_extra`). fastp reports are included in MultiQC.
- **Coverage QC report** (`report/hsmetrics_qc.pdf`, `report/hsmetrics_qc.tsv`): CollectHsMetrics of all samples in one table and plot, with a PASS / FLAG / FAIL class per sample based on the % of target bases covered > 30x (`qc_pass_30x`, `qc_flag_30x`).
- **Variant report** (`report/maf_summary.pdf`, `report/oncoplot.pdf`, `report/all_samples_filtered.maf`): the per-sample MAFs are filtered and merged, then summarised with maftools. Filters: covered regions of the `bed` file (plus `report_bed_padding`, default 10 bp, to keep splice sites), depth, alt reads, VAF, gnomAD population frequency (`POPAF`) and, optionally, features without a gene symbol. Samples with an excluded coverage QC class (default FAIL) are left out. Oncoplot with the top `report_top_genes` genes, a colour-blind-friendly palette and an automatic gene-name margin.
- **Cohort overview** (`report/cohort_flow.pdf`, `report/cohort_summary.tsv`): CONSORT-style diagram (samples sequenced, excluded for coverage with IDs, analysed, without somatic variants passing the filters) and a per-sample table with QC class, status and variant counts after each filter.
- **Conda environment** with all required tools, including Snakemake and R/maftools: `bash install_env.sh` (see `envs/avantseq.yaml`).

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
- README and docs updated (fastp, BQSR, orientation model, reports, new outputs and config keys), new "Notes and limitations" sections (tumor-only calling, merged PoN rationale, germline/PoN sites in the VCF, contamination on small panels).
- Workflow overview diagram (now including the report) and rule graphs updated.

## Upgrading from v1.0.x
- Update your `config/config.yaml` from the template: add `known_sites` (BQSR), and optionally the `fastp_*`, `qc_*` and `report_*` settings (defaults are used if missing).
- Install fastp and R with maftools, data.table, ggplot2 and ggrepel (or use `bash install_env.sh`); Atropos is no longer needed.
- Results differ from v1.0.x (BQSR, adapter trimming, orientation filtering), so re-run both pipelines, starting with CreatePoN.smk.
