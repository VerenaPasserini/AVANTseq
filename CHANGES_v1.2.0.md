# AVANTseq v1.2.0 – Reports, conda environment and germline filter fix

## New
- **Coverage QC report** (`report/hsmetrics_qc.pdf`, `report/hsmetrics_qc.tsv`): CollectHsMetrics of all samples in one table and plot, with a PASS / FLAG / FAIL class per sample based on the % of target bases covered > 30x (`qc_pass_30x`, `qc_flag_30x`).
- **Variant report** (`report/maf_summary.pdf`, `report/oncoplot.pdf`, `report/all_samples_filtered.maf`): the per-sample MAFs are filtered and merged, then summarised with maftools. Filters: covered regions of the `bed` file (plus `report_bed_padding`, default 10 bp, to keep splice sites), depth, alt reads, VAF, gnomAD population frequency (`POPAF`) and, optionally, features without a gene symbol. Samples with an excluded coverage QC class (default FAIL) are left out. Oncoplot with the top `report_top_genes` genes, a colour-blind-friendly palette and an automatic gene-name margin.
- **Cohort overview** (`report/cohort_flow.pdf`, `report/cohort_summary.tsv`): CONSORT-style diagram (samples sequenced, excluded for coverage with IDs, analysed, without reportable variants) and a per-sample table with QC class, status and variant counts after each filter.
- **Conda environment** with all required tools, including Snakemake and R/maftools: `bash install_env.sh` (see `envs/avantseq.yaml`).

## Fixed
- **Germline filter:** GetPileupSummaries again uses all common sites covered by reads (on- and off-target). Restricting it to the target regions in v1.1.0 produced a coarser tumor segmentation, so some heterozygous germline SNPs passed FilterMutectCalls.

## Upgrading from v1.1.0
- Add the report settings to your `config/config.yaml` (see the template): `qc_pass_30x`, `qc_flag_30x`, `report_exclude_qc`, `report_min_depth`, `report_min_alt`, `report_min_vaf`, `report_max_pop_af`, `report_drop_unnamed`, `report_top_genes`, `report_gene_margin`.
- R with maftools, data.table, ggplot2 and ggrepel must be available (included in the conda environment).
- Re-run `AVANTseq.smk`: Snakemake re-runs only GetPileupSummaries and the steps after it, then creates the reports. Alignment, BQSR and Mutect2 results are reused.
