# Snakemake Pipeline: Automated Somatic Variant Calling Workflow
# --------------------------------------------------------------
# This Snakemake pipeline automates somatic variant calling (SNVs and indels)
# in tumor-only mode from NGS data generated using hybrid capture-based
# targeted sequencing, using GATK Mutect2 with a custom Panel of Normals.
# --------------------------------------------------------------
# Author: Verena Passerini
# GitHub: https://github.com/VerenaPasserini/AVANTseq
# Last updated: October 2026 (v1.1.0)
# Snakemake version: 8.29.3
# --------------------------------------------------------------

# Configuration files (paths are relative to the directory snakemake is run from,
# i.e. the repository root). Override with --configfile if needed.
configfile: "config/config.yaml"
configfile: "config/samples_tumor.yaml"

from os.path import join
import re

# Restrict the {sample} wildcard to the sample names listed in the sample file,
# so that e.g. "{sample}.bam" cannot match "tumor1.dedup.bam"
wildcard_constraints:
    sample="|".join(re.escape(s) for s in config["samples"])

# Include rules
include: "rules/trim.smk"
include: "rules/align.smk"
include: "rules/qc.smk"
include: "rules/variants.smk"

# Pipeline output files
rule all:
    input:
        expand(join(config["work_dir"], "alignment/bams/{sample}.bam.bai"), sample=config["samples"]),
        expand(join(config["work_dir"], "alignment/qc/fastqc/{sample}_{read}_fastqc.html"), sample=config["samples"], read=["R1", "R2"]),
        join(config["work_dir"], "alignment/qc/multiqc_report.html"),
        join(config["work_dir"], "variants/qc/CoverageSummary.txt"),
        expand(join(config["work_dir"], "variants/qc/HsMetrics/{sample}_metrics.txt"), sample=config["samples"]),
        expand(join(config["work_dir"], "variants/mutect2/{sample}.vcf.gz"), sample=config["samples"]),
        expand(join(config["work_dir"], "variants/mutect2/{sample}.getpileupsummaries.table"), sample=config["samples"]),
        expand(join(config["work_dir"], "variants/mutect2/{sample}.calculatecontamination.table"), sample=config["samples"]),
        expand(join(config["work_dir"], "variants/mutect2/{sample}.segments.table"), sample=config["samples"]),
        expand(join(config["work_dir"], "variants/filtered/{sample}_filtered.vcf.gz"), sample=config["samples"]),
        expand(join(config["work_dir"], "variants/filtered/{sample}_filtered_norm_dec.vcf.gz"), sample=config["samples"]),
        expand(join(config["work_dir"], "variants/maf/{sample}.maf"), sample=config["samples"]),
