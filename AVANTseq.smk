
# Snakemake Pipeline: Automated Variant Calling Analysis Workflow
# --------------------------------------------------------------
# This Snakemake pipeline automates the process of variant calling 
# and genome-wide copy number variation (CNV) analysis from NGS data 
# generated using hybrid capture-based targeted sequencing.
# --------------------------------------------------------------
# Author: Verena Passerini
# GitHub: https://github.com/VerenaPass/AVANTseq
# Last updated: March 2025
# Snakemake version: 8.29.3
# --------------------------------------------------------------

# Define the config file with samples list and relevant paths

configfile: "config.yaml"
configfile: "samples.yaml"

from os.path import join
from pathlib import Path
import subprocess

# Include rules

include: "rules/trim.smk"
include: "rules/align.smk"
include: "rules/qc.smk"
include: "rules/variants.smk"

# Pipeline output files
rule all:
    input:
        expand(join(config["work_dir"], "alignment/bams/{sample}.bam.bai"), sample=config["samples"]),
        expand(join(config["work_dir"], "alignment/qc/fastqc/{sample}_fastqc.html"), sample=config["samples"]),
        expand(join(config["work_dir"], "alignment/qc/multiqc_report.html"), sample=config["samples"]),
        expand(join(config["work_dir"], "variants/qc/CoverageSummary.txt"), sample=config["samples"]),
        expand(join(config["work_dir"], "variants/qc/HsMetrics/{sample}_metrics.txt"), sample=config["samples"]),
        expand(join(config["work_dir"], "variants/mutect2/{sample}.vcf.gz"), sample=config["samples"]),
        expand(join(config["work_dir"], "variants/mutect2/{sample}.getpileupsummaries.table"), sample=config["samples"]),
        expand(join(config["work_dir"], "variants/mutect2/{sample}.calculatecontamination.table"), sample=config["samples"]),
        expand(join(config["work_dir"], "variants/mutect2/{sample}.segments.table"), sample=config["samples"]),
        expand(join(config["work_dir"], "variants/filtered/{sample}_filtered.vcf.gz"), sample=config["samples"]),
        expand(join(config["work_dir"], "variants/filtered/{sample}_filtered_norm_dec.vcf.gz"), sample=config["samples"]),
        expand(join(config["work_dir"], "variants/maf/{sample}.maf"), sample=config["samples"]),

