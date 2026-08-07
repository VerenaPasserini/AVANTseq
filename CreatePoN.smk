# Snakemake Pipeline: Automated Panel of Normal generation Workflow
# --------------------------------------------------------------
# This Snakemake pipeline automates the creation of a Panel of Normals (PoN) 
# for somatic variant calling using GATK’s Mutect2.
# --------------------------------------------------------------
# Author: Verena Passerini
# GitHub: https://github.com/VerenaPasserini/AVANTseq
# Last updated: May 2025
# Snakemake version: 8.29.3
# --------------------------------------------------------------

# Define the config file with samples list and relevant paths

configfile: "config_pon.yaml"
configfile: "samples_pon.yaml"

from os.path import join
from pathlib import Path
import subprocess

# Include the rules
include: "rules/trim.smk"
include: "rules/align.smk"
include: "rules/qc.smk"
include: "rules/pon.smk"

# Pipeline output files
rule all:
    input:
        expand(join(config["work_dir"], "alignment/bams/{sample}.bam.bai"), sample=config["samples"]),
        expand(join(config["work_dir"], "alignment/qc/fastqc/{sample}_fastqc.html"), sample=config["samples"]),
        expand(join(config["work_dir"], "alignment/qc/multiqc_report.html"), sample=config["samples"]),
        expand(join(config["work_dir"], "variants/qc/CoverageSummary.txt"), sample=config["samples"]),
        expand(join(config["work_dir"], "variants/qc/HsMetrics/{sample}_metrics.txt"), sample=config["samples"]),
        expand(join(config["work_dir"], "variants/mutect2/{sample}.vcf.gz"), sample=config["samples"]),
        join(config["work_dir"], "variants/mutect2/custom_pon.vcf.gz"),
        config["merged_pon"]

