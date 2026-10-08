# Snakemake Pipeline: Automated Panel of Normals generation Workflow
# --------------------------------------------------------------
# This Snakemake pipeline automates the creation of a Panel of Normals (PoN)
# for somatic variant calling using GATK's Mutect2.
# --------------------------------------------------------------
# Author: Verena Passerini
# GitHub: https://github.com/VerenaPasserini/AVANTseq
# Last updated: October 2026 (v1.1.0)
# Snakemake version: 8.29.3
# --------------------------------------------------------------

# Configuration files (paths are relative to the directory snakemake is run from,
# i.e. the repository root). Override with --configfile if needed.
configfile: "config/config.yaml"
configfile: "config/samples_normal.yaml"

from os.path import join
import re

# Restrict the {sample} wildcard to the sample names listed in the sample file,
# so that e.g. "{sample}.bam" cannot match "tumor1.dedup.bam"
wildcard_constraints:
    sample="|".join(re.escape(s) for s in config["samples"])

# Include the rules
include: "rules/trim.smk"
include: "rules/align.smk"
include: "rules/qc.smk"
include: "rules/pon.smk"

# Pipeline output files
rule all:
    input:
        expand(join(config["work_dir"], "alignment/bams/{sample}.bam.bai"), sample=config["samples"]),
        expand(join(config["work_dir"], "alignment/qc/fastqc/{sample}_{read}_fastqc.html"), sample=config["samples"], read=["R1", "R2"]),
        join(config["work_dir"], "alignment/qc/multiqc_report.html"),
        join(config["work_dir"], "variants/qc/CoverageSummary.txt"),
        expand(join(config["work_dir"], "variants/qc/HsMetrics/{sample}_metrics.txt"), sample=config["samples"]),
        expand(join(config["work_dir"], "variants/mutect2/pon/{sample}.vcf.gz"), sample=config["samples"]),
        join(config["work_dir"], "variants/mutect2/pon/custom_pon.vcf.gz"),
        config["merged_pon"]
