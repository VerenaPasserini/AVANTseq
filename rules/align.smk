# ----------------------------------------------------------------------------------------
# This module contains Snakemake rules for aligning paired-end trimmed reads to a reference
# genome using BWA-MEM, followed by duplicate marking with GATK MarkDuplicates, base quality
# score recalibration (GATK BQSR) and indexing with samtools.
# These steps generate analysis-ready BAM files for downstream processing,
# following the GATK data pre-processing best practices.
# ----------------------------------------------------------------------------------------

# Align trimmed reads
rule align:
    input:
        R1_trimmed=join(config["work_dir"], "alignment/trimmed/{sample}_R1-trimmed.fq.gz"),
        R2_trimmed=join(config["work_dir"], "alignment/trimmed/{sample}_R2-trimmed.fq.gz")
    output:
        temp(join(config["work_dir"], "alignment/bams/{sample}.temp.bam"))
    params:
        ref=config["ref_bwa"]
    threads: 16
    log:
        join(config["work_dir"], "alignment/log/{sample}_align.log")
    message:
        "Aligning {input.R1_trimmed} and {input.R2_trimmed}"
    shell:
        # Both bwa and samtools messages are written to the log file
        "(bwa mem -t {threads} -R '@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tLB:{wildcards.sample}\\tPL:ILLUMINA' "
        "{params.ref} {input.R1_trimmed} {input.R2_trimmed} "
        "| samtools sort -o {output}) "
        "> {log} 2>&1"


# Mark duplicates in aligned reads
rule mark_duplicates:
    input:
        bam=join(config["work_dir"], "alignment/bams/{sample}.temp.bam")  # Input BAM file from align rule
    output:
        bam_dedup=temp(join(config["work_dir"], "alignment/bams/{sample}.dedup.bam")),  # BAM file after marking duplicates
        bai_dedup=temp(join(config["work_dir"], "alignment/bams/{sample}.dedup.bai")),  # Index (needed by BaseRecalibrator)
        metrics=join(config["work_dir"], "alignment/qc/markdup/{sample}_dedup_metrics.txt")  # Duplicate metrics (collected by MultiQC)
    log:
        join(config["work_dir"], "alignment/log/{sample}_markdup.log")  # Log file for duplicate marking process
    message:
        "Marking duplicates in {input.bam}"
    shell:
        "gatk MarkDuplicates -I {input.bam} -O {output.bam_dedup} "
        "-M {output.metrics} --REMOVE_DUPLICATES false --CREATE_INDEX true "
        "> {log} 2>&1"


# Model systematic base quality errors (BQSR step 1)
rule base_recalibrator:
    input:
        bam=join(config["work_dir"], "alignment/bams/{sample}.dedup.bam"),
        bai=join(config["work_dir"], "alignment/bams/{sample}.dedup.bai"),
        ref=config["ref_fa"],
        known_sites=config["known_sites"],
        targets=config["targets"]
    output:
        join(config["work_dir"], "alignment/qc/bqsr/{sample}_recal.table")
    params:
        known=lambda wildcards, input: " ".join(f"--known-sites {k}" for k in input.known_sites)
    log:
        join(config["work_dir"], "alignment/log/{sample}_baserecalibrator.log")
    message:
        "Running BaseRecalibrator on {input.bam}"
    shell:
        "gatk BaseRecalibrator -R {input.ref} -I {input.bam} "
        "{params.known} "
        "-L {input.targets} --interval-padding 50 "
        "-O {output} > {log} 2>&1"


# Apply the recalibration model to produce the analysis-ready BAM (BQSR step 2)
rule apply_bqsr:
    input:
        bam=join(config["work_dir"], "alignment/bams/{sample}.dedup.bam"),
        bai=join(config["work_dir"], "alignment/bams/{sample}.dedup.bai"),
        ref=config["ref_fa"],
        table=join(config["work_dir"], "alignment/qc/bqsr/{sample}_recal.table")
    output:
        join(config["work_dir"], "alignment/bams/{sample}.bam")
    log:
        join(config["work_dir"], "alignment/log/{sample}_applybqsr.log")
    message:
        "Applying BQSR to {input.bam}"
    shell:
        "gatk ApplyBQSR -R {input.ref} -I {input.bam} "
        "--bqsr-recal-file {input.table} "
        "--create-output-bam-index false "
        "-O {output} > {log} 2>&1"


# Index analysis-ready BAM files
rule index:
    input:
        join(config["work_dir"], "alignment/bams/{sample}.bam")
    output:
        join(config["work_dir"], "alignment/bams/{sample}.bam.bai")
    log:
        join(config["work_dir"], "alignment/log/{sample}_index.log")
    message:
        "Indexing {input}"
    shell:
        "samtools index {input} > {log} 2>&1"
