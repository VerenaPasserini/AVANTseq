# ----------------------------------------------------------------------------------------
# This module contains Snakemake rules for aligning paired-end trimmed reads to a reference
# genome using BWA-MEM, followed by duplicate marking with GATK and indexing with samtools.
# These steps generate analysis-ready BAM files for downstream processing.
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
        "(bwa mem -t {threads} -R '@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tPL:ILLUMINA' "
        "{params.ref} {input.R1_trimmed} {input.R2_trimmed} "
        "| samtools sort -o {output}) "
        "> {log} 2>&1"


# Mark duplicates in aligned reads
rule mark_duplicates:
    input:
        bam=join(config["work_dir"], "alignment/bams/{sample}.temp.bam")  # Input BAM file from align rule
    output:
        bam_dedup=join(config["work_dir"], "alignment/bams/{sample}.bam"),  # Output BAM file after marking duplicates
        metrics=join(config["work_dir"], "alignment/log/{sample}_dedup_metrics.txt")  # Metrics file for duplicate marking
    log:
        join(config["work_dir"], "alignment/log/{sample}_markdup.log")  # Log file for duplicate marking process
    message:
        "Marking duplicates in {input.bam}"
    shell:
        "gatk MarkDuplicates -I {input.bam} -O {output.bam_dedup} "
        "-M {output.metrics} --REMOVE_DUPLICATES false "
        "> {log} 2>&1"

# Index bam files
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
