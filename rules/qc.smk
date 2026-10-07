# ---------------------------------------------------------------------------------------
# Quality Control Rules for raw reads and aligned BAM files
# ---------------------------------------------------------------------------------------
# This module includes QC steps to evaluate the quality and coverage of sequencing data.
#
# Rules included:
# - fastqc: Runs FastQC on the raw FASTQ files to assess sequencing quality.
# - qc_stats: Generates summary statistics and alignment counts using samtools stats and idxstats.
# - multiqc: Aggregates all QC reports (FastQC and samtools) into a single MultiQC report.
# - multibamsummary: Computes coverage summary across captured regions using multiBamSummary (deepTools).
# - collecthsmetrics: Uses GATK CollectHsMetrics to report hybrid selection metrics like coverage, on-target rate, and duplication.
#
# These steps ensure reliable input for downstream variant calling and allow detection of
# technical issues such as poor capture efficiency or low coverage.
# ---------------------------------------------------------------------------------------

import os

# Run FastQC on raw paired-end reads
rule fastqc:
    input:
        R1=join(config["work_dir"], "fastq/{sample}_R1.fastq.gz"),
        R2=join(config["work_dir"], "fastq/{sample}_R2.fastq.gz")
    output:
        html_R1=join(config["work_dir"], "alignment/qc/fastqc/{sample}_R1_fastqc.html"),
        html_R2=join(config["work_dir"], "alignment/qc/fastqc/{sample}_R2_fastqc.html"),
        zip_R1=join(config["work_dir"], "alignment/qc/fastqc/{sample}_R1_fastqc.zip"),
        zip_R2=join(config["work_dir"], "alignment/qc/fastqc/{sample}_R2_fastqc.zip")
    params:
        out_dir=lambda wildcards, output: os.path.dirname(output.html_R1)
    threads: 2
    log:
        join(config["work_dir"], "alignment/log/{sample}_fastqc.log")
    message:
        "Running FastQC on {input.R1} and {input.R2}"
    shell:
        "fastqc -t {threads} -o {params.out_dir} {input.R1} {input.R2} > {log} 2>&1"

# Run samtools stats and idxstats
rule qc_stats:
    input:
        bam=join(config["work_dir"], "alignment/bams/{sample}.bam"),
        idx=join(config["work_dir"], "alignment/bams/{sample}.bam.bai")
    output:
        stats=join(config["work_dir"], "alignment/qc/{sample}_stats.txt"),
        idxstats=join(config["work_dir"], "alignment/qc/{sample}_idxstats.txt")
    log:
        join(config["work_dir"], "alignment/log/{sample}_qc_stats.log")
    message:
        "Running samtools stats and idxstats on {input.bam}"
    shell:
        """
        samtools stats {input.bam} > {output.stats} 2> {log}
        samtools idxstats {input.bam} > {output.idxstats} 2>> {log}
        """

# Run MultiQC on FastQC and samtools stats
rule multiqc:
    input:
        expand(join(config["work_dir"], "alignment/qc/fastqc/{sample}_{read}_fastqc.zip"), sample=config["samples"], read=["R1", "R2"]),
        expand(join(config["work_dir"], "alignment/qc/{sample}_idxstats.txt"), sample=config["samples"]),
        expand(join(config["work_dir"], "alignment/qc/{sample}_stats.txt"), sample=config["samples"])
    output:
        join(config["work_dir"], "alignment/qc/multiqc_report.html")
    params:
        dir=lambda wildcards, output: os.path.dirname(output[0])
    log:
        join(config["work_dir"], "alignment/log/multiqc.log")
    message:
        "Running MultiQC on {params.dir}"
    shell:
        "multiqc -f -o {params.dir} {params.dir} > {log} 2>&1"

# Run multiBamSummary on captured regions
rule multibamsummary:
    input:
        bams=expand(join(config["work_dir"], "alignment/bams/{sample}.bam"), sample=config["samples"]),
        bais=expand(join(config["work_dir"], "alignment/bams/{sample}.bam.bai"), sample=config["samples"]),
        bed=config["bed"]
    output:
        counts=join(config["work_dir"], "variants/qc/CoverageSummary.txt"),
        npz=join(config["work_dir"], "variants/qc/CoverageSummary.npz")
    threads: 4
    log:
        join(config["work_dir"], "variants/log/multibamsummary.log")
    message:
        "Running multiBamSummary on {input.bams}"
    shell:
        "multiBamSummary BED-file --BED {input.bed} "
        "--bamfiles {input.bams} "
        "--numberOfProcessors {threads} "
        "--outFileName {output.npz} "
        "--outRawCounts {output.counts} "
        "> {log} 2>&1"

# Run CollectHsMetrics on captured regions
rule collecthsmetrics:
    input:
        bam=join(config["work_dir"], "alignment/bams/{sample}.bam"),
        baits=config["baits"],
        target=config["targets"]
    output:
        join(config["work_dir"], "variants/qc/HsMetrics/{sample}_metrics.txt")
    log:
        join(config["work_dir"], "variants/log/{sample}_metrics.log")
    message:
        "Running CollectHsMetrics on {input.bam}"
    shell:
        "gatk CollectHsMetrics -I {input.bam} -O {output} "
        "-BAIT_INTERVALS {input.baits} -TARGET_INTERVALS {input.target} "
        "-VALIDATION_STRINGENCY SILENT "
        "> {log} 2>&1"
