# ---------------------------------------------------------------------------------------
# Quality Control Rules for Aligned BAM Files
# ---------------------------------------------------------------------------------------
# This module includes QC steps to evaluate the quality and coverage of aligned sequencing data.
# 
# Rules included:
# - fastqc: Runs FastQC on BAM files to assess sequencing quality.
# - qc_stats: Generates summary statistics and alignment counts using samtools stats and idxstats.
# - multiqc: Aggregates all QC reports (FastQC and samtools) into a single MultiQC report.
# - multibamsummary: Computes coverage summary across captured regions using multiBamSummary (deepTools).
# - collecthsmetrics: Uses GATK CollectHsMetrics to report hybrid selection metrics like coverage, on-target rate, and duplication.
#
# These steps ensure reliable input for downstream variant calling and allow detection of 
# technical issues such as poor capture efficiency or low coverage.
# ---------------------------------------------------------------------------------------

# Run fastqc on aligned reads
rule fastqc:
    input:
        join(config["work_dir"], "alignment/bams/{sample}.bam")
    output:
        join(config["work_dir"], "alignment/qc/fastqc/{sample}_fastqc.html")
    log:
        join(config["work_dir"], "alignment/log/{sample}_fastqc.log")
    params:
        out_dir=join(config["work_dir"],"alignment/qc")
    message:
        "Running FastQC on {input}"
    shell:
        "fastqc -t 16 -o {params.out_dir} {input} > {log} 2>&1"

# Run samtools stats and idxstats
rule qc_stats:
    input:
        bam=join(config["work_dir"], "alignment/bams/{sample}.bam"),
        idx=join(config["work_dir"], "alignment/bams/{sample}.bam.bai")
    output:
        stats=join(config["work_dir"], "alignment/qc/{sample}_stats.txt"),
        idxstats=join(config["work_dir"], "alignment/qc/{sample}_idxstats.txt")
    message:
        "Running samtools stats and idxstats on {input.bam}"
    shell:
        """
        samtools stats {input.bam} > {output.stats}
        samtools idxstats {input.bam} > {output.idxstats}
        """

# Run multiqc on fastqc and samtools stats
rule multiqc:
    input:
        expand(join(config["work_dir"], "alignment/qc/{sample}_fastqc.html"), sample=config["samples"]),
        expand(join(config["work_dir"], "alignment/qc/{sample}_idxstats.txt"), sample=config["samples"]),
        expand(join(config["work_dir"], "alignment/qc/{sample}_stats.txt"), sample=config["samples"])
    output:
        join(config["work_dir"], "alignment/qc/multiqc_report.html")
    params:
        dir=join(config["work_dir"], "alignment/qc")
    log:
        join(config["work_dir"], "alignment/log/multiqc.log")
    message:
        "Running MultiQC on {params.dir}"
    shell:
        "multiqc -f -o {params.dir} {params.dir} > {log} 2>&1"

# Run multibam summary on captured regions
rule multibamsummary:
    input:
        expand(join(config["work_dir"], "alignment/bams/{sample}.bam"), sample=config["samples"]),
        expand(join(config["work_dir"], "alignment/bams/{sample}.bam.bai"), sample=config["samples"]),
        bed=config["bed"]
    output:
        join(config["work_dir"],"variants/qc/CoverageSummary.txt")
    params:
        bam_list=join(config["work_dir"],"alignment/bams/*.bam")
    log:
        join(config["work_dir"],"variants/log/multibamsummary.log")
    message:
        "Running multiBamSummary on {params.bam_list}"
    shell:
        "multiBamSummary BED-file --BED {input.bed} " 
        "--bamfiles {params.bam_list} " 
        "--outRawCounts {output} " 
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