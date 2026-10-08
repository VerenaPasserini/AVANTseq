# ---------------------------------------------------------------------------------------
# Snakemake rule for trimming paired-end FASTQ files using fastp.
# Automatically detects and removes adapter sequences (paired-end overlap analysis),
# trims low-quality 3' ends, and discards reads that become too short.
# The JSON report is collected by MultiQC.
# ---------------------------------------------------------------------------------------

# Trim raw FASTQ files
rule trim:
    input:
        R1=join(config["work_dir"], "fastq/{sample}_R1.fastq.gz"),
        R2=join(config["work_dir"], "fastq/{sample}_R2.fastq.gz")
    output:
        R1_trimmed=temp(join(config["work_dir"], "alignment/trimmed/{sample}_R1-trimmed.fq.gz")),
        R2_trimmed=temp(join(config["work_dir"], "alignment/trimmed/{sample}_R2-trimmed.fq.gz")),
        html=join(config["work_dir"], "alignment/qc/fastp/{sample}_fastp.html"),
        json=join(config["work_dir"], "alignment/qc/fastp/{sample}_fastp.json")
    params:
        cut_tail_quality=config.get("fastp_cut_tail_quality", 20),
        min_length=config.get("fastp_min_length", 25),
        extra=config.get("fastp_extra", "")
    threads: 8
    log:
        join(config["work_dir"], "alignment/log/{sample}_trim.log")
    message:
        "Trimming {input.R1} and {input.R2} with fastp"
    shell:
        "fastp --in1 {input.R1} --in2 {input.R2} "
        "--out1 {output.R1_trimmed} --out2 {output.R2_trimmed} "
        "--detect_adapter_for_pe "
        "--cut_tail --cut_tail_mean_quality {params.cut_tail_quality} "
        "--length_required {params.min_length} "
        "--thread {threads} {params.extra} "
        "--html {output.html} --json {output.json} "
        "> {log} 2>&1"
