# ---------------------------------------------------------------------------------------
# Snakemake rule for trimming paired-end FASTQ files using Atropos.
# Removes low-quality bases and adapters, producing cleaned reads for alignment.
# ---------------------------------------------------------------------------------------

# Trim fastqc files
rule trim:
    input:
        R1=join(config["work_dir"], "fastq/{sample}_R1.fastq.gz"),
        R2=join(config["work_dir"], "fastq/{sample}_R2.fastq.gz")
    output:
        R1_trimmed=temp(join(config["work_dir"], "alignment/trimmed/{sample}_R1-trimmed.fq.gz")),
        R2_trimmed=temp(join(config["work_dir"], "alignment/trimmed/{sample}_R2-trimmed.fq.gz"))
    log:
        join(config["work_dir"],"alignment/log/{sample}_trim.log")
    message:
        "Trimming {input.R1} and {input.R2}"
    shell:
        "atropos trim --threads 16 --quality-base 33 --format fastq --overlap 8 " 
        "--no-default-adapters --no-cache-adapters -pe1 {input.R1} -pe2 {input.R2} "
        "-o {output.R1_trimmed} -p {output.R2_trimmed} --quality-cutoff=5 --minimum-length=25 "
        "> {log} 2>&1"
