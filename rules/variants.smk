# ---------------------------------------------------------------------------------------
# Somatic Variant Calling and Annotation (Mutect2)
# ---------------------------------------------------------------------------------------
# This module contains rules for calling, filtering, normalizing, and annotating somatic variants
# using GATK's Mutect2 workflow, vt for normalization, and Funcotator for annotation.
#
# Rules included:
# - mutect2: Runs GATK Mutect2 on each sample to call potential somatic variants.
# - getpileupsummaries: Gathers pileup summary statistics used for contamination estimation.
# - calculatecontamination: Estimates contamination based on pileup summaries.
# - filtermutectcalls: Applies Mutect2-specific filters using contamination estimates.
# - vt_normalize_decompose: Normalizes and decomposes filtered VCFs for annotation.
# - funcotator: Annotates the final filtered variants in MAF format using Funcotator.
#
# This pipeline ensures robust somatic variant calling and annotation for tumor-only samples.
# ----------------------------------------------------------------------------------------

# Run Mutect2 on aligned reads
rule mutect2:
    input:
        bam=join(config["work_dir"], "alignment/bams/{sample}.bam"),
        bai=join(config["work_dir"], "alignment/bams/{sample}.bam.bai"),
        ref=config["ref_fa"],
        pon=config["pon"],
        germline_resource=config["germline_resource"],
        targets=config["targets"]
    output:
        join(config["work_dir"], "variants/mutect2/{sample}.vcf.gz")
    log:
        join(config["work_dir"], "variants/log/{sample}_mutect2.log")
    message:
        "Running Mutect2 on {input.bam}"
    shell:
        "gatk Mutect2 -R {input.ref} -I {input.bam} " 
        "-pon {input.pon} --germline-resource {input.germline_resource} " 
        "--intervals {input.targets} --genotype-germline-sites true " 
        "--genotype-pon-sites true --interval-padding 50 " 
        "-O {output} > {log} 2>&1"
    
# Run GetPileupSummaries on aligned reads
rule getpileupsummaries:
    input:
        bam=join(config["work_dir"], "alignment/bams/{sample}.bam"),
        bai=join(config["work_dir"], "alignment/bams/{sample}.bam.bai"),
        vcf=config["vcf_exac"]
    output:
        join(config["work_dir"], "variants/mutect2/{sample}.getpileupsummaries.table")
    log:
        join(config["work_dir"], "variants/log/{sample}_getpileupsummaries.log")
    message:
        "Running GetPileupSummaries on {input.bam}"
    shell:
        "gatk GetPileupSummaries -I {input.bam} " 
        "-V {input.vcf} -L {input.vcf} " 
        "-O {output} > {log} 2>&1"
    
# Run CalculateContamination on GetPileupSummaries
rule calculatecontamination:
    input:
        getpileupsummaries=join(config["work_dir"], "variants/mutect2/{sample}.getpileupsummaries.table")
    output:
        contamination_table=join(config["work_dir"], "variants/mutect2/{sample}.calculatecontamination.table"),
        segments=join(config["work_dir"], "variants/mutect2/{sample}.segments.table")
    log:
        join(config["work_dir"], "variants/log/{sample}_calculatecontamination.log")
    message:
        "Running CalculateContamination on {input.getpileupsummaries}"
    shell:
        "gatk CalculateContamination -I {input.getpileupsummaries} " 
        "-tumor-segmentation {output.segments} " 
        "-O {output.contamination_table} > {log} 2>&1"
    
# Run FilterMutectCalls on Mutect2 output
rule filtermutectcalls:
    input:
        ref=config["ref_fa"],
        vcf=join(config["work_dir"], "variants/mutect2/{sample}.vcf.gz"),
        contamination_table=join(config["work_dir"], "variants/mutect2/{sample}.calculatecontamination.table"),
        segments=join(config["work_dir"], "variants/mutect2/{sample}.segments.table")
    output:
        join(config["work_dir"], "variants/filtered/{sample}_filtered.vcf.gz")
    log:
        join(config["work_dir"], "variants/log/{sample}_filtermutectcalls.log")
    message:
        "Running FilterMutectCalls on {input.vcf}"
    shell:
        "gatk FilterMutectCalls -R {input.ref} -V {input.vcf} " 
        "--contamination-table {input.contamination_table} " 
        "--tumor-segmentation {input.segments} " 
        "-O {output} > {log} 2>&1"
    

# Normalize VCF files
rule vt_normalize_decompose:
    input:
        vcf=join(config["work_dir"], "variants/filtered/{sample}_filtered.vcf.gz"),
        ref=config["ref_fa"]
    output:
        vcf=join(config["work_dir"], "variants/filtered/{sample}_filtered_norm_dec.vcf.gz")
    params:
        temp=join(config["work_dir"], "variants/filtered/{sample}_temp.normalized.vcf")
    log:
        join(config["work_dir"], "variants/log/{sample}_vt_normalize_decompose.log")
    message:
        "Normalizing and decomposing {input.vcf}"
    shell:
        """
        vt normalize {input.vcf} -r {input.ref} -o {params.temp} 
        vt decompose {params.temp} -o {output.vcf}
        tabix -p vcf {output.vcf}
        rm {params.temp}
        """

# Run Funcotator on filtered variants
rule funcotator:
    input:
        vcf=join(config["work_dir"], "variants/filtered/{sample}_filtered_norm_dec.vcf.gz"),
        ref=config["ref_fa"],
        data=config["data_source"]
    output:
        join(config["work_dir"], "variants/maf/{sample}.maf")
    log:
        join(config["work_dir"], "variants/log/{sample}_funcotator.log")
    message:
        "Running Funcotator on {input.vcf}"
    shell:
        "gatk Funcotator --variant {input.vcf} --reference {input.ref} " 
        "--ref-version hg38 --data-sources-path {input.data} " 
        "--output {output} --output-file-format MAF " 
        "--remove-filtered-variants true " 
        "--transcript-selection-mode BEST_EFFECT > {log} 2>&1"
    
