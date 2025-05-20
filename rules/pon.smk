# ---------------------------------------------------------------------------------------
# Panel of Normals (PoN) Creation Workflow
# ---------------------------------------------------------------------------------------
# This set of rules creates a customized Panel of Normals (PoN) for somatic variant calling using Mutect2.
# The PoN helps filter out recurrent technical artifacts and germline variants in tumor-only samples.
#
# Steps included:
# - mutect2: Run Mutect2 on multiple normal samples in tumor only mode to generate raw variant calls.
# - genomicsdb_import: Import Mutect2 VCFs into a GenomicsDB workspace for joint processing.
# - create_pon: Generate a custom PoN VCF from the GenomicsDB workspace.
# - merge_pon_files: Merge the custom PoN with a public/reference PoN to produce the final PoN VCF used in variant calling.
#
# This workflow ensures improved specificity in somatic variant detection by leveraging both custom and public PoNs.
# ----------------------------------------------------------------------------------------

# Run Mutect2 on aligned reads
rule mutect2:
    input:
        bam=join(config["work_dir"], "alignment/bams/{sample}.bam"),
        bai=join(config["work_dir"], "alignment/bams/{sample}.bam.bai"),
        ref=config["ref_fa"],
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
        "--germline-resource {input.germline_resource} " 
        "--intervals {input.targets}" 
        "--interval-padding 50 -max-mnp-distance 0" 
        "-O {output} > {log} 2>&1"

# Run Genomics DB Import on mutect2 vcfs        
rule genomicsdb_import:
    input:
        vcfs=expand(join(config["work_dir"], "variants/mutect2/{sample}.vcf.gz"), sample=config["samples"]),
        ref=config["ref_fa"]
    output:
        join(config["work_dir"], "variants/mutect2/pon_db")
    params:
        interval_list=config["targets"],
        vcf_args=lambda wildcards, input: " ".join(f"-V {vcf}" for vcf in input.vcfs)
    message:
        "Running Genomics DB Import on mutect2 VCFs"
    shell:
        "gatk GenomicsDBImport "
        "-R {input.ref} "
        "-L {params.interval_list} "
        "{params.vcf_args} "
        "--genomicsdb-workspace-path {output} "

# Run CreateSomaticPanelOfNormals on mutect2 calls
rule create_pon:
    input:
        db=join(config["work_dir"], "variants/mutect2/pon_db"),
        ref=config["ref_fa"]
    output:
        join(config["work_dir"], "variants/mutect2/custom_pon.vcf.gz")
    log:
        join(config["work_dir"], "variants/log/createpon.log")
    message:
        "Running create panel of normal"
    shell:
        "gatk CreateSomaticPanelOfNormals " 
        "-R {input.ref}" 
        "-V gendb://{input.db}" 
        "-O {output} > {log} 2>&1"

# Merge custom pon with public pon to generate the final pon
rule merge_pon_files:
    input:
        pon=config["pon"],
        c_pon=join(config["work_dir"], "variants/mutect2/custom_pon.vcf.gz")
    output:
        config["merged_pon"]
    log:
        join(config["work_dir"], "variants/log/merge_pon.log")
    message:
        "Running merging pons"
    shell:
        "gatk MergeVcfs "
        "-I {input.pon} "
        "-I {input.c_pon}"
        "-O {output} > {log} 2>&1"
        
