# ---------------------------------------------------------------------------------------
# Report: clean and merge the per-sample MAFs (covered regions, depth, alt reads, VAF) and
# plot a maftools summary and an oncoplot (scripts/maf_report.R).
# ---------------------------------------------------------------------------------------

rule maf_report:
    input:
        mafs=expand(join(config["work_dir"], "variants/maf/{sample}.maf"), sample=config["samples"]),
        bed=config["bed"]  # covered regions: only variants inside them are reported
    output:
        maf=join(config["work_dir"], "report/all_samples_filtered.maf"),
        summary=join(config["work_dir"], "report/maf_summary.pdf"),
        oncoplot=join(config["work_dir"], "report/oncoplot.pdf")
    params:
        min_depth=config.get("report_min_depth", 20),  # minimum total depth (t_ref_count + t_alt_count)
        min_alt=config.get("report_min_alt", 4),       # minimum number of reads supporting the variant
        min_vaf=config.get("report_min_vaf", 0.10),    # minimum variant allele fraction
        max_pop_af=config.get("report_max_pop_af", 0.001)  # maximum gnomAD population allele frequency
    log:
        join(config["work_dir"], "report/maf_report.log")
    message:
        "Merging MAFs and plotting the maftools summary and oncoplot"
    script:
        "../scripts/maf_report.R"
