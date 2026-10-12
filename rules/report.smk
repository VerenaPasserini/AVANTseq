# ---------------------------------------------------------------------------------------
# QC report: collect CollectHsMetrics of all samples, classify them as PASS / FLAG / FAIL
# on the % of target bases covered > 30x, and plot coverage (scripts/qc_report.R).
# ---------------------------------------------------------------------------------------

rule qc_report:
    input:
        metrics=expand(join(config["work_dir"], "variants/qc/HsMetrics/{sample}_metrics.txt"), sample=config["samples"])
    output:
        table=join(config["work_dir"], "report/hsmetrics_qc.tsv"),
        plot=join(config["work_dir"], "report/hsmetrics_qc.pdf")
    params:
        pass_30x=config.get("qc_pass_30x", 0.8),  # PASS if > 80% of target bases are covered > 30x
        flag_30x=config.get("qc_flag_30x", 0.6)   # FLAG if > 60% (and <= 80%), FAIL otherwise
    log:
        join(config["work_dir"], "report/qc_report.log")
    message:
        "Collecting HsMetrics and plotting the coverage QC"
    script:
        "../scripts/qc_report.R"


# ---------------------------------------------------------------------------------------
# Report: clean and merge the per-sample MAFs (covered regions, depth, alt reads, VAF) and
# plot a maftools summary and an oncoplot (scripts/maf_report.R).
# ---------------------------------------------------------------------------------------

rule maf_report:
    input:
        mafs=expand(join(config["work_dir"], "variants/maf/{sample}.maf"), sample=config["samples"]),
        bed=config["bed"],  # covered regions: only variants inside them are reported
        qc=join(config["work_dir"], "report/hsmetrics_qc.tsv")  # samples with excluded QC status are skipped
    output:
        maf=join(config["work_dir"], "report/all_samples_filtered.maf"),
        summary=join(config["work_dir"], "report/maf_summary.pdf"),
        oncoplot=join(config["work_dir"], "report/oncoplot.pdf")
    params:
        min_depth=config.get("report_min_depth", 20),  # minimum total depth (t_ref_count + t_alt_count)
        min_alt=config.get("report_min_alt", 4),       # minimum number of reads supporting the variant
        min_vaf=config.get("report_min_vaf", 0.10),    # minimum variant allele fraction
        max_pop_af=config.get("report_max_pop_af", 0.001),  # maximum gnomAD population allele frequency
        exclude_qc=config.get("report_exclude_qc", ["FAIL"])  # QC classes left out of the report
    log:
        join(config["work_dir"], "report/maf_report.log")
    message:
        "Merging MAFs and plotting the maftools summary and oncoplot"
    script:
        "../scripts/maf_report.R"
