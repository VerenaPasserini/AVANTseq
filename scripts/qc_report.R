# -----------------------------------------------------------------------------
# AVANTseq QC report: collect GATK/Picard CollectHsMetrics results of all samples,
# classify each sample as PASS / FLAG / FAIL and plot target coverage.
# Run by Snakemake (rule qc_report in rules/report.smk).
# -----------------------------------------------------------------------------

log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggrepel) })

pass_30x <- as.numeric(snakemake@params[["pass_30x"]])   # PASS if PCT_TARGET_BASES_30X >  pass_30x
flag_30x <- as.numeric(snakemake@params[["flag_30x"]])   # FLAG if PCT_TARGET_BASES_30X >  flag_30x, else FAIL

# 1. Read the metrics table of each sample ------------------------------------
read_hs <- function(f) {
  lines <- readLines(f)
  i <- grep("^## METRICS CLASS", lines)[1]              # header + values follow this line
  x <- fread(text = lines[(i + 1):(i + 2)])
  x[, SAMPLE := sub("_metrics\\.txt$", "", basename(f))]
  x
}
hs <- rbindlist(lapply(snakemake@input[["metrics"]], read_hs), fill = TRUE)

hs[, QC := ifelse(PCT_TARGET_BASES_30X > pass_30x, "PASS",
           ifelse(PCT_TARGET_BASES_30X > flag_30x, "FLAG", "FAIL"))]
setcolorder(hs, c("SAMPLE", "QC"))
fwrite(hs, snakemake@output[["table"]], sep = "\t")
message("QC per sample (PCT_TARGET_BASES_30X):")
for (k in seq_len(nrow(hs)))
  message(sprintf("  %-20s %5.1f%%  %s", hs$SAMPLE[k], 100 * hs$PCT_TARGET_BASES_30X[k], hs$QC[k]))

# 2. Plot ------------------------------------------------------------------------
qc_colors <- c(FAIL = "red", FLAG = "orange", PASS = "green3")
p <- ggplot(hs, aes(x = PCT_TARGET_BASES_30X * 100, y = ZERO_CVG_TARGETS_PCT * 100, color = QC)) +
  geom_point(alpha = 0.8, size = 2) +
  scale_color_manual(values = qc_colors) +
  scale_x_continuous(breaks = seq(0, 100, 20), limits = c(0, 100)) +
  geom_vline(xintercept = flag_30x * 100, linetype = "dashed", color = "darkred", linewidth = 0.2) +
  geom_vline(xintercept = pass_30x * 100, linetype = "dotted", color = "darkgreen", linewidth = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "darkred", linewidth = 0.2) +
  geom_label_repel(aes(label = SAMPLE), box.padding = 0.35, point.padding = 0.5,
                   segment.color = "grey50", show.legend = FALSE) +
  labs(x = "% target bases with > 30x coverage", y = "% targets not covered", color = "QC") +
  theme_light()
ggsave(snakemake@output[["plot"]], p, width = 8, height = 6)
