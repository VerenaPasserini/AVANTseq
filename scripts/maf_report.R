# -----------------------------------------------------------------------------
# AVANTseq report: clean and merge the per-sample Funcotator MAFs, then plot a
# maftools summary and an oncoplot.
# Run by Snakemake (rule maf_report in rules/report.smk).
# -----------------------------------------------------------------------------

log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

suppressPackageStartupMessages({ library(maftools); library(data.table) })

bed_file  <- snakemake@input[["bed"]]
bed_pad   <- as.numeric(snakemake@params[["bed_padding"]])  # bp added on each side of the BED regions
min_depth <- as.numeric(snakemake@params[["min_depth"]])
min_alt   <- as.numeric(snakemake@params[["min_alt"]])
min_vaf   <- as.numeric(snakemake@params[["min_vaf"]])
max_popaf <- as.numeric(snakemake@params[["max_pop_af"]])   # max gnomAD allele frequency
top_genes <- as.numeric(snakemake@params[["top_genes"]])    # number of genes in the oncoplot
gene_mar_cfg <- as.numeric(snakemake@params[["gene_mar"]])  # 0 = automatic
drop_unnamed <- isTRUE(as.logical(snakemake@params[["drop_unnamed"]]))  # drop features without a gene symbol (ENSG IDs)
exclude_qc <- unlist(snakemake@params[["exclude_qc"]])      # QC classes to leave out (e.g. FAIL)

# Coverage QC from qc_report: leave out samples with an excluded QC class
qc <- fread(snakemake@input[["qc"]], select = c("SAMPLE", "QC", "PCT_TARGET_BASES_30X"),
            colClasses = list(character = "SAMPLE"))
excluded <- qc$SAMPLE[qc$QC %in% exclude_qc]
if (length(excluded) > 0)
  message("Samples excluded because of coverage QC (", paste(exclude_qc, collapse = "/"), "): ",
          paste(excluded, collapse = ", "), "\n")
maf_files <- snakemake@input[["mafs"]]
maf_files <- maf_files[!sub("\\.maf$", "", basename(maf_files)) %in% excluded]

# Covered regions (BED, 0-based start) -> 1-based intervals, merged per chromosome
bed_lines <- sub("\r$", "", readLines(bed_file, warn = FALSE))
bed_lines <- bed_lines[nzchar(bed_lines) & !grepl("^(#|track|browser)", bed_lines)]
f <- strsplit(bed_lines, "[[:space:]]+")
bed <- data.frame(chr   = vapply(f, `[`, "", 1),
                  start = suppressWarnings(as.numeric(vapply(f, `[`, "", 2))) + 1 - bed_pad,
                  end   = suppressWarnings(as.numeric(vapply(f, `[`, "", 3))) + bed_pad)
bed <- bed[!is.na(bed$start) & !is.na(bed$end), ]          # drops a possible header line
bed <- bed[order(bed$chr, bed$start), ]
bed_merged <- lapply(split(bed, bed$chr), function(b) {
  grp <- cumsum(c(TRUE, b$start[-1] > cummax(b$end)[-nrow(b)] + 1))
  data.frame(start = as.numeric(tapply(b$start, grp, min)), end = as.numeric(tapply(b$end, grp, max)))
})
message(sprintf("BED: %d regions on %d chromosomes, padded by %g bp on each side", nrow(bed), length(bed_merged), bed_pad))

# TRUE for variants overlapping a covered region
in_bed <- function(chr, s, e) {
  hit <- logical(length(chr))
  for (ch in intersect(unique(chr), names(bed_merged))) {
    m <- bed_merged[[ch]]
    i <- which(chr == ch)
    k <- findInterval(e[i], m$start)                  # last region starting before the variant end
    hit[i] <- k > 0 & m$end[pmax(k, 1)] >= s[i]
  }
  hit
}

# Funcotator classifications that maftools does not know -> MAF / maftools names
vc_map <- c(START_CODON_SNP         = "Translation_Start_Site",
            START_CODON_INS         = "Translation_Start_Site",
            START_CODON_DEL         = "Translation_Start_Site",
            DE_NOVO_START_IN_FRAME  = "5'UTR",
            DE_NOVO_START_OUT_FRAME = "5'UTR",
            COULD_NOT_DETERMINE     = "Unknown")
vc_known <- c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site",
              "Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del", "In_Frame_Ins",
              "Missense_Mutation", "Silent", "3'UTR", "5'UTR", "3'Flank", "5'Flank",
              "IGR", "Intron", "RNA", "lincRNA", "Splice_Region", "Targeted_Region", "Unknown")

# 1. Read, clean and label each sample's MAF -----------------------------------
counts <- data.frame(SAMPLE = character(), variants_all = integer(), variants_in_bed = integer(),
                     variants_depth_vaf = integer(), variants_reported = integer())  # per-sample counts
read_sample <- function(f) {
  sample <- sub("\\.maf$", "", basename(f))
  # skip the "#" header lines at the top only ("#" may also occur inside annotation fields);
  # Funcotator MAFs can have many thousands of header lines
  lines <- readLines(f, warn = FALSE)
  n_header <- match(FALSE, startsWith(lines, "#")) - 1
  df <- read.delim(text = lines[-seq_len(n_header)], comment.char = "", quote = "",
                   stringsAsFactors = FALSE)
  required <- c("Chromosome", "Start_Position", "End_Position", "t_alt_count", "t_ref_count")
  if (!all(required %in% names(df)))
    stop("unexpected MAF format in ", f, ": missing column(s) ", paste(setdiff(required, names(df)), collapse = ", "))
  n_all <- nrow(df)

  # keep variants overlapping the covered regions of the BED file
  keep <- in_bed(as.character(df$Chromosome), df$Start_Position, pmax(df$Start_Position, df$End_Position))
  dropped <- sort(table(df$Variant_Classification[!keep]), decreasing = TRUE)
  df <- df[keep, , drop = FALSE]
  # features without an official gene symbol (Funcotator reports the Ensembl ID, e.g. ENSG00000240040)
  if (drop_unnamed) df <- df[!grepl("^ENSG[0-9]+", df$Hugo_Symbol), , drop = FALSE]
  n_bed <- nrow(df)

  df$t_depth <- df$t_alt_count + df$t_ref_count
  df$VAF     <- df$t_alt_count / df$t_depth
  df <- df[df$t_depth >= min_depth & df$t_alt_count >= min_alt & df$VAF >= min_vaf, , drop = FALSE]
  n_reads <- nrow(df)

  # population frequency: POPAF (from Mutect2) is -log10 of the gnomAD allele frequency;
  # variants absent from gnomAD get a high POPAF and are kept
  popaf <- suppressWarnings(as.numeric(sub(",.*", "", df$POPAF)))
  df <- df[is.na(popaf) | popaf >= -log10(max_popaf), , drop = FALSE]

  message(sprintf("%-20s all: %5d   in BED (named genes): %5d   depth/alt/VAF: %5d   gnomAD AF <= %g: %5d",
                  sample, n_all, n_bed, n_reads, max_popaf, nrow(df)))
  if (length(dropped) > 0)
    message(sprintf("%-20s   outside BED: %s", "", paste(sprintf("%s %d", names(dropped), dropped), collapse = ", ")))
  counts <<- rbind(counts, data.frame(SAMPLE = sample, variants_all = n_all, variants_in_bed = n_bed,
                                      variants_depth_vaf = n_reads, variants_reported = nrow(df)))
  if (nrow(df) == 0) return(NULL)
  df$Tumor_Sample_Barcode <- sample
  df
}
message("Variants per sample (PASS variants from Funcotator):")
mut_df <- do.call(rbind, lapply(maf_files, read_sample))

# 1b. Cohort summary: table + CONSORT-style flow diagram -------------------------
all_samples <- sub("\\.maf$", "", basename(snakemake@input[["mafs"]]))
analysed    <- sub("\\.maf$", "", basename(maf_files))
no_variants <- counts$SAMPLE[counts$variants_reported == 0]
flagged     <- intersect(analysed, qc$SAMPLE[qc$QC == "FLAG"])
n_analysed  <- length(analysed)

cohort <- merge(data.frame(SAMPLE = all_samples), qc, by = "SAMPLE", all.x = TRUE)
cohort <- merge(cohort, counts, by = "SAMPLE", all.x = TRUE)
cohort$status <- ifelse(cohort$SAMPLE %in% excluded, "excluded (coverage QC)",
                 ifelse(cohort$SAMPLE %in% no_variants, "analysed, no somatic variant passing filters",
                        "analysed"))
write.table(cohort, snakemake@output[["cohort"]], sep = "\t", quote = FALSE, row.names = FALSE)
if (length(no_variants) > 0)
  message("WARNING: no somatic variant passed the filters in sample(s): ", paste(no_variants, collapse = ", "),
          " (they are counted in the cohort but do not appear in the plots)")

ids <- function(x) if (length(x) == 0) "" else paste(strwrap(paste(x, collapse = ", "), width = 45), collapse = "\n")
draw_flow <- function(file) {
  library(grid)
  pdf(file, width = 8.5, height = 7)
  grid.newpage()
  box <- function(x, y, w, h, title, body = "", fill = "#F2F5F8", border = "#4C78A8") {
    x <- unit(x, "npc"); y <- unit(y, "npc"); w <- unit(w, "npc"); h <- unit(h, "npc")
    grid.roundrect(x, y, w, h, r = unit(2, "mm"), gp = gpar(fill = fill, col = border, lwd = 1.5))
    grid.text(title, x, y + unit(if (nzchar(body)) 0.25 else 0, "in"), gp = gpar(fontsize = 12, fontface = "bold"))
    if (nzchar(body)) grid.text(body, x, y - unit(0.12, "in"), gp = gpar(fontsize = 9.5, col = "#333333"), vjust = 1)
  }
  arr <- function(x0, y0, x1, y1) grid.lines(c(x0, x1), c(y0, y1),
            arrow = arrow(length = unit(2.5, "mm"), type = "closed"), gp = gpar(col = "#555555", fill = "#555555"))
  grid.text("Cohort overview", 0.5, 0.95, gp = gpar(fontsize = 15, fontface = "bold"))
  # main column
  box(0.32, 0.80, 0.42, 0.11, sprintf("Samples sequenced: n = %d", length(all_samples)))
  box(0.32, 0.50, 0.42, 0.14, sprintf("Samples analysed: n = %d", n_analysed),
      if (length(flagged)) paste0("borderline coverage (FLAG): ", ids(flagged)) else "")
  box(0.32, 0.18, 0.42, 0.11, sprintf("Samples with somatic variants: n = %d", n_analysed - length(no_variants)))
  arr(0.32, 0.745, 0.32, 0.57); arr(0.32, 0.43, 0.32, 0.235)
  # exclusion boxes
  box(0.77, 0.65, 0.40, 0.17, sprintf("Excluded: insufficient coverage\n(%s): n = %d", paste(exclude_qc, collapse = "/"), length(excluded)),
      ids(excluded), fill = "#FDF0F0", border = "#E45756")
  box(0.77, 0.34, 0.40, 0.17, sprintf("Samples without somatic variants\npassing the filters: n = %d", length(no_variants)),
      ids(no_variants), fill = "#FFF7E8", border = "#F58518")
  arr(0.32, 0.65, 0.57, 0.65); arr(0.32, 0.34, 0.57, 0.34)
  grid.text(paste(sprintf("Coverage QC: PASS > %g%%, FLAG > %g%% of target bases covered > 30x.",
                          100 * as.numeric(snakemake@params[["pass_30x"]]), 100 * as.numeric(snakemake@params[["flag_30x"]])),
                  sprintf("Variant filters: covered regions (+%g bp), depth >= %g, alt reads >= %g, VAF >= %g, gnomAD AF <= %g.",
                          bed_pad, min_depth, min_alt, min_vaf, max_popaf), sep = "\n"),
            0.5, 0.045, gp = gpar(fontsize = 8, col = "#666666"))
  dev.off()
}
draw_flow(snakemake@output[["flow"]])

# 2. Harmonise Variant_Classification -------------------------------------------
if (!is.null(mut_df)) {
  vc <- mut_df$Variant_Classification
  vc[vc %in% names(vc_map)] <- vc_map[vc[vc %in% names(vc_map)]]
  unknown <- setdiff(unique(vc), vc_known)
  if (length(unknown) > 0) {
    message("WARNING: Variant_Classification values not recognised by maftools (counted as non-coding): ",
            paste(sprintf("%s (%d)", unknown, sapply(unknown, function(u) sum(vc == u))), collapse = ", "),
            "\n  -> add them to vc_map in scripts/maf_report.R if they should be mapped differently")
  }
  mut_df$Variant_Classification <- vc
}

# 3. Write the merged MAF and the plots -----------------------------------------
if (is.null(mut_df)) {
  message("no variants passed the filters in any sample")
  file.create(snakemake@output[["maf"]])
  for (p in c(snakemake@output[["summary"]], snakemake@output[["oncoplot"]])) {
    pdf(p); plot.new(); text(0.5, 0.5, "No variants passed the filters"); dev.off()
  }
} else {
  write.table(mut_df, snakemake@output[["maf"]], sep = "\t", quote = FALSE, row.names = FALSE)
  maf <- read.maf(mut_df)

  # Colour palette for the mutation types (colour-blind friendly, Tableau 10 based)
  vc_colors <- c(Missense_Mutation      = "#4C78A8",
                 Nonsense_Mutation      = "#E45756",
                 Frame_Shift_Del        = "#B279A2",
                 Frame_Shift_Ins        = "#FF9DA7",
                 In_Frame_Del           = "#72B7B2",
                 In_Frame_Ins           = "#54A24B",
                 Splice_Site            = "#F58518",
                 Translation_Start_Site = "#EECA3B",
                 Nonstop_Mutation       = "#BAB0AC",
                 Multi_Hit              = "#3B3B3B")

  pdf(snakemake@output[["summary"]], width = 11, height = 7)
  plotmafSummary(maf, addStat = "median", dashboard = TRUE, color = vc_colors)
  dev.off()

  # Oncoplot: top N genes, title with the actual number shown, margin sized to the longest gene name
  gene_summary <- getGeneSummary(maf)
  n_top    <- min(top_genes, nrow(gene_summary))
  if (n_top == 0) {
    pdf(snakemake@output[["oncoplot"]]); plot.new()
    text(0.5, 0.5, "No protein-altering variants to show"); dev.off()
    quit(save = "no")
  }
  genes    <- gene_summary$Hugo_Symbol[seq_len(n_top)]
  n_alt    <- length(unique(subsetMaf(maf, genes = genes, mafObj = FALSE)$Tumor_Sample_Barcode))
  title    <- sprintf("Top %d most frequently mutated genes - altered in %d of %d analysed samples (%.0f%%)",
                      n_top, n_alt, n_analysed, 100 * n_alt / n_analysed)
  # margin (in text lines) for the gene names: automatic from the longest name, or set in the config
  gene_mar <- if (gene_mar_cfg > 0) gene_mar_cfg else max(8, 1.0 * max(nchar(genes)) + 3)

  pdf(snakemake@output[["oncoplot"]], width = 11, height = max(5, 2.5 + 0.3 * n_top))
  tryCatch(oncoplot(maf, top = n_top, colors = vc_colors, titleText = title, cohortSize = n_analysed,
                    gene_mar = gene_mar, draw_titv = FALSE),
           error = function(e) { plot.new(); text(0.5, 0.5, paste("Oncoplot not possible:", conditionMessage(e))) })
  dev.off()
}
