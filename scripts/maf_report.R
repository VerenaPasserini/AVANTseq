# -----------------------------------------------------------------------------
# AVANTseq report: clean and merge the per-sample Funcotator MAFs, then plot a
# maftools summary and an oncoplot.
# Run by Snakemake (rule maf_report in rules/report.smk).
# -----------------------------------------------------------------------------

log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

suppressPackageStartupMessages({ library(maftools); library(data.table) })

bed_file  <- snakemake@input[["bed"]]
min_depth <- as.numeric(snakemake@params[["min_depth"]])
min_alt   <- as.numeric(snakemake@params[["min_alt"]])
min_vaf   <- as.numeric(snakemake@params[["min_vaf"]])
max_popaf <- as.numeric(snakemake@params[["max_pop_af"]])   # max gnomAD allele frequency
top_genes <- as.numeric(snakemake@params[["top_genes"]])    # number of genes in the oncoplot
drop_unnamed <- isTRUE(as.logical(snakemake@params[["drop_unnamed"]]))  # drop features without a gene symbol (ENSG IDs)
exclude_qc <- unlist(snakemake@params[["exclude_qc"]])      # QC classes to leave out (e.g. FAIL)

# Coverage QC from qc_report: leave out samples with an excluded QC class
qc <- fread(snakemake@input[["qc"]], select = c("SAMPLE", "QC"), colClasses = list(character = "SAMPLE"))
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
                  start = suppressWarnings(as.numeric(vapply(f, `[`, "", 2))) + 1,
                  end   = suppressWarnings(as.numeric(vapply(f, `[`, "", 3))))
bed <- bed[!is.na(bed$start) & !is.na(bed$end), ]          # drops a possible header line
bed <- bed[order(bed$chr, bed$start), ]
bed_merged <- lapply(split(bed, bed$chr), function(b) {
  grp <- cumsum(c(TRUE, b$start[-1] > cummax(b$end)[-nrow(b)] + 1))
  data.frame(start = as.numeric(tapply(b$start, grp, min)), end = as.numeric(tapply(b$end, grp, max)))
})
message(sprintf("BED: %d regions on %d chromosomes", nrow(bed), length(bed_merged)))

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
  df <- df[in_bed(as.character(df$Chromosome), df$Start_Position,
                  pmax(df$Start_Position, df$End_Position)), , drop = FALSE]
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
  if (nrow(df) == 0) return(NULL)
  df$Tumor_Sample_Barcode <- sample
  df
}
message("Variants per sample (PASS variants from Funcotator):")
mut_df <- do.call(rbind, lapply(maf_files, read_sample))

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
                 Frame_Shift_Ins        = "#9D755D",
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
  n_all    <- as.numeric(maf@summary$summary[maf@summary$ID == "Samples"])
  title    <- sprintf("Top %d most frequently mutated genes - altered in %d of %d samples (%.0f%%)",
                      n_top, n_alt, n_all, 100 * n_alt / n_all)
  gene_mar <- max(5, 0.6 * max(nchar(genes)) + 1)

  pdf(snakemake@output[["oncoplot"]], width = 11, height = max(5, 2.5 + 0.3 * n_top))
  tryCatch(oncoplot(maf, top = n_top, colors = vc_colors, titleText = title,
                    gene_mar = gene_mar, draw_titv = FALSE),
           error = function(e) { plot.new(); text(0.5, 0.5, paste("Oncoplot not possible:", conditionMessage(e))) })
  dev.off()
}
