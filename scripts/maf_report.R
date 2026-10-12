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

# Covered regions (BED, 0-based start) -> 1-based closed intervals
bed_lines <- grep("^(#|track|browser)", readLines(bed_file), value = TRUE, invert = TRUE)
bed <- fread(text = bed_lines, header = FALSE, select = 1:3, col.names = c("chr", "start", "end"))
bed[, `:=`(chr = as.character(chr), start = start + 1L)]
setkey(bed, chr, start, end)

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
  # skip the "#" header lines at the top only ("#" may also occur inside annotation fields)
  n_header <- sum(cumprod(startsWith(readLines(f, n = 2000), "#")))
  df <- read.delim(f, skip = n_header, comment.char = "", quote = "", stringsAsFactors = FALSE)
  n_all <- nrow(df)

  # keep variants overlapping the covered regions of the BED file
  pos <- data.table(chr = as.character(df$Chromosome), start = df$Start_Position,
                    end = pmax(df$Start_Position, df$End_Position))
  hit <- unique(foverlaps(pos, bed, type = "any", which = TRUE, nomatch = NULL)$xid)
  df <- df[sort(hit), , drop = FALSE]
  n_bed <- nrow(df)

  df$t_depth <- df$t_alt_count + df$t_ref_count
  df$VAF     <- df$t_alt_count / df$t_depth
  df <- df[df$t_depth >= min_depth & df$t_alt_count >= min_alt & df$VAF >= min_vaf, , drop = FALSE]
  n_reads <- nrow(df)

  # population frequency: POPAF (from Mutect2) is -log10 of the gnomAD allele frequency;
  # variants absent from gnomAD get a high POPAF and are kept
  popaf <- suppressWarnings(as.numeric(sub(",.*", "", df$POPAF)))
  df <- df[is.na(popaf) | popaf >= -log10(max_popaf), , drop = FALSE]

  message(sprintf("%-20s all: %5d   in BED: %5d   depth/alt/VAF: %5d   gnomAD AF <= %g: %5d",
                  sample, n_all, n_bed, n_reads, max_popaf, nrow(df)))
  if (nrow(df) == 0) return(NULL)
  df$Tumor_Sample_Barcode <- sample
  df
}
message("Variants per sample (PASS variants from Funcotator):")
mut_df <- do.call(rbind, lapply(snakemake@input[["mafs"]], read_sample))

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

  pdf(snakemake@output[["summary"]], width = 10, height = 7)
  plotmafSummary(maf, addStat = "median", dashboard = TRUE)
  dev.off()

  pdf(snakemake@output[["oncoplot"]], width = 10, height = 7)
  tryCatch(oncoplot(maf, draw_titv = FALSE),
           error = function(e) { plot.new(); text(0.5, 0.5, paste("Oncoplot not possible:", conditionMessage(e))) })
  dev.off()
}
