#' Generate Annotation Data Files for eRNAFinder
#'
#' Usage:
#' source("inst/scripts/generate_annotation_data.R")
#' generate_annotation_files(
#'   gencode_gtf_path = "/path/to/gencode.v49.basic.annotation.gtf",
#'   output_dir = "ext/"
#' )

generate_annotation_files <- function(gencode_gtf_path, output_dir = "extdata/") {

  if (!file.exists(gencode_gtf_path)) {
    stop("GENCODE GTF file not found: ", gencode_gtf_path)
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # Load and filter GENCODE
  gencode <- rtracklayer::import(gencode_gtf_path)
  gencode_tx <- gencode[gencode$type == "transcript"]

  # Exclusion transcripts (non-lncRNA)
  exclusion_transcripts <- gencode_tx[gencode_tx$gene_type != "lncRNA"]
  exclusion_transcripts$score <- 1000
  exclusion_output <- file.path(output_dir, "exclusion_transcripts.bed")
  rtracklayer::export(exclusion_transcripts, exclusion_output, format = "bed")

  # uaRNA (protein-coding promoters with flipped strand)
  protein_coding_tx <- gencode_tx[gencode_tx$gene_type == "protein_coding"]
  promoter_region <- GenomicRanges::promoters(protein_coding_tx, upstream = 300, downstream = 1000)
  GenomicRanges::strand(promoter_region) <- ifelse(
    GenomicRanges::strand(promoter_region) == "+", "-", "+"
  )
  promoter_region$score <- 1000
  uaRNA_output <- file.path(output_dir, "uaRNA.bed")
  rtracklayer::export(promoter_region, uaRNA_output, format = "bed")

  cat("Generated:\n")
  cat("  -", exclusion_output, "\n")
  cat("  -", uaRNA_output, "\n")
}
