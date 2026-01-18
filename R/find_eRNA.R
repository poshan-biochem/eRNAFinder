#' Find eRNA (Enhancer RNA) Transcripts
#'
#' @param stg_path Character. Path to StringTie GTF file
#' @param ignore.strand Logical. If FALSE (default), strand-specific matching
#'
#' @return GRanges object of eRNA transcripts
#'
#' @examples
#' \dontrun{
#'   eRNA <- find_eRNA("stringtie.gtf")
#' }
#'
#' @import GenomicRanges
#' @import rtracklayer
#' @import IRanges
#' @export
find_eRNA <- function(stg_path, ignore.strand = FALSE) {

  # Load StringTie
  stg <- rtracklayer::import(stg_path)
  stg <- stg[stg$source == "StringTie" & stg$type == "transcript"]

  # Load exclusion transcripts
  exclusion_path <- system.file("extdata", "exclusion_transcripts.bed", package = "eRNAFinder")
  if (!file.exists(exclusion_path)) {
    stop("Annotation files not found. Run: source(system.file('scripts/generate_annotation_data.R', package='eRNAFinder'))")
  }
  exclusion_transcripts <- rtracklayer::import(exclusion_path)

  # Load uaRNA
  uaRNA_path <- system.file("extdata", "uaRNA.bed", package = "eRNAFinder")
  uaRNA <- rtracklayer::import(uaRNA_path)

  # Find novel transcripts
  novel <- IRanges::subsetByOverlaps(stg, exclusion_transcripts, ignore.strand = ignore.strand, invert = TRUE)

  # Find eRNA
  eRNA <- IRanges::subsetByOverlaps(novel, uaRNA, ignore.strand = ignore.strand, invert = TRUE)

  return(eRNA)
}
