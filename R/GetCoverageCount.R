#' @export
GetCoverageCount <- function(snps_path, bigWig_path, coverage_cutoff) {
  #Load in bigWig file to get `coverage_count` and `filtered_snps_gr`.
  snps_gr <- readRDS(snps_path)
  cat("Loading in: ", bigWig_path, "\n")

  bigwig <- tryCatch(
    {
      import(bigWig_path, format = "bigwig")
    },
    error=function(cond) {
      message(paste("Error loading bigwig file: ", bigWig_path))
      message(cond)
      return(NA)
    },
    warning=function(cond) {
      message(paste("Warning loading bigwig file: ", bigWig_path))
      message(cond)
      return(NA)
    },
    finally={}
  )
  if(all(is.na(bigwig))) {
    return(NA)
  }
  overlap_loci <- findOverlaps(snps_gr, bigwig)
  bigwig_count <- bigwig$score[subjectHits(overlap_loci)]
  filter_idx <- which(bigwig_count >= coverage_cutoff)
  bigwig_count <- bigwig_count[filter_idx]
  filtered_snps_gr <- snps_gr[queryHits(overlap_loci)[filter_idx]]
  return(list(coverage_count = bigwig_count,
              filtered_snps_gr = filtered_snps_gr))
}


