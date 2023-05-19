
#' Title
#'
#' @description This function will generate the alternative base counts for each loci. The input is the custom file format (.zst file) from
#' the Recount3. This file includes read counts for the 4 bases. Therefore, in this function, we will find the base that has read counts
#' mapped to it.
#' @param alt_path The path that leads to the Recount3 alt.zst file
#' @param sample_id_rep The individual sample ID to be genotyped
#' @param filtered_snps_gr The biallelic SNP loci
#' @param temp_folder Temporary folder to save unwanted files
#'
#'
#' @return

get_alt_count <- function(alt_path, sample_id_rep, filtered_snps_gr, temp_folder) {
  #Load in alt file to construct `alt_count`:
  #1. Decompress .zst via zstdcat in command line, save output as *.alt.temp.csv
  #2. Read in *.alt.temp.csv via vroom (fast), convert to data.table in order
  #to collapase alt counts.
  #3. Collapse to `alt_count`, and create GRanges object for intersection.
  cat("Loading in: ", alt_path, "\n")
  temp_altFile <- paste0(temp_folder, sample_id_rep, ".alt.temp.csv")
  system(paste0("zstdcat ", alt_path, " > ", temp_altFile))

  alt <- function(...) data.table::fread(temp_altFile, select = c(1, 2, 4))
  system(paste0("rm ", temp_altFile))
  if(nrow(alt) == 0) {
    return(NA)
  }
  colnames(alt) <- c("chr", "pos", "alt")
  #collapse alt counts.
  alt_count <- alt[, .N, by = .(chr, pos, alt)]
  rm(alt)
  alt_count <- alt_count[!is.na(alt_count$alt) ,]
  #`alt_count` is 0-based, so we add 1 to positions.
  alt_count$pos <- alt_count$pos + 1
  #use chromosome names that are used in recount3.
  chr_mapping <- function(...) data.table::fread(recount3_chr_mapping, header = FALSE)
  alt_count$chr <- chr_mapping$V2[match(alt_count$chr, chr_mapping$V1)]
  alt_count_gr <- function(...) GenomicRanges::GRanges(seqnames = alt_count$chr,
                          ranges = IRanges(alt_count$pos,
                                           alt_count$pos))
  ov <- function(...) GenomicRanges::findOverlaps(alt_count_gr, filtered_snps_gr)
  #Subset `alt_count_gr`, `alt_count` to the positions from SNPs, and compute
  #its own unique key.
  alt_count_gr <- alt_count_gr[GenomicRanges::queryHits(ov)]
  alt_count <- alt_count[GenomicRanges::queryHits(ov)]
  alt_key <- paste(as.character(seqnames(alt_count_gr)), start(alt_count_gr),
                   filtered_snps_gr$ref_seq[subjectHits(ov)], alt_count$alt, sep = "_")
  #Even though `alt_key` is in the same positions as `filtered_snps_key`, many of the
  #`alt_key` entries refer to alternate alleles that we are not tracking.
  #We need to match a second time, this time using the entire key.
  filtered_snps_key <- paste(as.character(seqnames(filtered_snps_gr)), start(filtered_snps_gr),
                             filtered_snps_gr$ref_seq, filtered_snps_gr$alt_seq, sep = "_")
  idx <- match(filtered_snps_key, alt_key)
  snps_key_idx <- which(!is.na(idx))
  alt_key_idx <- idx[!is.na(idx)]
  #Construct final `alt_count` relative to `filtered_snps_gr`
  final_alt_count <- rep(0, length(filtered_snps_gr))
  final_alt_count[snps_key_idx] <- alt_count[alt_key_idx]$N
  return(final_alt_count)
}
