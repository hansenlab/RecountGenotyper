#' Get M and S from raw reads as well as the predicted accuracy for that SNP
#'
#' @name GetMandS
#' @description This function will generate the reference and alternative base counts for each loci. The alternative base input is the custom file format (.zst file) from
#' the Recount3. This file includes read counts for the 4 bases. Therefore, in this function, we will find the base that has read counts
#' mapped to it. The default coverage cutoff is set to 4. After the ref and alt counts are obtained, the M and S values will be calculated and returened.
#'
#' @param snps_gr The biallelic SNPs to be genotyped. This file should be a grange object containing seqnames,ranges,allele_freq. This file is not required.
#' @param bigWig_path Path to the bigWig file containing the total read counts for a single sample.
#' @param coverage_cutoff The minimum amount of read count mapped to a loci for that SNP to be included in the genotype calling. The default is 4 which was used in our model training and testing.
#' @param alt_path Path to alternative base read counts. This is outputed from the Recount3 pipeline as a .zst file.
#' @param sample_id_rep Single sample ID to be genotyped.
#' @param temp_folder Path to temporary folder.
#'
#' @importFrom S4Vectors subjectHits
#' @importFrom S4Vectors queryHits
#' @importFrom IRanges IRanges
#' @importFrom dplyr case_when
#'
#' @import GenomicRanges
#' @import rtracklayer
#' @import data.table
#' @import IRanges
#' @import dplyr
#' @import rms
#' @import tidyverse
#' @examples
#'
#' base_path<-system.file("extdata", package="RecountGenotyper")
#' alt_path=paste0(base_path,"/test.csv.zst")
#' bigWig_path=paste0(base_path,"/test.bw")
#'
#' sample_id_rep="test"
#' temp_folder= tempdir()
#'
#'  test_geno<-GetMandS(snps_gr=NULL, bigWig_path, coverage_cutoff=4, alt_path, sample_id_rep, temp_folder= tempdir())
#'
#' @export
GetMandS<-function(snps_gr=NULL, bigWig_path, coverage_cutoff=4,alt_path, sample_id_rep, temp_folder= tempdir()) {
 if (is.null(snps_gr)){
  #load in snp granges:
   cat("Loading in SNP granges")
   snps<- #This data is protected
   download.file(snps, destfile=paste0(tempdir(),"/biallelic_snp_gr.rds"), method="curl")
   snps_gr<-readRDS(paste0(tempdir(),"/biallelic_snp_gr.rds"))
  }
  #Load in bigWig file to get `coverage_count` and `filtered_snps_gr`.
  cat("Loading in: ", bigWig_path, "\n")

  bigwig <- tryCatch(
    {
     rtracklayer::import(bigWig_path, format = "bigwig")
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
  overlap_loci <- GenomicRanges::findOverlaps(snps_gr, bigwig)
  bigwig_count <- bigwig$score[subjectHits(overlap_loci)]
  filter_idx <- which(bigwig_count >= coverage_cutoff)
  bigwig_count <- bigwig_count[filter_idx]
  filtered_snps_gr <- snps_gr[queryHits(overlap_loci)[filter_idx]]



  #Load in alt file to construct `alt_count`:
  #1. Decompress .zst via zstdcat in command line, save output as *.alt.temp.csv
  #2. Read in *.alt.temp.csv via vroom (fast), convert to data.table in order
  #to collapase alt counts.
  #3. Collapse to `alt_count`, and create GRanges object for intersection.
  cat("Loading in: ", alt_path, "\n")
  temp_altFile <- paste0(temp_folder, sample_id_rep, ".alt.temp.csv")
  tryCatch({
    system(paste0("zstdcat ", alt_path, " > ", temp_altFile))
    },
           error = function(e){
             message("Error reading alt file:\n", e)
             message("'zstdcat' software is not installed")
           },
           warning = function(w){
             message("Warning reading alt file:\n", w)
             message("'zstdcat' software is not installed")
           },
           finally = {
           })


  alt <- data.table::fread(temp_altFile, select = c(1, 2, 4))
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
  alt_count$chr <- recount3_chr_mapping$V2[match(alt_count$chr, recount3_chr_mapping$V1)]
  alt_count_gr <- GenomicRanges::GRanges(seqnames = alt_count$chr,
                          ranges = IRanges(alt_count$pos,
                                           alt_count$pos))
  ov <- GenomicRanges::findOverlaps(alt_count_gr, filtered_snps_gr)
  #Subset `alt_count_gr`, `alt_count` to the positions from SNPs, and compute
  #its own unique key.
  alt_count_gr <- alt_count_gr[queryHits(ov)]
  alt_count <- alt_count[queryHits(ov)]
  alt_key <- paste(as.character(seqnames(alt_count_gr)), start(alt_count_gr),
                   filtered_snps_gr$ref_seq[subjectHits(ov)], alt_count$alt,  sep = "_")
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

  #Compute M and S values.
  ref_count <- bigwig_count - final_alt_count
  #We sometimes have cases where the alt counts > coverage counts (bigWig):
  #the alt reads were processed to keep overlapping pair-end reads
  #whereas the coverage counts (bigWig) did not keep overlapping pair-end reads.
  #this is an ad hoc way to deal with negative counts in `ref_count`.
  ref_count[ref_count < 0] <- 0
  #To prevent divde by 0 in M, S calculation, we add psuedocount of 1.
  M <- log2((ref_count + 1) / (final_alt_count + 1))
  S <- log2(sqrt((ref_count + 1) * (final_alt_count + 1)))



  ############
  expit <- function(x) {
    return(1/(1+exp(-x)))
  }


  eval_data <- data.frame(coverage = bigwig_count, AF = filtered_snps_gr$allele_freq)
  eval_data$major_AF <- case_when(eval_data$AF <= .5  ~ 1 - eval_data$AF,
                                  eval_data$AF > .5  ~ eval_data$AF)

  #prediction for low major AF
  id<-which(eval_data$major_AF < .95)
  eval_low_majorAF <- eval_data %>% filter(major_AF < .95)
  eval_low_majorAF$majorAF_bin <- cut(eval_low_majorAF$major_AF,
                                      breaks = c(0.500,0.632,0.908,0.839,0.749,0.950),
                                      include.lowest=TRUE)
  eval_low_majorAF <- left_join(eval_low_majorAF, accuracyModelLattice, by = c("coverage", "majorAF_bin"))
  eval_data$accuracy[id] <- unlist(eval_low_majorAF$predicted_accuracy)

  #prediction for high major AF
  id<-which(eval_data$major_AF >= .95)
  eval_high_majorAF <- eval_data %>% filter(major_AF >= .95)
  eval_high_majorAF$majorAF_bin <- "[0.95, 1]"
  eval_high_majorAF <- left_join(eval_high_majorAF, accuracyModelLattice, by = c("coverage", "majorAF_bin"))
  eval_data$accuracy[id] <- unlist(eval_high_majorAF$predicted_accuracy)
  #put the two together
  eval_data$accuracy[is.na(eval_data$accuracy)] <- 1 #for coverage > 100 outside of lattice, give it perfect prediction


return(data.table(chr = as.character(seqnames(filtered_snps_gr)),
                  start=start(filtered_snps_gr),
                  AF= filtered_snps_gr$allele_freq,
                  M=M,
                  S=S,
                  coverage=ref_count+final_alt_count,
                  predicted_geno_acc=eval_data$accuracy))
}
