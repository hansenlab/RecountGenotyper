#' @export
GetM_S <- function(coverage_count, alt_count){
#Compute M and S values.
ref_count <- coverage_count$coverage_count - alt_count
#We sometimes have cases where the alt counts > coverage counts (bigWig):
#the alt reads were processed to keep overlapping pair-end reads
#whereas the coverage counts (bigWig) did not keep overlapping pair-end reads.
#this is an ad hoc way to deal with negative counts in `ref_count`.
ref_count[ref_count < 0] <- 0
#To prevent divde by 0 in M, S calculation, we add psuedocount of 1.
M <- log2((ref_count + 1) / (alt_count + 1))
S <- log2(sqrt((ref_count + 1) * (alt_count + 1)))
return(list(M,S))
}
