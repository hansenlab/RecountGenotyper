
totalPtm <- proc.time()

recount3_chr_mapping <- "/dcl02/lieber/ajaffe/recount-pump/recount3.alts.chromosome_mappings.tsv"




#' @export
writeNullOutput = function(resultFile, latticeFile) {
	#Write a blank output when there is no data to genotype. We do not throw an error
	#because it allows the pipeline to keep running, and it will be apparent when the user
	#see that there's nothing that was genotyped.
	cat("Error loading in bigwig or no data to genotype! Writing blank results.\n")
	result <- data.table(chr = character(), start = numeric(), AF = numeric(), M = numeric(), S = numeric(), coverage = numeric(), pred_genotype = numeric())
	fwrite(result, file = resultFile)
	fwrite(result, file = latticeFile)
}


cat("Loading in SNPs and bigWig.\n")
ptm <- proc.time()
coverage_count_result <- get_coverage_count(snps_path = opt$SNPs,
											bigWig_path = opt$total,
											coverage_cutoff = as.numeric(opt$cutoff))

print(gc())
print(proc.time() - ptm)
cat("Finished loading in SNPs and bigWig.\n")

if(any(is.na(coverage_count_result)) | length(coverage_count_result$coverage_count) == 0) {
	writeNullOutput(opt$result, opt$outOfLattice)
	q()
}

cat("Loading in alt counts.\n")
ptm <- proc.time()

alt_count <- get_alt_count(alt_path = opt$alt,
						   sample_id_rep = opt$sample_id_rep,
						   filtered_snps_gr = coverage_count_result$filtered_snps_gr,
						   temp_folder = paste0(opt$tempFolder, "/"))

print(gc())
print(proc.time() - ptm)
cat("Finished loading in alt counts.\n")

if(all(is.na(alt_count))) {
	writeNullOutput(opt$result, opt$outOfLattice)
	q()
}

#Compute M and S values.
ref_count <- coverage_count_result$coverage_count - alt_count
#We sometimes have cases where the alt counts > coverage counts (bigWig):
#the alt reads were processed to keep overlapping pair-end reads
#whereas the coverage counts (bigWig) did not keep overlapping pair-end reads.
#this is an ad hoc way to deal with negative counts in `ref_count`.
ref_count[ref_count < 0] <- 0
#To prevent divde by 0 in M, S calculation, we add psuedocount of 1.
M <- log2((ref_count + 1) / (alt_count + 1))
S <- log2(sqrt((ref_count + 1) * (alt_count + 1)))

cat("Predicting genotype.\n")
ptm <- proc.time()

pred_genotype_result <- predict_genotype(model_path = opt$model,
										 model_lattice_path = opt$modelLattice,
										 prior_path = opt$prior,
										 M = M,
										 S = S)

print(gc())
print(proc.time() - ptm)
cat("Finished predicting genotype. Saving results now.\n")

#Predicted genotype result.
result <- data.table(chr = as.character(seqnames(coverage_count_result$filtered_snps_gr[pred_genotype_result$withinLattice_idx])),
					start = start(coverage_count_result$filtered_snps_gr[pred_genotype_result$withinLattice_idx]),
					AF = round(coverage_count_result$filtered_snps_gr$allele_freq[pred_genotype_result$withinLattice_idx], 4),
					M = M[pred_genotype_result$withinLattice_idx],
					S = S[pred_genotype_result$withinLattice_idx],
					coverage = ref_count[pred_genotype_result$withinLattice_idx] + alt_count[pred_genotype_result$withinLattice_idx],
					pred_genotype = pred_genotype_result$predicted_genotype)
fwrite(result, file = opt$result)

#M, S values that fall outside of our lattice, to be predicted with model file later.
outOfLattice_result <- data.table(chr = as.character(seqnames(coverage_count_result$filtered_snps_gr[pred_genotype_result$outsideLattice_idx])),
								start = start(coverage_count_result$filtered_snps_gr[pred_genotype_result$outsideLattice_idx]),
								AF = round(coverage_count_result$filtered_snps_gr$allele_freq[pred_genotype_result$outsideLattice_idx], 4),
								M = M[pred_genotype_result$outsideLattice_idx],
								S = S[pred_genotype_result$outsideLattice_idx],
								coverage = ref_count[pred_genotype_result$outsideLattice_idx] + alt_count[pred_genotype_result$outsideLattice_idx],
								pred_genotype = NA)
fwrite(outOfLattice_result, opt$outOfLattice)


cat("Total time elapsed:\n")
print(proc.time() - totalPtm)


