

base_path<-system.file("extdata", package="RecountGenotyper")
alt=paste0(base_path,"/test.csv.zst")
bw=paste0(base_path,"/test.bw")


test_that("M and S can be calculated", {
  expect_equal(c(snps_path, bigWig_path=bw, coverage_cutoff=4,alt_path=alt, sample_id_rep="test", temp_folder="~"),)
})
