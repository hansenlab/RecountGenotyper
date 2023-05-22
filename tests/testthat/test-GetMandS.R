test_that("Get M and S", {
  metadata<-read.csv("/dcs04/hansen/data/recount_genotype/pipeline/metadata/all_SRA.csv")
  snps_path<- "/dcl01/hansen/data/arazi/gtex_count/genotype/biallelic_SNP_gr.rda"
  bigWig_path<-"/dcl02/lieber/ajaffe/recount-pump/human_tranche_backups/sra_human_v3_6/66/DRP000366/97/DRR000897/sra_human_v3_61_in89_att0/DRR000897!DRP000366!hg38!sra.all.bw"
  alt_path<-metadata$alt[1]
  sample_id_rep<-metadata$sample_id[1]
  temp_folder<-"/dcs04/hansen/data/recount_genotype/pipeline/tmp/"
  recount3_chr_mapping<-"/dcl02/lieber/ajaffe/recount-pump/recount3.alts.chromosome_mappings.tsv"

  expect_equal(GetMandS(), 4)
})
