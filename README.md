
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RecountGenotyper

<!-- badges: start -->
<!-- badges: end -->

The goal of RecountGenotyper is to predict the genotype information from
the RNA-seq data in Recount3. The available data in Recount3 is in form
of total bigwig files (.bw) and alternative files (.zst).

## Installation

You can install the development version of RecountGenotyper from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("raziafrooz/RecountGenotyper")
```

## Genotype calling

This package has two functions:

1)  GetMandS()

- This function will take the SNP, bigwig, and alternative file path
  along with sample id and an tempurary folder.
- In this function, SNPs will be filtered if the total coverage is lower
  than 4. After the bigwig files and alternative files are loaded, the M
  and S values are calculated. The prediction accuracy based on the
  Allele frequency and total coverage will also be calculated.

``` r
#test<-GetMandS(snps_path, bigWig_path, coverage_cutoff=4,alt_path, sample_id_rep, temp_folder)
```

2)  GetGenotype()

- This function will take the prediciton model, M, and S values. M and S
  values are provided as columns in GetMandS() funciton.
- This function will output an array of the predicted genotypes. The
  order will be the same as GetMandS() funciton.

``` r
#test$predicted_genotype<-GetGenotype(model, test$M, test$S)
```
