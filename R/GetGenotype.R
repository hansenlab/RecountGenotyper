#' Title
#'
#' @param model_path
#' @param model_lattice_path
#' @param modelPrior
#' @param M
#' @param S
function(model_path, M, S) {
  lattice_max <- 8.5
  #Load in model information.
  if(all(is.na(modelPrior)) == T) {
    modelPrior = c(.94, .04, .02)
  }
  #Figure out which values of predictor S can be used for lattice lookup table.
  withinLattice_idx <- which(S <= lattice_max)
  outsideLattice_idx <- which(S > lattice_max)
  cat("Proportion of S values that can be inferred using lattice: ", length(withinLattice_idx)/length(S), "\n")
  cat("Number of S values to infer with model later: ", length(outsideLattice_idx), "\n")
  #Compute likelihood * modelPrior for each genotype.
  likelihood_times_modelPrior <- sapply(1:3, function(k){
    #Inference on mu_M and sd_M
    mu_M <- rep(NA, length(withinLattice_idx))
    sd_M <- rep(NA, length(withinLattice_idx))
    if(k == 1) {
      #modelLattice is a data.table with key set to its S values,
      #so we use the key to get lattice values.
      mu_M <- modelLattice[J(S[withinLattice_idx])]$mu_M_prediction_1
      sd_M <- modelLattice[J(S[withinLattice_idx])]$sd_M_prediction_1
    }else if(k == 2) {
      mu_M <- modelLattice[J(S[withinLattice_idx])]$mu_M_prediction_2
      sd_M <- modelLattice[J(S[withinLattice_idx])]$sd_M_prediction_2
    }else if(k == 3) {
      mu_M <- modelLattice[J(S[withinLattice_idx])]$mu_M_prediction_3
      sd_M <- modelLattice[J(S[withinLattice_idx])]$sd_M_prediction_3
    }
    sd_M[sd_M <= 0] = 0 #if we predict (extrapolate) any variance to be <= 0, set it to 0.
    return(modelPrior[k] * dnorm(x = M[withinLattice_idx], mean = mu_M, sd = sd_M, log = FALSE))
  })


  posterior <- sapply(1:3, function(k){
    if(is.null(dim(likelihood_times_modelPrior))) { #if we only have one SNP that needs to be genotyped
      return(likelihood_times_modelPrior[k] / sum(likelihood_times_modelPrior))
    }else { #else, we have many SNPs in a dataframe
      return(likelihood_times_modelPrior[, k] / rowSums(likelihood_times_modelPrior))
    }
  })

  #Pick genotype with the max posterior genotype as our prediction.
  if(is.null(dim(posterior))) { #if we only have one SNP that needs to be genotyped
    predicted_genotype <- which(posterior == max(posterior))
  }else { #else, we have many SNPs in a dataframe
    predicted_genotype <- apply(posterior, MARGIN = 1, FUN = function(x) which(x == max(x)))
  }


  return(list(predicted_genotype = predicted_genotype,
              withinLattice_idx = withinLattice_idx,
              outsideLattice_idx = outsideLattice_idx))
}