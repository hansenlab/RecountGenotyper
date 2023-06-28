#' Predict the sample genotype using M and S values
#'
#' @name GetGenotype
#' @description This function will use the M and S values generated from the GetMandS function to predict the genotype.
#'
#'
#'
#'
#' @param M M values calculated from GetMandS function
#' @param S S value calculated from GetMandS function
#' @param prior prior was calculated based on the 3 genotype distribution in our training set and default is c(0.93752717, 0.03951271, 0.02239503). This value is constant in our model.
#' @param model Genotyping model
#' @return array of the predicted genotype. The SNP order is the same as output from GetMandS function.
#'
#' @examples
#'
#' S=test$S
#' M=test$M
#'
#' test$predicted_genotype<-GetGenotype(model, M, S)
#'
#'
#' @export
GetGenotype<- function(model, prior=c(0.93752717, 0.03951271, 0.02239503), M, S) {
  lattice_max <- 8.5
  #Figure out which values of predictor S can be used for lattice lookup table.
  withinLattice_idx <- which(S <= lattice_max)
  outsideLattice_idx <- which(S > lattice_max)

  #Compute likelihood * prior for each genotype.
  likelihood_times_prior <- sapply(1:3, function(k){
    #Inference on mu_M and sd_M
    mu_M <- rep(NA, length(S))
    sd_M <- rep(NA, length(S))
    if(k == 1){
      #modelLattice is a data.table with key set to its S values,
      #so we use the key to get lattice values.
      mu_M[withinLattice_idx] <- modelLattice[J(S[withinLattice_idx])]$mu_M_prediction_1
      sd_M[withinLattice_idx] <- modelLattice[J(S[withinLattice_idx])]$sd_M_prediction_1
      #then, for SNPs out of the lattice, we use the model predictions.
      mu_M[outsideLattice_idx] <- predict(model[[1]][[1]], data.frame(mu_S = unlist(S[outsideLattice_idx])))
      sd_M[outsideLattice_idx] <- predict(model[[1]][[2]], data.frame(mu_S = unlist(S[outsideLattice_idx])))
    } else if(k == 2) {
      mu_M[withinLattice_idx] <- modelLattice[J(S[withinLattice_idx])]$mu_M_prediction_2
      sd_M[withinLattice_idx] <- modelLattice[J(S[withinLattice_idx])]$sd_M_prediction_2
      mu_M[outsideLattice_idx] <- predict(model[[2]][[1]], data.frame(mu_S = unlist(S[outsideLattice_idx])))
      sd_M[outsideLattice_idx] <- predict(model[[2]][[2]], data.frame(mu_S = unlist(S[outsideLattice_idx])))
    } else if(k == 3) {
      mu_M[withinLattice_idx] <- modelLattice[J(S[withinLattice_idx])]$mu_M_prediction_3
      sd_M[withinLattice_idx] <- modelLattice[J(S[withinLattice_idx])]$sd_M_prediction_3
      mu_M[outsideLattice_idx] <- predict(model[[3]][[1]], data.frame(mu_S = unlist(S[outsideLattice_idx])))
      sd_M[outsideLattice_idx] <- predict(model[[3]][[2]], data.frame(mu_S = unlist(S[outsideLattice_idx])))
    }
    sd_M[sd_M <= 0] <- 0 #if we predict (extrapolate) any variance to be <= 0, set it to 0.
    return(prior[k] * dnorm(x = M, mean = mu_M, sd = sd_M, log = FALSE))
  })


  posterior <- sapply(1:3, function(k){
    if(is.null(dim(likelihood_times_prior))) { #edge case: if we only have one SNP that needs to be genotyped
      return(likelihood_times_prior[k] / sum(likelihood_times_prior))
    }else { #else, we have many SNPs in a dataframe
      return(likelihood_times_prior[, k] / rowSums(likelihood_times_prior))
    }
  })

  #Pick genotype with the max posterior genotype as our prediction.
  if(is.null(dim(posterior))) { #edge case: if we only have one SNP that needs to be genotyped
    predicted_genotype <- which(posterior == max(posterior))
  }else { #else, we have many SNPs in a dataframe
    predicted_genotype <- apply(posterior, MARGIN = 1, FUN = function(x) which(x == max(x)))
  }

  return(predicted_genotype)
}
