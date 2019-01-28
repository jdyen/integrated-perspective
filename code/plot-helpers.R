# need a density dependence function
dens_scale <- function(x, param) {
  1 / (1 + x * param)
}

# need an update function
abund_est <- function(mat, density, inits, stage_age, detection, n_iter, n_burn = 0) {
  
  out <- matrix(0, ncol = n_iter + n_burn + 1, nrow = length(inits))
  out[, 1] <- inits
  for (i in seq_len(n_iter + n_burn)) {
    rescale <- dens_scale(sum(out[, i]), density)
    mat_tmp <- rescale * mat
    out[, i + 1] <- mat_tmp %*% out[, i]
  }
  
  out <- stage_age %*% out
  
  out <- detection * out
  
  out[, (n_burn + 2):ncol(out)]
  
}

# sample from the posterior and project abundances through time
sample_draws <- function(x, site, abund_data, chain = 1,
                         p = c(0.1, 0.5, 0.9),
                         nsim = 1000,
                         simplify = TRUE) {
  
  if (chain == "all") {
    draws <- do.call(rbind, draws)
  } else {
    draws <- draws[[chain]]
  }
  
  idx <- sample(seq_len(nrow(draws)),
                size = nsim, replace = FALSE)
  vars <- colnames(draws)
  
  
  x_sub <- draws[idx, ]
  
    
  mat_vals <- x_sub[, grep("mat_vec", vars)]
  inits <- x_sub[, grep("inits_mean", vars)]
  dens_param <- x_sub[, grep("density_depend", vars)]
  detection <- x_sub[, grep("detection_mean", vars)]
  stage_age_conv <- x_sub[, grep("stage_age", vars)]
  
  nyear <- length(abund_data)
    
  abund_out <- matrix(NA, nrow = nsim, ncol = nyear)
  for (i in seq_len(nsim)) {
    
    mat_full <- array(mat_vals[i, ], c(5, 5, 6))
    stage_age_tmp <- matrix(stage_age_conv[i, ], nrow = 8)
    abund_tmp <- abund_est(mat_full[, , site], dens_param[i], inits[i, ], stage_age_tmp, detection[i], nyear)
    abund_out[i, ] <- apply(abund_tmp, 2, sum)
    
  }
  

  if (simplify)
    abund_out <- apply(abund_out, 2, quantile, p = p)
  
  abund_out

}



