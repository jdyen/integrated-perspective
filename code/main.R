# optional: set working directory
# setwd("PATH/TO/DIR")

# install packages:
# install.packages("devtools")
# devtools::install_github("greta-dev/greta@master")
# install.packages("tensorflow")

# load packages:
# Use the "vary" branch from https://github.com/jdyen/greta.dynamics
# To install:
# devtools::install_github("jdyen/greta.dynamics@vary")
library(greta.dynamics)

# we need to use a few tensorflow fns
tf <- tensorflow::tf

# load pre-compiled data
pop_data <- readRDS("data/compiled-pop-data.rds")
stage_age_data <- readRDS("data/compiled-stage-age-data.rds")
catch_binary <- readRDS("data/compiled-binary-recapture-data.rds")

# setup matrix settings
classes <- c(ncol(stage_age_data), nrow(stage_age_data))

# setup matrix indices to fill Leslie matrix
idx <- row(diag(classes[1])) == col(diag(classes[1])) + 1 |
  row(diag(classes[1])) == classes[1] & col(diag(classes[1])) == classes[1]
idy <- row(diag(classes[1])) == 1 & col(diag(classes[1])) > 4

# set shared priors for Leslie matrix and related parameters
transition_alpha <- normal(0, 2)
transition_beta <- normal(0, 2)
transition_mean <- transition_alpha + transition_beta * seq_len(sum(idx))
transition_sd <- normal(0, 0.5, dim = sum(idx), truncation = c(0, Inf))
fecundity_mean <- normal(0, 3, dim = sum(idy), truncation = c(0, Inf))
fecundity_sd <- normal(0, 3, dim = sum(idy), truncation = c(0, Inf))
inits_mean <- normal(0, 50, dim = c(classes[1], 1L), truncation = c(0, Inf))
inits_sd <- normal(0, 50, dim = c(classes[1], 1L), truncation = c(0, Inf))
detection_mean <- beta(10, 10)

# setup density dependence function
beverton <- function(parameters) {
  
  # tf function to calculate density dependence
  form <- function(x, params) {
    tf$ones(1, dtype = x$dtype) / tf$ones(1, dtype = x$dtype) +
      x * params
  }
  
  # return outputs
  list(form = form, parameters = parameters, name = "beverton")
  
}

# setup Leslie matrices for each site
site_names <- as.character(unique(pop_data$site))
mat <- inits <- vector("list", length = length(site_names))
for (i in seq_along(site_names)) {
  
  # construct population matrix
  mat1 <- mat2 <- zeros(classes[1], classes[1])
  
  # vectorise matrix fill: transition
  mat1[idx] <- ilogit(normal(transition_mean, transition_sd))

  # vectorise matrix fill: fecundity
  mat2[idy] <- normal(fecundity_mean, fecundity_sd, truncation = c(0, Inf))

  # initial conditions prior
  inits[[i]] <- normal(inits_mean, inits_sd, truncation = c(0, Inf))

  # compile matrix
  mat[[i]] <- mat1 + mat2
  
}

# density dependence prior
density_depend <- uniform(1e-5, 0.2, dim = length(site_names))

# set prior for stage-to-age conversion
stage_eps <- 500000 / rep(1, nrow(stage_age_data))
alpha_set <- matrix(1, classes[2], classes[1])
stage_lower <- seq(0, -10, length = classes[2])
stage_upper <- seq(10, 0, length = classes[2])
for (i in seq_len(classes[2]))
  alpha_set[i, ] <- stage_eps[i] * dnorm(seq(stage_lower[i], stage_upper[i], length = classes[1]))
stage_age <- dirichlet(alpha = alpha_set)

# add likelihoods for observed population abundances
mu <- obs <- vector("list", length = length(site_names))
for (i in seq_along(site_names)) {
  
  data_tmp <- stage_data[, grep(site_names[i], stage_data_info)]
  n_iter <- ncol(data_tmp)
  iterated_tmp <- iterate_matrix(matrix = mat[[i]],
                                 initial_state = inits[[i]],
                                 niter = n_iter,
                                 density = beverton(density_depend[i]))
  iterated_states <- iterated_tmp$all_states
  modelled_states <- stage_age %*% iterated_states
  mu[[i]] <- detection_mean * c(modelled_states)
  obs[[i]] <- as_data(c(data_tmp))
  distribution(obs[[i]]) <- poisson(mu[[i]])
  
}

# add likelihood for age-stage conversion
size_sums <- apply(stage_age_data, 1, sum)
distribution(stage_age_data) <- multinomial(size = size_sums, p = stage_age)

# extract summary info from recapture data
first_obs <- apply(catch_binary, 1, function(x) min(which(x > 0)))
final_obs <- apply(catch_binary, 1, function(x) max(which(x > 0)))
obs_id <- apply(catch_binary, 1, function(x) seq(min(which(x > 0)), max(which(x > 0)), by = 1)[-1])
obs_id <- unlist(obs_id)
capture_vec <- apply(catch_binary, 1, function(x) x[min(which(x > 0)):max(which(x > 0))][-1])
capture_vec <- unlist(capture_vec)
n_time <- ncol(catch_binary)

# set prior for survival in CJS model
phi <- beta(1, 1, dim = n_time)

# derived parameter for CJS model
chi <- ones(n_time)
for (i in seq_len(n_time - 1)) {
  tn <- n_time - i
  chi[tn] <- (1 - phi[tn]) + phi[tn] * (1 - detection_mean) * chi[tn + 1]
}

# dummy variables used to set CJS model likelihood
alive_data <- ones(length(obs_id))               # definitely alive
not_seen_last <- final_obs != ncol(catch_binary) # ignore observations in last timestep
final_observation <- ones(sum(not_seen_last))    # final observation

# set likelihoods for CJS model
distribution(alive_data) <- bernoulli(phi[obs_id - 1])
distribution(capture_vec) <- bernoulli(detection_mean)
distribution(final_observation) <- bernoulli(chi[final_obs[not_seen_last]])  

# compile greta model
mat_vec <- do.call(c, mat)
inits_vec <- do.call(c, inits)
mu_vec <- do.call(c, mu)
mod <- model(transition_alpha, transition_beta,
             fecundity_mean, fecundity_sd,
             density_depend,
             detection_mean,
             phi,
             inits_vec,
             inits_mean, inits_sd,
             stage_age,
             mat_vec, inits_vec)

# sample from greta model; use adam() optimiser to set initial values
n_chains <- 4
n_samples <- 100
n_warmup <- 1000
n_thin <- 10
init_set <- initials(transition_alpha = 0,
                     transition_beta = 0.2,
                     fecundity_mean = rep(1, length(fecundity_mean)),
                     fecundity_sd = rep(1, length(fecundity_sd)),
                     inits_mean = rep(10, length(inits_mean)),
                     inits_sd = rep(1, length(inits_sd)),
                     detection_mean = 0.5,
                     density_depend = rep(1e-3, length(site_names)))
opt_est <- init_list <- vector("list", length = n_chains)
for (i in seq_len(n_chains)) {
  opt_est[[i]] <- opt(mod, optimiser = adam(),
                 initial_values = init_set,
                 max_iterations = 1000,
                 tolerance = 1e-8)
  init_list[[i]] <- initials(transition_alpha = opt_est[[i]]$par$transition_alpha,
                             transition_beta = opt_est[[i]]$par$transition_beta,
                             fecundity_mean = opt_est[[i]]$par$fecundity_mean,
                             fecundity_sd = opt_est[[i]]$par$fecundity_sd,
                             inits_mean = opt_est[[i]]$par$inits_mean,
                             inits_sd = opt_est[[i]]$par$inits_sd,
                             density_depend = opt_est[[i]]$par$density_depend,
                             detection_mean = opt_est[[i]]$par$detection_mean)
}

# sample short chain from posterior
draws <- mcmc(mod,
              sampler = rwmh(),
              n_samples = n_samples, warmup = n_warmup,
              chains = n_chains, thin = n_thin,
              initial_values = init_list)

# sample longer chain from posterior (use extra_samples to avoid slowdown in tf samplers)
draws <- extra_samples(draws, 
                       n_samples = 2000000,
                       thin = n_thin)

# reduce the draws (remove burn-in samples)
draws_reduced <- lapply(draws, function(x) coda::mcmc(x, start = 190011, end = 200010, thin = 1))
draws_reduced <- coda::as.mcmc.list(draws_reduced)

# save outputs
saveRDS(draws_reduced, file = "outputs/mod-mcmc-draws.rds")
mod_summary_reduced <- summary(draws_reduced)
saveRDS(mod_summary_reduced, file = "outputs/fitted-model-reduced.rds")

# save some extra information for model plots and summaries
plot_info <- list(sites = site_names, classes = classes, stage_data_info = stage_data_info,
                  stage_data = stage_data)
saveRDS(plot_info, file = "outputs/plot-info.rds")
