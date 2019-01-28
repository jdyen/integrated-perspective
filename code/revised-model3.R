# clear workspace
rm(list = ls())

# set working directory
setwd("~/Dropbox/research/integrated-perspective/")

# load packages
library(greta.dynamics)
library(lubridate)

# this needs a few tensorflow fns
tf <- tensorflow::tf

# source helper files
source("../catch-curves/code/length_to_mass_calculations.R")
mc_coefs <- length_weight_conversion$murraycod$coef
rm(length_weight_conversion)

# model settings
breaks <- c(0, 200, 500, 1000, 2000, 5000, 10000, 20000, 60000)

# load population data
pop_data <- readRDS("data/revised-pop-data.rds")
pop_data <- pop_data[pop_data$species == "murraycod", ]
pop_data <- pop_data[!is.na(pop_data$weight), ]
flow_data <- pop_data[, c(12:ncol(pop_data))]
pop_data <- data.frame(time = pop_data$year,
                       site = pop_data$system,
                       size = pop_data$weight)

# clean up population data
hist_fn <- function(x, breaks) {
  hist(x, breaks = breaks, plot = FALSE)$counts
}
stage_data <- tapply(pop_data$size, list(pop_data$site, pop_data$time), hist_fn, breaks = breaks)
stage_data_info <- paste(rep(rownames(stage_data), each = ncol(stage_data)),
                         rep(colnames(stage_data), times = nrow(stage_data)),
                         sep = "_")
stage_data_info <- stage_data_info[sapply(c(t(stage_data)), function(x) !is.null(x))]
stage_data <- do.call(cbind, t(stage_data))

# load growth data
growth_data <- read.csv("data/growth-data.csv")
growth_data <- data.frame(growth = growth_data$Otolith_growth,
                          id = growth_data$ID,
                          year = growth_data$Growth_year,
                          length = growth_data$Length,
                          age = growth_data$Age)
otolith_conversion <- data.frame(otolith_size = tapply(growth_data$growth, growth_data$id, sum),
                                 length = tapply(growth_data$length, growth_data$id, unique),
                                 age = tapply(growth_data$age, growth_data$id, max))
otolith_conversion$weight <- exp(mc_coefs[1] + mc_coefs[2] * log(otolith_conversion$length))
otolith_conversion$stage <- cut(otolith_conversion$weight, breaks, labels = FALSE)
otolith_conversion$age[otolith_conversion$age > 5] <- 5
stage_age_data <- matrix(0, nrow = length(breaks) - 1, ncol = max(otolith_conversion$age))
rownames(stage_age_data) <- seq_len(length(breaks) - 1)
colnames(stage_age_data) <- seq_len(max(otolith_conversion$age))
stage_age_tmp <- table(otolith_conversion$stage, otolith_conversion$age)
stage_age_data[rownames(stage_age_tmp), colnames(stage_age_tmp)] <-
  stage_age_data[rownames(stage_age_tmp), colnames(stage_age_tmp)] + 
  stage_age_tmp

# load CMR data
cmr_data <- read.csv("data/cmr-data.csv",
                     row.names = 1, stringsAsFactors = FALSE)
cmr_data <- cmr_data[cmr_data$code == "MC", ]
cmr_data <- data.frame(id = cmr_data$idfish,
                       time = cmr_data$year,
                       size = cmr_data$weight)
cmr_data <- cmr_data[!is.na(cmr_data$size), ]
cmr_data <- cmr_data[cmr_data$id != 0, ]

# clean CMR data
catch_size <- tapply(cmr_data$size, list(cmr_data$id, cmr_data$time), mean)
catch_size <- ifelse(is.na(catch_size), 0, catch_size)
observed <- apply(catch_size, 1, sum) > 0
catch_size <- catch_size[observed, ]

# setup matrix settings
classes <- c(ncol(stage_age_data), nrow(stage_age_data))

# setup matrix indices to fill
idx <- row(diag(classes[1])) == col(diag(classes[1])) + 1 |
  row(diag(classes[1])) == classes[1] & col(diag(classes[1])) == classes[1]
idy <- row(diag(classes[1])) == 1 & col(diag(classes[1])) > 4

# overall parameter priors
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

# setup matrices for each site
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

# stage-to-age conversion
stage_eps <- 500000 / rep(1, nrow(stage_age_data)) # / (apply(stage_age_data, 1, sum) + 1)
alpha_set <- matrix(1, classes[2], classes[1])
stage_lower <- seq(0, -10, length = classes[2])
stage_upper <- seq(10, 0, length = classes[2])
for (i in seq_len(classes[2]))
  alpha_set[i, ] <- stage_eps[i] * dnorm(seq(stage_lower[i], stage_upper[i], length = classes[1]))
stage_age <- dirichlet(alpha = alpha_set)

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

# add age-stage conversion
size_sums <- apply(stage_age_data, 1, sum)
distribution(stage_age_data) <- multinomial(size = size_sums, p = stage_age)

# recapture data
catch_size <- catch_size[apply(catch_size, 1, function(x) sum(x > 0)) > 1, ]
catch_binary <- ifelse(catch_size > 0, 1, 0)

# extract summary info from data
first_obs <- apply(catch_binary, 1, function(x) min(which(x > 0)))
final_obs <- apply(catch_binary, 1, function(x) max(which(x > 0)))
obs_id <- apply(catch_binary, 1, function(x) seq(min(which(x > 0)), max(which(x > 0)), by = 1)[-1])
obs_id <- unlist(obs_id)
capture_vec <- apply(catch_binary, 1, function(x) x[min(which(x > 0)):max(which(x > 0))][-1])
capture_vec <- unlist(capture_vec)
n_time <- ncol(catch_binary)

# priors
phi <- beta(1, 1, dim = n_time)

# derived parameter
chi <- ones(n_time)
for (i in seq_len(n_time - 1)) {
  tn <- n_time - i
  chi[tn] <- (1 - phi[tn]) + phi[tn] * (1 - detection_mean) * chi[tn + 1]
}

# dummy variables
alive_data <- ones(length(obs_id))            # definitely alive
not_seen_last <- final_obs != ncol(catch_binary) # ignore observations in last timestep
final_observation <- ones(sum(not_seen_last)) # final observation

# set likelihoods
distribution(alive_data) <- bernoulli(phi[obs_id - 1])
distribution(capture_vec) <- bernoulli(detection_mean)
distribution(final_observation) <- bernoulli(chi[final_obs[not_seen_last]])  

# compile greta model and sample
mat_vec <- do.call(c, mat)
inits_vec <- do.call(c, inits)
mu_vec <- do.call(c, mu)
mod <- model(transition_alpha, transition_beta,
             fecundity_mean, fecundity_sd,
             density_depend,
             detection_mean,
             phi,
             inits_mean, inits_sd,
             stage_age,
             mat_vec, inits_vec)

# sample from greta model
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

# sample from posterior
draws <- mcmc(mod,
              sampler = rwmh(),
              n_samples = n_samples, warmup = n_warmup,
              chains = n_chains, thin = n_thin,
              initial_values = init_list)

# sample from posterior
draws <- extra_samples(draws, 
                       n_samples = 2000000,
                       thin = n_thin)

# summarise fitted
mod_summary <- summary(draws, quantiles = c(0.025, 0.1, 0.5, 0.9, 0.975))

# save fitted model
saveRDS(mod_summary, file = "outputs/fitted-model.rds")

# save some extra information for model summaries
plot_info <- list(sites = site_names, classes = classes, stage_data_info = stage_data_info,
                  stage_data = stage_data)
saveRDS(plot_info, file = "outputs/plot-info.rds")

# maybe save draws and plot from those?
