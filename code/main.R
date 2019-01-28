# fit integrated model to data on individual growth, population size distributions,
#   and recapture histories

# load packages
library(integrated)
library(greta)
library(greta.dynamics)
library(future)

# model settings
breaks <- c(0, 200, 500, 1000, 2000, 5000, 10000, 20000, 60000)
classes <- length(breaks) - 1
replicates <- 1

# load population data
pop_data <- read.csv('./data/pop-data.csv',
                     row.names = 1, stringsAsFactors = FALSE)
pop_data <- pop_data[pop_data$species == 'murraycod', ]
pop_data <- pop_data[!is.na(pop_data$weight), ]
pop_data <- data.frame(time = pop_data$year,
                       site = pop_data$system,
                       size = pop_data$weight,
                       age = pop_data$age)

# load growth data
growth_data <- read.csv('./data/growth-data.csv')
growth_data <- data.frame(growth = growth_data$Otolith_growth,
                          id = growth_data$ID,
                          year = growth_data$Growth_year)

# load CMR data
cmr_data <- read.csv('./data/cmr-data.csv',
                     row.names = 1, stringsAsFactors = FALSE)
cmr_data <- cmr_data[cmr_data$code == 'MC', ]
cmr_data <- data.frame(id = cmr_data$idfish,
                       time = cmr_data$year,
                       size = cmr_data$weight)
cmr_data <- cmr_data[!is.na(cmr_data$size), ]
cmr_data <- cmr_data[cmr_data$id != 0, ]

# setup parallel run
plan(cluster)

# need somewhere to save fitted models
samples <- vector('list', length = length(unique(pop_data$site)))
test_data <- samples

# run through each river separately
for (i in seq_along(samples)) {
  
  # subset data to a single river
  pop_data_sub <- pop_data[pop_data$site == unique(pop_data$site)[i], ]
  
  pop_data_sub <- pop_data_sub[pop_data_sub$time < 2017, ]
  test_data[[i]] <- pop_data_sub[pop_data_sub$time >= 2017, ]
  
  # build process model
  mpm_process <- integrated_process(type = 'MPM',
                                    classes = classes,
                                    structure = 'stage',
                                    density_dependence = 'bh',
                                    replicates = replicates,
                                    params = list(fec_lower = 1, fec_upper = 1000,
                                                  surv_params1 = rep(1, classes),
                                                  surv_params2 = rep(1, classes),
                                                  density_lower = 0, density_upper = 0.5))
  
  # add population data
  abund_data_module <- integrated_data(data = pop_data_sub,
                                       integrated_process = mpm_process,
                                       process_link = 'stage_abundance',
                                       settings = list(breaks = breaks))
  
  # add cmr data
  cmr_data_module <- integrated_data(data = cmr_data,
                                     integrated_process = mpm_process,
                                     process_link = 'stage_recapture',
                                     settings = list(breaks = breaks))

  # add growth data
  size_data_module <- integrated_data(data = growth_data,
                                      integrated_process = mpm_process,
                                      process_link = 'individual_growth',
                                      settings = list(nbreaks = (classes + 1)))
  
  # combine data modules to create a joint likelihood object and return model parameters
  params <- integrated_model(integrated_process = mpm_process,
                             abund_data_module,
                             size_data_module,
                             cmr_data_module)
  
  # compile the model
  mod <- greta::model(params)
  
  # sample from the compiled model using futures::plan() defined above
  samples[[i]] <- greta::mcmc(mod, n_samples = 20000, warmup = 10000,
                              thin = 3, chains = 3)
  
}

# summarise model outputs for each river
mean_params <- vector('list', length = length(samples))
mu <- vector('list', length = length(samples))
mat_est <- vector('list', length = length(samples))
surv_vec <- vector('list', length = length(samples))
dens_param <- rep(NA, length(samples))
dens_all <- vector('list', length = length(samples))
for (i in seq_along(samples)) {
  
  mean_params_tmp <- matrix(NA, nrow = length(samples[[i]]), ncol = ncol(samples[[i]][[1]]))
  for (j in seq_along(samples[[i]])) {
    mean_params_tmp[j, ] <- apply(samples[[i]][[j]], 2, mean)
  }
  mean_params[[i]] <- apply(mean_params_tmp, 2, mean)
  
  surv_vec[[i]] <- mean_params[[i]][(2 * classes * classes + 1):(2 * classes * classes + classes)]
  
  mat_est[[i]] <- sweep(matrix(mean_params[[i]][1:(classes * classes)],
                               ncol = classes, byrow = FALSE),
                        2, surv_vec[[i]], '*') +
    matrix(mean_params[[i]][(classes * classes + (1:(classes * classes)))],
           ncol = classes, byrow = FALSE)
  
  dens_param[i] <- mean_params[[i]][((2 * classes + 3) * classes + 1)]
  
  dens_tmp <- matrix(NA, nrow = length(samples[[i]]), ncol = nrow(samples[[i]][[1]]))
  for (j in seq_along(samples[[i]])) {
    dens_tmp[j, ] <- samples[[i]][[j]][, ((2 * classes + 3) * classes + 1)]
  }
  dens_all[[i]] <- apply(dens_tmp, 2, mean)
  
  mu[[i]] <- mean_params[[i]][((2 * classes + 3) * classes + 2):(length(mean_params[[i]]))]

}

# plot the estimated transition matrices for each river
pdf(file = '.outputs/MS-Fig1.pdf', height = 7, width = 6)
system_lab <- c('Campaspe', 'Loddon', 'Broken', 'Goulburn',
                'Ovens', 'Murray')
size_vec <- breaks[c(1, 3, 5, 7, 9)] + round(diff(breaks)[c(1, 3, 5, 7, 9)] / 2)
layout(matrix(c(1, 2, 3, 4, 5, 6, 7, 7), nrow = 4, byrow = TRUE),
       heights = c(rep(1, 3), 0.2))
par(mar = c(3.8, 5.5, 3, 1))
for (i in seq_along(mat_est)) {
  
  image(t(log(round(mat_est[[i]], 2) + 0.001)),
        xaxt = 'n', yaxt = 'n',
        col = viridis::inferno(256), las = 1,
        zlim = log(range(sapply(mat_est, round, 2)) + 0.001), cex.axis = 1.25)
  axis(1, at = seq(0, 1, length = length(size_vec)), labels = size_vec)
  axis(2, at = seq(0, 1, length = length(size_vec)), labels = size_vec, las = 1)
  mtext('Size (g)', side = 1, adj = 0.5, line = 2.6, cex = 1)
  mtext('Size (g)', side = 2, adj = 0.5, line = 3.8, cex = 1)
  mtext(system_lab[i], side = 3, line = 0.7, adj = 1, cex = 1.2)
  
} 
par(mar = c(0.3, 5.5, 0, 1))
plot(c(0, 0.2) ~ c(0, 1), bty = 'n', type = 'n', xlab = '', ylab = '',
     xaxt = 'n', yaxt = 'n')
zrange <- log(range(sapply(mat_est, round, 2)) + 0.001)
oceanmap::set.colorbar(cbx = c(0, 1), cby = c(0.09, 0.14),
                       pal = viridis::inferno(256),
                       ticks = round(exp(seq(zrange[1], zrange[2],
                                         length = 6)), 2),
                       labels = c(0, 0.01, 0.1, 1, 10, 100),
                       cex.cb.ticks = 1.25,
                       cb.ticks.ypos = 0.04)
dev.off()

# projection helper function
bh_fun <- function(mat, n, dens_param, init = NULL) {
  
  out <- matrix(NA, nrow = nrow(mat), ncol = n)
  out[, 1] <- init
  for (i in seq_len(n)[-1]) {
    ntm1 <- sum(out[, (i - 1)])
    scale_fac <- 1 / (1 + dens_param * ntm1)
    out[, i] <- rpois(nrow(out), lambda = round((scale_fac * mat) %*% out[, (i - 1)]))
  }
  
  out
  
}

# calculate fitted and observed values and project the populations forward in time
nsim <- 100
ntime <- 10
out <- vector('list', length = length(mat_est))
real_abund <- vector('list', length = length(mat_est))
fitted_abund <- vector('list', length = length(mat_est))
mu_all <- vector('list', length = length(mat_est))
bayes_r2 <- matrix(NA, nrow = length(mat_est), ncol = nsim)
for (i in seq_along(mat_est)) {
  pop_data_sub <- pop_data[pop_data$site == unique(pop_data$site)[i], ]
  abund_data_tmp <- integrated_data(data = pop_data_sub,
                                    integrated_process = mpm_process,
                                    process_link = 'stage_abundance',
                                    settings = list(breaks = breaks))
  real_abund[[i]] <- abund_data_tmp$data[[1]]
  fitted_abund[[i]] <- apply(matrix(mu[[i]], nrow = nrow(real_abund[[i]])),
                             2, sum)
  
  out[[i]] <- vector('list', length = nsim)
  mu_all[[i]] <- vector('list', length = nsim)
  for (j in seq_len(nsim)) {
    
    params_tmp <- samples[[i]][[1]][sample(seq_len(nrow(samples[[i]][[1]])),
                                           size = 1), ]
    
    surv_tmp <- params_tmp[(2 * classes * classes + 1):(2 * classes * classes + classes)]
    
    mat_tmp <- sweep(matrix(params_tmp[1:(classes * classes)],
                                 ncol = classes, byrow = FALSE),
                          2, surv_tmp, '*') +
      matrix(params_tmp[(classes * classes + (1:(classes * classes)))],
             ncol = classes, byrow = FALSE)
    
    dens_tmp <- params_tmp[((2 * classes + 3) * classes + 1)]
    mu_all[[i]][[j]] <- params_tmp[((2 * classes + 3) * classes + 2):(length(params_tmp))]
    mu_mat_tmp <- matrix(mu[[i]], nrow = nrow(real_abund[[i]]))
    
    out[[i]][[j]] <- bh_fun(mat_tmp, n = ntime, dens_param = (dens_tmp),
                            init = mu_mat_tmp[, ncol(mu_mat_tmp)])
    
    bayes_r2[i, j] <- var(mu_all[[i]][[j]]) /
      (var(mu_all[[i]][[j]]) + var(mu_all[[i]][[j]] - c(real_abund[[i]])))

  }
}

# save Bayesian R2 values
write.csv(bayes_r2, file = './outputs/bayes_r2_values.csv')

# plot the projected population trajectories
pdf(file = './outputs/MS-Fig2.pdf', height = 7, width = 7)
par(mfrow = c(3, 2), mar = c(4.1, 5.1, 2.8, 1.5))
for (i in seq_along(mat_est)) {
  
  pop_data_sub <- pop_data[pop_data$site == unique(pop_data$site)[i], ]
  x_vals <- c(unique(pop_data_sub$time), (max(pop_data_sub$time) + 1):(max(pop_data_sub$time) + 10))
  sim_vals <- sapply(out[[i]], function(x) apply(x, 2, sum))
  plot_mean <- c(apply(real_abund[[i]], 2, sum),
                 apply(sim_vals, 1, mean))
  sim_mean <- apply(do.call('rbind', mu_all[[i]]),
                    2, mean)
  sim_upper <- apply(do.call('rbind', mu_all[[i]]),
                    2, quantile, p = 0.9)
  sim_lower <- apply(do.call('rbind', mu_all[[i]]),
                    2, quantile, p = 0.1)
  sim_mean <- apply(matrix(sim_mean, nrow = nrow(real_abund[[i]])),
                    2, sum)
  sim_upper <- apply(matrix(sim_upper, nrow = nrow(real_abund[[i]])),
                    2, sum)
  sim_lower <- apply(matrix(sim_lower, nrow = nrow(real_abund[[i]])),
                    2, sum)
  
  plot(plot_mean[-length(x_vals)] ~ x_vals[-length(x_vals)],
       type = "n", bty = "l", las = 1,
       xaxt = "n", yaxt = "n",
       ylim = range(c(plot_mean[-length(x_vals)], c(sim_mean, sim_upper, sim_lower),
                    c(sim_vals))),
       xlim = range(x_vals),
       xlab = "", ylab = "")
  
  # set colours
  col_pal_main <- viridis::inferno(256, alpha = 1)[c(20, seq(70, 200,
                                                             length = nsim))]
  col_pal_sub <- viridis::inferno(256, alpha = 0.4)[c(20, seq(70, 200,
                                                              length = nsim))]
  polygon(c(x_vals[1:ncol(real_abund[[i]])], rev(x_vals[1:ncol(real_abund[[i]])])),
          c(sim_upper, rev(sim_lower)),
          border = NA, col = ggplot2::alpha(col_pal_main[1], 0.5))
  lines(plot_mean[-length(x_vals)] ~ x_vals[-length(x_vals)],
        lwd = 1.2, col = col_pal_main[1])
  lines(sim_mean ~ x_vals[1:ncol(real_abund[[i]])],
        lwd = 1.2, col = col_pal_main[1])
  for (j in seq_len(nsim)) {
    lines(c(sim_mean[length(sim_mean)], sim_vals[, j]) ~ x_vals[(ncol(real_abund[[i]])):(length(x_vals))],
          lwd = 1.2, col = col_pal_main[j])
  }
  
  axis(1, at = seq(min(x_vals), max(x_vals), by = 3), las = 1)
  mtext("Year", side = 1, adj = 0.5, line = 2.5)
  axis(2, las = 1)
  mtext("Abundance", side = 2, adj = 0.5, line = 3.5)
  mtext(system_lab[i], side = 3, line = 0.7, adj = 1, cex = 1.2)

}
dev.off()

# plot the fitted and observed values
pdf(file = './outputs/MS-FigS1.pdf', height = 7, width = 7)
par(mfrow = c(3, 2), mar = c(4.1, 5.1, 2.8, 1.5))
for (i in seq_along(real_abund)) {
  real_tmp <- apply(real_abund[[i]], 2, sum)
  fitted_mean <- apply(matrix(apply(do.call('rbind', mu_all[[i]]),
                                    2, mean), nrow = nrow(real_abund[[i]])),
                       2, sum)
  fitted_lower <- apply(matrix(apply(do.call('rbind', mu_all[[i]]),
                                     2, quantile, p = 0.1), nrow = nrow(real_abund[[i]])),
                        2, sum)
  fitted_upper <- apply(matrix(apply(do.call('rbind', mu_all[[i]]),
                                     2, quantile, p = 0.9), nrow = nrow(real_abund[[i]])),
                        2, sum)
  pop_data_sub <- pop_data[pop_data$site == unique(pop_data$site)[i], ]
  x_tmp <- unique(pop_data_sub$time)
  plot(real_tmp ~ x_tmp,
       type = "n", bty = "l", las = 1,
       xaxt = "n", yaxt = "n",
       ylim = range(c(real_tmp, fitted_mean, fitted_lower, fitted_upper)),
       xlab = "", ylab = "")
  polygon(c(x_tmp, rev(x_tmp)), c(fitted_upper, rev(fitted_lower)),
          border = NA, col = col_pal_sub[90])
  lines(real_tmp ~ x_tmp, lwd = 2, col = col_pal_main[1])
  lines(fitted_mean ~ x_tmp, lwd = 2, col = col_pal_main[90])
  
  axis(1, at = x_tmp[seq(1, length(x_tmp), by = 1)], las = 1)
  mtext("Year", side = 1, adj = 0.5, line = 2.5)
  axis(2, las = 1)
  mtext("Abundance", side = 2, adj = 0.5, line = 3.5)
  mtext(system_lab[i], side = 3, line = 0.7, adj = 1, cex = 1.2)
  
}
dev.off()

# plot the effects of density dependence in each river
pdf(file = './outputs/MS-FigS2.pdf', width = 7, height = 7)
par(mfrow = c(3, 2), mar = c(5.1, 4.8, 2.8, 0.5))
x <- seq(0, 2000, by = 1)
plot_vals <- vector('list', length = length(dens_all))
for (i in seq_along(dens_param)) {
  dens_sum <- c(mean(dens_all[[i]]), quantile(dens_all[[i]], p = c(0.1, 0.5, 0.9)))
  plot_vals[[i]] <- matrix(NA, nrow = length(dens_sum), ncol = length(x))
  for (j in seq_along(dens_sum)) {
    plot_vals[[i]][j, ] <- 1 / (1 + dens_sum[j] * x)
  }
}
for (i in seq_along(dens_all)) {
  plot(plot_vals[[i]][1, ] ~ x, 
       type = "n", bty = "l", las = 1,
       ylim = range(unlist(plot_vals)),
       ylab = "Density parameter", xlab = "Abundance")
  polygon(c(x, rev(x)),
          c(plot_vals[[i]][2, ], rev(plot_vals[[i]][4, ])),
          border = NA, col = col_pal_sub[1])
  lines(plot_vals[[i]][1, ] ~ x,
        lwd = 2, col = col_pal_main[1])
  mtext(system_lab[i], side = 3, line = 0.7, adj = 1, cex = 1.2)
}
dev.off()
