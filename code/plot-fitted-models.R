# optional: set working directory
# setwd("PATH/TO/DIR")

# load some helpers
source("code/plot-helpers.R")

# load fitted model
mod_fitted <- readRDS("outputs/fitted-model-reduced.rds")
plot_info <- readRDS("outputs/plot-info.rds")
draws_red <- readRDS("outputs/mod-mcmc-draws.rds")

# pull out transitions and fecundity from summarised values
mat_est <- mod_fitted$quantiles[grep("mat_vec", rownames(mod_fitted$quantiles)), ]
mat_plot <- vector("list", length = ncol(mat_est))
idx <- row(diag(plot_info$classes[1])) == col(diag(plot_info$classes[1])) + 1 |
  row(diag(plot_info$classes[1])) == plot_info$classes[1] &
  col(diag(plot_info$classes[1])) == plot_info$classes[1]
idy <- row(diag(plot_info$classes[1])) == 1 &
  col(diag(plot_info$classes[1])) > 4
for (i in seq_along(mat_plot)) {
  mat_tmp <- array(mat_est[, i], dim = c(plot_info$classes[1],
                                         plot_info$classes[1],
                                         length(plot_info$sites)))
  mat_plot[[i]] <- apply(mat_tmp, 3, function(x) c(x[idx], x[idy]))
}

# what are the values in Yen et al. 2013?
lit_mean <- c(0.45, 0.60, 0.70, 0.75, 0.80,
              8500 * 0.5 * 0.0012)
lit_sd <- c(0.094, 0.12, 0.14, 0.11, 0.12,
            3000 * 0.075 * 0.001)

# plot settings for labels
system_names <- c("Literature",
                  rev(c("Broken", "Campaspe", "Goulburn", "Ovens", "Murray", "Loddon")))
parameter_lab <- c("Age 1 survival", "Age 2 survival", "Age 3 survival", 
                   "Age 4 survival", "Age 5 survival", "Fecundity")
system_lab <- c("Broken", "Campaspe", "Goulburn", "Ovens", "Murray", "Loddon")

# plot Fig. 2 from MS: fitted vital rates
pdf(file = "outputs/Fig2.pdf", height = 7, width = 7)
par(mfrow = c(3, 2), mar = c(4.1, 5.3, 2.8, 1.5))
xlow <- rep(0, 6)
xhigh <- c(rep(0, 5), 10)
for (i in seq_len(nrow(mat_plot[[1]]))) {
  
  # pull out values for all sites for parameter i
  plot_med <- mat_plot[[3]][i, ]
  plot_vl <- mat_plot[[2]][i, ]
  plot_vh <- mat_plot[[4]][i, ]
  
  # pull out literature value
  lit_val_med <- lit_mean[i]
  lit_val_lo <- lit_mean[i] - 1.96 * lit_sd[i]
  lit_val_hi <- lit_mean[i] + 1.96 * lit_sd[i]
  
  # plot the values
  plot(rev(c(plot_med, lit_val_med)),
       seq_len(length(plot_med) + 1),
       type = "n",
       las = 1,
       bty = "l",
       xlab = "", ylab = "",
       yaxt = "n", xaxt = "n",
       xlim = range(c(0, 1,
                      plot_med, plot_vl, plot_vh,
                      lit_val_med, lit_val_lo, lit_val_hi)),
       ylim = c(0.5, 7.5))
  for (j in seq_along(plot_med)) {
    lines(c(plot_vl[j], plot_vh[j]), c(7 - j + 1, 7 - j + 1), lwd = 3)
    points(plot_med[j], 7 - j + 1, pch = 16, cex = 1.25)
  }
  lines(c(lit_val_lo, lit_val_hi), c(1, 1), lwd = 3, col = "gray50")
  points(lit_val_med, 1, pch = 16, col = "gray50", cex = 1.25)
  
  # add some axis labels
  axis(1, las = 1)
  mtext("Parameter estimate", side = 1, adj = 0.5, line = 2.5, cex = 0.9)
  axis(2, at = seq_len(length(plot_med) + 1),
       labels = system_names, las = 1, cex = 0.9)
  mtext(parameter_lab[i], side = 3, line = 0.7, adj = 1, cex = 1.0)
  
}
dev.off()

# set a colour palette for subsequent plots
col_pal_main <- viridis::inferno(256, alpha = 1)[c(20, seq(70, 200,
                                                           length = 100))]
col_pal_sub <- viridis::inferno(256, alpha = 0.4)[c(20, seq(70, 200,
                                                            length = 100))]

# plot the fitted and observed values
pdf(file = "outputs/MS-FigS1.pdf", height = 7, width = 7)
par(mfrow = c(3, 2), mar = c(4.1, 5.1, 2.8, 1.5))
for (i in seq_along(plot_info$sites)) {
  
  real_abunds <- plot_info$stage_data[, grep(plot_info$sites[i], plot_info$stage_data_info)]
  real_tmp <- apply(real_abunds, 2, sum)
  
  fitted <- sample_draws(draws_red, site = i,
                         abund_data = real_tmp,
                         chain = 2,
                         nsim = 1000,
                         p = c(0.1, 0.5, 0.9))
  fitted_median <- fitted[2, ]
  fitted_lower <- fitted[1, ]
  fitted_upper <- fitted[3, ]
  
  x_tmp <- plot_info$stage_data_info[grep(plot_info$sites[i], plot_info$stage_data_info)]
  x_tmp <- as.numeric(sapply(strsplit(x_tmp, "_"), function(x) x[2]))
  plot(real_tmp ~ x_tmp,
       type = "n", bty = "l", las = 1,
       xaxt = "n", yaxt = "n",
       ylim = range(c(real_tmp, fitted_median, fitted_lower, fitted_upper)),
       xlab = "", ylab = "")
  polygon(c(x_tmp, rev(x_tmp)), c(fitted_upper, rev(fitted_lower)),
          border = NA, col = col_pal_sub[90])
  lines(real_tmp ~ x_tmp, lwd = 2, col = col_pal_main[1])
  lines(fitted_median ~ x_tmp, lwd = 2, col = col_pal_main[90])
  
  axis(1, at = x_tmp[seq(1, length(x_tmp), by = 1)], las = 1)
  mtext("Year", side = 1, adj = 0.5, line = 2.5)
  axis(2, las = 1)
  mtext("Abundance", side = 2, adj = 0.5, line = 3.5)
  mtext(system_lab[i], side = 3, line = 0.7, adj = 1, cex = 1.2)
  
}
dev.off()

# plot the effects of density dependence in each river
pdf(file = "outputs/MS-FigS2.pdf", width = 7, height = 7)
par(mfrow = c(3, 2), mar = c(5.1, 4.8, 2.8, 0.5))
dens_param <- draws[[2]][, grep("density_depend", colnames(draws_red[[2]]))]
dens_param <- t(apply(dens_param, 2, quantile, p = c(0.1, 0.5, 0.9)))
x <- seq(0, 2000, by = 1)
plot_vals <- vector("list", length = nrow(dens_param))
for (i in seq_len(nrow(dens_param))) {
  dens_sum <- dens_param[i, c(2, 1, 3)]
  plot_vals[[i]] <- matrix(NA, nrow = length(dens_sum), ncol = length(x))
  for (j in seq_along(dens_sum)) {
    plot_vals[[i]][j, ] <- 1 / (1 + dens_sum[j] * x)
  }
}
for (i in seq_len(nrow(dens_param))) {
  plot(plot_vals[[i]][1, ] ~ x, 
       type = "n", bty = "l", las = 1,
       ylim = range(unlist(plot_vals)),
       ylab = "Density parameter", xlab = "Abundance")
  polygon(c(x, rev(x)),
          c(plot_vals[[i]][2, ], rev(plot_vals[[i]][3, ])),
          border = NA, col = col_pal_sub[1])
  lines(plot_vals[[i]][1, ] ~ x,
        lwd = 2, col = col_pal_main[1])
  mtext(system_lab[i], side = 3, line = 0.7, adj = 1, cex = 1.2)
}
dev.off()

# calculate fitted and observed values and project the populations forward in time
bayes_r2 <- matrix(NA, ncol = length(plot_info$sites), nrow = 1000)
for (i in seq_along(plot_info$sites)) {

  real_abunds <- plot_info$stage_data[, grep(plot_info$sites[i], plot_info$stage_data_info)]
  real_tmp <- apply(real_abunds, 2, sum)
  fitted <- sample_draws(draws_red, site = i,
                         abund_data = real_tmp,
                         chain = 2,
                         p = c(0.1, 0.5, 0.9),
                         simplify = FALSE)

  for (j in seq_len(nrow(fitted))) {
    bayes_r2[j, i] <- var(fitted[j, ]) /
      (var(fitted[j, ]) + var(fitted[j, ] - c(real_abunds)))
  }
  
}
r2_out <- apply(bayes_r2, 2, quantile, p = c(0.1, 0.5, 0.9))
