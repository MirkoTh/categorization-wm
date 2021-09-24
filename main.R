library(tidyverse)
library(docstring)
library(grid)
library(gridExtra)
library(wiqid)

dirs_local <- c("R/utils.R", "R/rmc-thm.R")
walk(dirs_local, source)


# Settings ----------------------------------------------------------------

stepsize <- 1
max_x <- 2
x1 <- seq(1, max_x, by = stepsize)
x2 <- seq(1, max_x, by = stepsize)
x3 <- seq(1, max_x, by = stepsize)
tbl_features <- crossing(x1, x2, x3)


# Create Categories -------------------------------------------------------

# n_cat_per_feat <- c(3)
# stepsize_cat <- max_x / n_cat_per_feat
# l_tbl_categories <- map(n_cat_per_feat, create_categories, tbl = tbl_features)
# pl <- map2(l_tbl_categories, stepsize_cat, plot_clustered_grid)
# plot_arrangement(pl)
tbl_stimuli <- create_shepard_categories(tbl_features, "I", "x2")
l_tbl_categories <- list()
l_tbl_categories[[1]] <- tbl_stimuli

# Fit RMC to simulated Data -----------------------------------------------

idx_used <- 1
idxs <- 1:nrow(l_tbl_categories[[idx_used]])
n <- 96
# kind of badham, sanborn, & maylor
idx_shuffle <- sample(rep(rep(idxs, 2), 6), 96, replace = FALSE)
# actual goal of 2d feature space
#idx_shuffle <- sample(idxs, n, replace = TRUE)
tbl_used <- l_tbl_categories[[idx_used]][idx_shuffle, ]
params_init <- c(.31, 1.3)

optim(
  params_init, wrap_rmc, tbl_data = tbl_used,
  method = "L-BFGS-B",
  lower = c(.01, .01), upper = c(.6, 2),
  control = list(trace = 1)
)

library(future)
library(furrr)
coupling <- seq(.2, .35, by = .025)
phi <- seq(1, 2.5, by = .25)
tbl_params <- crossing(coupling, phi)

plan(multisession = parallel::detectCores())

l_out <- future_pmap(
  tbl_params,
  wrap_rmc,
  tbl_data = tbl_used
)

tbl_params$neg_ll <- map_dbl(l_out, 1)
tbl_params$n_cluster <- map_int(l_out, 2)

ggplot(tbl_params, aes(coupling, phi, fill = neg_ll)) +
  geom_raster(aes(fill = neg_ll)) +
  geom_text(aes(label = n_cluster), color = "white") +
  theme_bw() +
  scale_fill_viridis_c() +
  labs(
    title = "Neg. LL & Nr. Clusters",
    x = "Coupling Probability",
    y = "Phi"
  )



# To Play Around ----------------------------------------------------------


## With Feedback ----------------------------------------------------------

stimuli <- tbl_used[, c("x1", "x2", "x3")]
features_cat <- c("x1", "x2", "x3") 
features_cont <- c() 
n_values_cat <- 2
feedback <- NULL
print_posterior <- FALSE
previous_learning <- NULL
feedback <- tbl_used$category
params <- list(
  "coupling" = .5044,
  "phi" = .7738, #1,
  "salience_f" = .6888,#1,
  "salience_l" = .1615,#1,
  "a_0" = 2,
  "lambda_0" = 1,
  "sigma_sq_0" = (max(tbl_used[, c("x1", "x2")]) / 4) ^ 2,
  "mu_0" = (min(tbl_used[, c("x1", "x2")]) + max(tbl_used[, c("x1", "x2")])) / 2
)
n_categories <- length(unique(tbl_used$category))

l_pred <- predict_rmc_continuous(
  stimuli = tbl_used[, c("x1", "x2", "x3")], 
  features_cat =c("x1", "x2", "x3"),
  features_cont = c(),
  n_values_cat = 2,
  n_categories = n_categories,
  feedback = tbl_used$category,
  params = params,
  previous_learning = NULL, 
  print_posterior = FALSE
)

tbl_used$preds <- apply(
  l_pred$cat_probs, 1, FUN = function(x) which.max(x)
)
tbl_used$n_training <- n
tbl_used$n_categories <- length(unique(tbl_used$category))

l_summary <- summarize_blocks(tbl_used, 16)
plot_block_summary(l_summary)


## Without Feedback -------------------------------------------------------

previous_learning <- list(
  "stimuli" = l_pred[["stimuli"]],
  "feature_counts" = l_pred[["feature_counts"]],
  "cluster_counts" = l_pred[["cluster_counts"]],
  "assignments" = l_pred[["assignments"]]
)

l_pred_trained <- predict_rmc_continuous(
  stimuli = stimuli, 
  features_cat = features_cat,
  features_cont = features_cont,
  n_values_cat = n_values_cat,
  n_categories = n_categories,
  feedback = feedback,
  params = params,
  previous_learning = previous_learning, 
  print_posterior = print_posterior
)



# Explore Chi Square Distribution -----------------------------------------

pdf_chisq <- function(x, k) {
  (1 / (2 ^ (k / 2) * gamma(k / 2))) * ((x ^ (k / 2 - 1)) * exp(-x / 2))
}
pdf_chisq_inv <- function(x, v) {
  ((2 ^ -(v / 2)) / (gamma(v / 2))) * x ^ (-(v / 2) - 1) * exp(-1 / (2 * x))
}

to_tbl <- function(v, x) {
  tbl <- as_tibble(v)
  tbl$x <- x
  return(tbl)
}

plot_density <- function(tbl, y_lo, y_hi) {
  ggplot(tbl, aes(x, value)) +
    geom_line() +
    geom_point(color = "white", size = 2) +
    geom_point() +
    coord_cartesian(ylim = c(y_lo, y_hi)) +
    theme_bw()
}

plot_hist <- function(x, x_min, x_max) {
  ggplot(x %>% as_tibble(), aes(value)) + 
    geom_histogram() +
    coord_cartesian(xlim = c(x_min, x_max)) +
    theme_bw()
}

sample_norm <- function(lambda) {
  rnorm(1000, 0, 10/sqrt(lambda))
}

k <- seq(1, 9, by = 1)
x <- seq(0, 10, by = .01)
l_chisq <- map(k, pdf_chisq, x = x)
l_chisq_inv <- map(k, pdf_chisq_inv, x = x[x <= 1])
l_chisq <- pmap(list(l_chisq, list(x)), to_tbl)
l_chisq_inv <- pmap(list(l_chisq_inv, list(x[x <= 1])), to_tbl)

l_pl_chisq <- map(l_chisq, plot_density, y_lo = 0, y_hi = .5)
plot_arrangement(l_pl_chisq, n_cols = 3)

l_pl_chisq_inv <- map(l_chisq_inv, plot_density, y_lo = 0, y_hi = 5)
plot_arrangement(l_pl_chisq_inv, n_cols = 3)


l_norm <- map(k, sample_norm)

l_pl_norm <- map(l_norm, plot_hist, x_min = -30, x_max = 30)
plot_arrangement(l_pl_norm)


ggplot(
  dt(seq(-5, 100, by = .1), 10, 100) %>%
    as_tibble() %>%
    mutate(idx = seq(-5, 100, by = .1)),
  aes(idx, value)
) + 
  geom_line() +
  coord_cartesian(xlim = c(-100, 200))


pdf_dt_gen <- function(x, v, sd) {
  d <- dt2(x, 50, sd, v)
  d %>% as_tibble() %>% mutate(x = x, scale = sd, v = v)
  
}
x <- seq(0, 100, by = .5)
vs <- seq(3, 10, by = 1)
sds <- seq(1, 4, by = 1)
tbl_v_sd <- crossing(vs, sds) %>% rename(v = vs, sd = sds)
l_dt <- pmap(tbl_v_sd, pdf_dt_gen, x = x)

l_pl_dt <- map(l_dt, plot_density, y_lo = 0, y_hi = .2)
plot_arrangement(l_pl_dt, n_cols = 8)

