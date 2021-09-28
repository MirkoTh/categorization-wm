library(tidyverse)
library(docstring)
library(grid)
library(gridExtra)
library(wiqid)
library(future)
library(furrr)

dirs_local <- c("R/utils.R", "R/rmc-thm.R")
walk(dirs_local, source)


# Settings ----------------------------------------------------------------

stepsize <- 1
max_x <- 8
x1 <- seq(1, max_x, by = stepsize)
x2 <- seq(1, max_x, by = stepsize)
tbl_features <- crossing(x1, x2)


# Create Categories -------------------------------------------------------

n_cat_per_feat <- c(2, 4, 8)
stepsize_cat <- max_x / n_cat_per_feat
l_tbl_categories <- map(n_cat_per_feat, create_categories, tbl = tbl_features)
pl <- map2(l_tbl_categories, stepsize_cat, plot_clustered_grid)
plot_arrangement(pl)

# Fit RMC to simulated Data -----------------------------------------------

## for testing purposes
idx_used <- 1
ns <- c(50, 200)

n <- ns[[1]]
idxs <- 1:nrow(l_tbl_categories[[idx_used]])
idxs_shuffle <- map(rep(ns, each = length(l_tbl_categories)), sample, replace = TRUE, x = idxs)
l_tbl_used <- map2(idxs_shuffle, rep(l_tbl_categories, length(ns)), ~ .y[.x, ])
tbl_used <- l_tbl_used[[idx_used]]
n_categories <- rep(n_cat_per_feat^2, length(ns))


# order: coupling, phi, salience_l, a_0, lambda_0
# optimal for 100 trials and 9 categories: 0.24191989 2.62343018 0.01000000 0.04399989 0.01000000
params_init <- c(.1, 3, .01, .02, .01)
bounds_lower <- c(.01, .01, .00001, .00001, .00001)
bounds_upper <- c(.5, 10, 3, 1, 1)

plan(multisession = min((parallel::detectCores() - 1), length(l_tbl_used)))


wrap_optim <- function(tbl, n_categories, ...) {
  r <- optim(tbl_data = tbl, n_categories = n_categories, ...)
  beepr::beep()
  return(r)
}

# for comparing different training lengths
l_results <- future_map2(
  l_tbl_used,
  n_categories,
  safely(wrap_optim),
  par = params_init,
  fn = wrap_rmc,
  method = "L-BFGS-B",
  lower = bounds_lower,
  upper = bounds_upper,
  control = list(factr = .001)#.001
)


saveRDS(
  l_results, 
  file = str_c(
    "data/", lubridate::today(), "-optimization_results.Rda"
  )
)

ls_pred <- future_pmap(list(l_tbl_used, l_results, n_categories), predict_given_fit)
l_blocks <- pmap(list(ls_pred, l_tbl_used, n_categories), summarize_cat_probs, n_trials = 25)
tbl <- l_blocks %>%
  reduce(rbind) %>%
  mutate(length_training = as.factor(length_training))

ggplot(tbl, aes(block_nr, probability_mn, group = length_training)) +
  geom_line(aes(color = length_training)) +
  geom_point(color = "white", size = 3) +
  geom_point(aes(color = length_training)) +
  ggrepel::geom_label_repel(
    data = tbl %>% group_by(length_training) %>% filter(block_nr == max(block_nr)),
    aes(block_nr, probability_mn, label = str_c("Nr. Clusters = ", n_clusters))
  ) +
  facet_wrap(~ n_categories) +
  theme_bw() +
  scale_color_brewer(name = "Training Length\nNr. Trials", palette = "Set1") +
  scale_x_continuous(breaks = seq(1, max(tbl$block_nr), by = 1)) +
  labs(
    x = "Block Nr. (Block-Length = 25)",
    y = "Response Probability (Correct)"
  )

# To Play Around ----------------------------------------------------------


## With Feedback ----------------------------------------------------------

# using max probs instead of prob values
# tbl_used$preds <- apply(
#   l_pred$cat_probs, 1, FUN = function(x) which.max(x)
# )
# tbl_used$n_training <- n
# tbl_used$n_categories <- length(unique(tbl_used$category))
# l <- summarize_blocks(tbl_used, 17, FALSE)
# plot_block_summary(l)

stimuli <- tbl_used[, c("x1", "x2")]
features_cat <- c("x1", "x2") 
features_cont <- c() 
n_values_cat <- 2
feedback <- NULL
print_posterior <- FALSE
previous_learning <- NULL
feedback <- tbl_used$category
n_categories <- 100

params <- list(
  "coupling" = .5044,
  "phi" = .7738, #1,
  "salience_f" = .6888,#1,
  "salience_l" = .1615,#1
  "a_0" = 2,
  "lambda_0" = 1,
  "sigma_sq_0" = (max(tbl_used[, c("x1", "x2")]) / 5) ^ 2,
  "sigma_sq_0" = (max(tbl_used[, c("x1", "x2")]) / 4) ^ 2,
  "mu_0" = (min(tbl_used[, c("x1", "x2")]) + max(tbl_used[, c("x1", "x2")])) / 2
)

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

