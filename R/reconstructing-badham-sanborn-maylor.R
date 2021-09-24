library(tidyverse)
library(docstring)
library(grid)
library(gridExtra)
library(wiqid)

dirs_local <- c("R/utils.R", "R/rmc-thm.R", "badham-sanborn-maylor-2017/rmc.R")
walk(dirs_local, source)


# Settings ----------------------------------------------------------------

stepsize <- 1
max_x <- 2
x1 <- seq(1, max_x, by = stepsize)
x2 <- seq(1, max_x, by = stepsize)
x3 <- seq(1, max_x, by = stepsize)
tbl_features <- crossing(x1, x2, x3)


# Create Categories -------------------------------------------------------

n_categories <- 2
tbl_stimuli <- create_shepard_categories(tbl_features, "I", "x1")
l_tbl_categories <- list()
l_tbl_categories[[1]] <- tbl_stimuli

# Fit RMC to simulated Data -----------------------------------------------

# in the first two blocks, they presented all stimuli once in the first part 
# of the block, and once in the second part of the block

length_block <- 16
n_blocks <- 6
idx_used <- 1
idxs <- 1:nrow(l_tbl_categories[[idx_used]])
n <- n_blocks * length_block
n_b_12 <- 2 * length_block
n_b_after <- n - n_b_12


idx_b_after <- map(rep(length_block, 4), sample, x = rep(idxs, 2), replace = FALSE) %>% unlist()
idx_b_12 <- rep(rep(idxs, 2), 2)
idx_shuffle <- c(idx_b_12, idx_b_after)

tbl_used <- l_tbl_categories[[idx_used]][idx_shuffle, ]


# Using Their Function ----------------------------------------------------


stimuli <- tbl_used %>% select(x1, x2, x3) %>% as.matrix() - 1
feedback <- tbl_used %>% select(category) %>% as_vector() %>% unname() - 1

# these are the parameters from the paper
physical.salience <- 0.6888
label.salience <- 0.1615
coupling <- 0.5044
phi <- 0.7738
l_pred_orig <- rmc.predictions(
  stimuli, feedback, physical.salience, 
  label.salience, coupling, 
  print.posterior=FALSE, phi=phi
  )
pl_orig <- plot_resp_prob_by_block(l_pred_orig, feedback, length_block)


# Using Home-Grown Functions ----------------------------------------------

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

l_pred_thm <- predict_rmc_continuous(
  stimuli = tbl_used[, c("x1", "x2", "x3")], 
  features_cat = c("x1", "x2", "x3"),
  features_cont = c(),
  n_values_cat = 2,
  n_categories = n_categories,
  feedback = tbl_used$category,
  params = params,
  previous_learning = NULL, 
  print_posterior = FALSE
)

pl_thm <- plot_resp_prob_by_block(l_pred_thm, feedback, length_block)



# Compare Two Functions ---------------------------------------------------

plot_arrangement(list(pl_orig, pl_thm), 1)
