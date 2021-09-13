library(tidyverse)
library(docstring)
library(grid)
library(gridExtra)

dirs_local <- c("R/utils.R", "R/rmc-thm.R")
walk(dirs_local, source)


# Settings ----------------------------------------------------------------

stepsize <- 1
max_x <- 6
x1 <- seq(1, max_x, by = stepsize)
x2 <- seq(1, max_x, by = stepsize)
tbl_features <- crossing(x1, x2)


# Create Categories -------------------------------------------------------

n_cat_per_feat <- c(2, 3)
stepsize_cat <- max_x / n_cat_per_feat
l_tbl_categories <- map(n_cat_per_feat, create_categories, tbl = tbl_features)
pl <- map2(l_tbl_categories, stepsize_cat, plot_clustered_grid)
plot_all_clusterings(pl)


# Predict from RMC ------------------------------------------------------

## rmc parameters
salience_f <- 1
salience_l <- 1
coupling <- .2
max_clusters <- 100
phi <- 1

## experimental conditions
n <- list(100, 200, 500, 1000)
ns <- rep(n, each = length(l_tbl_categories))
l_tbl_categories_prep <- rep(l_tbl_categories, length(n))

## nr of experiments
n_experiments <- 10
l_tbl_categories_prep <- rep(l_tbl_categories_prep, n_experiments)
ns <- rep(ns, n_experiments)

## run experiments
l_tbl_categories <- map2(
  l_tbl_categories_prep,
  ns,
  predict_rmc_n,
  max_clusters = max_clusters,
  salience_f = salience_f,
  salience_l = salience_l,
  coupling = coupling,
  phi = phi
)

# Summarize Results -------------------------------------------------------

n_trials_per_block <- 50 # n trials per blocks for summarizing and plotting

l_summary <- map(
  l_tbl_categories, summarize_blocks, n_trials_per_block = n_trials_per_block, 
  show_plot = FALSE
)

plot_block_summary(l_summary)


