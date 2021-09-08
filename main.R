library(tidyverse)
library(docstring)
library(grid)
library(gridExtra)

source("R/utils.R")


# Settings ----------------------------------------------------------------

stepsize <- 1
max_x <- 100
x1 <- seq(1, max_x, by = stepsize)
x2 <- seq(1, max_x, by = stepsize)
tbl_features <- crossing(x1, x2)


# Create Categories -------------------------------------------------------

n_cat_per_feat <- c(2, 4, 8, 10)
stepsize_cat <- max_x / n_cat_per_feat
l_tbl_categories <- map(n_cat_per_feat, create_categories, tbl = tbl_features)
pl <- map2(l_tbl_categories, stepsize_cat, plot_clustered_grid)
plot_all_clusterings(pl)

