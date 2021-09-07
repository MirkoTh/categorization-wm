library(tidyverse)
library(docstring)

source("R/utils.R")


# Settings ----------------------------------------------------------------

stepsize <- 1
max_x <- 100
x1 <- seq(1, max_x, by = stepsize)
x2 <- seq(1, max_x, by = stepsize)
tbl_features <- crossing(x1, x2)


# Create Categories -------------------------------------------------------

n_cat_per_feat <- 5
stepsize_cat <- max_x / n_cat_per_feat
tbl_categories <- create_categories(tbl_features, n_cat_per_feat)

ggplot(tbl_categories, aes(x1, x2, group = category)) +
  geom_raster(aes(fill = category)) +
  theme_bw() + 
  theme(
    panel.background = element_rect(fill = NA),
    panel.ontop = TRUE,
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size = .1),
  ) +
  scale_fill_gradient(name = "Category\n",
                      low = "#FFFFFF",
                      high = "#012345") +
  scale_x_continuous(breaks = seq(0, max_x, stepsize_cat)) +
  scale_y_continuous(breaks = seq(0, max_x, stepsize_cat)) +
  labs(
    x = expression(X["1"]),
    y = expression(X["2"])
  )
