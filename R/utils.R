library(tidyverse)

create_categories <- function(tbl, n_cat_per_feat) {
  #' create categories from feature space
  #' 
  #' @description cut two-dimensional feature space with equidistant spacing
  #' into equally large categories
  #' columns starting with x are assumed to be features
  #' @param tbl \code{tibble} containing each of the two features in a column
  #' @param n_cat_per_feat \code{integer} stating how many categories to create
  #' @return the \code{tibble} with an added column stating the category
  #' 
  col_nms <- names(tbl)[starts_with("x", ignore.case = FALSE, vars = names(tbl))]
  mappings <- map(tbl[, col_nms], cut, n_cat_per_feat, labels = FALSE) %>% as_tibble()
  names(mappings) <- str_c(names(mappings), "_cat")
  tbl <- cbind(tbl, mappings) %>% as_tibble()
  tbl$category <- interaction(tbl$x1_cat, tbl$x2_cat)
  levels(tbl$category) <- 1:nlevels(tbl$category)
  tbl$category <- as.numeric(tbl$category)
  return(tbl)
}

plot_clustered_grid <- function(tbl, stepsize_cat) {
  #' plot clusters in grid
  #' 
  #' @description plot categories in a 2d space with equidistanct spacing
  #' @param tbl \code{tibble} containing the two features and the category as columns
  #' @param stepsize_cat \code{double} the stepsize to be shown on the axes and on the grid
  #' @return the ggplot2 object
  #' 
  ggplot(tbl, aes(x1, x2, group = category)) +
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
}

plot_all_clusterings <- function(pl, n_cols = 2) {
  #' plot all different clusterings on a page
  #' 
  #' @param pl all the ggplots
  #' @param n_cols nr columns of the page layout
  n_plots <- length(pl)
  n_cols <- 2
  n_rows <- ceiling(n_plots / n_cols)
  marrangeGrob(pl, nrow = n_rows, ncol = n_cols)
}


