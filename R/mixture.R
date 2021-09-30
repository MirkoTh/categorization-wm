library(MASS)
library(tidyverse)

dirs_local <- c("R/utils.R", "R/rmc-thm.R", "R/mixture-utils.R")
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


# Simulate Data From a Model ----------------------------------------------

n_trials <- 400
## Stimuli
tbl_categories_prep <- l_tbl_categories[[1]]
tbl_categories <- tbl_categories_prep[
  sample(1:nrow(tbl_categories_prep), n_trials, replace = TRUE), 
  ]


## Parameters

params <- list(
  pm = .5,
  pcont = .5,
  sd_error_1 = .1,
  sd_error_2 = .1,
  sd_m_cont_1 = .1,
  sd_m_cont_2 = .1, # same for simplicity; at fitting, they could differ
  sd_m_cat_1 = .1,
  sd_m_cat_2 = .1
)


## Sample from Model
tbl_resp_sample <- sample_from_model(
  params, n_trials, tbl_categories, tbl_categories_prep
  )

## combine stimuli and simulated responses
tbl_results <- tibble(cbind(tbl_categories, tbl_resp_sample))
tbl_results <- tbl_results %>%
  mutate(
    deviation_x1 = x1 - x1_resp,
    deviation_x2 = x2 - x2_resp
  )


# Visualize Responses -----------------------------------------------------


tbl_results %>%
  dplyr::select(bern_pm, bern_pcont, deviation_x1, deviation_x2) %>%
  mutate(
    bern_pm = str_c("Is in Memory = ", bern_pm),
    bern_pcont = str_c("Is Continuous = ", bern_pcont),
  ) %>%
  ggplot(aes(deviation_x1, deviation_x2)) +
  geom_hex(bins = 15) +
  #geom_density2d() +
  facet_wrap(bern_pcont ~ bern_pm) +
  theme_bw() +
  scale_fill_continuous(name = "Nr. Obser-\nvations") +
  labs(
    x = "Deviation X1",
    y ="Deviation X2"
  )

tbl_plt <- tibble(
  x = c(quo(x1_resp), quo(x2_resp)),
  y = c(quo(x1), quo(x2)),
  title = c("X1", "X2")
)

pl <- pmap(tbl_plt, plot_stim_resp, tbl = tbl_results)
plot_arrangement(pl, n_cols = 1)