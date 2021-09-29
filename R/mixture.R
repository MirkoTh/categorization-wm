library(MASS)
library(tidyverse)



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



# Simulate Data From a Model ----------------------------------------------

## Stimuli
tbl_categories <- l_tbl_categories[[1]]
tbl_categories <- tbl_categories[sample(1:nrow(tbl_categories), 200, replace = TRUE), ]

n_trials <- nrow(tbl_categories)

## Parameters
sd_error_1 <- .1
sd_error_2 <- .1
pm <- .5
sd_m_cont_1 <- .1
sd_m_cont_2 <- sd_m_cont_1 # same for simplicity; at fitting, they could differ


## Sample Data from Model
l_pm <- rbernoulli(n_trials, pm) %>% as.numeric()
lm1 <- rnorm(n_trials, tbl_categories$x1, sd_m_cont_1)
lg1 <- runif(n_trials, min(tbl_categories$x1), max(tbl_categories$x1))
lg2 <- runif(n_trials, min(tbl_categories$x2), max(tbl_categories$x2))
mu1 <- l_pm * lm1 + (1 - l_pm) * lg1
lm2 <- rnorm(n_trials, tbl_categories$x2, sd_m_cont_2)
mu2 <- l_pm * lm2 + (1 - l_pm) * lg2


mu <- cbind(mu1, mu2)
rownames(mu) <- 1:nrow(mu)

l_mu <- split(mu, sort(as.numeric(rownames(mu))))

m_resp_sample <- map(
  l_mu, mvrnorm, 
  n = 1, Sigma = matrix(c(sd_error_1, 0, 0, sd_error_2), nrow = 2, byrow = TRUE)
) %>%
  unlist() %>%
  matrix(ncol = ncol(mu), byrow = TRUE)
tbl_resp_sample <- as_tibble(m_resp_sample)
names(tbl_resp_sample) <- c("x1_resp", "x2_resp")

tbl_resp_sample$l_pm <- l_pm

tbl_results <- tibble(cbind(tbl_categories, tbl_resp_sample))
tbl_results <- tbl_results %>%
  mutate(
    deviation_x1 = x1 - x1_resp,
    deviation_x2 = x2 - x2_resp
  )


## Visualize Responses
tbl_results %>%
  dplyr::select(l_pm, deviation_x1, deviation_x2) %>%
  mutate(l_pm = str_c("Is in Memory = ", l_pm)) %>%
  ggplot(aes(deviation_x1, deviation_x2)) +
  geom_hex(bins = 15) +
  #geom_density2d() +
  facet_wrap(~ l_pm) +
  theme_bw() +
  scale_fill_continuous(name = "Nr. Obser-\nvations") +
  labs(
    x = "Deviation X1",
    y ="Deviation X2"
  )


plot_stim_resp <- function(tbl, x, y, title) {
  x <- enquo(x)
  y <- enquo(y)
  ggplot(tbl_results, aes(!!x, !!y)) +
    geom_point() +
    facet_wrap(~ l_pm) +
    theme_bw() +
    labs(
      x = "Response Value",
      y = "Stimulus Value",
      title = title
    )
}

tbl_plt <- tibble(
  x = c(quo(x1_resp), quo(x2_resp)),
  y = c(quo(x1), quo(x2)),
  title = c("X1", "X2")
)

pl <- pmap(tbl_plt, plot_stim_resp, tbl = tbl_results)
plot_arrangement(pl, n_cols = 1)




