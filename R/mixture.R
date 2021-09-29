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
pl <- map2(l_tbl_categories, stepsize_cat, plot_clustered_grid)
plot_arrangement(pl)


# Simulate Data From a Model ----------------------------------------------

## Stimuli
tbl_categories_prep <- l_tbl_categories[[1]]
tbl_categories <- tbl_categories_prep[sample(1:nrow(tbl_categories_prep), 200, replace = TRUE), ]

n_trials <- nrow(tbl_categories)

## Parameters
### pms
pm <- .5
pcont <- .5

### sds
sd_error_1 <- .1
sd_error_2 <- .1
sd_m_cont_1 <- .1
sd_m_cont_2 <- sd_m_cont_1 # same for simplicity; at fitting, they could differ
sd_m_cat_1 <- .1
sd_m_cat_2 <- sd_m_cat_1


## Sample Data from Model
bern_pm <- rbernoulli(n_trials, pm) %>% as.numeric()
bern_pcont <- rbernoulli(n_trials, pcont) %>% as.numeric()


tbl_cats <- tbl_categories_prep %>% 
  group_by(category) %>% 
  summarize(x1_avg = mean(x1), x2_avg = mean(x2))

category_deviation <- function(tbl_cats, x1, x2) {
  sqrt((x1 - tbl_cats$x1_avg) ^ 2 + (x2 - tbl_cats$x2_avg) ^ 2)
}

l_cat_devs <- tbl_categories %>% 
  dplyr::select(x1, x2) %>%
  pmap(category_deviation, tbl_cats = tbl_cats)
m_cat_devs <- l_cat_devs %>% unlist() %>% matrix(ncol = nrow(tbl_cats), byrow = TRUE)
m_cat_probs <- 1 - pnorm(m_cat_devs - mean(m_cat_devs))
cats_selected <- apply(m_cat_probs, 1, FUN = function(x) sample(1:4, 1, prob = x))



pm_mixture <- function(
  tbl_categories, tbl_cats, 
  sd_m_cont, sd_m_cat, bern_pcont,
  n_trials, var
) {
  lcont <- rnorm(n_trials, as_vector(tbl_categories[, var]), sd_m_cont)
  lcat <- rnorm(n_trials, as_vector(tbl_cats[, str_c(var, "_avg")])[cats_selected], sd_m_cat)
  lm <- bern_pcont * lcont + (1 - bern_pcont) * lcat
  return(lm)
}

tbl_tmp <- tibble(
  sd_m_cont = c(sd_m_cont_1, sd_m_cont_2),
  sd_m_cat = c(sd_m_cat_1, sd_m_cat_2),
  var = c("x1", "x2")
)

lms <- pmap(
  tbl_tmp, pm_mixture, 
  tbl_categories = tbl_categories,
  tbl_cats = tbl_cats,
  bern_pcont = bern_pcont,
  n_trials = n_trials
)

lm1 <- lms[[1]]
lm2 <- lms[[2]]
lg1 <- runif(n_trials, min(tbl_categories$x1), max(tbl_categories$x1))
lg2 <- runif(n_trials, min(tbl_categories$x2), max(tbl_categories$x2))
mu1 <- bern_pm * lm1 + (1 - bern_pm) * lg1
mu2 <- bern_pm * lm2 + (1 - bern_pm) * lg2


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

tbl_resp_sample$bern_pm <- bern_pm
tbl_resp_sample$bern_pcont <- bern_pcont

tbl_results <- tibble(cbind(tbl_categories, tbl_resp_sample))
tbl_results <- tbl_results %>%
  mutate(
    deviation_x1 = x1 - x1_resp,
    deviation_x2 = x2 - x2_resp
  )


## Visualize Responses
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


plot_stim_resp <- function(tbl, x, y, title) {
  x <- enquo(x)
  y <- enquo(y)
  ggplot(tbl_results, aes(!!x, !!y)) +
    geom_point() +
    facet_wrap(bern_pcont ~ bern_pm) +
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