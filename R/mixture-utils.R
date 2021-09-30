category_deviation <- function(tbl_cats, x1, x2) {
  #' calculate deviation from stimuli to category centers
  #' 
  #' @param tbl_cats \code{tibble} containing info about category centers 
  #' @param x1 x1 coordinates of presented stimuli
  #' @param x2 x2 coordinates of preswented stimuli
  #' @return deviation
  #' 
  sqrt((x1 - tbl_cats$x1_avg) ^ 2 + (x2 - tbl_cats$x2_avg) ^ 2)
}


pm_mixture <- function(
  tbl_categories, tbl_cats, cats_selected,
  sd_m_cont, sd_m_cat, bern_pcont,
  n_trials, var
) {
  #' mix categorical and continuous memory responses according to pcont
  #' 
  #' @param tbl_categories \code{tibble} with presented stimuli
  #' @param tbl_cats \code{tibble} with category centers
  #' @param cats_selected \code{vector} with selected categories
  #' @param sd_m_cont continuous sd
  #' @param sd_m_cat categorical sd
  #' @param bern_pcont \code{vector} stating whether continuous memory response
  #' @param n_trials nr of trials
  #' @param var what variable should be calculated?
  lcont <- rnorm(n_trials, as_vector(tbl_categories[, var]), sd_m_cont)
  lcat <- rnorm(n_trials, as_vector(tbl_cats[, str_c(var, "_avg")])[cats_selected], sd_m_cat)
  lm <- bern_pcont * lcont + (1 - bern_pcont) * lcat
  return(lm)
}


plot_stim_resp <- function(tbl, x, y, title) {
  #' plot univariate responses from different mpt nodes
  #' 
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

sample_from_model <- function(params, n_trials, tbl_categories, tbl_categories_prep) {
  #' sample from mixture model
  #' 
  #' @param params model parameters
  #' @param n_trials nr of trials
  #' @param tbl_categories \code{tibble} with presented stimuli
  #' @param tbl_categories_prep \code{tibble} with all possible stimuli 
  ## sample for each trials whether guessing, pcont, or pcat
  bern_pm <- rbernoulli(n_trials, params[["pm"]]) %>% as.numeric()
  bern_pcont <- rep(0, n_trials)
  pcont_cond <- rbernoulli(sum(bern_pm), params[["pcont"]]) %>% as.numeric()
  bern_pcont[bern_pm == 1] <- pcont_cond
  
  ## these are the category prototypes
  tbl_cats <- tbl_categories_prep %>% 
    group_by(category) %>% 
    summarize(x1_avg = mean(x1), x2_avg = mean(x2))
  
  ## map deviations from the category prototypes into choice probs
  l_cat_devs <- tbl_categories %>% 
    dplyr::select(x1, x2) %>%
    pmap(category_deviation, tbl_cats = tbl_cats)
  m_cat_devs <- l_cat_devs %>% unlist() %>% matrix(ncol = nrow(tbl_cats), byrow = TRUE)
  m_cat_probs <- 1 - pnorm(m_cat_devs - mean(m_cat_devs))
  cats_selected <- apply(m_cat_probs, 1, FUN = function(x) sample(1:4, 1, prob = x))
  
  ## mix continuous and categorical memory
  tbl_tmp <- tibble(
    sd_m_cont = c(params[["sd_m_cont_1"]], params[["sd_m_cont_2"]]),
    sd_m_cat = c(params[["sd_m_cat_1"]], params[["sd_m_cat_2"]]),
    var = c("x1", "x2")
  )
  lms <- pmap(
    tbl_tmp, pm_mixture, 
    tbl_categories = tbl_categories,
    cats_selected = cats_selected,
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
    n = 1, Sigma = matrix(
      c(
        params[["sd_error_1"]], 0, 0, params[["sd_error_2"]]
      ), nrow = 2, byrow = TRUE
    )
  ) %>%
    unlist() %>%
    matrix(ncol = ncol(mu), byrow = TRUE)
  tbl_resp_sample <- as_tibble(m_resp_sample)
  names(tbl_resp_sample) <- c("x1_resp", "x2_resp")
  
  tbl_resp_sample$bern_pm <- bern_pm
  tbl_resp_sample$bern_pcont <- bern_pcont
  
  return(tbl_resp_sample)
}