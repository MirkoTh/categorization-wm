# original code from Badham, Sanborn, & Maylor (2017)
# https://osf.io/m7tck/files/
# date of dump: 2021-09-08
# implementation of the Rational Model of Categorization (Anderson, 1991) for the 
# Shepard, Hovland, and Jenkins (1961) task

# predictions from the RMC
predict_rmc <- function(
  #' predict from Anderson's RMC
  #' 
  #' @description predict from Anderson's rational model of categorization (1991)
  #' using the original infererence algorithm for the posterior (aka local MAP)
  #' @param stimuli \code{matrix} containing the feature values of one stimulus per row
  #' @param n_values \code{integer} nr of values for the features assuming all 
  #' features have the same nr
  #' @param feedback \code{vector} containing category labels as integers
  #' @param salience_f \code{integer} feature-salience prior
  #' @param salience_l \code{integer} label-salience prior
  #' @param coupling \code{integer} coupling probability c
  #' @param phi \code{numeric} scaling parameter for response probabilities
  #' @param max_clusters \code{integer} max nr of clusters to use in the code
  #' @param assignments \code{vector} of integers stating the category
  #' assignments. Defaults ot NULL such that inferred categories are saved
  #' @param print_posterior {logical} stating whether the posterior should
  #' be printed while code is running
  #' @return the predicted category probabilities and category assignments
  #' as a \code{list}
  #' 
  stimuli,
  n_values, # nr of values for the features assuming all features have the same nr
  feedback, # the category assignment
  salience_f, # prior on obtaining *F*eatures with different values
  salience_l, # prior on obtaining different category *L*abels
  coupling,
  phi = 1, 
  max_clusters = 100, 
  assignments = NULL, 
  print_posterior = FALSE
) {
  n_stimuli <- nrow(stimuli)
  n_features <- ncol(stimuli) + 1 # including label feature
  n_categories <- length(unique(feedback))
  
  if(is.null(assignments)){
    assignments <- rep(0, n_stimuli) # assignment of stimuli to clusters
    update_assignments <- TRUE
  }else{
    update_assignments <- FALSE # assignments given to model in advance
  }
  cluster_counts <- rep(0, max_clusters) # counts of stimuli in each cluster
  salience <- c(rep(salience_f, times = ncol(stimuli)), salience_l)
  # feature counts is a cluster by feature by value array 
  # and includes pseudocounts from prior
  feature_counts <- array(rep(salience, each = max_clusters),
                          dim = c(max_clusters, n_features, max(n_values, n_categories)))
  n_clusters <- 1 # position of the lowest currently empty cluster
  
  out <- c()
  out$cat_probs <- matrix(nrow = nrow(stimuli), ncol = n_categories)
  for(i in 1:n_stimuli){
    # calculate prior, likelihoods, and posterior
    # clusters that already have stimuli
    log_prior <- log(coupling) + log(cluster_counts[1:n_clusters])
    # new cluster prior (in first empty position)
    log_prior[n_clusters] <- log(1 - coupling)
    log_prior <- log_prior - log(
      1 - coupling * (1 - sum(cluster_counts[1:n_clusters]))
    )
    # add labels to end of stimuli
    possible_stimuli <- cbind(
      matrix(rep(stimuli[i, ], n_categories), nrow = n_categories, byrow = TRUE),
      seq(1, (n_categories), by = 1)
    )
    # cols: cluster nr, feature nr (including cat label), values
    feature_index <- cbind(
      rep(1:n_clusters, times = n_features * n_categories),
      rep(1:n_features, times = n_categories, each = n_clusters),
      rep(t(possible_stimuli), each = n_clusters)
    )
    this_num <- array(
      feature_counts[feature_index], 
      dim = c(n_clusters, n_features, n_categories)
    )
    # what goes into den(ominator) to achieve a uniform prior?
    multiply_salience <- c(
      rep(n_values, (n_features - 1)), 
      n_categories
    )
    this_den <- array(
      outer(cluster_counts[1:n_clusters], multiply_salience * salience, FUN = "+"),
      dim=c(n_clusters, n_features, n_categories)
    )
    log_likelihood <- apply(
      log(this_num) - log(this_den), MARGIN = c(1, 3), FUN = sum
    )
    log_posterior <- matrix(
      log_prior, nrow = n_clusters, ncol = n_categories
    ) + log_likelihood
    
    if(print_posterior){
      print(exp(log_posterior[, feedback[i]]))
    }
    
    # compute prediction
    label_posterior <- colSums(exp(log_posterior - max(log_posterior)))
    out$cat_probs[i, ] <- label_posterior^phi / sum(label_posterior^phi)
    
    # update cluster assignment and count variables
    if(update_assignments){
      # using Anderson (1991) update rule
      assignments[i] <- which.max(log_posterior[, feedback[i]])
    }
    cluster_counts[assignments[i]] <- cluster_counts[assignments[i]] + 1
    feature_index_update <- cbind(rep(assignments[i], times=n_features), 
                                  1:n_features, 
                                  c(stimuli[i, ], feedback[i]))
    feature_counts[feature_index_update] <- (
      feature_counts[feature_index_update] + 1
    )
    n_clusters <- max(assignments) + 1
  }
  out$assignments <- assignments
  return(out)
}

predict_rmc_continuous <- function(
  #' predict from Anderson's RMC using continuous features
  #' 
  #' @description predict from Anderson's rational model of categorization (1991)
  #' using the original infererence algorithm for the posterior (aka local MAP)
  #' @param stimuli \code{matrix} containing the feature values of one stimulus per row
  #' @param features_cat \code{vector} with names of categorical features as strings
  #' @param features_cont \code{vector} with names of continuous features as strings
  #' @param n_values_cat \code{vector} nr of values for the categorical features as strings
  #' features have the same nr
  #' @param feedback \code{integer vector} containing category labels as integers
  #' @param a0 \code{integer} confidence in prior variance
  #' @param lambda0 \code{integer} confidence in prior mean
  #' @param coupling \code{integer} coupling probability c
  #' @param phi \code{numeric} scaling parameter for response probabilities
  #' @param max_clusters \code{integer} max nr of clusters to use in the code
  #' @param assignments \code{vector} of integers stating the category
  #' assignments. Defaults ot NULL such that inferred categories are saved
  #' @param print_posterior {logical} stating whether the posterior should
  #' be printed while code is running
  #' @return the predicted category probabilities and category assignments
  #' as a \code{list}
  #' 
  stimuli,
  features_cat,
  features_cont,
  n_values_cat,
  feedback,
  salience_f,
  salience_l,
  coupling,
  phi = 1, 
  max_clusters = 100, 
  assignments = NULL, 
  print_posterior = FALSE
) {
  n_stimuli <- nrow(stimuli)
  n_features_cat <- length(features_cat) + 1 # including label feature
  n_features_cont <- length(features_cont)
  n_categories <- length(unique(feedback))
  
  if(is.null(assignments)){
    assignments <- rep(0, n_stimuli) # assignment of stimuli to clusters
    update_assignments <- TRUE
  }else{
    update_assignments <- FALSE # assignments given to model in advance
  }
  cluster_counts <- rep(0, max_clusters) # counts of stimuli in each cluster
  salience <- c(rep(salience_f, times = ncol(stimuli)), salience_l)
  # feature counts is a cluster by feature by value array 
  # and includes pseudocounts from prior
  feature_counts <- array(rep(salience, each = max_clusters),
                          dim = c(max_clusters, n_features_cat, max(n_values_cat, n_categories)))
  n_clusters <- 1 # position of the lowest currently empty cluster
  
  out <- c()
  out$cat_probs <- matrix(nrow = nrow(stimuli), ncol = n_categories)
  
  mu_0 <- (min(stimuli) + max(stimuli)) / 2
  sigma_sq_0 <- (max(stimuli) / 5) ^ 2
  a_0 <- 2
  lambda_0 <- 1
  for(i in 1:n_stimuli){
    # calculate prior, likelihoods, and posterior
    # clusters that already have stimuli
    log_prior <- log(coupling) + log(cluster_counts[1:n_clusters])
    # new cluster prior (in first empty position)
    log_prior[n_clusters] <- log(1 - coupling)
    log_prior <- log_prior - log(
      1 - coupling * (1 - sum(cluster_counts[1:n_clusters]))
    )
    
    # likelihood for categorical features
    pdfs_cat_log <- pdf_cat_log(
      stimuli, i, features_cat, feature_counts, 
      n_values_cat, n_categories, n_clusters, n_features_cat
    )
    
    # likelihood for continuous features
    pdfs_cont_log <- pdf_cont_log(
      stimuli, assignments, i, features_cont, mu_0, lambda_0, a_0, sigma_sq_0
    )
    pdfs_cont_log <- array(
      rep(pdfs_cont_log, n_categories),
      dim = c(n_features_cont,n_clusters,  n_categories)
    ) %>% aperm(c(2, 1, 3))
    
    log_likelihood_prep <- abind::abind(
      pdfs_cat_log,
      pdfs_cont_log,
      along = 2
    )
    
    log_likelihood <- apply(
      log_likelihood_prep, MARGIN = c(1, 3), FUN = sum
    )
    log_posterior <- matrix(
      log_prior, nrow = n_clusters, ncol = n_categories
    ) + log_likelihood
    
    if(print_posterior){
      print(exp(log_posterior[, feedback[i]]))
    }
    
    # compute prediction
    label_posterior <- colSums(exp(log_posterior - max(log_posterior)))
    out$cat_probs[i, ] <- label_posterior^phi / sum(label_posterior^phi)
    
    # update cluster assignment and count variables
    if(update_assignments){
      # using Anderson (1991) update rule
      assignments[i] <- which.max(log_posterior[, feedback[i]])
    }
    cluster_counts[assignments[i]] <- cluster_counts[assignments[i]] + 1
    feature_index_update <- cbind(
      rep(assignments[i], times=n_features_cat), 
      1:n_features_cat, 
      c(stimuli[i, features_cat] %>% as_vector(), feedback[i])
    )
    feature_counts[feature_index_update] <- (
      feature_counts[feature_index_update] + 1
    )
    n_clusters <- max(assignments) + 1
  }
  out$assignments <- assignments
  return(out)
}

pdf_cont_log <- function(
  stimuli, assignments, i, features_cont, 
  mu_0, lambda_0, a_0, sigma_sq_0
) {
  tbl_part_1 <- tibble(stimuli[1:(i-1), ], assignments = assignments[1:(i-1)])
  tbl_part_2 <- crossing(
    stimuli[which(assignments[1:i] == 0), ],
    assignments = unique(assignments)
  )
  tbl_part_2$assignments[tbl_part_2$assignments == 0] <- max(assignments) + 1
  tbl_relevant <- rbind(tbl_part_1, tbl_part_2)
  n_clusters <- length(unique(tbl_relevant$assignments))
  tbl_summary <- tbl_relevant %>%
    select(all_of(features_cont), assignments) %>%
    group_by(assignments) %>%
    summarize_if(
      is.numeric, list(length, mean, var)
    )
  tbl_summary[is.na(tbl_summary)] <- 0
  names(tbl_summary) <- str_replace(names(tbl_summary), "fn1$", "n") %>% 
    str_replace("fn2$", "mean") %>%
    str_replace("fn3$", "var")
  tbl_summary$lambda_i <- tbl_summary[, 2] %>% as_vector() + lambda_0
  tbl_summary$a_i <- tbl_summary[, 2] %>% as_vector() + a_0
  col_longer <- names(tbl_summary)[str_detect(names(tbl_summary), "n$|mean$|var$")]
  tbl_summary <- tbl_summary %>%
    pivot_longer(all_of(col_longer)) %>%
    separate(name, c("variable", "aggregation"), sep = "_") %>%
    pivot_wider(names_from = variable, values_from = value)
  agg_stats <- tbl_summary %>%
    split(~ aggregation) %>%
    map(~ select(.x, all_of(features_cont)))
  mu_i <- (
    (lambda_0 * mu_0 + agg_stats[["n"]] * agg_stats[["mean"]]) /
      (lambda_0 + agg_stats[["n"]])
  )
  above_1 <- a_0 * sigma_sq_0
  above_2 <- (agg_stats[["n"]] - 1) * agg_stats[["var"]]
  above_3 <- (
    (lambda_0 * agg_stats[["n"]]) / 
      (lambda_0 + agg_stats[["n"]])
  ) * (mu_0 - agg_stats[["mean"]]) ^ 2
  below <- a_0 + agg_stats[["n"]]
  sigma_sq_i <- (above_1 + above_2 + above_3) / below
  update_i <- tbl_summary %>% 
    group_by(assignments) %>%
    summarize(
      lambda_i = min(lambda_i),
      a_i = min(a_i)
    ) %>% ungroup() %>% select(-assignments)
  scale <- sqrt(sigma_sq_i) * sqrt(1 + 1 / update_i$lambda_i)
  
  v_x <- rep(stimuli[i, features_cont] %>% as_vector(), n_clusters)
  v_mu_i <- as.vector(t(mu_i))
  v_scale <- as.vector(t(scale))
  v_ai <- rep(update_i$a_i, each = n_features_cont)
  
  v_out <- suppressWarnings(
    pmap_dbl(
      list(v_x, v_mu_i, v_scale, v_ai), dt2)
  ) %>% log()
  
  #  %>% 
  # matrix(ncol = n_features_cont, byrow = TRUE) %>%
  #   as_tibble()
  # names(tbl_out) <- features_cont
  return(v_out)
}


pdf_cat_log <- function(
  stimuli, i, features_cat, feature_counts, 
  n_values_cat, n_categories, n_clusters,
  n_features_cat
) {
  # add labels to end of stimuli
  possible_stimuli <- cbind(
    matrix(rep(stimuli[i, features_cat], n_categories), 
           nrow = n_categories, byrow = TRUE),
    seq(1, (n_categories), by = 1)
  )
  
  # cols: cluster nr, feature nr (including cat label), values
  feature_index <- cbind(
    rep(1:n_clusters, times = n_features_cat * n_categories),
    rep(1:n_features_cat, times = n_categories, each = n_clusters),
    rep(t(possible_stimuli) %>% as_vector(), each = n_clusters)
  )
  this_num <- array(
    feature_counts[feature_index], 
    dim = c(n_clusters, n_features_cat, n_categories)
  )
  # what goes into den(ominator) to achieve a uniform prior?
  multiply_salience <- c(
    rep(n_values_cat, (n_features_cat - 1)), 
    n_categories
  )
  this_den <- array(
    outer(cluster_counts[1:n_clusters], multiply_salience * salience, FUN = "+"),
    dim=c(n_clusters, n_features_cat, n_categories)
  )
  return(log(this_num) - log(this_den))
}



predict_rmc_n <- function(
  #' predict n randomly sampled stimuli from stimulus space using Anderson's RMC
  #' 
  #' @param tbl \code{tibble} with feature values and category labels as columns
  #' @param n number of trials to sample with replacement and predict
  #' @param max_clusters \code{integer} max nr of clusters to use in the code
  #' @param salience_f \code{integer} feature-salience prior
  #' @param salience_l \code{integer} label-salience prior
  #' @param coupling \code{integer} coupling probability c
  #' @param phi \code{numeric} scaling parameter for response probabilities
  #' @param assignments \code{vector} of integers stating the category
  #' assignments. Defaults to NULL such that inferred categories are saved
  #' @param print_posterior {logical} stating whether the posterior should
  #' be printed while code is running
  #' @return the predicted category probabilities and category assignments
  #' as a \code{list}
  #'
  tbl,
  n,
  max_clusters,
  salience_f = 1,
  salience_l = 1,
  coupling = .5,
  phi = 1
) {
  
  idxs <- 1:nrow(tbl)
  idx_shuffle <- sample(idxs, n, replace = TRUE)
  tbl_used <- tbl[idx_shuffle, ]
  
  stimuli <- as.matrix(tbl_used[, c("x1", "x2")])
  feedback <- tbl_used$category
  
  l_pred <- predict_rmc(
    stimuli = stimuli,
    n_values = length(unique(tbl_used$x1)), # assuming all features same
    feedback = feedback,
    salience_f = salience_f,
    salience_l = salience_l,
    coupling = coupling,
    max_clusters = max_clusters
  )
  tbl_used$preds <- apply(
    l_pred$cat_probs, 1, FUN = function(x) which.max(x)
  )
  tbl_used$n_training <- n
  tbl_used$n_categories <- length(unique(tbl_used$category))
  
  return(tbl_used)
}


summarize_blocks <- function(
  #' summarize categorized stimuli into n_blocks and 
  #' show plot summarized by block if required
  #' 
  #' @param tbl \code{tibble} with category labels and category predictions as columns
  #' @param n_trials_per_block number of traisl per block to summarize
  #' @param show_plot \code{logical} should a raster plot be shown by block
  #' to see categorization predictions along with true categories? defaults to TRUE
  
  #' @return a \code{list} with two tibbles; (a) the summarized results and
  #' (b) the assignments in the originally handed over tibble
  #'
  tbl, 
  n_trials_per_block, 
  show_plot = TRUE
) {
  n_blocks <- ceiling(nrow(tbl)/n_trials_per_block)
  tbl$block_nr <- rep(
    seq(1, n_blocks), 
    each = ceiling(nrow(tbl)/n_blocks)
  )[1:nrow(tbl)]
  # summarize accuracy per block
  tbl_results <- tbl %>%
    group_by(n_categories, n_training, block_nr) %>%
    summarize(
      accuracy = mean(category == preds)
    ) %>% ungroup()
  # prepare tbl for grouped plot
  tbl_long <- tbl %>%
    pivot_longer(
      cols = c(category, preds),
      names_to = "Label",
      values_to = "Value"
    )
  # plots: accuracy over blocks and raster
  pl <- ggplot(tbl_long, aes(x1, x2, group = Value)) +
    geom_raster(aes(fill = Value)) +
    facet_wrap(block_nr ~ Label, ncol = ceiling(n_blocks / 2) * 2)
  if (show_plot) {
    grid.draw(pl)
  }
  return(list(
    tbl_results, # summarized results
    tbl # category assignments
  ))
  
}


plot_block_summary <- function(l) {
  l_blocks <- map(l, 1)
  tbl_blocks <- reduce(l_blocks, rbind)
  tbl_blocks$n_training <- as.factor(tbl_blocks$n_training)
  tbl_summary <- tbl_blocks %>%
    group_by(n_categories, n_training, block_nr) %>%
    summarize(
      accuracy = mean(accuracy)
    )
  
  pl <- ggplot(tbl_summary, aes(block_nr, accuracy, group = n_training)) +
    geom_line(aes(color = n_training)) +
    geom_point(color = "white", size = 3) +
    geom_point(aes(color = n_training)) +
    facet_wrap(~ n_categories, ncol = 3) +
    theme_bw() +
    scale_color_brewer(name = "Nr. Training\nTrials", palette = "Set1") +
    scale_x_continuous(breaks = seq(1, max(tbl_blocks$block_nr))) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      x = "Block Nr.",
      y = "Proportion Correct"
    )
  grid.draw(pl)
}


