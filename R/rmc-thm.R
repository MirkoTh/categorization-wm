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
  #' assignments. Defaults ot NULL such that inferred categories are saved
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


