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
  #' @return the predicted category probabilities and category assignments as a {list}
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
                          dim = c(max_clusters, n_features, n_values))
  n_clusters <- 1 # position of the lowest currently empty cluster
  
  out <- c()
  out$cat_probs <- matrix(nrow = nrow(stimuli), ncol = n_values)
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
      matrix(rep(stimuli[i, ], n_values), nrow = n_values, byrow = TRUE),
      seq(0, (n_values - 1), by = 1)
    )
    # cols: cluster nr, feature nr (including cat label), values
    feature_index <- cbind(
      rep(1:n_clusters, times = n_features*n_values),
      rep(1:n_features, times = n_values, each = n_clusters),
      rep(t(possible_stimuli + 1), each = n_clusters)
    )
    this_num <- array(
      feature_counts[feature_index], 
      dim = c(n_clusters, n_features, n_values)
    )
    # what goes into den(ominator) to achieve a uniform prior?
    multiply_salience <- c(
      rep(n_values, (n_features - 1)), 
      length(unique(feedback))
      )
    this_den <- array(
      outer(cluster_counts[1:n_clusters], multiply_salience * salience, FUN = "+"),
      dim=c(n_clusters, n_features, n_values)
    )
    log_likelihood <- apply(
      log(this_num) - log(this_den), MARGIN = c(1, 3), FUN = sum
    )
    log_posterior <- matrix(
      log_prior, nrow = n_clusters, ncol = n_values
    ) + log_likelihood
    
    if(print_posterior){
      print(exp(log_posterior[, feedback[i] + 1]))
    }
    
    # compute prediction
    label_posterior <- colSums(exp(log_posterior - max(log_posterior)))
    out$cat_probs[i, ] <- label_posterior^phi / sum(label_posterior^phi)
    
    # update cluster assignment and count variables
    if(update_assignments){
      # using Anderson (1991) update rule
      assignments[i] <- which.max(log_posterior[, feedback[i] + 1])
    }
    cluster_counts[assignments[i]] <- cluster_counts[assignments[i]] + 1
    feature_index_update <- cbind(rep(assignments[i], times=n_features), 
                                  1:n_features, 
                                  c(stimuli[i, ], feedback[i]) + 1)
    feature_counts[feature_index_update] <- (
      feature_counts[feature_index_update] + 1
    )
    n_clusters <- max(assignments) + 1
  }
  out$assignments <- assignments
  return(out)
}
