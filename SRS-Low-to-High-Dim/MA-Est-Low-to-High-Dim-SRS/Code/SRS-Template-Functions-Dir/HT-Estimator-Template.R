# prototype function to mitigate upkeep between models
fixed_size_ht_est = function(model_data,
                                       n_sizes = c(100, 250, 500, 1000),cores = 3,
                                       m=1000,
                                       randomizer = NA){
  
  # Number of simulations per sample size controlled by m
  
  x_input_retrieve = function(dat){ # Function to determine relevant auxiliary variables
    x_columns = which(
      is.na(match(colnames(dat), # Retrieve any positions that are not named X0, y, y_star
                  c("X0", "y", "y_star"))))
    return(colnames(dat)[x_columns]) # Returns results
  }
  # Similar to function above, but gets the position of y_star in matrix or dataframe
  y_star_input_retrieve = function(dat){
    y_columns = which(is.na(match(colnames(dat), 
                                  c("y_star"))) == FALSE)
    return(c(colnames(dat)[y_columns]))
  }
  
  # sample sizes controlled by n_sizes
  # Cluster for running models in parallel
  cluster_1 = makeCluster(cores)
  registerDoSNOW(cluster_1)
  # Where the results will be stored. lasso_calculations_cvg is a vector
  fixed_size_ht_calc = foreach(samp_size = 1:length(n_sizes), .combine = "rbind",
                               .inorder = FALSE) %:%
    foreach( iteration = 1:m, .combine = "rbind",.inorder = FALSE) %dopar% {
      # setting seed for sampling
      set.seed(iteration + (1000 * (samp_size-1)))
      # Using SRS as I noticed the high bias in the SRS additive randomization
      drawing = sample(nrow(model_data), size= n_sizes[samp_size],replace = FALSE)
      # Getting the sample
      input_sampling = model_data[drawing,] 
      
      y_star = numeric(length(input_sampling[,"y"]))
      
      
      y_star = randomizer(input_sampling[,"y"])
      
      
      y_star_mu = mean(y_star)
      y_star_mu_df = data.frame(n_sizes[samp_size],y_star_mu)
      colnames(y_star_mu_df) <- c("n_size", "y_star_ht")
      return(y_star_mu_df)
    }
  stopCluster(cluster_1)
  return(fixed_size_ht_calc)
}

# randomizer needs to be a function that returns y_star as well a the means and
# variances of the scrambling distributions
# prototype function to mitigate upkeep between models
fixed_size_ht_est_variance = function(model_data,
                             n_sizes = c(100, 250, 500, 1000),cores = 3,
                             m=1000,
                             randomizer = NA){
  
  # Number of simulations per sample size controlled by m
  
  x_input_retrieve = function(dat){ # Function to determine relevant auxiliary variables
    x_columns = which(
      is.na(match(colnames(dat), # Retrieve any positions that are not named X0, y, y_star
                  c("X0", "y", "y_star"))))
    return(colnames(dat)[x_columns]) # Returns results
  }
  # Similar to function above, but gets the position of y_star in matrix or dataframe
  y_star_input_retrieve = function(dat){
    y_columns = which(is.na(match(colnames(dat), 
                                  c("y_star"))) == FALSE)
    return(c(colnames(dat)[y_columns]))
  }
  
  # sample sizes controlled by n_sizes
  # Cluster for running models in parallel
  cluster_1 = makeCluster(cores)
  registerDoSNOW(cluster_1)
  # Where the results will be stored. lasso_calculations_cvg is a vector
  fixed_size_ht_calc = foreach(samp_size = 1:length(n_sizes), .combine = "rbind",
                               .inorder = FALSE) %:%
    foreach( iteration = 1:m, .combine = "rbind",.inorder = FALSE) %dopar% {
      # setting seed for sampling
      set.seed(iteration + (1000 * (samp_size-1)))
      # Using SRS as I noticed the high bias in the SRS additive randomization
      drawing = sample(nrow(model_data), size= n_sizes[samp_size],replace = FALSE)
      # Getting the sample
      input_sampling = model_data[drawing,] 
      
      y_rand = randomizer(input_sampling[,"y"])
      y_star = y_rand$values
      
      N_pop = nrow(model_data)
      
      first_order_pi = (n_sizes[samp_size] / N_pop)
      second_order_pi = (n_sizes[samp_size] / N_pop) * 
        ((n_sizes[samp_size] - 1) / (N_pop - 1))
      
      
      y_star_mu = mean(y_star)
      
      alpha =  y_rand$alpha
      eta = y_rand$eta
      gamma = y_rand$gamma
      beta = y_rand$beta
      
      ht_variance = 0
      # temp = 0
      
      # for(item in 1:length(y_star)){
      #   temp = temp + (1 - first_order_pi) * 
      #     (y_star[item])^2
      #   temp = temp + ((second_order_pi - first_order_pi^2) / (second_order_pi)) * 
      #     sum((y_star[item] / first_order_pi) * (y_star[-item] / first_order_pi))
      # }
      ht_variance = ((1 - (n_sizes[samp_size] / 
                             nrow(model_data))) / n_sizes[samp_size]) * 
        var(y_star)#(1 / N_pop^2) * temp 
      
      
      rrt_variance = (1 / N_pop^2) * 
        sum( first_order_pi^(-2) * ((gamma/ beta^2) * 
              ((y_star^2 - (eta / beta^2)) *
              ((gamma / beta^2) + 1)^(-1)) + 
              (eta / beta^2)))
      
      variance =  rrt_variance + ht_variance
      
      # get upper and lower for confidence interval
      actual_mu = mean(model_data$y)
      upper = y_star_mu + 1.96 * sqrt(abs(variance))
      # Get lower bound 
      lower = y_star_mu - 1.96 * sqrt(abs(variance))
      
      covered_mu = (lower < actual_mu && actual_mu < upper)
     
      y_star_mu_df = data.frame(n_sizes[samp_size], y_star_mu, variance,
                                upper, lower, covered_mu)
      colnames(y_star_mu_df) <- c("n_size","y_star_mu", "variance",
                                  "upperCI", "lowerCI", "covered_mu")
      return(y_star_mu_df)
    }
  stopCluster(cluster_1)
  return(fixed_size_ht_calc)
}
