# prototype function to mitigate upkeep between models
linear_reg_model_results_fixed_size = function(model_data,
                                               n_sizes = c(100, 250, 500, 1000), cores = 3,
                                               m=1000,
                                               randomizer = NA,
                                               file_name = "linear_results"){
  
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
  
  all_x_input_retrieve = function(dat){
    x_columns = is.na(match(colnames(dat), # Retrieve any positions that are y or y_star
                            c("y", "y_star")))
    return(c(colnames(dat)[x_columns]))
  }
  
  # Function to find variable with highest correlation to y
  locate_high_corr <- function(dat){
    # Get a vector of the correlations of input variables with y
    corr_vector = cor(dat[,c(2:ncol(dat))])[1:(ncol(dat) -2),"y"]
    # Get the index of the variable with the highest correlation to y
    var_index = which(corr_vector == max(corr_vector))
    return(var_index)
  }
  
  # Construct a matrix of the population data
  pop_x_col = x_input_retrieve(model_data) # get column positions for auxiliary variables of population
  pop_all_x_col = all_x_input_retrieve(model_data)
  # Get the population matrix
  pop_with_out_constant_mat = as.matrix(model_data[,pop_x_col])
  
  pop_with_out_constant_mat = cbind(X0 = rep(1,nrow(pop_with_out_constant_mat)), 
                                    pop_with_out_constant_mat)
  
  # sample sizes controlled by n_sizes
  # Cluster for running models in parallel
  cluster_1 = makeCluster(cores)
  registerDoSNOW(cluster_1)
  # Where the results will be stored. lasso_calculations_cvg is a vector
  linear_calculations_cvg = foreach(samp_size = 1:length(n_sizes), 
                                    .combine = "rbind", 
                                    .inorder = FALSE) %:%
    foreach( iteration = 1:m, .combine = "rbind",.inorder = FALSE) %dopar% {
      
      # setting seed for sampling
      set.seed(iteration + 1000 * (samp_size-1))
      
      
      drawing = sample(nrow(model_data), size= n_sizes[samp_size],replace = FALSE)
      # Getting the sample
      input_sampling = model_data[drawing,] 
      
      y_star = numeric(length(input_sampling[,"y"]))
      
      # Getting the sample
      input_sampling = model_data[drawing,] 
      
      
      
      y_star = randomizer(input_sampling[,"y"])
      
      #print("y_star")
      #print(length(y_star))
      # Get the cross validated model for glmnet
      # Get the auxiliary variables column positions
      
      input_sampling = cbind(input_sampling, y_star) # Bind y_star column
      
      x_input_col = x_input_retrieve(input_sampling)
      
      # Get the y_star position
      y_input_col = y_star_input_retrieve(input_sampling)
      #print(colnames(input_sampling))
      
      # Get all the predictors and the response y
      input_sampling_without_const = as.matrix(input_sampling[,c(x_input_col,
                                                                 y_input_col)])
      
      input_sampling_without_const = cbind(X0 = rep(1,nrow(input_sampling_without_const)), 
                                           input_sampling_without_const)
      
      x_input_col = all_x_input_retrieve(input_sampling_without_const)
      
      # Get the y_star position
      y_input_col = y_star_input_retrieve(input_sampling_without_const)
      
      # Building the LASSO model
      #print(colnames(input_sampling_without_const[,x_input_col]))
      #print(colnames(input_sampling_without_const[,x_input_col]))
      
      betas_lin =  solve(t(input_sampling_without_const[,x_input_col])  %*%
                           input_sampling_without_const[,x_input_col], tol = NULL) %*% t(input_sampling_without_const[,x_input_col]) %*%
        y_star
      
      
      x_cols = all_x_input_retrieve(input_sampling_without_const)
      
      y_star_lin_sample =  input_sampling_without_const[,x_cols] %*% betas_lin
      
      
      y_star_lin_pop = pop_with_out_constant_mat %*% betas_lin
      
      
      #print(paste("population", y_star_lin_pop))
      y_star_lin_estimate = mean((y_star - y_star_lin_sample)) + 
        mean(y_star_lin_pop) # Difference estimator result
      #print(dim(y_star_lin_sample))
      #print(dim(y_star_lin_pop))
      # print(paste("estimate", y_star_lin_estimate))
      
      
      beta_values = as.vector(t(betas_lin[,1]))
      #print(beta_values)
      #print(length(beta_values))
      #beta_values = matrix(beta_values, nrow = 1, ncol = length(beta_values))
      
      linear_df = cbind(c(n_sizes[samp_size]),c(y_star_lin_estimate), 
                        matrix(beta_values,nrow = 1))
      # print(linear_df)
      #print(dim(linear_df))
      # Getting the names of the betas
      beta_names = paste0("beta",0:(nrow(betas_lin)-1))
      # Adding the column names
      print(length(c("n_size","y_star_linear",beta_names)))
      colnames(linear_df) <- c("n_size","y_star_linear",beta_names)
      # Returning row of the dataframe
      #print("data frame diensions")
      #print(dim(lasso_df))
      #print(colnames(lasso_df))
      return(as.data.frame(linear_df))
    }
  stopCluster(cluster_1)
  return(linear_calculations_cvg)
}

linear_reg_fixed_var = function(model_data,
                                               n_sizes = c(100, 250, 500, 1000), cores = 3,
                                               m=1000,
                                               randomizer = NA,
                                               file_name = "linear_results"){
  
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
  
  all_x_input_retrieve = function(dat){
    x_columns = is.na(match(colnames(dat), # Retrieve any positions that are y or y_star
                            c("y", "y_star")))
    return(c(colnames(dat)[x_columns]))
  }
  
  # Function to find variable with highest correlation to y
  locate_high_corr <- function(dat){
    # Get a vector of the correlations of input variables with y
    corr_vector = cor(dat[,c(2:ncol(dat))])[1:(ncol(dat) -2),"y"]
    # Get the index of the variable with the highest correlation to y
    var_index = which(corr_vector == max(corr_vector))
    return(var_index)
  }
  
  # Variance for model-assisted estimator, plugin
  var_calculation_ma = function(act_values, est_values, pop_N,
                                first_order_pi, second_order_pi)
  {
    temp = 0
    for(item in 1:length(act_values))
    {
      temp = temp + (1 - first_order_pi) * 
        ((act_values[item] - est_values[item]) / first_order_pi)^2
      
      temp = temp + ((second_order_pi - first_order_pi^2) / (second_order_pi)) * 
        sum(((act_values[item] - est_values[item]) / first_order_pi) * 
              ((act_values[-item] - est_values[-item]) / first_order_pi))
    }
    var_ma = (1 / pop_N^2) * temp
    return(var_ma)
  }
  
  # Construct a matrix of the population data
  pop_x_col = x_input_retrieve(model_data) # get column positions for auxiliary variables of population
  pop_all_x_col = all_x_input_retrieve(model_data)
  # Get the population matrix
  pop_with_out_constant_mat = as.matrix(model_data[,pop_x_col])
  
  pop_with_out_constant_mat = cbind(X0 = rep(1,nrow(pop_with_out_constant_mat)), 
                                    pop_with_out_constant_mat)
  
  # sample sizes controlled by n_sizes
  # Cluster for running models in parallel
  cluster_1 = makeCluster(cores)
  registerDoSNOW(cluster_1)
  # Where the results will be stored. lasso_calculations_cvg is a vector
  linear_calculations_cvg = foreach(samp_size = 1:length(n_sizes), 
                                    .combine = "rbind", 
                                    .inorder = FALSE) %:%
    foreach( iteration = 1:m, .combine = "rbind",.inorder = FALSE) %dopar% {
      
      # setting seed for sampling
      set.seed(iteration + 1000 * (samp_size-1))
      
      
      drawing = sample(nrow(model_data), size= n_sizes[samp_size],replace = FALSE)
      # Getting the sample
      input_sampling = model_data[drawing,] 
      
      y_rand = randomizer(input_sampling[,"y"])
      y_star = y_rand$values
      
      # Getting the sample
      input_sampling = model_data[drawing,] 
      
      
      
      # Populations size
      N_pop = nrow(model_data)
      # First and second order inclusion probs for SRS
      first_order_pi = (n_sizes[samp_size] / N_pop)
      second_order_pi = (n_sizes[samp_size] / N_pop) * 
        ((n_sizes[samp_size] - 1) / (N_pop - 1))
      
      # Randomization components
      alpha =  y_rand$alpha
      eta = y_rand$eta
      gamma = y_rand$gamma
      beta = y_rand$beta
      
      
      input_sampling = cbind(input_sampling, y_star) # Bind y_star column
      
      x_input_col = x_input_retrieve(input_sampling)
      
      # Get the y_star position
      y_input_col = y_star_input_retrieve(input_sampling)
     
      
      # Get all the predictors and the response y
      input_sampling_without_const = as.matrix(input_sampling[,c(x_input_col,
                                                                 y_input_col)])
      
      input_sampling_without_const = cbind(X0 = rep(1,nrow(input_sampling_without_const)), 
                                           input_sampling_without_const)
      
      x_input_col = all_x_input_retrieve(input_sampling_without_const)
      
      # Get the y_star position
      y_input_col = y_star_input_retrieve(input_sampling_without_const)
      
      
      betas_lin =  solve(t(input_sampling_without_const[,x_input_col])  %*%
                           input_sampling_without_const[,x_input_col], tol = NULL) %*% t(input_sampling_without_const[,x_input_col]) %*%
        y_star
      
      
      x_cols = all_x_input_retrieve(input_sampling_without_const)
      
      y_star_lin_sample =  input_sampling_without_const[,x_cols] %*% betas_lin
      
      
      y_star_lin_pop = pop_with_out_constant_mat %*% betas_lin
      
      
      y_star_lin_estimate = mean((y_star - y_star_lin_sample)) + 
        mean(y_star_lin_pop) # Difference estimator result
      
      
      
      linear_var_est = var_calculation_ma(y_star, y_star_lin_sample, 
                                          pop_N = N_pop,
                                          first_order_pi = first_order_pi, 
                                          second_order_pi = second_order_pi)
      
      # Variance of the randomized response method, linear method
      rrt_variance = (1 / N_pop^2) * 
        sum( first_order_pi^(-2) * (gamma/ beta^2) * 
               ((y_star^2 - (eta / beta^2)) *
                  ((gamma / beta^2) + 1)^(-1)) + 
               (eta / beta^2))
      
      # Total variance
      variance =  rrt_variance + linear_var_est
      
      # get upper and lower for confidence interval
      actual_mu = mean(model_data$y)
      upper = y_star_lin_estimate  + 1.96 * sqrt(variance)
      # Get lower bound 
      lower = y_star_lin_estimate - 1.96 * sqrt(variance)
      
      covered_mu = (lower < actual_mu && actual_mu < upper)
      
      
      beta_values = as.vector(t(betas_lin[,1]))
   
      
      linear_df = cbind(c(n_sizes[samp_size]),c(y_star_lin_estimate), 
                        upper, lower, covered_mu, variance, rrt_variance,
                        matrix(beta_values,nrow = 1))
      
      # Getting the names of the betas
      beta_names = paste0("beta",0:(nrow(betas_lin)-1))
      # Adding the column names
      colnames(linear_df) <- c("n_size","y_star_linear", "upper.CI", "lower.CI",
                               "covered_mu","variance","rrt_variance",beta_names)
      # Returning row of the dataframe
      
      return(as.data.frame(linear_df))
    }
  stopCluster(cluster_1)
  return(linear_calculations_cvg)
}

## R lm() function ---------------
linear_reg_model_plain = function(model_data,
         n_sizes = c(100, 250, 500, 1000), cores = 3,
         m=1000,
         randomizer = NA,
         file_name = "linear_results"){
  
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
  
  all_x_input_retrieve = function(dat){
    x_columns = is.na(match(colnames(dat), # Retrieve any positions that are y or y_star
                            c("y", "y_star")))
    return(c(colnames(dat)[x_columns]))
  }
  
  # Variance for model-assisted estimator, plugin
  var_calculation_ma = function(act_values, est_values, pop_N,
                                first_order_pi, second_order_pi)
  {
    temp = 0
    for(item in 1:length(act_values))
    {
      temp = temp + (1 - first_order_pi) * 
        ((act_values[item] - est_values[item]) / first_order_pi)^2
      
      temp = temp + ((second_order_pi - first_order_pi^2) / (second_order_pi)) * 
        sum(((act_values[item] - est_values[item]) / first_order_pi) * 
              ((act_values[-item] - est_values[-item]) / first_order_pi))
    }
    var_ma = (1 / pop_N^2) * temp
    return(var_ma)
  }

  
  # Construct a matrix of the population data
  pop_x_col = x_input_retrieve(model_data) # get column positions for auxiliary variables of population
  pop_all_x_col = all_x_input_retrieve(model_data)
  # Get the population matrix
  pop_with_out_constant_mat = as.data.frame(model_data[,pop_x_col])
  
  pop_with_out_constant_mat = pop_with_out_constant_mat
  
  # sample sizes controlled by n_sizes
  # Cluster for running models in parallel
  cluster_1 = makeCluster(cores)
  registerDoSNOW(cluster_1)
  # Where the results will be stored. lasso_calculations_cvg is a vector
  linear_calculations_cvg = foreach(samp_size = 1:length(n_sizes), 
                                    .combine = "rbind",
                                    .packages = c("sampling"), 
                                    .inorder = FALSE) %:%
    foreach( iteration = 1:m, .combine = "rbind",.inorder = FALSE) %dopar% {
      
      # setting seed for sampling
      set.seed(iteration + 1000 * (samp_size-1))
      drawing = sample(nrow(model_data), size= n_sizes[samp_size],replace = FALSE)
      
      # Getting the sample
      input_sampling = model_data[drawing,] 
      
      y_star = numeric(length(input_sampling[,"y"]))
      
      # Getting the sample
      input_sampling = model_data[drawing,] 
      
      y_rand = randomizer(input_sampling[,"y"])
      y_star = y_rand$values
      
      # Populations size
      N_pop = nrow(model_data)
      # First and second order inclusion probs for SRS
      first_order_pi = (n_sizes[samp_size] / N_pop)
      second_order_pi = (n_sizes[samp_size] / N_pop) * 
        ((n_sizes[samp_size] - 1) / (N_pop - 1))
      
      # Randomization components
      alpha =  y_rand$alpha
      eta = y_rand$eta
      gamma = y_rand$gamma
      beta = y_rand$beta
      
      
      # Get the cross validated model for glmnet
      # Get the auxiliary variables column positions
      
      input_sampling = cbind(input_sampling, y_star) # Bind y_star column
      
      x_input_col = x_input_retrieve(input_sampling)
      
      # Get the y_star position
      y_input_col = y_star_input_retrieve(input_sampling)
      
      
      input_sampling_without_const = 
        data.frame(input_sampling[,c(x_input_col, y_input_col)])
      
      linear_model = lm(y_star ~ .,
                        data = input_sampling_without_const)
      
      
      y_star_lin_sample =  predict(linear_model,input_sampling_without_const)
      
      y_star_lin_pop =  predict(linear_model,pop_with_out_constant_mat)
      
      y_star_lin_est = mean(y_star - y_star_lin_sample) + y_star_lin_pop
      
      beta_values = as.vector(coef(linear_model))
      
      
      linear_df = cbind(c(n_sizes[samp_size]),c(y_star_lin_est), 
                        matrix(beta_values,nrow = 1))
      
      # Getting the names of the betas
      beta_names = paste0("beta",0:(length(beta_values)-1))
      # Adding the column names
      colnames(linear_df) <- c("n_size","y_star_linear",beta_names)
      # Returning row of the dataframe
      return(as.data.frame(linear_df))
    }
  stopCluster(cluster_1)
  return(linear_calculations_cvg)
}


## R lm() function ---------------------------------------------------------------------------------------------------------------
# Get the variance using the RRT and the MA variance
linear_reg_model_plain_var = function(model_data,
                                  n_sizes = c(100, 250, 500, 1000), cores = 3,
                                  m=1000,
                                  randomizer = NA,
                                  file_name = "linear_results"){
  
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
  
  all_x_input_retrieve = function(dat){
    x_columns = is.na(match(colnames(dat), # Retrieve any positions that are y or y_star
                            c("y", "y_star")))
    return(c(colnames(dat)[x_columns]))
  }
  
  
  # Variance for model-assisted estimator, plugin
  var_calculation_ma = function(act_values, est_values, pop_N,
                                first_order_pi, second_order_pi)
  {
    temp = 0
    for(item in 1:length(act_values))
    {
      temp = temp + (1 - first_order_pi) * 
        ((act_values[item] - est_values[item]) / first_order_pi)^2
      
      temp = temp + ((second_order_pi - first_order_pi^2) / (second_order_pi)) * 
        sum(((act_values[item] - est_values[item]) / first_order_pi) * 
              ((act_values[-item] - est_values[-item]) / first_order_pi))
    }
    var_ma = (1 / pop_N^2) * temp
    return(var_ma)
  }
  
  # Construct a matrix of the population data
  pop_x_col = x_input_retrieve(model_data) # get column positions for auxiliary variables of population
  pop_all_x_col = all_x_input_retrieve(model_data)
  # Get the population matrix
  pop_with_out_constant_mat = as.data.frame(model_data[,pop_x_col])
  
  
  # sample sizes controlled by n_sizes
  # Cluster for running models in parallel
  cluster_1 = makeCluster(cores)
  registerDoSNOW(cluster_1)
  # Where the results will be stored. lasso_calculations_cvg is a vector
  linear_calculations_cvg = foreach(samp_size = 1:length(n_sizes), 
                                    .combine = "rbind",
                                    .packages = c("sampling"), 
                                    .inorder = FALSE) %:%
    foreach( iteration = 1:m, .combine = "rbind",.inorder = FALSE) %dopar% {
      
      # setting seed for sampling
      set.seed(iteration + 1000 * (samp_size-1))
      drawing = sample(nrow(model_data), size= n_sizes[samp_size],replace = FALSE)
      
      # Getting the sample
      input_sampling = model_data[drawing,] 
      
      y_star = numeric(length(input_sampling[,"y"]))
      
      # Getting the sample
      input_sampling = model_data[drawing,] 
      
      y_rand = randomizer(input_sampling[,"y"])
      y_star = y_rand$values
      
      # Populations size
      N_pop = nrow(model_data)
      # First and second order inclusion probs for SRS
      first_order_pi = (n_sizes[samp_size] / N_pop)
      second_order_pi = (n_sizes[samp_size] / N_pop) * 
        ((n_sizes[samp_size] - 1) / (N_pop - 1))
      
      # Randomization components
      alpha =  y_rand$alpha
      eta = y_rand$eta
      gamma = y_rand$gamma
      beta = y_rand$beta
      
      
      # Get the cross validated model for glmnet
      # Get the auxiliary variables column positions
      
      input_sampling = cbind(input_sampling, y_star) # Bind y_star column
      
      x_input_col = x_input_retrieve(input_sampling)
      
      # Get the y_star position
      y_input_col = y_star_input_retrieve(input_sampling)
      
      
      input_sampling_without_const = 
        data.frame(input_sampling[,c(x_input_col, y_input_col)])
      
     
      linear_model = lm(y_star ~ .,
                                data = input_sampling_without_const)
      
      
      y_star_lin_sample =  predict(linear_model,input_sampling_without_const)
      
      # Get predictions on the population data
      y_star_linear_pop = predict(linear_model,
                                 pop_with_out_constant_mat)
      # # get the result
      y_star_linear_est = mean(y_star - y_star_lin_sample) + mean(y_star_linear_pop)
      
      linear_var_est = var_calculation_ma(y_star, y_star_lin_sample, 
                                                pop_N = N_pop,
                                                first_order_pi = first_order_pi, 
                                                second_order_pi = second_order_pi)
      
      # Variance of the randomized response method, linear method
      rrt_variance = (1 / N_pop^2) * 
        sum( first_order_pi^(-2) * (gamma/ beta^2) * 
               ((y_star^2 - (eta / beta^2)) *
                  ((gamma / beta^2) + 1)^(-1)) + 
               (eta / beta^2))
      
      # Total variance
      variance =  rrt_variance + linear_var_est
      
      # get upper and lower for confidence interval
      actual_mu = mean(model_data$y)
      upper = y_star_linear_est + 1.96 * sqrt(variance)
      # Get lower bound 
      lower = y_star_linear_est - 1.96 * sqrt(variance)
      
      covered_mu = (lower < actual_mu && actual_mu < upper)
      
      beta_values = as.vector(coef(linear_model))
      
      
      linear_df = cbind(c(n_sizes[samp_size]),c(y_star_linear_est),c(variance), 
                        c(upper), c(lower), c(covered_mu), matrix(beta_values,nrow = 1))
      
      # Getting the names of the betas
      beta_names = paste0("beta",0:(length(beta_values)-1))
      # Adding the column names
      colnames(linear_df) <- c("n_size","y_star_linear","variance","upperCI",
                               "lowerCI","covered_mu",beta_names)
      # Returning row of the dataframe
      return(as.data.frame(linear_df))
    }
  stopCluster(cluster_1)
  return(linear_calculations_cvg)
}
