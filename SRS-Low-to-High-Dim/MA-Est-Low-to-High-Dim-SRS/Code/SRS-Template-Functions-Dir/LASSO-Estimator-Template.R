
lasso_reg_est_base = function(model_data,
                         n_sizes = c(100, 250, 500, 1000), cores = 3,
                         m=1000,
                         numb_lambda = 100,
                         randomizer = NA,
                         boot_yes = FALSE){
  
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
  # Construct a matrix of the population data
  pop_x_col = x_input_retrieve(model_data) # get column positions for auxiliary variables of population
  # Get the population matrix
  pop_with_out_constant_mat = as.matrix(model_data[,pop_x_col])
  
  # sample sizes controlled by n_sizes
  # Cluster for running models in parallel
  cluster_1 = makeCluster(cores)
  registerDoSNOW(cluster_1)
  # Where the results will be stored. lasso_calculations_cvg is a vector
  lasso_calculations_cvg = foreach(samp_size = 1:length(n_sizes), .combine = "rbind",
                                   .packages = ("glmnet"), .inorder = FALSE) %:%
    foreach( iteration = 1:m, .combine = "rbind",.inorder = FALSE) %dopar% {
      # setting seed for sampling
      set.seed(iteration + (1000 * (samp_size-1)))
      
      drawing = sample(nrow(model_data), size= n_sizes[samp_size],replace = FALSE)
      # Getting the sample
      input_sampling = model_data[drawing,] 
      
      y_star = numeric(length(input_sampling[,"y"]))
      
      
      y_star = randomizer(input_sampling[,"y"])
      
      
      # Get the cross validated model for glmnet
      # Get the auxiliary variables column positions
      input_sampling = cbind(input_sampling, y_star) # Bind y_star column
      x_input_col = x_input_retrieve(input_sampling)
      
      # Get the y_star position
      y_input_col = y_star_input_retrieve(input_sampling)
      
      # Get all the predictors and the response y
      input_sampling_without_const = as.matrix(input_sampling[,c(x_input_col,
                                                                 y_input_col)]) 
      x_input_col = x_input_retrieve(input_sampling_without_const)
      
      # Get the y_star position
      y_input_col = y_star_input_retrieve(input_sampling_without_const)
      
      # Building the LASSO model
      #print(colnames(input_sampling_without_const[,x_input_col]))
      #print(colnames(input_sampling_without_const[,x_input_col]))
      lasso_model = cv.glmnet(x = input_sampling_without_const[,x_input_col],
                              y = input_sampling_without_const[,"y_star"],
                              family = "gaussian", alpha = 1,
                              nlambda = numb_lambda,
                              maxit = 10^9,
                              thresh = 1e-8,
                              type.measure = "mse",
                              nfolds = 10)
      # Get the predictions using the best lambda found on sample
      
      y_star_lasso = predict(lasso_model,
                             input_sampling_without_const[,x_input_col],
                             type = "response",
                             s = lasso_model$lambda.min)
      # Get predictions on the population data
      
      y_star_lasso_pop = predict(lasso_model,
                                 pop_with_out_constant_mat,
                                 type = "response",
                                 s = lasso_model$lambda.min)
      # # get the result
      y_star_lasso_est = mean(y_star - y_star_lasso) + mean(y_star_lasso_pop)
      
      # get the number of non-zero coefficients
      lasso_coef = as.vector(coef(lasso_model, s = lasso_model$lambda.min))
      # Organzing to bind as to the dataframe
      lasso_coef = matrix(lasso_coef, nrow = 1, ncol = length(lasso_coef))
      # Constructing the dataframe
      lasso_df = cbind(n_sizes[samp_size],y_star_lasso_est, lasso_coef, 
                       lasso_model$lambda.min)
      # Getting the names of the betas
      beta_names = paste0("beta",0:(length(lasso_coef)-1))
      # Adding the column names
      colnames(lasso_df) <- c("n_size","y_star_lasso",beta_names,"lambda")
      # Returning row of the dataframe
      return(as.data.frame(lasso_df))
    }
  stopCluster(cluster_1)
  
  return(lasso_calculations_cvg)
}


## LASSO regression function with Bootstrap estimation ---------------------------------------------------------------------------------------------------------------
lasso_reg_est = function(model_data,
                                   n_sizes = c(100, 250, 500, 1000), cores = 3,
                                   m=1000,
                                   numb_lambda = 100,
                                   randomizer = NA,
                         boot_yes = FALSE, 
                         num_replicates = 250, # Dr. Mostafa Advised 200 or 500 to cut down on computation time
                         lasso_reg_boot_var_function = lasso_reg_bootstrap.var_est){
  
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
  # Construct a matrix of the population data
  pop_x_col = x_input_retrieve(model_data) # get column positions for auxiliary variables of population
  # Get the population matrix
  pop_with_out_constant_mat = as.matrix(model_data[,pop_x_col])
  
  # sample sizes controlled by n_sizes
  # Cluster for running models in parallel
  cluster_1 = makeCluster(cores)
  registerDoSNOW(cluster_1)
  # Where the results will be stored. lasso_calculations_cvg is a vector
  lasso_calculations_cvg = foreach(samp_size = 1:length(n_sizes), .combine = "rbind",
                                   .packages = ("glmnet"), .inorder = FALSE) %:%
    foreach( iteration = 1:m, .combine = "rbind",.inorder = FALSE) %dopar% {
      # setting seed for sampling
      set.seed(iteration + (1000 * (samp_size-1)))
      
      run_seed = iteration + (1000 * (samp_size-1))
      drawing = sample(nrow(model_data), size= n_sizes[samp_size],replace = FALSE)
      # Getting the sample
      input_sampling = model_data[drawing,] 
      
      
      y_star = numeric(length(input_sampling[,"y"]))
      
      
      y_star = randomizer(input_sampling[,"y"])$values
      
      
      # Get the cross validated model for glmnet
      # Get the auxiliary variables column positions
      input_sampling = cbind(input_sampling, y_star) # Bind y_star column
      x_input_col = x_input_retrieve(input_sampling)
      
      # Get the y_star position
      y_input_col = y_star_input_retrieve(input_sampling)
      
      # Get all the predictors and the response y
      input_sampling_without_const = as.matrix(input_sampling[,c(x_input_col,
                                                                 y_input_col)]) 
      x_input_col = x_input_retrieve(input_sampling_without_const)
      
      # Get the y_star position
      y_input_col = y_star_input_retrieve(input_sampling_without_const)
     
      # Building the LASSO model
      #print(colnames(input_sampling_without_const[,x_input_col]))
      #print(colnames(input_sampling_without_const[,x_input_col]))
      lasso_model = cv.glmnet(x = input_sampling_without_const[,x_input_col],
                              y = input_sampling_without_const[,"y_star"],
                              family = "gaussian", alpha = 1,
                              nlambda = numb_lambda,
                              maxit = 10^9,
                              thresh = 1e-8,
                              type.measure = "mse",
                              nfolds = 10)
      # Get the predictions using the best lambda found on sample
      
      y_star_lasso_samp = predict(lasso_model,
                             input_sampling_without_const[,x_input_col],
                             type = "response",
                             s = lasso_model$lambda.min)
      # Get predictions on the population data
      
      y_star_lasso_pop = predict(lasso_model,
                                 pop_with_out_constant_mat,
                                 type = "response",
                                 s = lasso_model$lambda.min)
      # # get the result
      y_star_lasso_est = mean(y_star - y_star_lasso_samp) + mean(y_star_lasso_pop)
      
      
      
      
      # get the number of non-zero coefficients
      lasso_coef = as.vector(coef(lasso_model, s = lasso_model$lambda.min))
      # Organzing to bind as to the dataframe
      lasso_coef = matrix(lasso_coef, nrow = 1, ncol = length(lasso_coef))
      # Constructing the dataframe
      lasso_df = cbind(n_sizes[samp_size],y_star_lasso_est, lasso_coef, 
                       lasso_model$lambda.min, run_seed)
      # Getting the names of the betas
      beta_names = paste0("beta",0:(length(lasso_coef)-1))
      # Adding the column names
      colnames(lasso_df) <- c("n_size","y_star_lasso",beta_names,"lambda", "run_seed")
      # Returning row of the dataframe
      
      lasso_df = as.data.frame(lasso_df)
      
      ## Get the bootstrapping results
      
      if(boot_yes == TRUE)
      {
        boot_results = lasso_reg_boot_var_function(boot_model_data = model_data,
                                                    samp_data = model_data[drawing,],
                                                   cores = 1,
                                                   boot_m = num_replicates,
                                                   numb_lambda = 100,
                                                    randomizer = randomizer)
        actual_mu = mean(model_data[,"y"])
        boot_results$actual_mu = actual_mu
        boot_results = boot_results %>% # 
          dplyr::select(n_size, boot_var, lower.CI, upper.CI, actual_mu, y_star_boot_est) %>% distinct() %>%
          rename(lower.CI.percent = lower.CI, upper.CI.percent = upper.CI) # percentile bootstrap confidence intervals
        # Looking at lasso_df, thus the y_star_lasso is used from lasso_df
        lasso_df = left_join(lasso_df,
                                           boot_results,
                                           by = "n_size") %>%
          rowwise() %>%
          mutate(mu_star_covered = ifelse(y_star_lasso > lower.CI.percent && # Check that population mean estimate is in the bootstrap 95% conf. interval
                                            y_star_lasso < upper.CI.percent, TRUE, FALSE),
                 mu_covered = ifelse(actual_mu < y_star_lasso  + 1.96 * sqrt(boot_var) && # 95% bootstrap confidence interval coverage (non-percentile) for theactual population mean
                                       actual_mu > y_star_lasso  - 1.96 * sqrt(boot_var), TRUE, FALSE),
                 upper.ci = y_star_lasso  + 1.96 * sqrt(boot_var),
                 lower.ci = y_star_lasso  - 1.96 * sqrt(boot_var)) %>%
          ungroup()
      }
      
      return(lasso_df)
    }
  stopCluster(cluster_1)
  
  return(lasso_calculations_cvg)
}



## LASSO estimator with bootstrapping. It can be broken inton batches.
lasso_reg_est.start.stop = function(model_data,
                         n_sizes = c(100, 250, 500, 1000), cores = 3,
                         m_start = 1,
                         m=1000,
                         numb_lambda = 100,
                         randomizer = NA,
                         boot_yes = FALSE, 
                         num_replicates = 250, # Dr. Mostafa Advised 200 or 500 to cut down on computation time
                         lasso_reg_boot_var_function = lasso_reg_bootstrap.var_est){
  
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
  # Construct a matrix of the population data
  pop_x_col = x_input_retrieve(model_data) # get column positions for auxiliary variables of population
  # Get the population matrix
  pop_with_out_constant_mat = as.matrix(model_data[,pop_x_col])
  
  # sample sizes controlled by n_sizes
  # Cluster for running models in parallel
  cluster_1 = makeCluster(cores)
  registerDoSNOW(cluster_1)
  # Where the results will be stored. lasso_calculations_cvg is a vector
  lasso_calculations_cvg = foreach(samp_size = 1:length(n_sizes), .combine = "rbind",
                                   .packages = ("glmnet"), .inorder = FALSE) %:%
    foreach( iteration = m_start:m, .combine = "rbind",.inorder = FALSE) %dopar% {
      # setting seed for sampling
      set.seed(iteration + (1000 * (samp_size-1)))
      
      run_seed = iteration + (1000 * (samp_size-1))
      drawing = sample(nrow(model_data), size= n_sizes[samp_size],replace = FALSE)
      # Getting the sample
      input_sampling = model_data[drawing,] 
      
      
      y_star = numeric(length(input_sampling[,"y"]))
      
      
      y_star = randomizer(input_sampling[,"y"])$values
      
      
      # Get the cross validated model for glmnet
      # Get the auxiliary variables column positions
      input_sampling = cbind(input_sampling, y_star) # Bind y_star column
      x_input_col = x_input_retrieve(input_sampling)
      
      # Get the y_star position
      y_input_col = y_star_input_retrieve(input_sampling)
      
      # Get all the predictors and the response y
      input_sampling_without_const = as.matrix(input_sampling[,c(x_input_col,
                                                                 y_input_col)]) 
      x_input_col = x_input_retrieve(input_sampling_without_const)
      
      # Get the y_star position
      y_input_col = y_star_input_retrieve(input_sampling_without_const)
      
      # Building the LASSO model
      #print(colnames(input_sampling_without_const[,x_input_col]))
      #print(colnames(input_sampling_without_const[,x_input_col]))
      lasso_model = cv.glmnet(x = input_sampling_without_const[,x_input_col],
                              y = input_sampling_without_const[,"y_star"],
                              family = "gaussian", alpha = 1,
                              nlambda = numb_lambda,
                              maxit = 10^9,
                              thresh = 1e-8,
                              type.measure = "mse",
                              nfolds = 10)
      # Get the predictions using the best lambda found on sample
      
      y_star_lasso_samp = predict(lasso_model,
                                  input_sampling_without_const[,x_input_col],
                                  type = "response",
                                  s = lasso_model$lambda.min)
      # Get predictions on the population data
      
      y_star_lasso_pop = predict(lasso_model,
                                 pop_with_out_constant_mat,
                                 type = "response",
                                 s = lasso_model$lambda.min)
      # # get the result
      y_star_lasso_est = mean(y_star - y_star_lasso_samp) + mean(y_star_lasso_pop)
      
      
      
      
      # get the number of non-zero coefficients
      lasso_coef = as.vector(coef(lasso_model, s = lasso_model$lambda.min))
      # Organzing to bind as to the dataframe
      lasso_coef = matrix(lasso_coef, nrow = 1, ncol = length(lasso_coef))
      # Constructing the dataframe
      lasso_df = cbind(n_sizes[samp_size],y_star_lasso_est, lasso_coef, 
                       lasso_model$lambda.min, run_seed)
      # Getting the names of the betas
      beta_names = paste0("beta",0:(length(lasso_coef)-1))
      # Adding the column names
      colnames(lasso_df) <- c("n_size","y_star_lasso",beta_names,"lambda", "run_seed")
      # Returning row of the dataframe
      
      lasso_df = as.data.frame(lasso_df)
      
      ## Get the bootstrapping results
      
      if(boot_yes == TRUE)
      {
        boot_results = lasso_reg_boot_var_function(boot_model_data = model_data,
                                                   samp_data = model_data[drawing,],
                                                   cores = 1,
                                                   boot_m = num_replicates,
                                                   numb_lambda = 100,
                                                   randomizer = randomizer)
        actual_mu = mean(model_data[,"y"])
        boot_results$actual_mu = actual_mu
        boot_results = boot_results %>% # 
          dplyr::select(n_size, boot_var, lower.CI, upper.CI, actual_mu, y_star_boot_est) %>% distinct() %>%
          rename(lower.CI.percent = lower.CI, upper.CI.percent = upper.CI) # percentile bootstrap confidence intervals
        # Looking at lasso_df, thus the y_star_lasso is used from lasso_df
        lasso_df = left_join(lasso_df,
                             boot_results,
                             by = "n_size") %>%
          rowwise() %>%
          mutate(mu_star_covered = ifelse(y_star_lasso > lower.CI.percent && # Check that population mean estimate is in the bootstrap 95% conf. interval
                                            y_star_lasso < upper.CI.percent, TRUE, FALSE),
                 mu_covered = ifelse(actual_mu < y_star_lasso  + 1.96 * sqrt(boot_var) && # 95% bootstrap confidence interval coverage (non-percentile) for theactual population mean
                                       actual_mu > y_star_lasso  - 1.96 * sqrt(boot_var), TRUE, FALSE),
                 upper.ci = y_star_lasso  + 1.96 * sqrt(boot_var),
                 lower.ci = y_star_lasso  - 1.96 * sqrt(boot_var)) %>%
          ungroup()
      }
      
      return(lasso_df)
    }
  stopCluster(cluster_1)
  
  return(lasso_calculations_cvg)
}





# LASSO Variance Estimation for Plug-in Estimator -------------------------------------
lasso_reg_var_est = function(model_data,
                         n_sizes = c(100, 250, 500, 1000), cores = 3,
                         m=1000,
                         numb_lambda = 100,
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
  
  # Variance for model-assisted estimator, plugin
  var_calculation_ma = function(act_values, est_values, pop_N,
                                first_order_pi, second_order_pi)
  {
    temp = 0
    for(item in 1:length(act_values))
    {
      temp = temp + (1 - first_order_pi) * (1 / first_order_pi^2) * 
        ((act_values[item] - est_values[item]))^2
      
      temp = temp + ((second_order_pi - first_order_pi^2) / (second_order_pi)) * 
        sum(((act_values[item] - est_values[item]) / first_order_pi) * 
              ((act_values[-item] - est_values[-item]) / first_order_pi))
    }
    var_ma = (1 / pop_N^2) * temp
    return(var_ma)
  }
  
  # Construct a matrix of the population data
  pop_x_col = x_input_retrieve(model_data) # get column positions for auxiliary variables of population
  # Get the population matrix
  pop_with_out_constant_mat = as.matrix(model_data[,pop_x_col])
  
  # sample sizes controlled by n_sizes
  # Cluster for running models in parallel
  cluster_1 = makeCluster(cores)
  registerDoSNOW(cluster_1)
  # Where the results will be stored. lasso_calculations_cvg is a vector
  lasso_calculations_cvg = foreach(samp_size = 1:length(n_sizes), .combine = "rbind",
                                   .packages = ("glmnet"), .inorder = FALSE) %:%
    foreach( iteration = 1:m, .combine = "rbind",.inorder = FALSE) %dopar% {
      # setting seed for sampling
      set.seed(iteration + (1000 * (samp_size-1)))
      run_seed = iteration + (1000 * (samp_size-1))
      
      drawing = sample(nrow(model_data), size= n_sizes[samp_size],replace = FALSE)
      # Getting the sample
      input_sampling = model_data[drawing,] 
      
      y_star = numeric(length(input_sampling[,"y"]))
      
      # Randomization for y_star and the components of randomization
      y_rand = randomizer(input_sampling[,"y"])
      y_star = y_rand$values
      
      # Populations size
      N_pop = nrow(model_data)
      # First and second order inclusion probs for SRS
      first_order_pi = (n_sizes[samp_size] / N_pop)
      second_order_pi = (n_sizes[samp_size] / N_pop) * 
        ((n_sizes[samp_size] - 1) / (N_pop - 1))
      
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
      
      # Get all the predictors and the response y
      input_sampling_without_const = as.matrix(input_sampling[,c(x_input_col,
                                                                 y_input_col)]) 
      x_input_col = x_input_retrieve(input_sampling_without_const)
      
      # Get the y_star position
      y_input_col = y_star_input_retrieve(input_sampling_without_const)
      
      # Building the LASSO model
      
      lasso_model = cv.glmnet(x = input_sampling_without_const[,x_input_col],
                              y = input_sampling_without_const[,"y_star"],
                              family = "gaussian", alpha = 1,
                              nlambda = numb_lambda,
                              maxit = 10^9,
                              thresh = 1e-8,
                              type.measure = "mse",
                              nfolds = 10)
      # Get the predictions using the best lambda found on sample
      
      y_star_lasso = predict(lasso_model,
                             input_sampling_without_const[,x_input_col],
                             type = "response",
                             s = lasso_model$lambda.min)
      
      # Get predictions on the population data
      y_star_lasso_pop = predict(lasso_model,
                                 pop_with_out_constant_mat,
                                 type = "response",
                                 s = lasso_model$lambda.min)
      # # get the result
      y_star_lasso_est = mean(y_star - y_star_lasso) + mean(y_star_lasso_pop)
      
      
      # # get the result
      y_star_lasso_var_est = 0
      
      y_star_lasso_var_est = var_calculation_ma(y_star, y_star_lasso, 
                                                pop_N = N_pop,
                                    first_order_pi = first_order_pi, 
                                    second_order_pi = second_order_pi)
      
      # Variance of the randomized response method, linear method
      rrt_variance = (1 / N_pop^2) * 
        sum( first_order_pi^(-2) * ((gamma/ beta^2) * 
               ((y_star^2 - (eta / beta^2)) *
                  ((gamma / beta^2) + 1)^(-1)) + 
               (eta / beta^2)))
      
      # Total variance
      variance =  rrt_variance + y_star_lasso_var_est
      
      # get upper and lower for confidence interval
      actual_mu = mean(model_data$y)
      upper = y_star_lasso_est + 1.96 * sqrt(variance)
      # Get lower bound 
      lower = y_star_lasso_est - 1.96 * sqrt(variance)
      
      covered_mu = (lower < actual_mu && actual_mu < upper)
      
      # get the number of non-zero coefficients
      lasso_coef = as.vector(coef(lasso_model, s = lasso_model$lambda.min))
      # Organzing to bind as to the dataframe
      lasso_coef = matrix(lasso_coef, nrow = 1, ncol = length(lasso_coef))
      # Constructing the dataframe
      lasso_df = cbind(n_sizes[samp_size], y_star_lasso_est,
                       rrt_variance, variance, lower, upper, covered_mu, 
                       lasso_coef, 
                       lasso_model$lambda.min, run_seed)
      # Getting the names of the betas
      beta_names = paste0("beta",0:(length(lasso_coef)-1))
      # Adding the column names
      colnames(lasso_df) <- c("n_size","y_star_lasso","rrt_var","variance",
                              "lower.CI", "upper.CI", "covered_mu",
                              beta_names,"lambda", "run_seed")
      
      
      # Returning row of the dataframe
      return(as.data.frame(lasso_df))
    }
  stopCluster(cluster_1)
  return(lasso_calculations_cvg)
}




## Function for Bootstrap Variance Estiamtion ------------------------------
lasso_reg_bootstrap.var_est = function(boot_model_data,
                                       samp_data,
                             cores = 1,
                             boot_m=1000,
                             numb_lambda = 100,
                             randomizer = NA){
  library(foreach)
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
  pop_x_col = x_input_retrieve(boot_model_data) # get column positions for auxiliary variables of population
  # Get the population matrix
  pop_with_out_constant_mat = as.matrix(boot_model_data[,pop_x_col])
  
  # sample sizes controlled by n_sizes
  # Cluster for running models in parallel
 
  # Where the results will be stored. lasso_calculations_cvg is a vector
  lasso_calculations_cvg =
    foreach( iteration = 1:boot_m, 
             .packages = c("glmnet","tidyverse"), .combine = "rbind",.inorder = FALSE) %do% {
      # setting seed for sampling
      set.seed(iteration + (1000 * (nrow(samp_data)-1)))
      
      # Bootstrap variance estimation 
      drawing = sample(nrow(samp_data), size= nrow(samp_data),replace = TRUE)
      # Getting the sample
      input_sampling = samp_data[drawing,] 
      
      y_star = numeric(length(input_sampling[,"y"]))
      
      # Randomization for y_star and the components of randomization
      y_rand = randomizer(input_sampling[,"y"])
      y_star = y_rand$values
      
      # Populations size
      N_pop = nrow(boot_model_data)
      
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
      
      # Get all the predictors and the response y
      input_sampling_without_const = as.matrix(input_sampling[,c(x_input_col,
                                                                 y_input_col)]) 
      x_input_col = x_input_retrieve(input_sampling_without_const)
      
      # Get the y_star position
      y_input_col = y_star_input_retrieve(input_sampling_without_const)
      
      # Building the LASSO model
      
      lasso_model = cv.glmnet(x = input_sampling_without_const[,x_input_col],
                              y = input_sampling_without_const[,"y_star"],
                              family = "gaussian", alpha = 1,
                              nlambda = numb_lambda,
                              maxit = 10^9,
                              thresh = 1e-8,
                              type.measure = "mse",
                              nfolds = 10)
      # Get the predictions using the best lambda found on sample
      
      y_star_lasso = predict(lasso_model,
                             input_sampling_without_const[,x_input_col],
                             type = "response",
                             s = lasso_model$lambda.min)
      
      # Get predictions on the population data
      y_star_lasso_pop = predict(lasso_model,
                                 pop_with_out_constant_mat,
                                 type = "response",
                                 s = lasso_model$lambda.min)
      # # get the result
      y_star_lasso_est = mean(y_star - y_star_lasso) + mean(y_star_lasso_pop)
      
      
      # get the number of non-zero coefficients
      lasso_coef = as.vector(coef(lasso_model, s = lasso_model$lambda.min))
      # Organzing to bind as to the dataframe
      lasso_coef = matrix(lasso_coef, nrow = 1, ncol = length(lasso_coef))
      # Constructing the dataframe
      actual_mu = mean(boot_model_data[,"y"])
      lasso_df = cbind(nrow(samp_data), y_star_lasso_est,
                       lasso_coef, 
                       lasso_model$lambda.min, actual_mu)
      # Getting the names of the betas
      beta_names = paste0("beta",0:(length(lasso_coef)-1))
      # Adding the column names
      
      colnames(lasso_df) <- c("n_size","y_star_lasso_boot",beta_names,"lambda",
                              "actual_mu")
      
      
      # Returning row of the dataframe
      return(as.data.frame(lasso_df))
    }

  lasso_calculations_cvg_bootstrap = lasso_calculations_cvg %>% 
    group_by(n_size) %>%
    summarise(boot_var = var(y_star_lasso_boot),
              lower.CI = quantile(y_star_lasso_boot, probs = 0.025),
              upper.CI = quantile(y_star_lasso_boot, probs = 0.975),
              y_star_boot_est = mean(y_star_lasso_boot)) %>%
    dplyr::select(n_size, boot_var, lower.CI, upper.CI, y_star_boot_est)
  lasso_calculations_cvg = left_join(lasso_calculations_cvg,
                                     lasso_calculations_cvg_bootstrap,
                                     by = "n_size")
  lasso_calculations_cvg = lasso_calculations_cvg
  return(lasso_calculations_cvg)
}

### Mash Variance -------
# Mashregi implementation of the bootstrap variance estimator
lasso_reg_bootstrap.var_est.mash = function(boot_model_data,
                                       samp_data,
                                       cores = 1,
                                       boot_m=1000,
                                       numb_lambda = 100,
                                       randomizer = NA){
  library(foreach)
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
  pop_x_col = x_input_retrieve(boot_model_data) # get column positions for auxiliary variables of population
  # Get the population matrix
  pop_with_out_constant_mat = as.matrix(boot_model_data[,pop_x_col])
  
  # sample sizes controlled by n_sizes
  # Cluster for running models in parallel
  
  # Where the results will be stored. lasso_calculations_cvg is a vector
  lasso_calculations_cvg =
    foreach( iteration = 1:boot_m, 
             .packages = c("glmnet","tidyverse"), .combine = "rbind",.inorder = FALSE) %do% {
               # setting seed for sampling
               set.seed(iteration + (1000 * (nrow(samp_data)-1)))
               
               # Store the pseudo population indices
               pseudo_pop_indices = NA
               # Part from duplicating the sample
               if((nrow(pop_with_out_constant_mat) / nrow(samp_data) %% 2) == 0 )
               {
                 pseudo_pop_indices = rep(1:nrow(samp_data), floor(nrow(pop_with_out_constant_mat) / nrow(samp_data)))
               }
               else
               {
                 pseudo_U_1 = rep(1:nrow(samp_data), floor(nrow(pop_with_out_constant_mat) / nrow(samp_data)))
                 
                 pseudo_U_2 = sample(1:nrow(samp_data), 
                                     size = nrow(pop_with_out_constant_mat) - nrow(samp_data) * 
                                       floor(nrow(pop_with_out_constant_mat) / nrow(samp_data)), replace = FALSE )
                 
                 pseudo_pop_indices = c(pseudo_U_1, pseudo_U_2)
               }
               
               pseudo_pop = samp_data[pseudo_pop_indices,]
               
               # Bootstrap variance estimation 
               drawing = sample(nrow(pseudo_pop), size= nrow(samp_data),replace = TRUE)
               # Getting the sample
               input_sampling = pseudo_pop[drawing,] 
               
               y_star = numeric(length(input_sampling[,"y"]))
               
               # Randomization for y_star and the components of randomization
               y_rand = randomizer(input_sampling[,"y"])
               y_star = y_rand$values
               
               # Populations size
               N_pop = nrow(boot_model_data)
               
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
               
               # Get all the predictors and the response y
               input_sampling_without_const = as.matrix(input_sampling[,c(x_input_col,
                                                                          y_input_col)]) 
               x_input_col = x_input_retrieve(input_sampling_without_const)
               
               # Get the y_star position
               y_input_col = y_star_input_retrieve(input_sampling_without_const)
               
               # Building the LASSO model
               
               lasso_model = cv.glmnet(x = input_sampling_without_const[,x_input_col],
                                       y = input_sampling_without_const[,"y_star"],
                                       family = "gaussian", alpha = 1,
                                       nlambda = numb_lambda,
                                       maxit = 10^9,
                                       thresh = 1e-8,
                                       type.measure = "mse",
                                       nfolds = 10)
               # Get the predictions using the best lambda found on sample
               
               y_star_lasso = predict(lasso_model,
                                      input_sampling_without_const[,x_input_col],
                                      type = "response",
                                      s = lasso_model$lambda.min)
               
               # Get predictions on the population data
               y_star_lasso_pop = predict(lasso_model,
                                          pop_with_out_constant_mat,
                                          type = "response",
                                          s = lasso_model$lambda.min)
               # # get the result
               y_star_lasso_est = mean(y_star - y_star_lasso) + mean(y_star_lasso_pop)
               
               
               # get the number of non-zero coefficients
               lasso_coef = as.vector(coef(lasso_model, s = lasso_model$lambda.min))
               # Organzing to bind as to the dataframe
               lasso_coef = matrix(lasso_coef, nrow = 1, ncol = length(lasso_coef))
               # Constructing the dataframe
               actual_mu = mean(boot_model_data[,"y"])
               lasso_df = cbind(nrow(samp_data), y_star_lasso_est,
                                lasso_coef, 
                                lasso_model$lambda.min, actual_mu)
               # Getting the names of the betas
               beta_names = paste0("beta",0:(length(lasso_coef)-1))
               # Adding the column names
               
               colnames(lasso_df) <- c("n_size","y_star_lasso_boot",beta_names,"lambda",
                                       "actual_mu")
               
               
               # Returning row of the dataframe
               return(as.data.frame(lasso_df))
             }
  
  lasso_calculations_cvg_bootstrap = lasso_calculations_cvg %>% 
    group_by(n_size) %>%
    summarise(boot_var = var(y_star_lasso_boot),
              lower.CI = quantile(y_star_lasso_boot, probs = 0.025),
              upper.CI = quantile(y_star_lasso_boot, probs = 0.975),
              y_star_boot_est = mean(y_star_lasso_boot)) %>%
    dplyr::select(n_size, boot_var, lower.CI, upper.CI, y_star_boot_est)
  lasso_calculations_cvg = left_join(lasso_calculations_cvg,
                                     lasso_calculations_cvg_bootstrap,
                                     by = "n_size")
  lasso_calculations_cvg = lasso_calculations_cvg
  return(lasso_calculations_cvg)
}

