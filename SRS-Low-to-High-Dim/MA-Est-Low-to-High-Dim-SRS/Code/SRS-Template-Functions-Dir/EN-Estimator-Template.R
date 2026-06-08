elastic_net_reg_model_results = function(model_data, n_sizes = c(100,250,500,1000),
                                         cores = 3, linear_cores = 1,
                                         m = 1000,
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
  # Construct a matrix of the population data
  pop_x_col = x_input_retrieve(model_data) # get column positions for auxiliary variables of population
  # Get the population matrix
  pop_with_out_constant = as.matrix(model_data[,pop_x_col])
  
  # Caret function for elastic net
  elastic_train = function(input_dat, sampling_weights = NA)
  {
    # Train control for elastic net
    elastic_control = trainControl(method="cv",
                                   number = 5)
    # The parameter grid
    elastic_grid = expand.grid(alpha = seq(0.1,0.9,0.1), # Avoid pure ridge or lasso regression
                               lambda = seq(4,0,length.out = 100))
    elastic_net_model = NA
    # Building model
    if(is.na(sampling_weights)==TRUE)
    {
      elastic_net_model = train(y_star ~.,
                                data = input_dat,
                                method = "glmnet",
                                trControl = elastic_control,
                                tuneGrid = elastic_grid,
                                verbose=FALSE)
    }
    else{
      elastic_net_model = train(y_star ~.,
                                data = input_dat,
                                method = "glmnet",
                                trControl = elastic_control,
                                tuneGrid = elastic_grid,
                                weights = sampling_weights,
                                verbose=FALSE)
    }
    return(elastic_net_model)
  }
  
  # sample sizes controlled by n_sizes
  # Cluster for running models in parallel
  cluster_1 = makeCluster(cores)#, 
  #outfile = "./Code/output/test_new_elastic_framework.txt")
  registerDoSNOW(cluster_1)
  # Where the results will be stored. lasso_calculations_cvg is a vector
  elastic_calculations = foreach(samp_size = 1:length(n_sizes), .combine = "rbind",
                                 .packages = c("glmnet", "caret"), .inorder = FALSE) %:%
    foreach( iteration = 1:m, .combine = "rbind",.inorder = FALSE) %dopar% {
      # setting seed for sampling
      set.seed(iteration + (1000 * (samp_size-1)))
      
      drawing = sample(nrow(model_data), size= n_sizes[samp_size],replace = FALSE)
      # Getting the sample
      input_sampling = model_data[drawing,] 
      
      y_star = numeric(length(input_sampling[,"y"]))
      
      
      y_rand = randomizer(input_sampling[,"y"]) # Utilizaing randomizer
      y_star = y_rand$values
      
      
      
      # Randomization components
      alpha =  y_rand$alpha
      eta = y_rand$eta
      gamma = y_rand$gamma
      beta = y_rand$beta
      
      
      # Get the auxiliary variables column positions
      input_sampling = cbind(input_sampling, y_star) # Bind y_star column
      x_input_col = x_input_retrieve(input_sampling)
      
      # Get the y_star position
      y_input_col = y_star_input_retrieve(input_sampling)
      
      # Get all the predictors and the response y
      # Keeping as a dataframe
      input_sampling_without_const = as.matrix(input_sampling[,c(x_input_col,
                                                                 y_input_col)])
      x_input_col = x_input_retrieve(input_sampling_without_const)
      
      # Get the y_star position
      y_input_col = y_star_input_retrieve(input_sampling_without_const)
      
      
      # Building the Elastic Net model
      elastic_model = elastic_train(input_sampling_without_const)
      
      # Elastic net predictions for the sample
      y_star_elastic_sample = predict.train(elastic_model,
                                            newdata = input_sampling_without_const)

      
      # Elastic net predictions for the population
      y_star_elastic_pop = predict.train(elastic_model,
                                         newdata = pop_with_out_constant)
      # 
      y_star_elastic_estimate =  mean(y_star - y_star_elastic_sample) + 
        mean(y_star_elastic_pop)
      
      
      
      
      # Constructing elastic net data frame
      elastic_net_df = data.frame(n_sizes[samp_size], y_star_elastic_estimate,
                                  elastic_model$bestTune$alpha,
                                  elastic_model$bestTune$lambda)
      
      # Adding the column names
      colnames(elastic_net_df) <- c("n_size","y_star_elastic",
                                    "alpha", "lambda")
      # Returning row of the dataframe
      return(elastic_net_df)
    }
  stopCluster(cluster_1)
  return(elastic_calculations)
}

## Variance Calculations ------------------

elastic_net_var_results = function(model_data, n_sizes = c(100,250,500,1000),
                                         cores = 3, linear_cores = 1,
                                         m = 1000,
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
  # Construct a matrix of the population data
  pop_x_col = x_input_retrieve(model_data) # get column positions for auxiliary variables of population
  # Get the population matrix
  pop_with_out_constant = as.matrix(model_data[,pop_x_col])
  
  # Caret function for elastic net
  elastic_train = function(input_dat, sampling_weights = NA)
  {
    # Train control for elastic net
    elastic_control = trainControl(method="cv",
                                   number = 5)
    # The parameter grid
    elastic_grid = expand.grid(alpha = seq(0.1,0.9,0.1), # Avoid pure ridge or lasso regression
                               lambda = seq(4,0,length.out = 100))
    elastic_net_model = NA
    # Building model
    if(is.na(sampling_weights)==TRUE)
    {
      elastic_net_model = train(y_star ~.,
                                data = input_dat,
                                method = "glmnet",
                                trControl = elastic_control,
                                tuneGrid = elastic_grid,
                                verbose=FALSE)
    }
    else{
      elastic_net_model = train(y_star ~.,
                                data = input_dat,
                                method = "glmnet",
                                trControl = elastic_control,
                                tuneGrid = elastic_grid,
                                weights = sampling_weights,
                                verbose=FALSE)
    }
    return(elastic_net_model)
  }
  
  # Variance for model-assisted estimator, plugin
  var_calculation_ma = function(act_values, est_values, pop_N,
                                first_order_pi, second_order_pi)
  {
    temp = 0
    for(item in 1:length(act_values))
    {
      temp = temp + (1 - first_order_pi) * (1 / first_order_pi)^2 * 
        ((act_values[item] - est_values[item]))^2
      
      temp = temp + ((second_order_pi - first_order_pi^2) / (second_order_pi)) * 
        sum(((act_values[item] - est_values[item]) / first_order_pi) * 
              ((act_values[-item] - est_values[-item]) / first_order_pi))
    }
    var_ma = (1 / pop_N^2) * temp
    return(var_ma)
  }
  
  
  # sample sizes controlled by n_sizes
  # Cluster for running models in parallel
  cluster_1 = makeCluster(cores)#, 
  #outfile = "./Code/output/test_new_elastic_framework.txt")
  registerDoSNOW(cluster_1)
  # Where the results will be stored. lasso_calculations_cvg is a vector
  elastic_calculations = foreach(samp_size = 1:length(n_sizes), .combine = "rbind",
                                 .packages = c("glmnet", "caret"), .inorder = FALSE) %:%
    foreach( iteration = 1:m, .combine = "rbind",.inorder = FALSE) %dopar% {
      # setting seed for sampling
      set.seed(iteration + (1000 * (samp_size-1)))
      
      drawing = sample(nrow(model_data), size= n_sizes[samp_size],replace = FALSE)
      # Getting the sample
      input_sampling = model_data[drawing,] 
      
      y_star = numeric(length(input_sampling[,"y"]))
      
      
      y_rand = randomizer(input_sampling[,"y"]) # Utilizaing randomizer
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
      
      
      # Get the auxiliary variables column positions
      input_sampling = cbind(input_sampling, y_star) # Bind y_star column
      x_input_col = x_input_retrieve(input_sampling)
      
      # Get the y_star position
      y_input_col = y_star_input_retrieve(input_sampling)
      
      # Get all the predictors and the response y
      # Keeping as a dataframe
      input_sampling_without_const = as.matrix(input_sampling[,c(x_input_col,
                                                                 y_input_col)])
      x_input_col = x_input_retrieve(input_sampling_without_const)
      
      # Get the y_star position
      y_input_col = y_star_input_retrieve(input_sampling_without_const)
      
      
      # Building the Elastic Net model
      elastic_model = elastic_train(input_sampling_without_const)
      
      # Elastic net predictions for the sample
      y_star_elastic_sample = predict.train(elastic_model,
                                            newdata = input_sampling_without_const)
      
      # get the result
      elastic_var_est = var_calculation_ma(y_star, y_star_elastic_sample, 
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
      variance =  rrt_variance + elastic_var_est
      
      
      
      # Elastic net predicitons for the population
      y_star_elastic_pop = predict.train(elastic_model,
                                         newdata = pop_with_out_constant)
      # 
      y_star_elastic_estimate =  mean(y_star - y_star_elastic_sample) + 
        mean(y_star_elastic_pop)
      
      
      # get upper and lower for confidence interval
      actual_mu = mean(model_data$y)
      upper = y_star_elastic_estimate + 1.96 * sqrt(variance)
      # Get lower bound 
      lower = y_star_elastic_estimate - 1.96 * sqrt(variance)
      
      covered_mu = (lower < actual_mu && actual_mu < upper)
      
      # Constructing elastic net data frame
      elastic_net_df = data.frame(n_sizes[samp_size], y_star_elastic_estimate,
                                  elastic_model$bestTune$alpha,
                                  elastic_model$bestTune$lambda,
                                  variance, rrt_variance, upper, lower, covered_mu)
      
      # Adding the column names
      colnames(elastic_net_df) <- c("n_size","y_star_elastic", "alpha",
                                    "lambda",
                                    "variance", "rrt_variance",
                                    "upperCI", "lowerCI", "covered_mu")
      # Returning row of the dataframe
      return(elastic_net_df)
    }
  stopCluster(cluster_1)
  return(elastic_calculations)
}


### glmnet only for elastic net CV

elastic_net_var_results_gnet = function(model_data, n_sizes = c(100,250,500,1000),
                                        cores = 3, linear_cores = 1,
                                        m = 1000,
                                        randomizer = NA,
                                        alpha_vec = seq(0.1,0.9,by = 0.1)){
  
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
  pop_with_out_constant = as.matrix(model_data[,pop_x_col])
  
  # Caret function for elastic net
  # elastic_train = function(input_dat, sampling_weights = NA)
  # {
  #   # Train control for elastic net
  #   elastic_control = trainControl(method="cv",
  #                                  number = 5)
  #   # The parameter grid
  #   elastic_grid = expand.grid(alpha = seq(0.1,0.9,0.1), # Avoid pure ridge or lasso regression
  #                              lambda = seq(4,0,length.out = 100))
  #   elastic_net_model = NA
  #   # Building model
  #   if(is.na(sampling_weights)==TRUE)
  #   {
  #     elastic_net_model = train(y_star ~.,
  #                               data = input_dat,
  #                               method = "glmnet",
  #                               trControl = elastic_control,
  #                               tuneGrid = elastic_grid,
  #                               verbose=FALSE)
  #   }
  #   else{
  #     elastic_net_model = train(y_star ~.,
  #                               data = input_dat,
  #                               method = "glmnet",
  #                               trControl = elastic_control,
  #                               tuneGrid = elastic_grid,
  #                               weights = sampling_weights,
  #                               verbose=FALSE)
  #   }
  #   return(elastic_net_model)
  # }
  
  elastic_cv_ver_2 = function(alphas = seq(0.1,0.9,by = 0.1),
                              input_x, input_y, weights = NULL, num_folds=10){
    input_x = as.matrix(input_x)
    input_y = as.matrix(input_y)
    elastic_mod = list()
    # mse = c()
    mse_2 = c()
    for(hyper_param in 1:length(alphas)){
      fold_id = 1:nrow(input_x) %% num_folds + 1
      elastic_mod[[hyper_param]] = cv.glmnet(x=input_x, y = input_y, alpha = alphas[hyper_param],
                                             foldid = fold_id, weights = weights)
      # mse[hyper_param] = 1 / length(as.vector(y)) *
      #   sum(input_y - predict(elastic_mod[[hyper_param]],
      #                         input_x, s= "lambda.min"))^2
      min_index = which.min(elastic_mod[[hyper_param]]$cvm)
      mse_2[hyper_param] = elastic_mod[[hyper_param]]$cvm[min_index]
    }
    best_model_index = which.min(mse_2)
    # print(alphas[best_model_index])
    # print(mse)
    # print(mse_2)
    return(list("model" = elastic_mod[[best_model_index]],
                "alpha" = alphas[best_model_index]))
  }
  
  # Variance for model-assisted estimator, plugin
  var_calculation_ma = function(act_values, est_values, pop_N,
                                first_order_pi, second_order_pi)
  {
    temp = 0
    for(item in 1:length(act_values))
    {
      temp = temp + (1 - first_order_pi) * (1 / first_order_pi)^2 *
        ((act_values[item] - est_values[item]))^2
      
      temp = temp + ((second_order_pi - first_order_pi^2) / (second_order_pi)) *
        sum(((act_values[item] - est_values[item]) / first_order_pi) *
              ((act_values[-item] - est_values[-item]) / first_order_pi))
    }
    var_ma = (1 / pop_N^2) * temp
    return(var_ma)
  }
  
  
  # sample sizes controlled by n_sizes
  # Cluster for running models in parallel
  cluster_1 = makeCluster(cores)#,
  #outfile = "./Code/output/test_new_elastic_framework.txt")
  registerDoSNOW(cluster_1)
  # Where the results will be stored. lasso_calculations_cvg is a vector
  elastic_calculations = foreach(samp_size = 1:length(n_sizes), .combine = "rbind",
                                 .packages = c("glmnet", "caret"), .inorder = FALSE) %:%
    foreach( iteration = 1:m, .combine = "rbind",.inorder = FALSE) %dopar% {
      # setting seed for sampling
      run_seed = iteration + (1000 * (samp_size-1))
      set.seed(run_seed)
      
      drawing = sample(nrow(model_data), size= n_sizes[samp_size],replace = FALSE)
      # Getting the sample
      input_sampling = model_data[drawing,]
      
      y_star = numeric(length(input_sampling[,"y"]))
      
      
      y_rand = randomizer(input_sampling[,"y"]) # Utilizaing randomizer
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
      
      
      # Get the auxiliary variables column positions
      input_sampling = cbind(input_sampling, y_star) # Bind y_star column
      x_input_col = x_input_retrieve(input_sampling)
      
      # Get the y_star position
      y_input_col = y_star_input_retrieve(input_sampling)
      
      # Get all the predictors and the response y
      # Keeping as a dataframe
      input_sampling_without_const = as.matrix(input_sampling[,c(x_input_col,
                                                                 y_input_col)])
      x_input_col = x_input_retrieve(input_sampling_without_const)
      
      # Get the y_star position
      y_input_col = y_star_input_retrieve(input_sampling_without_const)
      
      
      # Building the Elastic Net model
      elastic_model_obj = elastic_cv_ver_2(input_x = input_sampling_without_const[,x_input_col],
                                           input_y = input_sampling_without_const[,y_input_col],
                                           alphas = alpha_vec)
      
      elastic_model = elastic_model_obj$model
      
      # Elastic net predictions for the sample
      y_star_elastic_sample = predict(elastic_model,
                                      input_sampling_without_const[,x_input_col],
                                      type = "response",
                                      s = elastic_model$lambda.min)
      
      # get the result
      elastic_var_est = var_calculation_ma(y_star, y_star_elastic_sample,
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
      variance =  rrt_variance + elastic_var_est
      
      
      
      # Elastic net predicitons for the population
      y_star_elastic_pop = predict(elastic_model,
                                   pop_with_out_constant[,x_input_col],
                                   type = "response",
                                   s = elastic_model$lambda.min)
      #
      y_star_elastic_estimate =  mean(y_star - y_star_elastic_sample) +
        mean(y_star_elastic_pop)
      
      
      # get upper and lower for confidence interval
      actual_mu = mean(model_data$y)
      upper = y_star_elastic_estimate + 1.96 * sqrt(variance)
      # Get lower bound
      lower = y_star_elastic_estimate - 1.96 * sqrt(variance)
      
      covered_mu = (lower < actual_mu && actual_mu < upper)
      
      # Constructing elastic net data frame
      elastic_net_df = data.frame(run_seed, n_sizes[samp_size], y_star_elastic_estimate,
                                  elastic_model_obj$alpha,
                                  elastic_model$lambda.min,
                                  variance, rrt_variance, upper, lower, covered_mu)
      
      # Adding the column names
      colnames(elastic_net_df) <- c("run_seed","n_size","y_star_elastic", "alpha",
                                    "lambda",
                                    "variance", "rrt_variance",
                                    "upperCI", "lowerCI", "covered_mu")
      # Returning row of the dataframe
      return(elastic_net_df)
    }
  stopCluster(cluster_1)
  return(elastic_calculations)
}
