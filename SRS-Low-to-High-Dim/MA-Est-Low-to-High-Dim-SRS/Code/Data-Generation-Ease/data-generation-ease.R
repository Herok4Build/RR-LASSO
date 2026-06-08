#' Thomas Johnson III
#' 09-29-2024
#' Data Generation functions for targets pipeline update.
#' Mainly for cleaning the targets pipeline for ease of reproducibility versus running a bunch of scripts.
#' 


## ----message = FALSE, warning=FALSE-----------------------------------------------------------------------------
library(Matrix);library("doSNOW"); library("MASS"); library("tidyverse");

# Model 1 --------------

# Generate Model 1 with different predictors code---------------------------------------------------------------------------------------------------------------
# Adapted from Jeremy's code
correlate_unifrom_m.1_p <- function(N, rho, max_pred = 20)
{
  base_pred = 8
  num_pred = max_pred - base_pred # extra predictors that are noise
  # not genrating sigma squared equal to 0.16 without adding to diagonal
  sigma.mat = diag(max_pred) + ifelse(diag(max_pred) == 0, rho, 0)
  x_mat_model.1 = cbind(rep(1,N),apply(cbind(
    mvrnorm(N, mu = rep(0,max_pred), Sigma = sigma.mat, tol = 1)), 
    2, pnorm))
  colnames(x_mat_model.1) <- paste0("X", 0:max_pred)
  # Account for extra predictors
  coeff.mat = matrix(c(1, 0, 1, 0, 1.5, 0, 0, 0, 1, rep(0,num_pred)), 
                     nrow =max_pred+1, ncol = 1, byrow = TRUE)
  
  # error for the model
  error_m.1 = rnorm(N, 0, 0.1)
  model_1.rho = x_mat_model.1 %*% coeff.mat
  y <- model_1.rho + error_m.1
  
  m.1_model_rho = data.frame(x_mat_model.1, y)
  return(m.1_model_rho)
}


## ---------------------------------------------------------------------------------------------------------------
# For seed
data.seedling_1 = 12
N = 5000


model_1_data_list_create_25_50_pred = function(N = 5000, rho = c(0.1, 0.3, 0.5, 0.7), 
                                    data.seedling = data.seedling_1)
{
  set.seed(data.seedling_1)
  
  m.1_reg_data_rho.1_50_pred = correlate_unifrom_m.1_p(N, .1, max_pred = 50)
  m.1_reg_data_rho.3_50_pred = correlate_unifrom_m.1_p(N, .3, max_pred = 50)
  m.1_reg_data_rho.5_50_pred = correlate_unifrom_m.1_p(N, .5, max_pred = 50)
  m.1_reg_data_rho.7_50_pred = correlate_unifrom_m.1_p(N, .7, max_pred = 50)
  
  m.1_reg_data_rho.1_25_pred = correlate_unifrom_m.1_p(N, .1, max_pred = 25)
  m.1_reg_data_rho.3_25_pred = correlate_unifrom_m.1_p(N, .3, max_pred = 25)
  m.1_reg_data_rho.5_25_pred = correlate_unifrom_m.1_p(N, .5, max_pred = 25)
  m.1_reg_data_rho.7_25_pred = correlate_unifrom_m.1_p(N, .7, max_pred = 25)
  
  model_1_25_50_pred = list("25" = list("rho_0.1" = m.1_reg_data_rho.1_25_pred,
                                        "rho_0.3" = m.1_reg_data_rho.3_25_pred,
                                        "rho_0.5" = m.1_reg_data_rho.5_25_pred,
                                        "rho_0.7" = m.1_reg_data_rho.7_25_pred),
                            "50" = list("rho_0.1" = m.1_reg_data_rho.1_50_pred,
                                        "rho_0.3" = m.1_reg_data_rho.3_50_pred,
                                        "rho_0.5" = m.1_reg_data_rho.5_50_pred,
                                        "rho_0.7" = m.1_reg_data_rho.7_50_pred))
  return(model_1_25_50_pred)
}

# Model 1

model_1_data_list_create = function(N = 5000, data.seedling = data.seedling_1)
{
  set.seed(data.seedling)
  
  m.1_reg_data_rho.1_50_pred = correlate_unifrom_m.1_p(N, .1, max_pred = 50)
  m.1_reg_data_rho.3_50_pred = correlate_unifrom_m.1_p(N, .3, max_pred = 50)
  m.1_reg_data_rho.5_50_pred = correlate_unifrom_m.1_p(N, .5, max_pred = 50)
  m.1_reg_data_rho.7_50_pred = correlate_unifrom_m.1_p(N, .7, max_pred = 50)
  
  m.1_reg_data_rho.1_25_pred = correlate_unifrom_m.1_p(N, .1, max_pred = 25)
  m.1_reg_data_rho.3_25_pred = correlate_unifrom_m.1_p(N, .3, max_pred = 25)
  m.1_reg_data_rho.5_25_pred = correlate_unifrom_m.1_p(N, .5, max_pred = 25)
  m.1_reg_data_rho.7_25_pred = correlate_unifrom_m.1_p(N, .7, max_pred = 25)
  
  set.seed(data.seedling)
  
  m.1_reg_data_rho.1_100_pred = correlate_unifrom_m.1_p(N, .1, max_pred = 100)
  m.1_reg_data_rho.3_100_pred = correlate_unifrom_m.1_p(N, .3, max_pred = 100)
  m.1_reg_data_rho.5_100_pred = correlate_unifrom_m.1_p(N, .5, max_pred = 100)
  m.1_reg_data_rho.7_100_pred = correlate_unifrom_m.1_p(N, .7, max_pred = 100)
  
  model_1_data_100_pred = list("rho_0.1" =m.1_reg_data_rho.1_100_pred,
                               "rho_0.3" = m.1_reg_data_rho.3_100_pred,
                               "rho_0.5" = m.1_reg_data_rho.5_100_pred,
                               "rho_0.7" = m.1_reg_data_rho.7_100_pred)
  
  set.seed(data.seedling)
  
  m.1_reg_data_rho.1_200_pred = correlate_unifrom_m.1_p(N, .1, max_pred = 200)
  m.1_reg_data_rho.3_200_pred = correlate_unifrom_m.1_p(N, .3, max_pred = 200)
  m.1_reg_data_rho.5_200_pred = correlate_unifrom_m.1_p(N, .5, max_pred = 200)
  m.1_reg_data_rho.7_200_pred = correlate_unifrom_m.1_p(N, .7, max_pred = 200)
  
  model_1_data_200_pred = list("rho_0.1" =m.1_reg_data_rho.1_200_pred,
                               "rho_0.3" = m.1_reg_data_rho.3_200_pred,
                               "rho_0.5" = m.1_reg_data_rho.5_200_pred,
                               "rho_0.7" = m.1_reg_data_rho.7_200_pred)
  
  model_1_data_list = list("25" = list("rho_0.1" = m.1_reg_data_rho.1_25_pred,
                                        "rho_0.3" = m.1_reg_data_rho.3_25_pred,
                                        "rho_0.5" = m.1_reg_data_rho.5_25_pred,
                                        "rho_0.7" = m.1_reg_data_rho.7_25_pred),
                            "50" = list("rho_0.1" = m.1_reg_data_rho.1_50_pred,
                                        "rho_0.3" = m.1_reg_data_rho.3_50_pred,
                                        "rho_0.5" = m.1_reg_data_rho.5_50_pred,
                                        "rho_0.7" = m.1_reg_data_rho.7_50_pred),
                           "100" = model_1_data_100_pred,
                           "200" = model_1_data_200_pred)
  
  
  return(model_1_data_list)
}

model_1_data_list_create_high_dim_300_600 = function(N = 5000, data.seedling = data.seedling_1)
{
  
  set.seed(data.seedling)
  
  m.1_reg_data_rho.1_300_pred = correlate_unifrom_m.1_p(N, .1, max_pred = 300)
  m.1_reg_data_rho.3_300_pred = correlate_unifrom_m.1_p(N, .3, max_pred = 300)
  m.1_reg_data_rho.5_300_pred = correlate_unifrom_m.1_p(N, .5, max_pred = 300)
  m.1_reg_data_rho.7_300_pred = correlate_unifrom_m.1_p(N, .7, max_pred = 300)
  
  model_1_data_300_pred = list("rho_0.1" =m.1_reg_data_rho.1_300_pred,
                               "rho_0.3" = m.1_reg_data_rho.3_300_pred,
                               "rho_0.5" = m.1_reg_data_rho.5_300_pred,
                               "rho_0.7" = m.1_reg_data_rho.7_300_pred)
  
  set.seed(data.seedling)
  
  m.1_reg_data_rho.1_600_pred = correlate_unifrom_m.1_p(N, .1, max_pred = 600)
  m.1_reg_data_rho.3_600_pred = correlate_unifrom_m.1_p(N, .3, max_pred = 600)
  m.1_reg_data_rho.5_600_pred = correlate_unifrom_m.1_p(N, .5, max_pred = 600)
  m.1_reg_data_rho.7_600_pred = correlate_unifrom_m.1_p(N, .7, max_pred = 600)
  
  model_1_data_600_pred = list("rho_0.1" =m.1_reg_data_rho.1_600_pred,
                               "rho_0.3" = m.1_reg_data_rho.3_600_pred,
                               "rho_0.5" = m.1_reg_data_rho.5_600_pred,
                               "rho_0.7" = m.1_reg_data_rho.7_600_pred)
  
  model_1_data_list = list("300" = model_1_data_300_pred,
                           "600" = model_1_data_600_pred)
  
  
  return(model_1_data_list)
}

# Model 2 -------------------------------

# Incorrect version
# The variance is too small and will not trigger the indicator functions
# mod.2_generate = function(N = 5000, rho_value = 0.1, max_pred = 25, dat.seed = 13)
# {
#   # need to rerun as correlation was not correctly captured
#   set.seed(dat.seed)
#   # Generate the correlation-covariance matrix
#   matrix.cov.sigma = diag(rep(150^2,max_pred)) + 
#     ifelse(diag(rep(150^2,max_pred)) == 0, (rho_value * 150^2), 0)
#   # Generate data for Model
#   mat_data = mvrnorm(N, mu = rep(0,max_pred), 
#                      Sigma = matrix.cov.sigma)
#   # mat_data = apply(mat_data, 2, pnorm, mean = 0, sd = 1)
#   # 
#   # mat_data = apply(mat_data, 2, dnorm, mean = 150, sd = 10)
#   
#   mod_data = cbind(rep(500,N),
#                    mat_data)
#   
#   
#   model_x_result = mod_data[,1] + 2 * mod_data[,5] + 
#     400 * (mod_data[,6] > 156) - 400 * (mod_data[,6] <= 156) + 
#     1000 * (mod_data[,3] > 190) + 300 * (mod_data[,6] > 200)
#   epsilon = rnorm(N,0, 0.1) # Changed to be consistent with other models and to avoid issues with initial randomization (like too much noise)
#   y = model_x_result + epsilon
#   
#   mod_4_frame = as.data.frame(cbind(mod_data, y))
#   # Name the columns
#   colnames(mod_4_frame) <- c(paste0("X", (0:max_pred)), "y")
#   return(mod_4_frame)
# }

# Model 2 to run
mod.2_generate = function(N = 5000, rho_value = 0.1, max_pred = 25, dat.seed = 13)
{
  # need to rerun as correlation was not correctly captured
  set.seed(dat.seed)
  # Generate the correlation-covariance matrix
  matrix.cov.sigma = diag(rep(150^2,max_pred)) + 
    ifelse(diag(rep(150^2,max_pred)) == 0, (rho_value * 150^2), 0)
  # Generate data for Model
  mat_data = mvrnorm(N, mu = rep(0,max_pred), 
                     Sigma = matrix.cov.sigma)
  # mat_data = apply(mat_data, 2, pnorm, mean = 0, sd = 1)
  # 
  # mat_data = apply(mat_data, 2, dnorm, mean = 150, sd = 10)
  
  mod_data = cbind(rep(500,N),
                   mat_data)
  
  
  model_x_result = mod_data[,1] + 2 * mod_data[,5] + 
    400 * (mod_data[,6] > 156) - 400 * (mod_data[,6] <= 156) + 
    1000 * (mod_data[,3] > 190) + 300 * (mod_data[,6] > 200)
  epsilon = rnorm(N,0, 0.1) # Changed to be consistent with other models and to avoid issues with initial randomization (like too much noise)
  y = model_x_result + epsilon
  
  mod_4_frame = as.data.frame(cbind(mod_data, y))
  # Name the columns
  colnames(mod_4_frame) <- c(paste0("X", (0:max_pred)), "y")
  return(mod_4_frame)
}

model_2_data_list_create = function(N = 5000, data.seedling = data.seedling_1)
{
  
  set.seed(data.seedling)
  
  model_2_data_25_pred =
    list("rho_0.1" = mod.2_generate(N = 5000, rho_value = 0.1, max_pred = 25, dat.seed = data.seedling),
         "rho_0.3" = mod.2_generate(N = 5000, rho_value = 0.3, max_pred = 25, dat.seed = data.seedling),
         "rho_0.5" = mod.2_generate(N = 5000, rho_value = 0.5, max_pred = 25, dat.seed = data.seedling),
         "rho_0.7" = mod.2_generate(N = 5000, rho_value = 0.7, max_pred = 25, dat.seed = data.seedling))
  
  set.seed(data.seedling)
  
  model_2_data_50_pred =
    list("rho_0.1" = mod.2_generate(N = 5000, rho_value = 0.1, max_pred = 50, dat.seed = data.seedling),
         "rho_0.3" = mod.2_generate(N = 5000, rho_value = 0.3, max_pred = 50, dat.seed = data.seedling),
         "rho_0.5" = mod.2_generate(N = 5000, rho_value = 0.5, max_pred = 50, dat.seed = data.seedling),
         "rho_0.7" = mod.2_generate(N = 5000, rho_value = 0.7, max_pred = 50, dat.seed = data.seedling))
  
  set.seed(data.seedling)
  
  model_2_data_100_pred =
    list("rho_0.1" = mod.2_generate(N = 5000, rho_value = 0.1, max_pred = 100, dat.seed = data.seedling),
         "rho_0.3" = mod.2_generate(N = 5000, rho_value = 0.3, max_pred = 100, dat.seed = data.seedling),
         "rho_0.5" = mod.2_generate(N = 5000, rho_value = 0.5, max_pred = 100, dat.seed = data.seedling),
         "rho_0.7" = mod.2_generate(N = 5000, rho_value = 0.7, max_pred = 100, dat.seed = data.seedling))
  
  set.seed(data.seedling)
  
  model_2_data_200_pred =
    list("rho_0.1" = mod.2_generate(N = 5000, rho_value = 0.1, max_pred = 200, dat.seed = data.seedling),
         "rho_0.3" = mod.2_generate(N = 5000, rho_value = 0.3, max_pred = 200, dat.seed = data.seedling),
         "rho_0.5" = mod.2_generate(N = 5000, rho_value = 0.5, max_pred = 200, dat.seed = data.seedling),
         "rho_0.7" = mod.2_generate(N = 5000, rho_value = 0.7, max_pred = 200, dat.seed = data.seedling))
  
  model_2_data_list = list("25" = model_2_data_25_pred,
                           "50" = model_2_data_50_pred,
                           "100" = model_2_data_100_pred,
                           "200" = model_2_data_200_pred)
  
  
  return(model_2_data_list )
}

model_2_data_list_create_high_dim = function(N = 5000, data.seedling = data.seedling_1)
{
  
  
  set.seed(data.seedling)
  
  model_2_data_300_pred =
    list("rho_0.1" = mod.2_generate(N = 5000, rho_value = 0.1, max_pred = 300, dat.seed = data.seedling),
         "rho_0.3" = mod.2_generate(N = 5000, rho_value = 0.3, max_pred = 300, dat.seed = data.seedling),
         "rho_0.5" = mod.2_generate(N = 5000, rho_value = 0.5, max_pred = 300, dat.seed = data.seedling),
         "rho_0.7" = mod.2_generate(N = 5000, rho_value = 0.7, max_pred = 300, dat.seed = data.seedling))
  
  set.seed(data.seedling)
  
  model_2_data_600_pred =
    list("rho_0.1" = mod.2_generate(N = 5000, rho_value = 0.1, max_pred = 600, dat.seed = data.seedling),
         "rho_0.3" = mod.2_generate(N = 5000, rho_value = 0.3, max_pred = 600, dat.seed = data.seedling),
         "rho_0.5" = mod.2_generate(N = 5000, rho_value = 0.5, max_pred = 600, dat.seed = data.seedling),
         "rho_0.7" = mod.2_generate(N = 5000, rho_value = 0.7, max_pred = 600, dat.seed = data.seedling))
  
  model_2_data_list = list("300" = model_2_data_300_pred,
                           "600" = model_2_data_600_pred)
  
  
  return(model_2_data_list )
}



# Model 3 ---------------------------

m.3_correlated_uniform <- function(N, rho_value, max_pred, dat.seed = 12)
{
  set.seed(dat.seed)
  base_pred = 20
  num_pred = abs(max_pred - 20)
  matrix.sigma.cov = diag(max_pred) + ifelse(diag(max_pred) == 0, rho_value, 0)
  X_model.3 = apply(cbind(rep(1,N), 
                          mvrnorm(N, mu = rep(0,max_pred), Sigma = matrix.sigma.cov)), 
                    2, pnorm)
  colnames(X_model.3) <- paste0("X",0:max_pred)
  fresh_coef = c(0, 0.15, 0.25, 0, 0.34, 0, 0.4, 0, 0, 0.19)
  # Number of non-significant coefficients
  non.sig_coeffs = rep(0,(max_pred - length(fresh_coef)))
  coeff.mat <- matrix(c(1,c(fresh_coef, non.sig_coeffs)), 
                      nrow =max_pred +1, ncol = 1, 
                      byrow =TRUE)
  model.3.func = X_model.3 %*% coeff.mat
  print("The coefficients")
  print(which(as.vector(t(coeff.mat))!= 0))
  print("The values")
  print(coeff.mat[which(as.vector(t(coeff.mat))!= 0)])
  error_model.3 = rnorm(N, 0, 0.1)
  y <- model.3.func + error_model.3
  
  model.3.data_total = data.frame(X_model.3,y)
  return(model.3.data_total)
}

model_3_data_list_create = function(N = 5000, data.seedling = data.seedling_1)
{
  
  set.seed(data.seedling)
  
  model_3_data_25_pred =
    list("rho_0.1" = m.3_correlated_uniform(N = 5000, rho_value = 0.1, max_pred = 25, dat.seed = data.seedling),
         "rho_0.3" = m.3_correlated_uniform(N = 5000, rho_value = 0.3, max_pred = 25, dat.seed = data.seedling),
         "rho_0.5" = m.3_correlated_uniform(N = 5000, rho_value = 0.5, max_pred = 25, dat.seed = data.seedling),
         "rho_0.7" = m.3_correlated_uniform(N = 5000, rho_value = 0.7, max_pred = 25, dat.seed = data.seedling))
  
  model_3_data_50_pred =
    list("rho_0.1" = m.3_correlated_uniform(N = 5000, rho_value = 0.1, max_pred = 50, dat.seed = data.seedling),
         "rho_0.3" = m.3_correlated_uniform(N = 5000, rho_value = 0.3, max_pred = 50, dat.seed = data.seedling),
         "rho_0.5" = m.3_correlated_uniform(N = 5000, rho_value = 0.5, max_pred = 50, dat.seed = data.seedling),
         "rho_0.7" = m.3_correlated_uniform(N = 5000, rho_value = 0.7, max_pred = 50, dat.seed = data.seedling))
  
  set.seed(data.seedling)
  
  model_3_data_100_pred =
    list("rho_0.1" = m.3_correlated_uniform(N = 5000, rho_value = 0.1, max_pred = 100, dat.seed = data.seedling),
         "rho_0.3" = m.3_correlated_uniform(N = 5000, rho_value = 0.3, max_pred = 100, dat.seed = data.seedling),
         "rho_0.5" = m.3_correlated_uniform(N = 5000, rho_value = 0.5, max_pred = 100, dat.seed = data.seedling),
         "rho_0.7" = m.3_correlated_uniform(N = 5000, rho_value = 0.7, max_pred = 100, dat.seed = data.seedling))
  
  set.seed(data.seedling)
  
  model_3_data_200_pred =
    list("rho_0.1" = m.3_correlated_uniform(N = 5000, rho_value = 0.1, max_pred = 200, dat.seed = data.seedling),
         "rho_0.3" = m.3_correlated_uniform(N = 5000, rho_value = 0.3, max_pred = 200, dat.seed = data.seedling),
         "rho_0.5" = m.3_correlated_uniform(N = 5000, rho_value = 0.5, max_pred = 200, dat.seed = data.seedling),
         "rho_0.7" = m.3_correlated_uniform(N = 5000, rho_value = 0.7, max_pred = 200, dat.seed = data.seedling))
  
  model_3_data_list = list("25" = model_3_data_25_pred,
                           "50" = model_3_data_50_pred,
                           "100" = model_3_data_100_pred,
                           "200" = model_3_data_200_pred)
  
  
  return(model_3_data_list )
}



model_3_data_high_dim_list_create = function(N = 5000, data.seedling = data.seedling_1)
{
  
  set.seed(data.seedling)
  
  model_3_data_300_pred =
    list("rho_0.1" = m.3_correlated_uniform(N = N, rho_value = 0.1, max_pred = 300, dat.seed = data.seedling),
         "rho_0.3" = m.3_correlated_uniform(N = N, rho_value = 0.3, max_pred = 300, dat.seed = data.seedling),
         "rho_0.5" = m.3_correlated_uniform(N = N, rho_value = 0.5, max_pred = 300, dat.seed = data.seedling),
         "rho_0.7" = m.3_correlated_uniform(N = N, rho_value = 0.7, max_pred = 300, dat.seed = data.seedling))
  
  model_3_data_600_pred =
    list("rho_0.1" = m.3_correlated_uniform(N = N, rho_value = 0.1, max_pred = 600, dat.seed = data.seedling),
         "rho_0.3" = m.3_correlated_uniform(N = N, rho_value = 0.3, max_pred = 600, dat.seed = data.seedling),
         "rho_0.5" = m.3_correlated_uniform(N = N, rho_value = 0.5, max_pred = 600, dat.seed = data.seedling),
         "rho_0.7" = m.3_correlated_uniform(N = N, rho_value = 0.7, max_pred = 600, dat.seed = data.seedling))
  
  
  model_3_data_list = list("300" = model_3_data_300_pred,
                           "600" = model_3_data_600_pred)
  
  
  return(model_3_data_list )
}

# Model 4 ----------------------

# From Meidi, Model 8
m_4_regression_func = function(N = 5000, rho_value, #coef_vec = c(-1,1,1,-1),
                               max_pred = 50, dat.seed = 12){
  set.seed(dat.seed)
  scale_vars = function(old_var, lower_bound=-1, upper_bound=1){
    new_var = scales::rescale(old_var, to = c(lower_bound, upper_bound))
    return(new_var)
  }
  v_vars_num = max_pred # Two categorical variables from multinomial
  sigma.mat = diag(rep(2,v_vars_num)) + ifelse(diag(rep(2,v_vars_num)) == 0, 2 * rho_value, 0)
  v_vars =
    apply(mvrnorm(N, mu = rep(0,v_vars_num), Sigma = sigma.mat),2,scale_vars)
  
  # Generate data for Model
  mod_data = cbind(rep(3, N), # intercept
                   v_vars)
  mod_data = as.data.frame(mod_data)
  # Name the columns
  colnames(mod_data) <- c(paste0("X", (0:max_pred)))
  
  model_x_result = mod_data[,"X0"] +
    mod_data[,"X1"] * mod_data[,"X2"] + mod_data[,"X3"]^2 +
    mod_data[,"X4"] * mod_data[,"X7"] +
    mod_data[,"X8"] * mod_data[,"X10"] -
    mod_data[,"X6"]
  
  epsilon = rnorm(N, 0, 0.1)
  y = model_x_result + epsilon
  
  mod_4_frame = as.data.frame(cbind(mod_data, y))
  
  colnames(mod_4_frame) <- c(paste0("X", (0:max_pred)), "y")
  return(mod_4_frame)
}


model_4_data_list_create = function(N = 5000, data.seedling = data.seedling_1)
{
  
  set.seed(data.seedling)
  
  model_4_data_25_pred =
    list("rho_0.1" = m_4_regression_func(N = 5000, rho_value = 0.1, max_pred = 25, dat.seed = data.seedling),
         "rho_0.3" = m_4_regression_func(N = 5000, rho_value = 0.3, max_pred = 25, dat.seed = data.seedling),
         "rho_0.5" = m_4_regression_func(N = 5000, rho_value = 0.5, max_pred = 25, dat.seed = data.seedling),
         "rho_0.7" = m_4_regression_func(N = 5000, rho_value = 0.7, max_pred = 25, dat.seed = data.seedling))
  
  model_4_data_50_pred =
    list("rho_0.1" = m_4_regression_func(N = 5000, rho_value = 0.1, max_pred = 50, dat.seed = data.seedling),
         "rho_0.3" = m_4_regression_func(N = 5000, rho_value = 0.3, max_pred = 50, dat.seed = data.seedling),
         "rho_0.5" = m_4_regression_func(N = 5000, rho_value = 0.5, max_pred = 50, dat.seed = data.seedling),
         "rho_0.7" = m_4_regression_func(N = 5000, rho_value = 0.7, max_pred = 50, dat.seed = data.seedling))
  
  set.seed(data.seedling)
  
  model_4_data_100_pred =
    list("rho_0.1" = m_4_regression_func(N = 5000, rho_value = 0.1, max_pred = 100, dat.seed = data.seedling),
         "rho_0.3" = m_4_regression_func(N = 5000, rho_value = 0.3, max_pred = 100, dat.seed = data.seedling),
         "rho_0.5" = m_4_regression_func(N = 5000, rho_value = 0.5, max_pred = 100, dat.seed = data.seedling),
         "rho_0.7" = m_4_regression_func(N = 5000, rho_value = 0.7, max_pred = 100, dat.seed = data.seedling))
  
  set.seed(data.seedling)
  
  model_4_data_200_pred =
    list("rho_0.1" = m_4_regression_func(N = 5000, rho_value = 0.1, max_pred = 200, dat.seed = data.seedling),
         "rho_0.3" = m_4_regression_func(N = 5000, rho_value = 0.3, max_pred = 200, dat.seed = data.seedling),
         "rho_0.5" = m_4_regression_func(N = 5000, rho_value = 0.5, max_pred = 200, dat.seed = data.seedling),
         "rho_0.7" = m_4_regression_func(N = 5000, rho_value = 0.7, max_pred = 200, dat.seed = data.seedling))
  
  model_4_data_list = list("25" = model_4_data_25_pred,
                           "50" = model_4_data_50_pred,
                           "100" = model_4_data_100_pred,
                           "200" = model_4_data_200_pred)
  
  
  return(model_4_data_list )
}

model_4_data_high_dim_list_create = function(N = 5000, data.seedling = data.seedling_1)
{
  
  set.seed(data.seedling)
  
  model_4_data_300_pred =
    list("rho_0.1" = m_4_regression_func(N = N, rho_value = 0.1, max_pred = 300, dat.seed = data.seedling),
         "rho_0.3" = m_4_regression_func(N = N, rho_value = 0.3, max_pred = 300, dat.seed = data.seedling),
         "rho_0.5" = m_4_regression_func(N = N, rho_value = 0.5, max_pred = 300, dat.seed = data.seedling),
         "rho_0.7" = m_4_regression_func(N = N, rho_value = 0.7, max_pred = 300, dat.seed = data.seedling))
  
  model_4_data_600_pred =
    list("rho_0.1" = m_4_regression_func(N = N, rho_value = 0.1, max_pred = 600, dat.seed = data.seedling),
         "rho_0.3" = m_4_regression_func(N = N, rho_value = 0.3, max_pred = 600, dat.seed = data.seedling),
         "rho_0.5" = m_4_regression_func(N = N, rho_value = 0.5, max_pred = 600, dat.seed = data.seedling),
         "rho_0.7" = m_4_regression_func(N = N, rho_value = 0.7, max_pred = 600, dat.seed = data.seedling))
  
  model_4_data_list = list("300" = model_4_data_300_pred,
                           "600" = model_4_data_600_pred)
  
  
  return(model_4_data_list )
}



