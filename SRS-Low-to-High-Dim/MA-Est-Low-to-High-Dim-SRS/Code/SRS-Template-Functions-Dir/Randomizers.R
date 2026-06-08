# Randomizers
# Thomas Johnson III
# The randomization funciton for y

rand_y_2norm_proc = function(y){
  
  alpha = 0 # Mean of the normal distribution used
  # Mran of the F distribution with df1 = 10, df2 = 10
  beta = 1.5
  gamma = (.2/1.5)^2
  eta = 0.1^2
  #print(rand_select)
  # From Chen, Zhao and Wang
  T_dist_val = rnorm(length(y),beta,.2/1.5)
  S_norm_dist_val = rnorm(length(y), alpha,0.1)
  r = y * T_dist_val + S_norm_dist_val
  y_star = (r - alpha) / beta
  rrt_obj = list("values" = y_star,
                 "beta" = beta,
                 "gamma" = gamma,
                 "eta" = eta,
                 "alpha" = alpha)
  return(rrt_obj)
}