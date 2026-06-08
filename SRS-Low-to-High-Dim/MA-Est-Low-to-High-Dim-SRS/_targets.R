# Targets R script for reproducible R code
# Thomas Johnson III
# 07-23-2023

# import targets
library(targets)


# Loading LASSO funcitons
# The initial test run will be with relaxed LASSO

srs_templates_path = paste0(here::here(), "/Code/SRS-Template-Functions-Dir/")
srs_target_ext_paths =  paste0(here::here(), "/Code/targets-est/")
srs_target_data_gen_path =  paste0(here::here(), "/Code/Data-Generation-Ease/data-generation-ease.R")


source(paste0(srs_templates_path, "LASSO-Estimator-Template.R"))
source(paste0(srs_templates_path, "Linear-Reg-Estimator-Template.R"))
source(paste0(srs_templates_path, "Randomizers.R"))
source(paste0(srs_templates_path, "Basic-Funcs.R"))
source(paste0(srs_templates_path, "HT-Estimator-Template.R"))
source(paste0(srs_templates_path, "EN-Estimator-Template.R"))
source(srs_target_data_gen_path)

# Other targets files
source(paste0(srs_target_ext_paths, "xgbtree-targets.R"))
source(paste0(srs_target_ext_paths, "Ranger-RF-Targets.R"))
source(paste0(srs_target_ext_paths, "LASSO-PEML-Bootstrap.R"))

# Setup data storage
srs_data_store_model_1 = paste0(here::here(), "/Data/SRS-Model-1-Dim/Variance/")

if(dir.exists(srs_data_store_model_1 ) == FALSE)
{
  dir.create(srs_data_store_model_1 )
}


srs_data_store_model_1_peml = paste0(here::here(), "/Data/SRS-Model-1-Dim/Variance/peml/")

if(dir.exists(srs_data_store_model_1_peml) == FALSE)
{
  dir.create(srs_data_store_model_1_peml)
}

srs_data_store_model_2_peml = paste0(here::here(), "/Data/SRS-Model-2-Dim/Variance/peml/")

if(dir.exists(srs_data_store_model_2_peml) == FALSE)
{
  dir.create(srs_data_store_model_2_peml)
}

srs_data_store_model_3_peml = paste0(here::here(), "/Data/SRS-Model-3-Static-Dim/Variance/peml/")

if(dir.exists(srs_data_store_model_3_peml) == FALSE)
{
  dir.create(srs_data_store_model_3_peml)
}

srs_data_store_model_4_peml = paste0(here::here(), "/Data/SRS-Model-4-Dim/Variance/peml/")

if(dir.exists(srs_data_store_model_4_peml) == FALSE)
{
  dir.create(srs_data_store_model_4_peml)
}

srs_data_store_model_2 = paste0(here::here(), "/Data/SRS-Model-2-Dim/Variance/")

srs_data_store_model_2_rf_term_node_9 = paste0(srs_data_store_model_2, "term_node_9_20")

if(dir.exists(srs_data_store_model_2_rf_term_node_9) == FALSE)
{
  dir.create(srs_data_store_model_2_rf_term_node_9)
}

srs_data_store_model_3 = paste0(here::here(), "/Data/SRS-Model-3-Static-Dim/Variance/")

srs_data_store_model_4 = paste0(here::here(), "/Data/SRS-Model-4-Dim/Variance/")

srs_model_3_data_dir = paste0(here::here(), "/Data/SRS-Model-3-Static-Dim/")

core_alloc = parallel::detectCores() - 2

tar_option_set(packages = c("tidyverse", "doSNOW"))

sample_sizes = c(  250, 500, 1000)

samp_sizes_resampling = c(250, 500)

list(
  # Model 1 --------------------------------------------------
  # 25 predictors
  tar_target(model_1_data_list,model_1_data_list_create()),
  tar_target(model_1_data_high_dim_list,model_1_data_list_create_high_dim_300_600()),
  
  tar_target(save_model_1_25_pred, saveRDS(model_1_data_list[["25"]],
                                           paste0(srs_data_store_model_1 ,
                                                  "model_1_data_25_pred.rds"))),
  tar_target(save_model_1_50_pred, saveRDS(model_1_data_list[["50"]],
                                           paste0(srs_data_store_model_1 ,
                                                  "model_1_data_50_pred.rds"))),
  tar_target(save_model_1_100_pred, saveRDS(model_1_data_list[["100"]],
                                           paste0(srs_data_store_model_1 ,
                                                  "model_1_data_100_pred.rds"))),
  tar_target(save_model_1_200_pred, saveRDS(model_1_data_list[["200"]],
                                           paste0(srs_data_store_model_1 ,
                                                  "model_1_data_200_pred.rds"))),
  
  tar_target(save_model_1_300_pred, saveRDS(model_1_data_high_dim_list[["300"]],
                                            paste0(srs_data_store_model_1 ,
                                                   "model_1_data_300_pred.rds"))),
  
  tar_target(save_model_1_600_pred, saveRDS(model_1_data_high_dim_list[["600"]],
                                            paste0(srs_data_store_model_1 ,
                                                   "model_1_data_600_pred.rds"))),
  tar_target(starting_m, seq(1,1000,by = 100)),
  tar_target(ending_m,  seq(100,1000,by = 100)),

  ## HT Estimator --------
  ### 25 Predictors -------
  tar_target(ht_estimator_model_1_25_pred,
             map_dfr(model_1_data_list[["25"]], fixed_size_ht_est_variance,
                     n_sizes = c(  250, 500, 1000),cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc, 
                     .id = "rho")),
  tar_target(ht_estimator_model_1_25_pred_results,
             produce_csv(ht_estimator_model_1_25_pred, 
                         srs_data_store_model_1, "ht-plugin-var-sim-model-1-25-pred.csv")),
  ### 50 Predictors ------
  tar_target(ht_estimator_model_1_50_pred,
             map_dfr(model_1_data_list[["50"]], fixed_size_ht_est_variance,
                     n_sizes = c(  250, 500, 1000),cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc, 
                     .id = "rho")),
  tar_target(ht_estimator_model_1_50_pred_results,
             produce_csv(ht_estimator_model_1_50_pred, 
                         srs_data_store_model_1, "ht-plugin-var-sim-model-1-50-pred.csv")),
  ### 100 Predictors ------
  tar_target(ht_estimator_model_1_100_pred,
             map_dfr(model_1_data_list[["100"]], fixed_size_ht_est_variance,
                     n_sizes = c(  250, 500, 1000),cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc, 
                     .id = "rho")),
  tar_target(ht_estimator_model_1_100_pred_results,
             produce_csv(ht_estimator_model_1_100_pred, 
                         srs_data_store_model_1, "ht-plugin-var-sim-model-1-100-pred.csv")),
  ### 200 Predictors ---------------------
  tar_target(ht_estimator_model_1_200_pred,
             map_dfr(model_1_data_list[["200"]], fixed_size_ht_est_variance,
                     n_sizes = c(  250, 500, 1000),cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc, 
                     .id = "rho")),
  tar_target(ht_estimator_model_1_200_pred_results,
             produce_csv(ht_estimator_model_1_200_pred, 
                         srs_data_store_model_1, "ht-plugin-var-sim-model-1-200-pred.csv")),
  
  ### 300 Predictors ---------------------
  tar_target(ht_estimator_model_1_300_pred,
             map_dfr(model_1_data_high_dim_list[["300"]], fixed_size_ht_est_variance,
                     n_sizes = c(  250, 500, 1000),cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc, 
                     .id = "rho")),
  tar_target(ht_estimator_model_1_300_pred_results,
             produce_csv(ht_estimator_model_1_300_pred, 
                         srs_data_store_model_1, "ht-plugin-var-sim-model-1-300-pred.csv")),
  
  ### 600 Predictors ---------------------
  tar_target(ht_estimator_model_1_600_pred,
             map_dfr(model_1_data_high_dim_list[["600"]], fixed_size_ht_est_variance,
                     n_sizes = c(  250, 500, 1000),cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc, 
                     .id = "rho")),
  tar_target(ht_estimator_model_1_600_pred_results,
             produce_csv(ht_estimator_model_1_600_pred, 
                         srs_data_store_model_1, "ht-plugin-var-sim-model-1-600-pred.csv")),
  ## GREG/Linear Regression --------
  ### 25 predictors -----------
  tar_target(linear_reg_estimator_model_1_25_pred,
             map_dfr(model_1_data_list[["25"]], linear_reg_fixed_var,
                     n_sizes = c(  250, 500, 1000), cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc,.id = "rho")),
  tar_target(linear_reg_estimator_model_1_25_pred_results,
             produce_csv(linear_reg_estimator_model_1_25_pred, 
                         srs_data_store_model_1, "linear-plugin-var-sim-model-1-25-pred.csv")),
  ### 50 predictors ----------
  tar_target(linear_reg_estimator_model_1_50_pred,
             map_dfr(model_1_data_list[["50"]], linear_reg_fixed_var,
                     n_sizes = c(  250, 500, 1000), cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc,.id = "rho")),
  tar_target(linear_reg_estimator_model_1_50_pred_results,
             produce_csv(linear_reg_estimator_model_1_50_pred, 
                         srs_data_store_model_1, "linear-plugin-var-sim-model-1-50-pred.csv")),
  
  ### 100 predictors -------
  tar_target(linear_reg_estimator_model_1_100_pred,
             map_dfr(model_1_data_list[["100"]], linear_reg_fixed_var,
                     n_sizes = c(  250, 500, 1000), cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc,.id = "rho")),
  tar_target(linear_reg_estimator_model_1_100_pred_results,
             produce_csv(linear_reg_estimator_model_1_100_pred, 
                         srs_data_store_model_1, "linear-plugin-var-sim-model-1-100-pred.csv")),
  
  ### 200 predictors ------
  tar_target(linear_reg_estimator_model_1_200_pred,
             map_dfr(model_1_data_list[["200"]], linear_reg_fixed_var,
                     n_sizes = c(  250, 500, 1000), cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc,.id = "rho")),
  tar_target(linear_reg_estimator_model_1_200_pred_results,
             produce_csv(linear_reg_estimator_model_1_200_pred, 
                         srs_data_store_model_1, "linear-plugin-var-sim-model-1-200-pred.csv")),
  
  ### 300 predictors ------
  tar_target(linear_reg_estimator_model_1_300_pred,
             map_dfr(model_1_data_high_dim_list[["300"]], linear_reg_fixed_var,
                     n_sizes = c( 500, 1000), cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc,.id = "rho")),
  tar_target(linear_reg_estimator_model_1_300_pred_results,
             produce_csv(linear_reg_estimator_model_1_300_pred, 
                         srs_data_store_model_1, "linear-plugin-var-sim-model-1-300-pred.csv")),
  
  ## LASSO Estimator ---------------------------
  ### 25 Predictors -------------
  tar_target(lasso_estimator_model_1_25_pred,
             map_dfr(model_1_data_list[["25"]], lasso_reg_var_est ,
                     n_sizes = c(  250, 500, 1000),cores = core_alloc,
                     m=1000, 
                     randomizer = rand_y_2norm_proc, .id = "rho")),
  tar_target(lasso_estimator_model_1_25_pred_results,
             produce_csv(lasso_estimator_model_1_25_pred, 
                         srs_data_store_model_1, "lasso-plugin-var-sim-model-1-25-pred.csv")),
  
  #### Bootstrap Simulations -----------------------
  tar_target(lasso_estimator_bootstrap_model_1_25_pred,
             map_dfr(model_1_data_list[["25"]], lasso_reg_est,
                     n_sizes = c(  250, 500),cores = core_alloc,
                     m=1000, boot_yes = TRUE,
                     randomizer = rand_y_2norm_proc, .id = "rho")),
  tar_target(lasso_estimator_bootstrap_model_1_25_pred_results,
             produce_csv(lasso_estimator_bootstrap_model_1_25_pred, 
                         srs_data_store_model_1, "lasso-boot-sim-model-1-25-pred.csv")),
  
  # Mashregi version
  tar_target(lasso_estimator_mash.boot_model_1_25_pred,
             map_dfr(model_1_data_list[["25"]], lasso_reg_est,
                     n_sizes = c(  250, 500),cores = core_alloc,
                     m=1000, boot_yes = TRUE,
                     num_replicates = 250,
                     lasso_reg_boot_var_function = lasso_reg_bootstrap.var_est.mash,
                     randomizer = rand_y_2norm_proc, .id = "rho")),
  tar_target(lasso_estimator_mash.boot_model_1_25_pred_results,
             produce_csv(lasso_estimator_mash.boot_model_1_25_pred, 
                         srs_data_store_model_1, "lasso-mash-boot-sim-model-1-25-pred.csv")),
  
  
  ### 50 Predictors -------------
  tar_target(lasso_estimator_model_1_50_pred,
             map_dfr(model_1_data_list[["50"]], lasso_reg_var_est ,
                     n_sizes = c(  250, 500, 1000),cores = core_alloc,
                     m=1000, 
                     randomizer = rand_y_2norm_proc, .id = "rho")),
  tar_target(lasso_estimator_model_1_50_pred_results,
             produce_csv(lasso_estimator_model_1_50_pred, 
                         srs_data_store_model_1, "lasso-plugin-var-sim-model-1-50-pred.csv")),
  
  #### Bootstrap Simulations -----------------------
  tar_target(lasso_estimator_bootstrap_model_1_50_pred,
             map_dfr(model_1_data_list[["50"]], lasso_reg_est,
                     n_sizes = c(  250, 500),cores = core_alloc,
                     m=1000, boot_yes = TRUE,
                     randomizer = rand_y_2norm_proc, .id = "rho")),
  tar_target(lasso_estimator_bootstrap_model_1_50_pred_results,
             produce_csv(lasso_estimator_bootstrap_model_1_50_pred,
                         srs_data_store_model_1, "lasso-boot-sim-model-1-50-pred.csv")),
  
  # Mashregi version
  tar_target(lasso_estimator_mash.boot_model_1_50_pred,
             map_dfr(model_1_data_list[["50"]], lasso_reg_est,
                     n_sizes = c(  250, 500),cores = core_alloc,
                     m=1000, boot_yes = TRUE,
                     num_replicates = 250,
                     lasso_reg_boot_var_function = lasso_reg_bootstrap.var_est.mash,
                     randomizer = rand_y_2norm_proc, .id = "rho")),
  tar_target(lasso_estimator_mash.boot_model_1_50_pred_results,
             produce_csv(lasso_estimator_mash.boot_model_1_50_pred, 
                         srs_data_store_model_1, "lasso-mash-boot-sim-model-1-50-pred.csv")),
  
  
  ### 100 Predictors ----------------
  tar_target(lasso_estimator_model_1_100_pred,
             map_dfr(model_1_data_list[["100"]], lasso_reg_var_est ,
                     n_sizes = c(  250, 500, 1000),cores = core_alloc,
                     m=1000, 
                     randomizer = rand_y_2norm_proc, .id = "rho")),
  tar_target(lasso_estimator_model_1_100_pred_results,
             produce_csv(lasso_estimator_model_1_100_pred, 
                         srs_data_store_model_1, "lasso-plugin-var-sim-model-1-100-pred.csv")),
  
  #### Bootstrap Simulations -----------------------
  tar_target(lasso_estimator_bootstrap_model_1_100_pred,
             map_dfr(model_1_data_list[["100"]], lasso_reg_est,
                     n_sizes = c(  250, 500),cores = core_alloc,
                     m=1000, boot_yes = TRUE,
                     randomizer = rand_y_2norm_proc, .id = "rho")),
  tar_target(lasso_estimator_bootstrap_model_1_100_pred_results,
             produce_csv(lasso_estimator_bootstrap_model_1_100_pred,
                         srs_data_store_model_1, "lasso-boot-sim-model-1-100-pred.csv")),
  
  # Mashregi version
  tar_target(lasso_estimator_mash.boot_model_1_100_pred,
             map_dfr(model_1_data_list[["100"]], lasso_reg_est,
                     n_sizes = c(  250, 500),cores = core_alloc,
                     m=1000, boot_yes = TRUE,
                     num_replicates = 250,
                     lasso_reg_boot_var_function = lasso_reg_bootstrap.var_est.mash,
                     randomizer = rand_y_2norm_proc, .id = "rho")),
  tar_target(lasso_estimator_mash.boot_model_1_100_pred_results,
             produce_csv(lasso_estimator_mash.boot_model_1_100_pred, 
                         srs_data_store_model_1, "lasso-mash-boot-sim-model-1-100-pred.csv")),
  
  ### 200 Predictors ----------------
  tar_target(lasso_estimator_model_1_200_pred,
             map_dfr(model_1_data_list[["200"]], lasso_reg_var_est ,
                     n_sizes = c(  250, 500, 1000),cores = core_alloc,
                     m=1000, 
                     randomizer = rand_y_2norm_proc, .id = "rho")),
  tar_target(lasso_estimator_model_1_200_pred_results,
             produce_csv(lasso_estimator_model_1_200_pred, 
                         srs_data_store_model_1, "lasso-plugin-var-sim-model-1-200-pred.csv")),
  
  #### Bootstrap Simulations -----------------------
  tar_target(lasso_estimator_bootstrap_model_1_200_pred,
             map_dfr(model_1_data_list[["200"]], lasso_reg_est,
                     n_sizes = c(  250, 500),cores = core_alloc,
                     m=1000, boot_yes = TRUE,
                     randomizer = rand_y_2norm_proc, .id = "rho")),
  tar_target(lasso_estimator_bootstrap_model_1_200_pred_results,
             produce_csv(lasso_estimator_bootstrap_model_1_200_pred,
                         srs_data_store_model_1, "lasso-boot-sim-model-1-200-pred.csv")),
  
  # Mashregi version
  tar_target(lasso_estimator_mash.boot_model_1_200_pred,
             map_dfr(model_1_data_list[["200"]], lasso_reg_est,
                     n_sizes = c(  250, 500),cores = core_alloc,
                     m=1000, boot_yes = TRUE,
                     num_replicates = 250,
                     lasso_reg_boot_var_function = lasso_reg_bootstrap.var_est.mash,
                     randomizer = rand_y_2norm_proc, .id = "rho")),
  tar_target(lasso_estimator_mash.boot_model_1_200_pred_results,
             produce_csv(lasso_estimator_mash.boot_model_1_200_pred, 
                         srs_data_store_model_1, "lasso-mash-boot-sim-model-1-200-pred.csv")),
  
  
  ### 300 Predictors ----------------
  tar_target(lasso_estimator_model_1_300_pred,
             map_dfr(model_1_data_high_dim_list[["300"]], lasso_reg_var_est ,
                     n_sizes = c(  250, 500, 1000),cores = core_alloc,
                     m=1000, 
                     randomizer = rand_y_2norm_proc, .id = "rho")),
  tar_target(lasso_estimator_model_1_300_pred_results,
             produce_csv(lasso_estimator_model_1_300_pred, 
                         srs_data_store_model_1, "lasso-plugin-var-sim-model-1-300-pred.csv")),
  # Mashregi version
  tar_target(lasso_estimator_mash.boot_model_1_300_pred,
             map_dfr(model_1_data_high_dim_list[["300"]], lasso_reg_est,
                     n_sizes = c(  250, 500),cores = core_alloc,
                     m=1000, boot_yes = TRUE,
                     num_replicates = 250,
                     lasso_reg_boot_var_function = lasso_reg_bootstrap.var_est.mash,
                     randomizer = rand_y_2norm_proc, .id = "rho")),
  tar_target(lasso_estimator_mash.boot_model_1_300_pred_results,
             produce_csv(lasso_estimator_mash.boot_model_1_300_pred, 
                         srs_data_store_model_1, "lasso-mash-boot-sim-model-1-300-pred.csv")),
  
  ### 600 Predictors ----------------
  tar_target(lasso_estimator_model_1_600_pred,
             map_dfr(model_1_data_high_dim_list[["600"]], lasso_reg_var_est ,
                     n_sizes = c(  250, 500),cores = core_alloc,
                     m=1000, 
                     randomizer = rand_y_2norm_proc, .id = "rho"),
             pattern = ),
  tar_target(lasso_estimator_model_1_600_pred_results,
             produce_csv(lasso_estimator_model_1_600_pred, 
                         srs_data_store_model_1, "lasso-plugin-var-sim-model-1-600-pred.csv")),
  
  # Mashregi version
  tar_target(lasso_estimator_mash.boot_model_1_600_pred,
             map_dfr(model_1_data_high_dim_list[["600"]], lasso_reg_est.start.stop,
                     n_sizes = c(  250, 500),cores = core_alloc,
                     m_start = starting_m,
                     m=ending_m, boot_yes = TRUE,
                     num_replicates = 250,
                     lasso_reg_boot_var_function = lasso_reg_bootstrap.var_est.mash,
                     randomizer = rand_y_2norm_proc, .id = "rho"),
             pattern=map(starting_m, ending_m)),
  tar_target(lasso_estimator_mash.boot_model_1_600_pred_results,
             produce_csv(lasso_estimator_mash.boot_model_1_600_pred, 
                         srs_data_store_model_1, "lasso-mash-boot-sim-model-1-600-pred.csv")),
  
  ## ENET -------------------
  ### 25 Predictors --------------
  tar_target(enet_estimator_model_1_25_pred,
             map_dfr(model_1_data_list[["25"]], elastic_net_var_results,
                     n_sizes = c(  250, 500, 1000), cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc,
                     .id = "rho")),
  tar_target(enet_estimator_model_1_25_pred_results,
             produce_csv(enet_estimator_model_1_25_pred, 
                         srs_data_store_model_1, "elastic-net-plugin-var-sim-model-1-25-pred.csv")),
  
  ### 50 Predictors --------------
  tar_target(enet_estimator_model_1_50_pred,
             map_dfr(model_1_data_list[["50"]], elastic_net_var_results,
                     n_sizes = c(  250, 500, 1000), cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc,
                     .id = "rho")),
  tar_target(enet_estimator_model_1_50_pred_results,
             produce_csv(enet_estimator_model_1_50_pred, 
                         srs_data_store_model_1, "elastic-net-plugin-var-sim-model-1-50-pred.csv")),
  
  ### 100 Predictors --------------
  tar_target(enet_estimator_model_1_100_pred,
             map_dfr(model_1_data_list[["100"]], elastic_net_var_results,
                     n_sizes = c(  250, 500, 1000), cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc,
                     .id = "rho")),
  tar_target(enet_estimator_model_1_100_pred_results,
             produce_csv(enet_estimator_model_1_100_pred, 
                         srs_data_store_model_1, "elastic-net-plugin-var-sim-model-1-100-pred.csv")),
  
  ### 200 Predictors --------------
  tar_target(enet_estimator_model_1_200_pred,
             map_dfr(model_1_data_list[["200"]], elastic_net_var_results,
                     n_sizes = c(  250, 500, 1000), cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc,
                     .id = "rho")),
  tar_target(enet_estimator_model_1_200_pred_results,
             produce_csv(enet_estimator_model_1_200_pred, 
                         srs_data_store_model_1, "elastic-net-plugin-var-sim-model-1-200-pred.csv")),
  
  ### 300 Predictors --------------
  tar_target(enet_estimator_model_1_300_pred,
             map_dfr(model_1_data_high_dim_list[["300"]], elastic_net_var_results,
                     n_sizes = c(  250, 500, 1000), cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc,
                     .id = "rho")),
  tar_target(enet_estimator_model_1_300_pred_results,
             produce_csv(enet_estimator_model_1_300_pred, 
                         srs_data_store_model_1, "elastic-net-plugin-var-sim-model-1-300-pred.csv")),
  
  ### 200 Predictors --------------
  tar_target(enet_estimator_model_1_600_pred,
             map_dfr(model_1_data_high_dim_list[["600"]], elastic_net_var_results,
                     n_sizes = c(  250, 500, 1000), cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc,
                     .id = "rho")),
  tar_target(enet_estimator_model_1_600_pred_results,
             produce_csv(enet_estimator_model_1_600_pred, 
                         srs_data_store_model_1, "elastic-net-plugin-var-sim-model-1-600-pred.csv")),
  
  # Model 2 -----------------------------------
  tar_target(model_2_data_list, model_2_data_list_create(N = 5000, data.seedling = 12)),
  tar_target(model_2_data_list_high_dim_list,model_2_data_list_create_high_dim(N = 5000, data.seedling = 12)),
  tar_target(save_model_2_25_pred, saveRDS(model_2_data_list[["25"]],
                                           paste0(srs_data_store_model_2 ,
                                                  "model_2_data_25_pred.rds"))),
  tar_target(save_model_2_50_pred, saveRDS(model_2_data_list[["50"]],
                                           paste0(srs_data_store_model_2 ,
                                                  "model_2_data_50_pred.rds"))),
  tar_target(save_model_2_100_pred, saveRDS(model_2_data_list[["100"]],
                                            paste0(srs_data_store_model_2 ,
                                                   "model_2_data_100_pred.rds"))),
  tar_target(save_model_2_200_pred, saveRDS(model_2_data_list[["200"]],
                                            paste0(srs_data_store_model_2 ,
                                                   "model_2_data_200_pred.rds"))),
  tar_target(save_model_2_300_pred, saveRDS(model_2_data_list_high_dim_list[["300"]],
                                            paste0(srs_data_store_model_2 ,
                                                   "model_2_data_300_pred.rds"))),
  tar_target(save_model_2_600_pred, saveRDS(model_2_data_list_high_dim_list[["600"]],
                                            paste0(srs_data_store_model_2 ,
                                                   "model_2_data_600_pred.rds"))),
  
  
  
  ## HT Estimator -----------------
  ### 25 Predictors -------
  tar_target(ht_estimator_model_2_25_pred,
             map_dfr(model_2_data_list[["25"]], fixed_size_ht_est_variance,
                     n_sizes = c(  250, 500, 1000),cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc, 
                     .id = "rho")),
  tar_target(ht_estimator_model_2_25_pred_results,
             produce_csv(ht_estimator_model_2_25_pred, 
                         srs_data_store_model_2, "ht-plugin-var-sim-model-2-25-pred.csv")),
  ### 50 Predictors ------
  tar_target(ht_estimator_model_2_50_pred,
             map_dfr(model_2_data_list[["50"]], fixed_size_ht_est_variance,
                     n_sizes = c(  250, 500, 1000),cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc, 
                     .id = "rho")),
  tar_target(ht_estimator_model_2_50_pred_results,
             produce_csv(ht_estimator_model_2_50_pred, 
                         srs_data_store_model_2, "ht-plugin-var-sim-model-2-50-pred.csv")),
  ### 100 Predictors ------
  tar_target(ht_estimator_model_2_100_pred,
             map_dfr(model_2_data_list[["100"]], fixed_size_ht_est_variance,
                     n_sizes = c(  250, 500, 1000),cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc, 
                     .id = "rho")),
  tar_target(ht_estimator_model_2_100_pred_results,
             produce_csv(ht_estimator_model_2_100_pred, 
                         srs_data_store_model_2, "ht-plugin-var-sim-model-2-100-pred.csv")),
  ### 200 Predictors ---------------------
  tar_target(ht_estimator_model_2_200_pred,
             map_dfr(model_2_data_list[["200"]], fixed_size_ht_est_variance,
                     n_sizes = c(  250, 500, 1000),cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc, 
                     .id = "rho")),
  tar_target(ht_estimator_model_2_200_pred_results,
             produce_csv(ht_estimator_model_2_200_pred, 
                         srs_data_store_model_2, "ht-plugin-var-sim-model-2-200-pred.csv")),
  
  ### 300 Predictors ---------------------
  tar_target(ht_estimator_model_2_300_pred,
             map_dfr(model_2_data_list_high_dim_list[["300"]], fixed_size_ht_est_variance,
                     n_sizes = c(  250, 500, 1000),cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc, 
                     .id = "rho")),
  tar_target(ht_estimator_model_2_300_pred_results,
             produce_csv(ht_estimator_model_2_300_pred, 
                         srs_data_store_model_2, "ht-plugin-var-sim-model-2-300-pred.csv")),
  
  ### 600 Predictors ---------------------
  tar_target(ht_estimator_model_2_600_pred,
             map_dfr(model_2_data_list_high_dim_list[["600"]], fixed_size_ht_est_variance,
                     n_sizes = c(  250, 500, 1000),cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc, 
                     .id = "rho")),
  tar_target(ht_estimator_model_2_600_pred_results,
             produce_csv(ht_estimator_model_2_600_pred, 
                         srs_data_store_model_2, "ht-plugin-var-sim-model-2-600-pred.csv")),
  
  ## GREG/Linear Regression --------
  ### 25 predictors -----------
  tar_target(linear_reg_estimator_model_2_25_pred,
             map_dfr(model_2_data_list[["25"]], linear_reg_fixed_var,
                     n_sizes = c(  250, 500, 1000), cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc,.id = "rho")),
  tar_target(linear_reg_estimator_model_2_25_pred_results,
             produce_csv(linear_reg_estimator_model_2_25_pred, 
                         srs_data_store_model_2, "linear-plugin-var-sim-model-2-25-pred.csv")),
  ### 50 predictors ----------
  tar_target(linear_reg_estimator_model_2_50_pred,
             map_dfr(model_2_data_list[["50"]], linear_reg_fixed_var,
                     n_sizes = c(  250, 500, 1000), cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc,.id = "rho")),
  tar_target(linear_reg_estimator_model_2_50_pred_results,
             produce_csv(linear_reg_estimator_model_2_50_pred, 
                         srs_data_store_model_2, "linear-plugin-var-sim-model-2-50-pred.csv")),
  
  ### 100 predictors -------
  tar_target(linear_reg_estimator_model_2_100_pred,
             map_dfr(model_2_data_list[["100"]], linear_reg_fixed_var,
                     n_sizes = c(  250, 500, 1000), cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc,.id = "rho")),
  tar_target(linear_reg_estimator_model_2_100_pred_results,
             produce_csv(linear_reg_estimator_model_2_100_pred, 
                         srs_data_store_model_2, "linear-plugin-var-sim-model-2-100-pred.csv")),
  
  ### 200 predictors ------
  tar_target(linear_reg_estimator_model_2_200_pred,
             map_dfr(model_2_data_list[["200"]], linear_reg_fixed_var,
                     n_sizes = c(  250, 500, 1000), cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc,.id = "rho")),
  tar_target(linear_reg_estimator_model_2_200_pred_results,
             produce_csv(linear_reg_estimator_model_2_200_pred, 
                         srs_data_store_model_2, "linear-plugin-var-sim-model-2-200-pred.csv")),
  
  ### 300 predictors ------
  tar_target(linear_reg_estimator_model_2_300_pred,
             map_dfr(model_2_data_list_high_dim_list[["300"]], linear_reg_fixed_var,
                     n_sizes = c(  500, 1000), cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc,.id = "rho")),
  tar_target(linear_reg_estimator_model_2_300_pred_results,
             produce_csv(linear_reg_estimator_model_2_300_pred, 
                         srs_data_store_model_2, "linear-plugin-var-sim-model-2-300-pred.csv")),
  
  
  ## LASSO Estimator ---------------------------
  ### 25 Predictors -------------
  tar_target(lasso_estimator_model_2_25_pred,
             map_dfr(model_2_data_list[["25"]], lasso_reg_var_est ,
                     n_sizes = c(  250, 500, 1000),cores = core_alloc,
                     m=1000, 
                     randomizer = rand_y_2norm_proc, .id = "rho")),
  tar_target(lasso_estimator_model_2_25_pred_results,
             produce_csv(lasso_estimator_model_2_25_pred, 
                         srs_data_store_model_2, "lasso-plugin-var-sim-model-2-25-pred.csv")),
  
  #### Bootstrap Simulations -----------------------
  tar_target(lasso_estimator_bootstrap_model_2_25_pred,
             map_dfr(model_2_data_list[["25"]], lasso_reg_est,
                     n_sizes = c(  250, 500),cores = core_alloc,
                     m=1000, boot_yes = TRUE,
                     randomizer = rand_y_2norm_proc, .id = "rho")),
  tar_target(lasso_estimator_bootstrap_model_2_25_pred_results,
             produce_csv(lasso_estimator_bootstrap_model_2_25_pred,
                         srs_data_store_model_2, "lasso-boot-sim-model-2-25-pred.csv")),
  
  # Mashregi version
  tar_target(lasso_estimator_mash_bootstrap_model_2_25_pred,
             map_dfr(model_2_data_list[["25"]], lasso_reg_est,
                     n_sizes = c(  250, 500),cores = core_alloc,
                     m=1000, boot_yes = TRUE,
                     lasso_reg_boot_var_function = lasso_reg_bootstrap.var_est.mash,
                     randomizer = rand_y_2norm_proc, .id = "rho")),
  tar_target(lasso_estimator_mash_bootstrap_model_2_25_pred_results,
             produce_csv(lasso_estimator_mash_bootstrap_model_2_25_pred,
                         srs_data_store_model_2, "lasso-mash-boot-sim-model-2-25-pred.csv")),
  
  
  ### 50 Predictors -------------
  tar_target(lasso_estimator_model_2_50_pred,
             map_dfr(model_2_data_list[["50"]], lasso_reg_var_est ,
                     n_sizes = c(  250, 500, 1000),cores = core_alloc,
                     m=1000, 
                     randomizer = rand_y_2norm_proc, .id = "rho")),
  tar_target(lasso_estimator_model_2_50_pred_results,
             produce_csv(lasso_estimator_model_2_50_pred, 
                         srs_data_store_model_2, "lasso-plugin-var-sim-model-2-50-pred.csv")),
  
  #### Bootstrap Simulations -----------------------
  tar_target(lasso_estimator_bootstrap_model_2_50_pred,
             map_dfr(model_2_data_list[["50"]], lasso_reg_est,
                     n_sizes = c(  250, 500),cores = core_alloc,
                     m=1000, boot_yes = TRUE,
                     randomizer = rand_y_2norm_proc, .id = "rho")),
  tar_target(lasso_estimator_bootstrap_model_2_50_pred_results,
             produce_csv(lasso_estimator_bootstrap_model_2_50_pred,
                         srs_data_store_model_2, "lasso-boot-sim-model-2-50-pred.csv")),
  
  # Mashregi version
  tar_target(lasso_estimator_mash_bootstrap_model_2_50_pred,
             map_dfr(model_2_data_list[["50"]], lasso_reg_est,
                     n_sizes = c(  250, 500),cores = core_alloc,
                     m=1000, boot_yes = TRUE,
                     lasso_reg_boot_var_function = lasso_reg_bootstrap.var_est.mash,
                     randomizer = rand_y_2norm_proc, .id = "rho")),
  tar_target(lasso_estimator_mash_bootstrap_model_2_50_pred_results,
             produce_csv(lasso_estimator_mash_bootstrap_model_2_50_pred,
                         srs_data_store_model_2, "lasso-mash-boot-sim-model-2-50-pred.csv")),

  
  ### 100 Predictors ----------------
  tar_target(lasso_estimator_model_2_100_pred,
             map_dfr(model_2_data_list[["100"]], lasso_reg_var_est ,
                     n_sizes = c(  250, 500, 1000),cores = core_alloc,
                     m=1000, 
                     randomizer = rand_y_2norm_proc, .id = "rho")),
  tar_target(lasso_estimator_model_2_100_pred_results,
             produce_csv(lasso_estimator_model_2_100_pred, 
                         srs_data_store_model_2, "lasso-plugin-var-sim-model-2-100-pred.csv")),
  
  #### Bootstrap Simulations -----------------------
  tar_target(lasso_estimator_bootstrap_model_2_100_pred,
             map_dfr(model_2_data_list[["100"]], lasso_reg_est,
                     n_sizes = c(  250, 500),cores = core_alloc,
                     m=1000, boot_yes = TRUE,
                     randomizer = rand_y_2norm_proc, .id = "rho")),
  tar_target(lasso_estimator_bootstrap_model_2_100_pred_results,
             produce_csv(lasso_estimator_bootstrap_model_2_100_pred,
                         srs_data_store_model_2, "lasso-boot-sim-model-2-100-pred.csv")),
  
  # Mashregi version
  tar_target(lasso_estimator_mash_bootstrap_model_2_100_pred,
             map_dfr(model_2_data_list[["100"]], lasso_reg_est,
                     n_sizes = c(  250, 500),cores = core_alloc,
                     m=1000, boot_yes = TRUE,
                     lasso_reg_boot_var_function = lasso_reg_bootstrap.var_est.mash,
                     randomizer = rand_y_2norm_proc, .id = "rho")),
  tar_target(lasso_estimator_mash_bootstrap_model_2_100_pred_results,
             produce_csv(lasso_estimator_mash_bootstrap_model_2_100_pred,
                         srs_data_store_model_2, "lasso-mash-boot-sim-model-2-100-pred.csv")),
  
  
  ### 200 Predictors ----------------
  tar_target(lasso_estimator_model_2_200_pred,
             map_dfr(model_2_data_list[["200"]], lasso_reg_var_est ,
                     n_sizes = c(  250, 500, 1000),cores = core_alloc,
                     m=1000, 
                     randomizer = rand_y_2norm_proc, .id = "rho")),
  tar_target(lasso_estimator_model_2_200_pred_results,
             produce_csv(lasso_estimator_model_2_200_pred, 
                         srs_data_store_model_2, "lasso-plugin-var-sim-model-2-200-pred.csv")),
  
  #### Bootstrap Simulations -----------------------
  tar_target(lasso_estimator_bootstrap_model_2_200_pred,
             map_dfr(model_2_data_list[["200"]], lasso_reg_est,
                     n_sizes = c(  250, 500),cores = core_alloc,
                     m=1000, boot_yes = TRUE,
                     randomizer = rand_y_2norm_proc, .id = "rho")),
  tar_target(lasso_estimator_bootstrap_model_2_200_pred_results,
             produce_csv(lasso_estimator_bootstrap_model_2_200_pred,
                         srs_data_store_model_2, "lasso-boot-sim-model-2-200-pred.csv")),
  
  # Mashregi bootstrap
  tar_target(lasso_estimator_mash_bootstrap_model_2_200_pred,
             map_dfr(model_2_data_list[["200"]], lasso_reg_est,
                     n_sizes = c(  250, 500),cores = core_alloc,
                     m=1000, boot_yes = TRUE,
                     lasso_reg_boot_var_function = lasso_reg_bootstrap.var_est.mash,
                     randomizer = rand_y_2norm_proc, .id = "rho")),
  tar_target(lasso_estimator_mash_bootstrap_model_2_200_pred_results,
             produce_csv(lasso_estimator_mash_bootstrap_model_2_200_pred,
                         srs_data_store_model_2, "lasso-mash-boot-sim-model-2-200-pred.csv")),
  ### 300 Predictors ----------------
  tar_target(lasso_estimator_model_2_300_pred,
             map_dfr(model_2_data_list_high_dim_list[["300"]], lasso_reg_var_est ,
                     n_sizes = c(  250, 500, 1000),cores = core_alloc,
                     m=1000, 
                     randomizer = rand_y_2norm_proc, .id = "rho")),
  tar_target(lasso_estimator_model_2_300_pred_results,
             produce_csv(lasso_estimator_model_2_300_pred, 
                         srs_data_store_model_2, "lasso-plugin-var-sim-model-2-300-pred.csv")),
  #### Bootstrap Simulations --------------------
  # Mashregi Bootstrap
  tar_target(lasso_estimator_mash_bootstrap_model_2_300_pred,
             map_dfr(model_2_data_list_high_dim_list[["300"]], lasso_reg_est.start.stop,
                     n_sizes = c(  250, 500),cores = core_alloc,
                     m_start = starting_m,
                     m=ending_m, boot_yes = TRUE,
                     lasso_reg_boot_var_function = lasso_reg_bootstrap.var_est.mash,
                     randomizer = rand_y_2norm_proc, .id = "rho"),
             pattern=map(starting_m, ending_m)),
  tar_target(lasso_estimator_mash_bootstrap_model_2_300_pred_results,
             produce_csv(lasso_estimator_mash_bootstrap_model_2_300_pred,
                         srs_data_store_model_2, "lasso-mash-boot-sim-model-2-300-pred.csv")),
  
  ### 600 Predictors ----------------
  tar_target(lasso_estimator_model_2_600_pred,
             map_dfr(model_2_data_list_high_dim_list[["600"]], lasso_reg_var_est ,
                     n_sizes = c(  250, 500, 1000),cores = core_alloc,
                     m=1000, 
                     randomizer = rand_y_2norm_proc, .id = "rho")),
  tar_target(lasso_estimator_model_2_600_pred_results,
             produce_csv(lasso_estimator_model_2_600_pred, 
                         srs_data_store_model_2, "lasso-plugin-var-sim-model-2-600-pred.csv")),
  #### Bootstrap Simulations --------------------
  # Mashregi Bootstrap
  tar_target(lasso_estimator_mash_bootstrap_model_2_600_pred,
             map_dfr(model_2_data_list_high_dim_list[["600"]], lasso_reg_est.start.stop,
                     n_sizes = c(  250, 500),cores = core_alloc,
                     m_start = starting_m,
                     m=ending_m, boot_yes = TRUE,
                     lasso_reg_boot_var_function = lasso_reg_bootstrap.var_est.mash,
                     randomizer = rand_y_2norm_proc, .id = "rho"),
             pattern=map(starting_m, ending_m)),
  tar_target(lasso_estimator_mash_bootstrap_model_2_600_pred_results,
             produce_csv(lasso_estimator_mash_bootstrap_model_2_600_pred,
                         srs_data_store_model_2, "lasso-mash-boot-sim-model-2-600-pred.csv")),
  
  
  
  
  ## ENET -------------------
  ### 25 Predictors --------------
  tar_target(enet_estimator_model_2_25_pred,
             map_dfr(model_2_data_list[["25"]], elastic_net_var_results,
                     n_sizes = c(  250, 500, 1000), cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc,
                     .id = "rho")),
  tar_target(enet_estimator_model_2_25_pred_results,
             produce_csv(enet_estimator_model_2_25_pred, 
                         srs_data_store_model_2, "elastic-net-plugin-var-sim-model-2-25-pred.csv")),
  
  ### 50 Predictors --------------
  tar_target(enet_estimator_model_2_50_pred,
             map_dfr(model_2_data_list[["50"]], elastic_net_var_results,
                     n_sizes = c(  250, 500, 1000), cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc,
                     .id = "rho")),
  tar_target(enet_estimator_model_2_50_pred_results,
             produce_csv(enet_estimator_model_2_50_pred, 
                         srs_data_store_model_2, "elastic-net-plugin-var-sim-model-2-50-pred.csv")),
  
  ### 100 Predictors --------------
  tar_target(enet_estimator_model_2_100_pred,
             map_dfr(model_2_data_list[["100"]], elastic_net_var_results,
                     n_sizes = c(  250, 500, 1000), cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc,
                     .id = "rho")),
  tar_target(enet_estimator_model_2_100_pred_results,
             produce_csv(enet_estimator_model_2_100_pred, 
                         srs_data_store_model_2, "elastic-net-plugin-var-sim-model-2-100-pred.csv")),
  
  ### 200 Predictors --------------
  tar_target(enet_estimator_model_2_200_pred,
             map_dfr(model_2_data_list[["200"]], elastic_net_var_results,
                     n_sizes = c(  250, 500, 1000), cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc,
                     .id = "rho")),
  tar_target(enet_estimator_model_2_200_pred_results,
             produce_csv(enet_estimator_model_2_200_pred, 
                         srs_data_store_model_2, "elastic-net-plugin-var-sim-model-2-200-pred.csv")),
  
  ### 300 Predictors --------------
  tar_target(enet_estimator_model_2_300_pred,
             map_dfr(model_2_data_list_high_dim_list[["300"]], elastic_net_var_results,
                     n_sizes = c(  250, 500, 1000), cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc,
                     .id = "rho")),
  tar_target(enet_estimator_model_2_300_pred_results,
             produce_csv(enet_estimator_model_2_300_pred, 
                         srs_data_store_model_2, "elastic-net-plugin-var-sim-model-2-300-pred.csv")),
  
  
  ### 600 Predictors --------------
  tar_target(enet_estimator_model_2_600_pred,
             map_dfr(model_2_data_list_high_dim_list[["600"]], elastic_net_var_results,
                     n_sizes = c(  250, 500, 1000), cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc,
                     .id = "rho")),
  tar_target(enet_estimator_model_2_600_pred_results,
             produce_csv(enet_estimator_model_2_600_pred, 
                         srs_data_store_model_2, "elastic-net-plugin-var-sim-model-2-600-pred.csv")),
  
  
  # Model 3 ----------------------
  tar_target(model_3_data_list, model_3_data_list_create(N = 5000, data.seedling = 12)),
  tar_target(model_3_data_high_dim_list, model_3_data_high_dim_list_create(N = 5000, data.seedling = 12)),
  tar_target(save_model_3_25_pred, saveRDS(model_3_data_list[["25"]],
                                           paste0(srs_data_store_model_3 ,
                                                  "model_3_data_25_pred.rds"))),
  tar_target(save_model_3_50_pred, saveRDS(model_3_data_list[["50"]],
                                           paste0(srs_data_store_model_3 ,
                                                  "model_3_data_50_pred.rds"))),
  tar_target(save_model_3_100_pred, saveRDS(model_3_data_list[["100"]],
                                            paste0(srs_data_store_model_3 ,
                                                   "model_3_data_100_pred.rds"))),
  tar_target(save_model_3_200_pred, saveRDS(model_3_data_list[["200"]],
                                            paste0(srs_data_store_model_3 ,
                                                   "model_3_data_200_pred.rds"))),
  tar_target(save_model_3_300_pred, saveRDS(model_3_data_high_dim_list[["300"]],
                                            paste0(srs_data_store_model_3 ,
                                                   "model_3_data_300_pred.rds"))),
  tar_target(save_model_3_600_pred, saveRDS(model_3_data_high_dim_list[["600"]],
                                            paste0(srs_data_store_model_3 ,
                                                   "model_3_data_600_pred.rds"))),
  

  
  
  
  ## HT Estimator --------
  ### 25 Predictors -------
  tar_target(ht_estimator_model_3_25_pred,
             map_dfr(model_3_data_list[["25"]], fixed_size_ht_est_variance,
                     n_sizes = c(  250, 500, 1000),cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc, 
                     .id = "rho")),
  tar_target(ht_estimator_model_3_25_pred_results,
             produce_csv(ht_estimator_model_3_25_pred, 
                         srs_data_store_model_3, "ht-plugin-var-sim-model-3-25-pred.csv")),
  ### 50 Predictors ------
  tar_target(ht_estimator_model_3_50_pred,
             map_dfr(model_3_data_list[["50"]], fixed_size_ht_est_variance,
                     n_sizes = c(  250, 500, 1000),cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc, 
                     .id = "rho")),
  tar_target(ht_estimator_model_3_50_pred_results,
             produce_csv(ht_estimator_model_3_50_pred, 
                         srs_data_store_model_3, "ht-plugin-var-sim-model-3-50-pred.csv")),
  ### 100 Predictors ------
  tar_target(ht_estimator_model_3_100_pred,
             map_dfr(model_3_data_list[["100"]], fixed_size_ht_est_variance,
                     n_sizes = c(  250, 500, 1000),cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc, 
                     .id = "rho")),
  tar_target(ht_estimator_model_3_100_pred_results,
             produce_csv(ht_estimator_model_3_100_pred, 
                         srs_data_store_model_3, "ht-plugin-var-sim-model-3-100-pred.csv")),
  ### 200 Predictors ---------------------
  tar_target(ht_estimator_model_3_200_pred,
             map_dfr(model_3_data_list[["200"]], fixed_size_ht_est_variance,
                     n_sizes = c(  250, 500, 1000),cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc, 
                     .id = "rho")),
  tar_target(ht_estimator_model_3_200_pred_results,
             produce_csv(ht_estimator_model_3_200_pred, 
                         srs_data_store_model_3, "ht-plugin-var-sim-model-3-200-pred.csv")),
  
  
  ### 300 Predictors ---------------------
  tar_target(ht_estimator_model_3_300_pred,
             map_dfr(model_3_data_high_dim_list[["300"]], fixed_size_ht_est_variance,
                     n_sizes = c(  250, 500, 1000),cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc, 
                     .id = "rho")),
  tar_target(ht_estimator_model_3_300_pred_results,
             produce_csv(ht_estimator_model_3_300_pred, 
                         srs_data_store_model_3, "ht-plugin-var-sim-model-3-300-pred.csv")),
  
  ### 600 Predictors ---------------------
  tar_target(ht_estimator_model_3_600_pred,
             map_dfr(model_3_data_high_dim_list[["600"]], fixed_size_ht_est_variance,
                     n_sizes = c(  250, 500, 1000),cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc, 
                     .id = "rho")),
  tar_target(ht_estimator_model_3_600_pred_results,
             produce_csv(ht_estimator_model_3_600_pred, 
                         srs_data_store_model_3, "ht-plugin-var-sim-model-3-600-pred.csv")),
  
  
  
  ## GREG/Linear Regression --------
  ### 25 predictors -----------
  tar_target(linear_reg_estimator_model_3_25_pred,
             map_dfr(model_3_data_list[["25"]], linear_reg_fixed_var,
                     n_sizes = c(  250, 500, 1000), cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc,.id = "rho")),
  tar_target(linear_reg_estimator_model_3_25_pred_results,
             produce_csv(linear_reg_estimator_model_3_25_pred, 
                         srs_data_store_model_3, "linear-plugin-var-sim-model-3-25-pred.csv")),
  ### 50 predictors ----------
  tar_target(linear_reg_estimator_model_3_50_pred,
             map_dfr(model_3_data_list[["50"]], linear_reg_fixed_var,
                     n_sizes = c(  250, 500, 1000), cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc,.id = "rho")),
  tar_target(linear_reg_estimator_model_3_50_pred_results,
             produce_csv(linear_reg_estimator_model_3_50_pred, 
                         srs_data_store_model_3, "linear-plugin-var-sim-model-3-50-pred.csv")),
  
  ### 100 predictors -------
  tar_target(linear_reg_estimator_model_3_100_pred,
             map_dfr(model_3_data_list[["100"]], linear_reg_fixed_var,
                     n_sizes = c(  250, 500, 1000), cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc,.id = "rho")),
  tar_target(linear_reg_estimator_model_3_100_pred_results,
             produce_csv(linear_reg_estimator_model_3_100_pred, 
                         srs_data_store_model_3, "linear-plugin-var-sim-model-3-100-pred.csv")),
  
  ### 200 predictors ------
  tar_target(linear_reg_estimator_model_3_200_pred,
             map_dfr(model_3_data_list[["200"]], linear_reg_fixed_var,
                     n_sizes = c(  250, 500, 1000), cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc,.id = "rho")),
  tar_target(linear_reg_estimator_model_3_200_pred_results,
             produce_csv(linear_reg_estimator_model_3_200_pred, 
                         srs_data_store_model_3, "linear-plugin-var-sim-model-3-200-pred.csv")),
  
  ### 300 predictors ------
  tar_target(linear_reg_estimator_model_3_300_pred,
             map_dfr(model_3_data_high_dim_list[["300"]], linear_reg_fixed_var,
                     n_sizes = c( 500, 1000), cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc,.id = "rho")),
  tar_target(linear_reg_estimator_model_3_300_pred_results,
             produce_csv(linear_reg_estimator_model_3_300_pred, 
                         srs_data_store_model_3, "linear-plugin-var-sim-model-3-300-pred.csv")),
  
  
  
  ## LASSO Estimator ---------------------------
  ### 25 Predictors -------------
  tar_target(lasso_estimator_model_3_25_pred,
             map_dfr(model_3_data_list[["25"]], lasso_reg_var_est ,
                     n_sizes = c(  250, 500, 1000),cores = core_alloc,
                     m=1000, 
                     randomizer = rand_y_2norm_proc, .id = "rho")),
  tar_target(lasso_estimator_model_3_25_pred_results,
             produce_csv(lasso_estimator_model_3_25_pred, 
                         srs_data_store_model_3, "lasso-plugin-var-sim-model-3-25-pred.csv")),
  
  #### Bootstrap Simulations -----------------------
  tar_target(lasso_estimator_bootstrap_model_3_25_pred,
             map_dfr(model_3_data_list[["25"]], lasso_reg_est,
                     n_sizes = c(  250, 500),cores = core_alloc,
                     m=1000, boot_yes = TRUE,
                     randomizer = rand_y_2norm_proc, .id = "rho")),
  tar_target(lasso_estimator_bootstrap_model_3_25_pred_results,
             produce_csv(lasso_estimator_bootstrap_model_3_25_pred,
                         srs_data_store_model_3, "lasso-boot-sim-model-3-25-pred.csv")),
  
  # Mashregi bootstrap
  tar_target(lasso_estimator_mash_bootstrap_model_3_25_pred,
             map_dfr(model_3_data_list[["25"]], lasso_reg_est,
                     n_sizes = c(  250, 500),cores = core_alloc,
                     m=1000, boot_yes = TRUE,
                     lasso_reg_boot_var_function = lasso_reg_bootstrap.var_est.mash,
                     randomizer = rand_y_2norm_proc, .id = "rho")),
  tar_target(lasso_estimator_mash_bootstrap_model_3_25_pred_results,
             produce_csv(lasso_estimator_mash_bootstrap_model_3_25_pred,
                         srs_data_store_model_3, "lasso-mash-boot-sim-model-3-25-pred.csv")),
  
  ### 50 Predictors -------------
  tar_target(lasso_estimator_model_3_50_pred,
             map_dfr(model_3_data_list[["50"]], lasso_reg_var_est ,
                     n_sizes = c(  250, 500, 1000),cores = core_alloc,
                     m=1000, 
                     randomizer = rand_y_2norm_proc, .id = "rho")),
  tar_target(lasso_estimator_model_3_50_pred_results,
             produce_csv(lasso_estimator_model_3_50_pred, 
                         srs_data_store_model_3, "lasso-plugin-var-sim-model-3-50-pred.csv")),
  
  #### Bootstrap Simulations -----------------------
  tar_target(lasso_estimator_bootstrap_model_3_50_pred,
             map_dfr(model_3_data_list[["50"]], lasso_reg_est,
                     n_sizes = c(  250, 500),cores = core_alloc,
                     m=1000, boot_yes = TRUE,
                     randomizer = rand_y_2norm_proc, .id = "rho")),
  tar_target(lasso_estimator_bootstrap_model_3_50_pred_results,
             produce_csv(lasso_estimator_bootstrap_model_3_50_pred,
                         srs_data_store_model_3, "lasso-boot-sim-model-3-50-pred.csv")),
  
  # Mashregi bootstrap
  tar_target(lasso_estimator_mash_bootstrap_model_3_50_pred,
             map_dfr(model_3_data_list[["50"]], lasso_reg_est,
                     n_sizes = c(  250, 500),cores = core_alloc,
                     m=1000, boot_yes = TRUE,
                     lasso_reg_boot_var_function = lasso_reg_bootstrap.var_est.mash,
                     randomizer = rand_y_2norm_proc, .id = "rho")),
  tar_target(lasso_estimator_mash_bootstrap_model_3_50_pred_results,
             produce_csv(lasso_estimator_mash_bootstrap_model_3_50_pred,
                         srs_data_store_model_3, "lasso-mash-boot-sim-model-3-50-pred.csv")),
  
  
  ### 100 Predictors ----------------
  tar_target(lasso_estimator_model_3_100_pred,
             map_dfr(model_3_data_list[["100"]], lasso_reg_var_est ,
                     n_sizes = c(  250, 500, 1000),cores = core_alloc,
                     m=1000, 
                     randomizer = rand_y_2norm_proc, .id = "rho")),
  tar_target(lasso_estimator_model_3_100_pred_results,
             produce_csv(lasso_estimator_model_3_100_pred, 
                         srs_data_store_model_3, "lasso-plugin-var-sim-model-3-100-pred.csv")),
  
  #### Bootstrap Simulations -----------------------
  tar_target(lasso_estimator_bootstrap_model_3_100_pred,
             map_dfr(model_3_data_list[["100"]], lasso_reg_est,
                     n_sizes = c(  250, 500),cores = core_alloc,
                     m=1000, boot_yes = TRUE,
                     randomizer = rand_y_2norm_proc, .id = "rho")),
  tar_target(lasso_estimator_bootstrap_model_3_100_pred_results,
             produce_csv(lasso_estimator_bootstrap_model_3_100_pred,
                         srs_data_store_model_3, "lasso-boot-sim-model-3-100-pred.csv")),
  
  # Mashregi bootstrap
  tar_target(lasso_estimator_mash_bootstrap_model_3_100_pred,
             map_dfr(model_3_data_list[["100"]], lasso_reg_est,
                     n_sizes = c(  250, 500),cores = core_alloc,
                     m=1000, boot_yes = TRUE,
                     lasso_reg_boot_var_function = lasso_reg_bootstrap.var_est.mash,
                     randomizer = rand_y_2norm_proc, .id = "rho")),
  tar_target(lasso_estimator_mash_bootstrap_model_3_100_pred_results,
             produce_csv(lasso_estimator_mash_bootstrap_model_3_100_pred,
                         srs_data_store_model_3, "lasso-mash-boot-sim-model-3-100-pred.csv")),
  
  ### 200 Predictors ----------------
  tar_target(lasso_estimator_model_3_200_pred,
             map_dfr(model_3_data_list[["200"]], lasso_reg_var_est ,
                     n_sizes = c(  250, 500, 1000),cores = core_alloc,
                     m=1000, 
                     randomizer = rand_y_2norm_proc, .id = "rho")),
  tar_target(lasso_estimator_model_3_200_pred_results,
             produce_csv(lasso_estimator_model_3_200_pred, 
                         srs_data_store_model_3, "lasso-plugin-var-sim-model-3-200-pred.csv")),
  
  #### Bootstrap Simulations -----------------------
  tar_target(lasso_estimator_bootstrap_model_3_200_pred,
             map_dfr(model_3_data_list[["200"]], lasso_reg_est,
                     n_sizes = c(  250, 500),cores = core_alloc,
                     m=1000, boot_yes = TRUE,
                     randomizer = rand_y_2norm_proc, .id = "rho")),
  tar_target(lasso_estimator_bootstrap_model_3_200_pred_results,
             produce_csv(lasso_estimator_bootstrap_model_3_200_pred,
                         srs_data_store_model_3, "lasso-boot-sim-model-3-200-pred.csv")),
  
  # Mashregi bootstrap
  tar_target(lasso_estimator_mash_bootstrap_model_3_200_pred,
             map_dfr(model_3_data_list[["200"]], lasso_reg_est,
                     n_sizes = c(  250, 500),cores = core_alloc,
                     m=1000, boot_yes = TRUE,
                     lasso_reg_boot_var_function = lasso_reg_bootstrap.var_est.mash,
                     randomizer = rand_y_2norm_proc, .id = "rho")),
  tar_target(lasso_estimator_mash_bootstrap_model_3_200_pred_results,
             produce_csv(lasso_estimator_mash_bootstrap_model_3_200_pred,
                         srs_data_store_model_3, "lasso-mash-boot-sim-model-3-200-pred.csv")),
  
  ### 300 Predictors ----------------
  tar_target(lasso_estimator_model_3_300_pred,
             map_dfr(model_3_data_high_dim_list[["300"]], lasso_reg_var_est ,
                     n_sizes = c(  250, 500, 1000),cores = core_alloc,
                     m=1000, 
                     randomizer = rand_y_2norm_proc, .id = "rho")),
  tar_target(lasso_estimator_model_3_300_pred_results,
             produce_csv(lasso_estimator_model_3_300_pred, 
                         srs_data_store_model_3, "lasso-plugin-var-sim-model-3-300-pred.csv")),
  
  #### Bootstrap Simulations -----------------------
  # Mashregi bootstrap
  tar_target(lasso_estimator_mash_bootstrap_model_3_300_pred,
             map_dfr(model_3_data_high_dim_list[["300"]], lasso_reg_est.start.stop,
                     n_sizes = c(  250, 500),cores = core_alloc,
                     m_start = starting_m,
                     m=ending_m, boot_yes = TRUE,
                     lasso_reg_boot_var_function = lasso_reg_bootstrap.var_est.mash,
                     randomizer = rand_y_2norm_proc, .id = "rho"),
             pattern=map(starting_m, ending_m)),
  tar_target(lasso_estimator_mash_bootstrap_model_3_300_pred_results,
             produce_csv(lasso_estimator_mash_bootstrap_model_3_300_pred,
                         srs_data_store_model_3, "lasso-mash-boot-sim-model-3-300-pred.csv")),
  
  ### 600 Predictors ----------------
  tar_target(lasso_estimator_model_3_600_pred,
             map_dfr(model_3_data_high_dim_list[["600"]], lasso_reg_var_est ,
                     n_sizes = c(  250, 500, 1000),cores = core_alloc,
                     m=1000, 
                     randomizer = rand_y_2norm_proc, .id = "rho")),
  tar_target(lasso_estimator_model_3_600_pred_results,
             produce_csv(lasso_estimator_model_3_600_pred, 
                         srs_data_store_model_3, "lasso-plugin-var-sim-model-3-600-pred.csv")),
  
  #### Bootstrap Simulations -----------------------
  # Mashregi bootstrap
  tar_target(lasso_estimator_mash_bootstrap_model_3_600_pred,
             map_dfr(model_3_data_high_dim_list[["600"]], lasso_reg_est.start.stop,
                     n_sizes = c(  250, 500),cores = core_alloc,
                     m_start = starting_m, m = ending_m, boot_yes = TRUE,
                     lasso_reg_boot_var_function = lasso_reg_bootstrap.var_est.mash,
                     randomizer = rand_y_2norm_proc, .id = "rho"),
             pattern=map(starting_m, ending_m)),
  tar_target(lasso_estimator_mash_bootstrap_model_3_600_pred_results,
             produce_csv(lasso_estimator_mash_bootstrap_model_3_600_pred,
                         srs_data_store_model_3, "lasso-mash-boot-sim-model-3-600-pred.csv")),
  
  
  
  
  
  
  
  ## ENET -------------------
  ### 25 Predictors --------------
  tar_target(enet_estimator_model_3_25_pred,
             map_dfr(model_3_data_list[["25"]], elastic_net_var_results,
                     n_sizes = c(  250, 500, 1000), cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc,
                     .id = "rho")),
  tar_target(enet_estimator_model_3_25_pred_results,
             produce_csv(enet_estimator_model_3_25_pred, 
                         srs_data_store_model_3, "elastic-net-plugin-var-sim-model-3-25-pred.csv")),
  
  ### 50 Predictors --------------
  tar_target(enet_estimator_model_3_50_pred,
             map_dfr(model_3_data_list[["50"]], elastic_net_var_results,
                     n_sizes = c(  250, 500, 1000), cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc,
                     .id = "rho")),
  tar_target(enet_estimator_model_3_50_pred_results,
             produce_csv(enet_estimator_model_3_50_pred, 
                         srs_data_store_model_3, "elastic-net-plugin-var-sim-model-3-50-pred.csv")),
  
  ### 100 Predictors --------------
  tar_target(enet_estimator_model_3_100_pred,
             map_dfr(model_3_data_list[["100"]], elastic_net_var_results,
                     n_sizes = c(  250, 500, 1000), cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc,
                     .id = "rho")),
  tar_target(enet_estimator_model_3_100_pred_results,
             produce_csv(enet_estimator_model_3_100_pred, 
                         srs_data_store_model_3, "elastic-net-plugin-var-sim-model-3-100-pred.csv")),
  
  ### 200 Predictors --------------
  tar_target(enet_estimator_model_3_200_pred,
             map_dfr(model_3_data_list[["200"]], elastic_net_var_results,
                     n_sizes = c(  250, 500, 1000), cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc,
                     .id = "rho")),
  tar_target(enet_estimator_model_3_200_pred_results,
             produce_csv(enet_estimator_model_3_200_pred, 
                         srs_data_store_model_3, "elastic-net-plugin-var-sim-model-3-200-pred.csv")),
  
  ### 300 Predictors --------------
  tar_target(enet_estimator_model_3_300_pred,
             map_dfr(model_3_data_high_dim_list[["300"]], elastic_net_var_results,
                     n_sizes = c(  250, 500, 1000), cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc,
                     .id = "rho")),
  tar_target(enet_estimator_model_3_300_pred_results,
             produce_csv(enet_estimator_model_3_300_pred, 
                         srs_data_store_model_3, "elastic-net-plugin-var-sim-model-3-300-pred.csv")),
  
  ### 600 Predictors --------------
  tar_target(enet_estimator_model_3_600_pred,
             map_dfr(model_3_data_high_dim_list[["600"]], elastic_net_var_results,
                     n_sizes = c(  250, 500, 1000), cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc,
                     .id = "rho")),
  tar_target(enet_estimator_model_3_600_pred_results,
             produce_csv(enet_estimator_model_3_600_pred, 
                         srs_data_store_model_3, "elastic-net-plugin-var-sim-model-3-600-pred.csv")),
  
  
  # Model 4 --------------------------
  tar_target(model_4_data_list, model_4_data_list_create(N = 5000, data.seedling = 12)),
  tar_target(model_4_data_high_dim_list,model_4_data_high_dim_list_create(N = 5000, data.seedling = 12)),
  tar_target(save_model_4_25_pred, saveRDS(model_4_data_list[["25"]],
                                           paste0(srs_data_store_model_4 ,
                                                  "model_4_data_25_pred.rds"))),
  tar_target(save_model_4_50_pred, saveRDS(model_4_data_list[["50"]],
                                           paste0(srs_data_store_model_4 ,
                                                  "model_4_data_50_pred.rds"))),
  tar_target(save_model_4_100_pred, saveRDS(model_4_data_list[["100"]],
                                            paste0(srs_data_store_model_4 ,
                                                   "model_4_data_100_pred.rds"))),
  tar_target(save_model_4_200_pred, saveRDS(model_4_data_list[["200"]],
                                            paste0(srs_data_store_model_4 ,
                                                   "model_4_data_200_pred.rds"))),
  tar_target(save_model_4_300_pred, saveRDS(model_4_data_high_dim_list[["300"]],
                                            paste0(srs_data_store_model_4 ,
                                                   "model_4_data_300_pred.rds"))),
  tar_target(save_model_4_600_pred, saveRDS(model_4_data_high_dim_list[["600"]],
                                            paste0(srs_data_store_model_4 ,
                                                   "model_4_data_600_pred.rds"))),
  
  
  
  
  ## HT Estimator --------
  ### 25 Predictors -------
  tar_target(ht_estimator_model_4_25_pred,
             map_dfr(model_4_data_list[["25"]], fixed_size_ht_est_variance,
                     n_sizes = c(  250, 500, 1000),cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc, 
                     .id = "rho")),
  tar_target(ht_estimator_model_4_25_pred_results,
             produce_csv(ht_estimator_model_4_25_pred, 
                         srs_data_store_model_4, "ht-plugin-var-sim-model-4-25-pred.csv")),
  ### 50 Predictors ------
  tar_target(ht_estimator_model_4_50_pred,
             map_dfr(model_4_data_list[["50"]], fixed_size_ht_est_variance,
                     n_sizes = c(  250, 500, 1000),cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc, 
                     .id = "rho")),
  tar_target(ht_estimator_model_4_50_pred_results,
             produce_csv(ht_estimator_model_4_50_pred, 
                         srs_data_store_model_4, "ht-plugin-var-sim-model-4-50-pred.csv")),
  ### 100 Predictors ------
  tar_target(ht_estimator_model_4_100_pred,
             map_dfr(model_4_data_list[["100"]], fixed_size_ht_est_variance,
                     n_sizes = c(  250, 500, 1000),cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc, 
                     .id = "rho")),
  tar_target(ht_estimator_model_4_100_pred_results,
             produce_csv(ht_estimator_model_4_100_pred, 
                         srs_data_store_model_4, "ht-plugin-var-sim-model-4-100-pred.csv")),
  ### 200 Predictors ---------------------
  tar_target(ht_estimator_model_4_200_pred,
             map_dfr(model_4_data_list[["200"]], fixed_size_ht_est_variance,
                     n_sizes = c(  250, 500, 1000),cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc, 
                     .id = "rho")),
  tar_target(ht_estimator_model_4_200_pred_results,
             produce_csv(ht_estimator_model_4_200_pred, 
                         srs_data_store_model_4, "ht-plugin-var-sim-model-4-200-pred.csv")),
  
  ### 300 Predictors ---------------------
  tar_target(ht_estimator_model_4_300_pred,
             map_dfr(model_4_data_high_dim_list[["300"]], fixed_size_ht_est_variance,
                     n_sizes = c(  250, 500, 1000),cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc, 
                     .id = "rho")),
  tar_target(ht_estimator_model_4_300_pred_results,
             produce_csv(ht_estimator_model_4_300_pred, 
                         srs_data_store_model_4, "ht-plugin-var-sim-model-4-300-pred.csv")),
  
  ### 600 Predictors ---------------------
  tar_target(ht_estimator_model_4_600_pred,
             map_dfr(model_4_data_high_dim_list[["600"]], fixed_size_ht_est_variance,
                     n_sizes = c(  250, 500, 1000),cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc, 
                     .id = "rho")),
  tar_target(ht_estimator_model_4_600_pred_results,
             produce_csv(ht_estimator_model_4_600_pred, 
                         srs_data_store_model_4, "ht-plugin-var-sim-model-4-600-pred.csv")),
  
  
  
  ## GREG/Linear Regression --------
  ### 25 predictors -----------
  tar_target(linear_reg_estimator_model_4_25_pred,
             map_dfr(model_4_data_list[["25"]], linear_reg_fixed_var,
                     n_sizes = c(  250, 500, 1000), cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc,.id = "rho")),
  tar_target(linear_reg_estimator_model_4_25_pred_results,
             produce_csv(linear_reg_estimator_model_4_25_pred, 
                         srs_data_store_model_4, "linear-plugin-var-sim-model-4-25-pred.csv")),
  ### 50 predictors ----------
  tar_target(linear_reg_estimator_model_4_50_pred,
             map_dfr(model_4_data_list[["50"]], linear_reg_fixed_var,
                     n_sizes = c(  250, 500, 1000), cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc,.id = "rho")),
  tar_target(linear_reg_estimator_model_4_50_pred_results,
             produce_csv(linear_reg_estimator_model_4_50_pred, 
                         srs_data_store_model_4, "linear-plugin-var-sim-model-4-50-pred.csv")),
  
  ### 100 predictors -------
  tar_target(linear_reg_estimator_model_4_100_pred,
             map_dfr(model_4_data_list[["100"]], linear_reg_fixed_var,
                     n_sizes = c(  250, 500, 1000), cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc,.id = "rho")),
  tar_target(linear_reg_estimator_model_4_100_pred_results,
             produce_csv(linear_reg_estimator_model_4_100_pred, 
                         srs_data_store_model_4, "linear-plugin-var-sim-model-4-100-pred.csv")),
  
  ### 200 predictors ------
  tar_target(linear_reg_estimator_model_4_200_pred,
             map_dfr(model_4_data_list[["200"]], linear_reg_fixed_var,
                     n_sizes = c(  250, 500, 1000), cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc,.id = "rho")),
  tar_target(linear_reg_estimator_model_4_200_pred_results,
             produce_csv(linear_reg_estimator_model_4_200_pred, 
                         srs_data_store_model_4, "linear-plugin-var-sim-model-4-200-pred.csv")),
  
  ### 300 predictors ------
  tar_target(linear_reg_estimator_model_4_300_pred,
             map_dfr(model_4_data_high_dim_list[["300"]], linear_reg_fixed_var,
                     n_sizes = c( 500, 1000), cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc,.id = "rho")),
  tar_target(linear_reg_estimator_model_4_300_pred_results,
             produce_csv(linear_reg_estimator_model_4_300_pred, 
                         srs_data_store_model_4, "linear-plugin-var-sim-model-4-300-pred.csv")),
  
  
  
  ## LASSO Estimator ---------------------------
  ### 25 Predictors -------------
  tar_target(lasso_estimator_model_4_25_pred,
             map_dfr(model_4_data_list[["25"]], lasso_reg_var_est ,
                     n_sizes = c(  250, 500, 1000),cores = core_alloc,
                     m=1000, 
                     randomizer = rand_y_2norm_proc, .id = "rho")),
  tar_target(lasso_estimator_model_4_25_pred_results,
             produce_csv(lasso_estimator_model_4_25_pred, 
                         srs_data_store_model_4, "lasso-plugin-var-sim-model-4-25-pred.csv")),
  
  #### Bootstrap Simulations -----------------------
  tar_target(lasso_estimator_bootstrap_model_4_25_pred,
             map_dfr(model_4_data_list[["25"]], lasso_reg_est,
                     n_sizes = c(  250, 500),cores = core_alloc,
                     m=1000, boot_yes = TRUE,
                     randomizer = rand_y_2norm_proc, .id = "rho")),
  tar_target(lasso_estimator_bootstrap_model_4_25_pred_results,
             produce_csv(lasso_estimator_bootstrap_model_4_25_pred,
                         srs_data_store_model_4, "lasso-boot-sim-model-4-25-pred.csv")),
  
  # Mashregi bootstrap
  tar_target(lasso_estimator_mash_bootstrap_model_4_25_pred,
             map_dfr(model_4_data_list[["25"]], lasso_reg_est,
                     n_sizes = c(  250, 500),cores = core_alloc,
                     m=1000, boot_yes = TRUE,
                     lasso_reg_boot_var_function = lasso_reg_bootstrap.var_est.mash,
                     randomizer = rand_y_2norm_proc, .id = "rho")),
  tar_target(lasso_estimator_mash_bootstrap_model_4_25_pred_results,
             produce_csv(lasso_estimator_mash_bootstrap_model_4_25_pred,
                         srs_data_store_model_4, "lasso-mash-boot-sim-model-4-25-pred.csv")),
  
  ### 50 Predictors -------------
  tar_target(lasso_estimator_model_4_50_pred,
             map_dfr(model_4_data_list[["50"]], lasso_reg_var_est ,
                     n_sizes = c(  250, 500, 1000),cores = core_alloc,
                     m=1000, 
                     randomizer = rand_y_2norm_proc, .id = "rho")),
  tar_target(lasso_estimator_model_4_50_pred_results,
             produce_csv(lasso_estimator_model_4_50_pred, 
                         srs_data_store_model_4, "lasso-plugin-var-sim-model-4-50-pred.csv")),
  
  #### Bootstrap Simulations -----------------------
  tar_target(lasso_estimator_bootstrap_model_4_50_pred,
             map_dfr(model_4_data_list[["50"]], lasso_reg_est,
                     n_sizes = c(  250, 500),cores = core_alloc,
                     m=1000, boot_yes = TRUE,
                     randomizer = rand_y_2norm_proc, .id = "rho")),
  tar_target(lasso_estimator_bootstrap_model_4_50_pred_results,
             produce_csv(lasso_estimator_bootstrap_model_4_50_pred,
                         srs_data_store_model_4, "lasso-boot-sim-model-4-50-pred.csv")),
  
  # Mashregi bootstrap
  tar_target(lasso_estimator_mash_bootstrap_model_4_50_pred,
             map_dfr(model_4_data_list[["50"]], lasso_reg_est,
                     n_sizes = c(  250, 500),cores = core_alloc,
                     m=1000, boot_yes = TRUE,
                     lasso_reg_boot_var_function = lasso_reg_bootstrap.var_est.mash,
                     randomizer = rand_y_2norm_proc, .id = "rho")),
  tar_target(lasso_estimator_mash_bootstrap_model_4_50_pred_results,
             produce_csv(lasso_estimator_mash_bootstrap_model_4_50_pred,
                         srs_data_store_model_4, "lasso-mash-boot-sim-model-4-50-pred.csv")),
  
  
  ### 100 Predictors ----------------
  tar_target(lasso_estimator_model_4_100_pred,
             map_dfr(model_4_data_list[["100"]], lasso_reg_var_est ,
                     n_sizes = c(  250, 500, 1000),cores = core_alloc,
                     m=1000, 
                     randomizer = rand_y_2norm_proc, .id = "rho")),
  tar_target(lasso_estimator_model_4_100_pred_results,
             produce_csv(lasso_estimator_model_4_100_pred, 
                         srs_data_store_model_4, "lasso-plugin-var-sim-model-4-100-pred.csv")),
  
  #### Bootstrap Simulations -----------------------
  tar_target(lasso_estimator_bootstrap_model_4_100_pred,
             map_dfr(model_4_data_list[["100"]], lasso_reg_est,
                     n_sizes = c(  250, 500),cores = core_alloc,
                     m=1000, boot_yes = TRUE,
                     randomizer = rand_y_2norm_proc, .id = "rho")),
  tar_target(lasso_estimator_bootstrap_model_4_100_pred_results,
             produce_csv(lasso_estimator_bootstrap_model_4_100_pred,
                         srs_data_store_model_4, "lasso-boot-sim-model-4-100-pred.csv")),
  # Mashregi bootstrap
  tar_target(lasso_estimator_mash_bootstrap_model_4_100_pred,
             map_dfr(model_4_data_list[["100"]], lasso_reg_est,
                     n_sizes = c(  250, 500),cores = core_alloc,
                     m=1000, boot_yes = TRUE,
                     lasso_reg_boot_var_function = lasso_reg_bootstrap.var_est.mash,
                     randomizer = rand_y_2norm_proc, .id = "rho")),
  tar_target(lasso_estimator_mash_bootstrap_model_4_100_pred_results,
             produce_csv(lasso_estimator_mash_bootstrap_model_4_100_pred,
                         srs_data_store_model_4, "lasso-mash-boot-sim-model-4-100-pred.csv")),
  
  
  ### 200 Predictors ----------------
  tar_target(lasso_estimator_model_4_200_pred,
             map_dfr(model_4_data_list[["200"]], lasso_reg_var_est ,
                     n_sizes = c(  250, 500, 1000),cores = core_alloc,
                     m=1000, 
                     randomizer = rand_y_2norm_proc, .id = "rho")),
  tar_target(lasso_estimator_model_4_200_pred_results,
             produce_csv(lasso_estimator_model_4_200_pred, 
                         srs_data_store_model_4, "lasso-plugin-var-sim-model-4-200-pred.csv")),
  
  #### Bootstrap Simulations -----------------------
  tar_target(lasso_estimator_bootstrap_model_4_200_pred,
             map_dfr(model_4_data_list[["200"]], lasso_reg_est,
                     n_sizes = c(  250, 500),cores = core_alloc,
                     m=1000, boot_yes = TRUE,
                     randomizer = rand_y_2norm_proc, .id = "rho")),
  tar_target(lasso_estimator_bootstrap_model_4_200_pred_results,
             produce_csv(lasso_estimator_bootstrap_model_4_200_pred,
                         srs_data_store_model_4, "lasso-boot-sim-model-4-200-pred.csv")),
  
  # Mashregi bootstrap
  tar_target(lasso_estimator_mash_bootstrap_model_4_200_pred,
             map_dfr(model_4_data_list[["200"]], lasso_reg_est,
                     n_sizes = c(  250, 500),cores = core_alloc,
                     m=1000, boot_yes = TRUE,
                     lasso_reg_boot_var_function = lasso_reg_bootstrap.var_est.mash,
                     randomizer = rand_y_2norm_proc, .id = "rho")),
  tar_target(lasso_estimator_mash_bootstrap_model_4_200_pred_results,
             produce_csv(lasso_estimator_mash_bootstrap_model_4_200_pred,
                         srs_data_store_model_4, "lasso-mash-boot-sim-model-4-200-pred.csv")),
  
  
  ### 300 Predictors ----------------
  tar_target(lasso_estimator_model_4_300_pred,
             map_dfr(model_4_data_high_dim_list[["300"]], lasso_reg_var_est ,
                     n_sizes = c(  250, 500, 1000),cores = core_alloc,
                     m=1000, 
                     randomizer = rand_y_2norm_proc, .id = "rho")),
  tar_target(lasso_estimator_model_4_300_pred_results,
             produce_csv(lasso_estimator_model_4_300_pred, 
                         srs_data_store_model_4, "lasso-plugin-var-sim-model-4-300-pred.csv")),
  
  #### Bootstrap Simulations -----------------------
  # Mashregi bootstrap
  tar_target(lasso_estimator_mash_bootstrap_model_4_300_pred,
             map_dfr(model_4_data_high_dim_list[["300"]], lasso_reg_est,
                     n_sizes = c(  250, 500),cores = core_alloc,
                     m=1000, boot_yes = TRUE,
                     lasso_reg_boot_var_function = lasso_reg_bootstrap.var_est.mash,
                     randomizer = rand_y_2norm_proc, .id = "rho")),
  tar_target(lasso_estimator_mash_bootstrap_model_4_300_pred_results,
             produce_csv(lasso_estimator_mash_bootstrap_model_4_300_pred,
                         srs_data_store_model_4, "lasso-mash-boot-sim-model-4-300-pred.csv")),
  
  ### 600 Predictors ----------------
  tar_target(lasso_estimator_model_4_600_pred,
             map_dfr(model_4_data_high_dim_list[["600"]], lasso_reg_var_est ,
                     n_sizes = c(  250, 500, 1000),cores = core_alloc,
                     m=1000, 
                     randomizer = rand_y_2norm_proc, .id = "rho")),
  tar_target(lasso_estimator_model_4_600_pred_results,
             produce_csv(lasso_estimator_model_4_600_pred, 
                         srs_data_store_model_4, "lasso-plugin-var-sim-model-4-600-pred.csv")),
  
  #### Bootstrap Simulations -----------------------
  # Mashregi bootstrap
  tar_target(lasso_estimator_mash_bootstrap_model_4_600_pred,
             map_dfr(model_4_data_high_dim_list[["600"]], lasso_reg_est,
                     n_sizes = c(  250, 500),cores = core_alloc,
                     m=1000, boot_yes = TRUE,
                     lasso_reg_boot_var_function = lasso_reg_bootstrap.var_est.mash,
                     randomizer = rand_y_2norm_proc, .id = "rho")),
  tar_target(lasso_estimator_mash_bootstrap_model_4_600_pred_results,
             produce_csv(lasso_estimator_mash_bootstrap_model_4_600_pred,
                         srs_data_store_model_4, "lasso-mash-boot-sim-model-4-600-pred.csv")),
  
  
  
  
  
  ## ENET -------------------
  ### 25 Predictors --------------
  tar_target(enet_estimator_model_4_25_pred,
             map_dfr(model_4_data_list[["25"]], elastic_net_var_results,
                     n_sizes = c(  250, 500, 1000), cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc,
                     .id = "rho")),
  tar_target(enet_estimator_model_4_25_pred_results,
             produce_csv(enet_estimator_model_4_25_pred, 
                         srs_data_store_model_4, "elastic-net-plugin-var-sim-model-4-25-pred.csv")),
  
  ### 50 Predictors --------------
  tar_target(enet_estimator_model_4_50_pred,
             map_dfr(model_4_data_list[["50"]], elastic_net_var_results,
                     n_sizes = c(  250, 500, 1000), cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc,
                     .id = "rho")),
  tar_target(enet_estimator_model_4_50_pred_results,
             produce_csv(enet_estimator_model_4_50_pred, 
                         srs_data_store_model_4, "elastic-net-plugin-var-sim-model-4-50-pred.csv")),
  
  ### 100 Predictors --------------
  tar_target(enet_estimator_model_4_100_pred,
             map_dfr(model_4_data_list[["100"]], elastic_net_var_results,
                     n_sizes = c(  250, 500, 1000), cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc,
                     .id = "rho")),
  tar_target(enet_estimator_model_4_100_pred_results,
             produce_csv(enet_estimator_model_4_100_pred, 
                         srs_data_store_model_4, "elastic-net-plugin-var-sim-model-4-100-pred.csv")),
  
  ### 200 Predictors --------------
  tar_target(enet_estimator_model_4_200_pred,
             map_dfr(model_4_data_list[["200"]], elastic_net_var_results,
                     n_sizes = c(  250, 500, 1000), cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc,
                     .id = "rho")),
  tar_target(enet_estimator_model_4_200_pred_results,
             produce_csv(enet_estimator_model_4_200_pred, 
                         srs_data_store_model_4, "elastic-net-plugin-var-sim-model-4-200-pred.csv")),
  
  ### 300 Predictors --------------
  tar_target(enet_estimator_model_4_300_pred,
             map_dfr(model_4_data_high_dim_list[["300"]], elastic_net_var_results,
                     n_sizes = c(  250, 500, 1000), cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc,
                     .id = "rho")),
  tar_target(enet_estimator_model_4_300_pred_results,
             produce_csv(enet_estimator_model_4_300_pred, 
                         srs_data_store_model_4, "elastic-net-plugin-var-sim-model-4-300-pred.csv")),
  
  ### 600 Predictors --------------
  tar_target(enet_estimator_model_4_600_pred,
             map_dfr(model_4_data_high_dim_list[["600"]], elastic_net_var_results,
                     n_sizes = c(  250, 500, 1000), cores = core_alloc,
                     m=1000,
                     randomizer = rand_y_2norm_proc,
                     .id = "rho")),
  tar_target(enet_estimator_model_4_600_pred_results,
             produce_csv(enet_estimator_model_4_600_pred, 
                         srs_data_store_model_4, "elastic-net-plugin-var-sim-model-4-600-pred.csv"))
  
  
  #xgbtree_target_obj,
  #ranger_rf_target_obj,
  #lasso_peml_obj
)


# targets::tar_make(linear_reg_estimator_model_1_300_pred)

