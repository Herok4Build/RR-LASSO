#' Thomas Johnson III
#' 09-26-2024
#' Relative bias table
#' 
#' 

library(tidyverse)
library(tables)

data_locations = c(paste0(here::here(), 
                          "/Data/SRS-Model-1-Dim/Rel-Bias/relative-bias-model-1.csv"),
                   paste0(here::here(), 
                          "/Data/SRS-Model-2-Dim/Rel-Bias/relative-bias-model-2.csv"),
                   paste0(here::here(), 
                          "/Data/SRS-Model-3-Dim/Rel-Bias/relative-bias-model-3.csv"),
                   paste0(here::here(), 
                          "/Data/SRS-Model-4-Dim/Rel-Bias/relative-bias-model-4.csv"))


model_1_bias = read_csv(data_locations[1],lazy = TRUE) %>%
  mutate(model = "Model 1", rho = as.character(rho))

model_2_bias = read_csv(data_locations[2],lazy = TRUE)%>%
  mutate(model = "Model 2", rho = as.character(rho))

model_3_bias = read_csv(data_locations[3],lazy = TRUE)%>%
  mutate(model = "Model 3", rho = as.character(rho))

model_4_bias = read_csv(data_locations[4],lazy = TRUE)%>%
  mutate(model = "Model 4", rho = as.character(rho))

srs_models_bias = bind_rows(list(model_1_bias,
                                 model_2_bias,
                                 model_3_bias,
                                 model_4_bias)) %>%
  mutate(rho = case_when( rho == "0.1" ~"rho_0.1",
                          rho == "0.3" ~"rho_0.3",
                          rho == "0.5" ~"rho_0.5",
                          rho == "0.7" ~"rho_0.7",
                          .default = as.character(rho)),
         model = case_when(model == "Model 2" ~ "Model 4",
                           model == "Model 3" ~ "Model 2",
                           model == "Model 4" ~ "Model 3",
                           .default = model))


construct_rel_bias_table = 
  function(relative_bias_results)
  {
    relative_bias_results = relative_bias_results %>% select(-mse) %>%
      filter(n_size != 100, pred != 75 ) %>%
      mutate(rel_bias = round(rel_bias,3)) %>%
      pivot_wider( values_from = "rel_bias", names_from = estimator)
    rel_bias_table = tabular(Factor(pred) * Factor(rho) * Factor(n_size)  ~ 
                               Factor(model) * (HT + GREG + LASSO + RF + ENET) * (identity),
                             data = relative_bias_results %>% filter(n_size %in% c(250, 500)))
    return(latex.tabular(rel_bias_table))
  }

construct_rel_bias_table(srs_models_bias)


construct_rel_bias_table(srs_models_bias %>% filter(model %in% c("Model 1", "Model 2")))


construct_rel_bias_table(srs_models_bias %>% filter(model %in% c("Model 3", "Model 4")))
