
library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(survextrap)
library(flexsurv)
library(rslurm)
library(rsimsum)
library(stringr)
library(readr)
library(rstan)
library(tibble)
library(rstpm2)
library(pracma)
library(posterior)
library(stats)
library(data.table)
library(scales)
library(ggpubr)
library(grid)
library(gridExtra) 

jobname <- "mix_weib_full1"
user <- Sys.info()["user"]
project_directory <- paste0("/projects/aa/", user, "/")
store_res <- paste0(project_directory, "simsurvextrap_slurm_", jobname, "/")

setwd(store_res)


performance_res <- readRDS("rmst_and_irmst_performance.rds")


scenarios <- readRDS("scenarios.rds")

scenarios <-  scenarios %>%
   filter(weibull_model_id == "weibull_mod1") %>%
  # filter based on analysis 1 .
  filter(design_id == "single_arm",
         stan_fit_method == "mcmc",
         mspline_df == 10,
         data_cut_model_id == "data_cut_mod2",
         add_knots %in% c("none", "extra_knots1"),
         backhaz == T) %>%
  mutate(new_model_id = paste0(external_bias_model_id, "_", add_knots)) %>%
  mutate(new_model_id_labels = 0) %>%
  slice(c(2,3,5,4,1,6,7))

################################################
# External data labels.
################################################

external_data_models <- readRDS("external_data_models.rds")
external_data_models_labels <- readRDS("external_data_models_labels.rds")

external_data_models_index <- external_data_models_labels$external_bias_model_id  
external_data_labels_index <- external_data_models_labels$external_data_label

external_data_models_labels <- 
  external_data_models %>%
  arrange(-loghaz_bias) %>%
  mutate(haz_bias = -100+100*exp(loghaz_bias),
         haz_bias_temp = round(haz_bias)) %>%
  mutate(haz_bias = as.character(haz_bias)) %>%
  mutate(haz_bias = if_else(haz_bias_temp > 0, paste0("External data with +",haz_bias_temp, "% bias"), haz_bias))  %>%
  mutate(haz_bias = if_else(haz_bias_temp < 0, paste0("External data with \u2212",abs(haz_bias_temp), "% bias"), haz_bias))  %>%
  mutate(haz_bias = if_else(haz_bias == "0", "Unbiased external data", haz_bias)) %>%
  select(-haz_bias_temp)  %>%
  rename(external_data_label = haz_bias) %>%
  # arrange(abs(loghaz_bias)) %>%
  select(external_bias_model_id, external_data_label) %>%
  add_row(external_bias_model_id = "none", external_data_label = "No external data",
          .before = 1)# %>%



################################################
# Extra knots labels.
################################################

extra_knots_settings <- readRDS("extra_knots_settings.rds")
extra_knots_models <- readRDS("extra_knots_models.rds")
extra_knots_models_index <- extra_knots_models$extra_knots_id
extra_knots_labels_index <- extra_knots_models$extra_knots_labels

estimand_labels <- readRDS("estimand_labels.rds")
extra_knots_settings <- readRDS("extra_knots_settings.rds")

extra_knots_models <- tibble(extra_knots_id = names(extra_knots_settings)) %>%
  mutate(extra_knots_labels = 0)

extra_knots_models$extra_knots_labels <- c(",\nno extra knots",
                                           ",\nwith extra knots at t=5,10,25",
                                           ",\nwith extra knots at t=5,10,25",
                                           0)



################################################
# New labels.
################################################

scenarios <- readRDS("scenarios.rds")

new_model <-  scenarios %>%
  mutate(add_knots = if_else(is.na(add_knots), "default", add_knots)) %>%
  mutate(include_external_data = if_else(include_external_data, "yes", "no")) %>%
  mutate(bsmooth = if_else(bsmooth == T, "Smoothed", "Standard")) %>%
  # mutate(add_knots %in% c("default", "extra_knots5")) %>%
  filter(weibull_model_id == "weibull_mod1") %>%
  # filter based on analysis 1 .
  filter(design_id == "single_arm",
         stan_fit_method == "mcmc",
         mspline_df == 10,
         data_cut_model_id == "data_cut_mod2",
         add_knots %in% c("none", "extra_knots1"),
         backhaz == T) %>%
  mutate(new_model_id = paste0(external_bias_model_id, "_", add_knots)) %>%
  mutate(new_model_id_labels = 0) %>%
  slice(c(2,3,5,4,1,6,7))


for(i in 1:nrow(new_model)){
  # i <- 2
  if(new_model$external_bias_model_id[i] == "none"){
    #new_model$add_knots
    label1 <- external_data_models_labels$external_data_label[
      external_data_models_labels$external_bias_model_id == new_model$external_bias_model_id[i]]
    label2 <-  extra_knots_models$extra_knots_labels[
      extra_knots_models$extra_knots_id == new_model$add_knots[i]]
    
    # label2 <- "*`,`"
    # label3 <- "~`with`~"
    # label4 <-  extra_knots_models$extra_knots_labels[
    #   extra_knots_models$extra_knots_id == new_model$add_knots[i]]
    
    new_model$new_model_id_labels[i] <- paste0(label1, label2)
    
  } else{ 
    
    new_model$new_model_id_labels[i] <- external_data_models_labels$external_data_label[
      external_data_models_labels$external_bias_model_id == new_model$external_bias_model_id[i]]
    
  }
}

# scenarios_test$new_model_id
new_model$new_model_id
new_model$new_model_id_labels
#new_model

labels <- new_model %>%
  select(scenario_fit_id, new_model_id, new_model_id_labels)

table_single_arm <- 
  new_model  %>%
  left_join(performance_res, by = "scenario_fit_id")  %>%
  filter(estimand_id == "rmst2_trt0") %>%
  mutate("Scenarios" = factor(new_model_id, levels = new_model_id,
                              labels = new_model_id_labels)) %>%
  mutate("External data" = factor(external_bias_model_id,
                                  levels = external_data_models_index,
                                  labels = external_data_labels_index)) %>%
  mutate("Extra knots" = factor(add_knots,
                                levels = extra_knots_models_index,
                                labels = extra_knots_labels_index)) %>%
  mutate("GPM rates" = if_else(backhaz, "Included", "Not included")) %>%
  mutate("All models" = paste0(external_bias_model_id, "_", add_knots)) %>%
  #arrange(`Extra knots`) %>%
  arrange(`External data`) %>%
  arrange(backhaz) %>%
  arrange(include_external_data) %>%
  filter(stat %in% c("bias", "cover", "empse",  "mse", "modelse", "mean", "true")) %>%
  mutate(est_round = sprintf("%.2f", round(est, 2)),
           mcse_round = sprintf("%.3f", round(mcse, 3)),
           est_round = as.character(est_round),
           mcse_round = as.character(mcse_round)) %>%                                                    
  mutate(est_str = paste0(est_round, " (", mcse_round, ")")) %>%
  mutate(est_str = 
           case_when(stat == "true" ~ est_round,
                     .default = est_str)
         ) %>%
  select(c(scenario_fit_id, stat, est_str)) %>%
  pivot_wider(names_from = "stat",
             values_from = "est_str" ) %>%
  left_join(labels, join_by("scenario_fit_id")) %>%
  select(new_model_id_labels, true, bias, mse, modelse, cover) #%>%
  #View()

#summary(as.factor(test$stat))
#summ

write_csv(table_single_arm, "plots/single_arm_performance_table.csv")

#########################################################################
# Two arm scenarios.
#########################################################################

scenarios <- readRDS("scenarios.rds")

################################################
# External data labels.
################################################

external_data_models <- readRDS("external_data_models.rds")
external_data_models_labels <- readRDS("external_data_models_labels.rds")

external_data_models_labels <- 
  external_data_models %>%
  arrange(-loghaz_bias) %>%
  mutate(haz_bias = -100+100*exp(loghaz_bias),
         haz_bias_temp = round(haz_bias)) %>%
  mutate(haz_bias = as.character(haz_bias)) %>%
  mutate(haz_bias = if_else(haz_bias_temp > 0, paste0("External data with +",haz_bias_temp, "% bias"), haz_bias))  %>%
  mutate(haz_bias = if_else(haz_bias_temp < 0, paste0("External data with \u2212",abs(haz_bias_temp), "% bias"), haz_bias))  %>%
  mutate(haz_bias = if_else(haz_bias == "0", "Unbiased external data", haz_bias)) %>%
  select(-haz_bias_temp)  %>%
  rename(external_data_label = haz_bias) %>%
  # arrange(abs(loghaz_bias)) %>%
  select(external_bias_model_id, external_data_label) %>%
  add_row(external_bias_model_id = "none", external_data_label = "No external data",
          .before = 1)# %>%

################################################
# Extra knots labels.
################################################

extra_knots_settings <- readRDS("extra_knots_settings.rds")
extra_knots_models <- readRDS("extra_knots_models.rds")
extra_knots_models_index <- extra_knots_models$extra_knots_id
extra_knots_labels_index <- extra_knots_models$extra_knots_labels

estimand_labels <- readRDS("estimand_labels.rds")
extra_knots_settings <- readRDS("extra_knots_settings.rds")

extra_knots_models <- tibble(extra_knots_id = names(extra_knots_settings)) %>%
  mutate(extra_knots_labels = 0)

extra_knots_models$extra_knots_labels <- c(",\nno extra knots",
                                           ",\nwith extra knots at t=5,10,25",
                                           ",\nwith extra knots at t=5,10,25",
                                           0)

################################################
# New labels.
################################################

new_model <-  scenarios %>%
  mutate(add_knots = if_else(is.na(add_knots), "default", add_knots)) %>%
  mutate(include_external_data = if_else(include_external_data, "yes", "no")) %>%
  rename(Model = model) %>%
  mutate(Model = case_when(Model == "ph" ~ 1,
                           Model == "nonph" ~ 2,
                           Model == "separate" ~ 3)) %>%
  mutate(Model = factor(Model, levels = c(1,2,3), 
                        labels = c("PH", 
                                   "Non-PH", 
                                   "Separate arms" ))) %>%
  filter(weibull_model_id == "weibull_mod1",
         design_id == "two_arm",
         stan_fit_method == "mcmc",
         mspline_df == 10,
         data_cut_model_id == "data_cut_mod2",
         add_knots %in% c("none", "extra_knots1"),
         backhaz == T,
         waning_model_id == "waning_mod1") %>%
  mutate(new_model_id = paste0(external_bias_model_id, "_", add_knots)) %>%
  mutate(new_model_id_labels = 0) #%>%
#slice(c(2,3,5,4,1,6,7))

summary(as.factor(new_model$new_model_id))
#View(new_model)

for(i in 1:nrow(new_model)){
  # i <- 2
  if(new_model$external_bias_model_id[i] == "none"){
    #new_model$add_knots
    label1 <- external_data_models_labels$external_data_label[
      external_data_models_labels$external_bias_model_id == new_model$external_bias_model_id[i]]
    label2 <-  extra_knots_models$extra_knots_labels[
      extra_knots_models$extra_knots_id == new_model$add_knots[i]]
    
    # label2 <- "*`,`"
    # label3 <- "~`with`~"
    # label4 <-  extra_knots_models$extra_knots_labels[
    #   extra_knots_models$extra_knots_id == new_model$add_knots[i]]
    
    new_model$new_model_id_labels[i] <- paste0(label1, label2)
    
  } else{ 
    
    new_model$new_model_id_labels[i] <- external_data_models_labels$external_data_label[
      external_data_models_labels$external_bias_model_id == new_model$external_bias_model_id[i]]
    
  }
}
new_model_levels <- unique(new_model$new_model_id)[c(2,3,5,4,1,6,7)]
new_model_labels <- unique(new_model$new_model_id_labels)[c(2,3,5,4,1,6,7)]
new_model_levels
new_model_labels    

# scenarios_test$new_model_id
new_model$new_model_id
new_model$new_model_id_labels
new_model  %>%
  #  slice(c(2,3,5,4,1,6,7)) %>%
  pull(new_model_id_labels)


####################################################
# Difference in RMST across arms.
####################################################

estimand_labels <- readRDS("estimand_labels.rds")

irmst_estimand_vec <- estimand_labels %>%
  filter(estimand == "irmst") %>%
  pull(estimand_id)
names(new_model)

labels <- new_model %>%
  mutate("Scenarios" = factor(new_model_id, levels = new_model_id,
                              labels = new_model_id_labels)) %>%
  mutate("External data" = factor(external_bias_model_id,
                                  levels = external_data_models_index,
                                  labels = external_data_labels_index)) %>%
  mutate("Extra knots" = factor(add_knots,
                                levels = extra_knots_models_index,
                                labels = extra_knots_labels_index)) %>%
  mutate("GPM rates" = if_else(backhaz, "Included", "Not included")) %>%
  mutate("All models" = paste0(external_bias_model_id, "_", add_knots)) %>%
  select(scenario_fit_id, trt_effect_model_id, Scenarios, new_model_id_labels, Model)

summary(as.factor(labels$Scenarios))

table_two_arm <- 
  new_model  %>%
  left_join(performance_res, by = "scenario_fit_id")  %>%
  filter(estimand_id == "irmst2") %>%
  filter(stat %in% c("bias", "cover", "empse",  "mse", "modelse", "mean", "true")) %>%
  mutate(est_round = sprintf("%.2f", round(est, 2)),
         mcse_round = sprintf("%.3f", round(mcse, 3)),
         est_round = as.character(est_round),
         mcse_round = as.character(mcse_round)) %>%                                                    
  mutate(est_str = paste0(est_round, " (", mcse_round, ")")) %>%
  mutate(est_str = 
           case_when(stat == "true" ~ est_round,
                     .default = est_str)
  ) %>%
  select(c(scenario_fit_id, stat, est_str)) %>%
  pivot_wider(names_from = "stat",
              values_from = "est_str" ) %>%
  left_join(labels, join_by("scenario_fit_id")) %>%
  select(trt_effect_model_id,Scenarios, Model, true, bias, mse, modelse, cover) %>%
  filter(Scenarios %in% c("Unbiased external data",
                          "No external data,\nno extra knots",
                          "No external data,\nwith extra knots at t=5,10,25")) %>%
  arrange(Model) %>%
  arrange(Scenarios) %>%
  arrange(trt_effect_model_id)

write_csv(table_two_arm, "plots/two_arm_performance_table.csv")
