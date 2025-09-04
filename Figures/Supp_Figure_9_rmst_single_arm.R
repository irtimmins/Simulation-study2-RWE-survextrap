
#########################################################
# Supp Figure #, Simulation study, rmst plot.
# Generate plots for each rmst time point.
#########################################################

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


setwd("/scratch/klvq491/simsurvextrap_slurm_mix_weib_full1")

################################################
# RMST for single arm.
################################################

scenarios <- readRDS("scenarios.rds")
scenarios_freq <- readRDS("scenarios_single_arm_freq.rds")

performance_res <- readRDS("rmst_and_irmst_performance.rds")
performance_freq_res <- readRDS("rmst_and_irmst_performance_single_arm_freq.rds")

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

new_model <-  scenarios %>%
  mutate(add_knots = if_else(is.na(add_knots), "default", add_knots)) %>%
  mutate(include_external_data = if_else(include_external_data, "yes", "no")) %>%
  mutate(bsmooth = if_else(bsmooth == T, "Smoothed", "Standard")) %>%
  # mutate(add_knots %in% c("default", "extra_knots5")) %>%
  filter(weibull_model_id == "weibull_mod1") %>%
  # filter based on analysis 1 .
  filter(weibull_model_id == "weibull_mod1",
         design_id == "single_arm",
         stan_fit_method == "mcmc",
         mspline_df == 10,
         data_cut_model_id == "data_cut_mod2",
         add_knots %in% c("none", "extra_knots1"),
         backhaz == T,
         waning_model_id == "waning_mod1") %>%
  mutate(new_model_id = paste0(external_bias_model_id, "_", add_knots)) %>%
  mutate(new_model_id_labels = 0) 
View(new_model)

for(i in 1:nrow(new_model)){
  
  if(new_model$external_bias_model_id[i] == "none"){
    #new_model$add_knots
    label1 <- external_data_models_labels$external_data_label[
      external_data_models_labels$external_bias_model_id == new_model$external_bias_model_id[i]]
    label2 <-  extra_knots_models$extra_knots_labels[
      extra_knots_models$extra_knots_id == new_model$add_knots[i]]
    
    new_model$new_model_id_labels[i] <- paste0(label1, label2)
    
  } else{ 
    
    new_model$new_model_id_labels[i] <- external_data_models_labels$external_data_label[
      external_data_models_labels$external_bias_model_id == new_model$external_bias_model_id[i]]
    
  }
}

new_model_levels <- unique(new_model$new_model_id)[c(2,3,5,4,1,6,7)]
new_model_labels <- unique(new_model$new_model_id_labels)[c(2,3,5,4,1,6,7)]
#new_model_levels
#new_model_labels 

# scenarios_test$new_model_id
new_model$new_model_id
new_model$new_model_id_labels
new_model

estimand_labels <- readRDS("estimand_labels.rds")

rmst_estimand_vec <- estimand_labels %>%
  filter(trt == 0) %>%
  filter(estimand == "rmst") %>%
  pull(estimand_id)

rmst_estimand_vec

for(rmst_estimand in rmst_estimand_vec){
  
  #rmst_estimand <- "rmst1_trt0"
  
  
  #j <- 1
  scen_survextrap_all_df <- 
    new_model  %>%
    filter(external_bias_model_id %in% c("none", paste0("external_bias_mod", c(1,3,5)))) %>%
    left_join(performance_res, by = "scenario_fit_id")  %>%
    filter(estimand_id == rmst_estimand) %>%
    mutate(Model = new_model_id) %>%
    select(estimand,  t, trt, stat,est,lower, upper, scenario_fit_id, Model, trt_effect_model_id)
  

  scen_freq_all_df <- scenarios_freq %>%
    #  filter(trt_effect_model_id == paste0("trt_effect_mod", j)) %>%
    left_join(performance_freq_res %>% 
                mutate(single_arm_freq_model_id = scenario_fit_id),
              by = "single_arm_freq_model_id")  %>%
    filter(estimand_id == rmst_estimand) %>%
    mutate(Model = parametric_model) %>%
    mutate(scenario_fit_id = single_arm_freq_model_id) %>%
    select(estimand, t, trt, stat, est, lower, upper, scenario_fit_id, Model, trt_effect_model_id)
  

      new_levels <- c("exp", 
                      "weibull", 
                      "spline",
                      "none_none", 
                      "none_extra_knots1" ,
                      "external_bias_mod3_extra_knots1",
                      "external_bias_mod1_extra_knots1",
                      "external_bias_mod5_extra_knots1")
      
      new_labels <- c( "Exponential model",
                       "Weibull model",
                       "Royston-Parmar spline model,\ndf = 3",
                       "survextrap model, no external data,\nno extra knots",
                       "survextrap model, no external data,\nwith extra knots at t=5,10,25",
                       "survextrap model,\nwith external data with +20% bias",
                       "survextrap model,\nwith unbiased external data",
                       "survextrap model,\nwith external data with −20% bias")

      scen_df <- bind_rows(scen_freq_all_df, 
                               scen_survextrap_all_df) %>%
        mutate(Model = factor(Model, levels = new_levels, labels = new_labels)) %>%
        arrange(Model) %>%
        group_by(stat) %>%
        mutate(y_height = -row_number()) %>%
        ungroup() 
    
    # summary(scen_df$Model)

    y_min <- scen_df %>%
      select(y_height) %>%
      pull() %>%
      min() - 1
    
    y_max <- scen_df %>%
      select(y_height) %>%
      pull() %>%
      max() + 1
    
    true_value <- scen_df %>%
      filter(stat == "true") %>%
      select(est) %>%
      distinct() %>%
      pull()
    
    x_min <- scen_df %>%
      filter(stat == "mean") %>% 
      select(lower) %>%
      filter(!is.na(lower)) %>%
      pull() %>%
      min() %>%
      min(true_value) - 0.3
    
    x_max <- scen_df %>%
      filter(stat == "mean") %>% 
      select(lower) %>%
      filter(!is.na(lower)) %>%
      pull() %>%
      max() %>%
      max(true_value) + 0.3

    x_axis_label <- scen_df %>%
      select(estimand, t) %>%
      distinct() %>%
      mutate(label = paste0("RMST", " at ", t, "-y")) %>%
      pull(label)
    
    
    colour_values <- c("gray40", "gray40", "gray40", hue_pal()(6)[1], hue_pal()(6)[c(1,2,4,6)])
    fill_values <- c( "white","white","white","white", hue_pal()(6)[1],  hue_pal()(6)[c(2,4,6)])
    shape_values <- c( 3, 4, 8, rep(21, 5))
    size_values <- c(1.2, 1.2, 1.2, rep(2, 5))
    
   forest_plot <- 
      scen_df %>%
      filter(stat == "mean") %>%
      mutate(true_value = true_value) %>%
      ggplot(aes(x = est, y = y_height)) + 
      theme_classic()+
      theme(panel.grid = element_blank(),
            panel.border = element_blank(),
            axis.line.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            axis.title.y=element_blank(),
            axis.text.x = element_text( size = 6), #, angle = 45, vjust = 0.5, hjust=0.5),
            axis.title.x = element_text( face="bold",size = 8),
            legend.background = element_blank(),
            legend.box.background = element_rect(colour = "black"),
            legend.key.spacing.y = unit(3, "pt"),
            legend.text = element_text(size=8)) + 
      geom_segment(x = true_value, 
                   xend = true_value, 
                   y = y_min-0.5, 
                   yend = y_max, 
                   colour = "gray50",
                   alpha = 0.3)+
      geom_point(aes(group = scenario_fit_id, 
                     colour = Model,
                     fill = Model,
                     shape = Model,
                     size = Model),
                 na.rm = FALSE,
                 alpha = 1,
                 stroke = 1) + 
      geom_segment(aes(x= true_value, xend = est, y = y_height, yend = y_height,
                       colour = Model),  alpha = 0.35)+
      # geom_point(aes(x=lower, y = y_height, colour = Scenarios), 
      #            shape = 91, size = 3, alpha = 0.5)+
      # geom_point(aes(x=upper, y = y_height, colour = Scenarios), 
      #            shape = 93, size = 3, alpha = 0.5)+
      scale_colour_manual("Models", values = colour_values)+
      scale_fill_manual("Models", values = fill_values)+
      scale_shape_manual("Models", values = shape_values)+
      scale_size_manual("Models", values = size_values)+
      scale_x_continuous(x_axis_label, limits = c(x_min, x_max))
  
    forest_plot
    saveRDS(forest_plot, paste0("plots/single_arm/forest_", rmst_estimand, "_cond2.rds"))

  tiff(file = paste0("plots/single_arm/forest_", rmst_estimand, "_cond2.tiff"),
       width = 4.6, 
       height = 3.2,
       units = 'in',  
       res = 300, 
       compression = "lzw")
  print( forest_plot)
  dev.off()
}

