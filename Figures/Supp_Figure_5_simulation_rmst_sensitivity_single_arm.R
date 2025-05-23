
###########################################
# Figure RMST Sensitivity.
###########################################
###########################################

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
library(forcats)

jobname <- "mix_weib_full1"
user <- Sys.info()["user"]
project_directory <- paste0("/projects/aa/", user, "/")
store_res <- paste0(project_directory, "simsurvextrap_slurm_", jobname, "/")

source("R/estimands_aim2.R")

setwd(store_res)

############################################################################
# Sensitivity on MCMC vs Laplace "Opt" approximation
############################################################################

scenarios <- readRDS("scenarios.rds")
performance_res <- readRDS("rmst_and_irmst_performance.rds")

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

scenarios <- readRDS("scenarios.rds")

new_model <-  
  scenarios %>%
  mutate(add_knots = if_else(is.na(add_knots), "default", add_knots)) %>%
  mutate(include_external_data = if_else(include_external_data, "yes", "no")) %>%
  mutate(bsmooth = if_else(bsmooth == T, "Smoothed", "Standard")) %>%
  # mutate(add_knots %in% c("default", "extra_knots5")) %>%
  filter(weibull_model_id == "weibull_mod1") %>%
  # filter based on analysis 1 .
  filter(design_id == "single_arm",
         mspline_df == 10,
         data_cut_model_id == "data_cut_mod2",
         add_knots %in% c("none", "extra_knots1"),
         backhaz == T) %>%
  mutate(new_model_id = paste0(external_bias_model_id, "_", add_knots)) %>%
  mutate(new_model_id_labels = 0)


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
reorder_vec <- c(2,3,5, 4,1,6,7)
new_model_levels <- unique(new_model$new_model_id)[reorder_vec]
new_model_labels <- unique(new_model$new_model_id_labels)[reorder_vec]
  
estimand_labels <- readRDS("estimand_labels.rds")

rmst_estimand_vec <- estimand_labels %>%
  filter(estimand == "rmst", trt == 0) %>%
  pull(estimand_id)


for(rmst_estimand in rmst_estimand_vec){
  
  # #rmst_estimand <- "rmst5_trt0"
  # 
  # for(i in 1:2){

    #extra_knots_models
    scen_df <-
      new_model  %>%
      mutate("Scenarios" = factor(new_model_id, 
                                  levels = new_model_levels,
                                  labels = new_model_labels)) %>%
      arrange(stan_fit_method) %>%
      arrange(Scenarios) %>%
      left_join(performance_res, by = "scenario_fit_id")  %>%
      filter(estimand_id == rmst_estimand) %>%
      mutate("External data" = factor(external_bias_model_id,
                                      levels = external_data_models_index,
                                      labels = external_data_labels_index)) %>%
      mutate("Extra knots" = factor(add_knots,
                                    levels = extra_knots_models_index,
                                    labels = extra_knots_labels_index)) %>%
      mutate("GPM rates" = if_else(backhaz, "Included", "Not included")) %>%
      mutate("All models" = paste0(external_bias_model_id, "_", add_knots)) %>%
      mutate("Stan fit method" = if_else(stan_fit_method == "mcmc", "MCMC", "Opt" )) %>%
      #arrange(`Extra knots`) %>%
      arrange(`External data`) %>%
      arrange(backhaz) %>%
      arrange(include_external_data) %>%
      group_by(stat) %>%
      mutate(y_height = -row_number()) %>%
      ungroup() 
    

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
      select(upper) %>%
      filter(!is.na(upper)) %>%
      pull() %>%
      max() %>%
      max(true_value) + 0.3
    
    y_min <- scen_df %>%
      select(y_height) %>%
      pull() %>%
      min() - 1
    
    y_max <- scen_df %>%
      select(y_height) %>%
      pull() %>%
      max() + 0.5
    
    x_axis_label <- scen_df %>%
      select(estimand, t) %>%
      distinct() %>%
      mutate(label = paste0(toupper(estimand), " at ", t, "-y")) %>%
      pull(label)
    
    #View(scen_df)
   
    colour_values <- c(hue_pal()(6)[1], hue_pal()(6))
    fill_values <-c("white",hue_pal()(6)[1],  hue_pal()(6)[2:6])
    
     
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
      geom_point(aes(colour = Scenarios,
                     fill = Scenarios,
                     shape = `Stan fit method`),
                 na.rm = FALSE,
                 alpha = 1,
                 stroke = 1)+ 
      geom_segment(aes(x= true_value, xend = est, y = y_height,
                       yend = y_height, colour = Scenarios),  alpha = 0.35)+
      # geom_point(aes(x=lower, y = y_height, colour = Scenarios), 
      #            shape = 91, size = 3, alpha = 0.5)+
      # geom_point(aes(x=upper, y = y_height, colour = Scenarios), 
      #            shape = 93, size = 3, alpha = 0.5)+
      scale_colour_manual("Settings", values = colour_values)+
      scale_fill_manual("Settings", values = fill_values)+            
      scale_shape_manual(values = c(21,24))+   
      scale_x_continuous(x_axis_label, limits = c(x_min, x_max))+
        guides(                              
        colour = guide_legend(override.aes=list(shape = 21,
                                               stroke = 1,
                                               size = 1.7,
                                               alpha = 1,
                                               colour = colour_values,
                                               fill = fill_values),
        order = 1),
        fill = guide_legend(override.aes=list(shape = 21,
                                                stroke = 1,
                                                size = 1.7,
                                                alpha = 1,
                                                colour = colour_values,
                                                fill = fill_values),
                              order = 1),
       shape = guide_legend(override.aes=list(colour = "gray40",
                                                   order = 2)))

    forest_plot  
    saveRDS(forest_plot, paste0("plots/single_arm/forest_",rmst_estimand,"_single_arm_plot_stan.rds"))
    
  tiff(file = paste0("plots/single_arm/forest_",rmst_estimand,"_single_arm_stan.tiff"),   
       width = 6.5, 
       height = 4.9,
       units = 'in',  
       res = 300, 
       compression = "lzw")
  print(forest_plot)
  dev.off()
  
  print(paste0("plots/single_arm/forest_",rmst_estimand,"_single_arm_stan.tiff"))
  
}


############################################################################
# Sensitivity on knot locations.
############################################################################

scenarios <- readRDS("scenarios.rds")
performance_res <- readRDS("rmst_and_irmst_performance.rds")

external_data_models <- readRDS("external_data_models.rds")
external_data_models_labels <- 
  external_data_models %>%
  arrange(-loghaz_bias) %>%
  mutate(haz_bias = -100+100*exp(loghaz_bias),
         haz_bias_temp = round(haz_bias)) %>%
  mutate(haz_bias = as.character(haz_bias)) %>%
  mutate(haz_bias = if_else(haz_bias_temp > 0, paste0("`+`*`",haz_bias_temp, "`*`%`~`bias`"), haz_bias))  %>%
  mutate(haz_bias = if_else(haz_bias_temp < 0, paste0("`\u2212`*`",abs(haz_bias_temp), "`*`%`~`bias`"), haz_bias))  %>%
  mutate(haz_bias = if_else(haz_bias == "0", "`Unbiased`~`external`~`data`~`in`~`model`", haz_bias)) %>%
  select(-haz_bias_temp)  %>%
  rename(external_data_label = haz_bias) %>%
  # arrange(abs(loghaz_bias)) %>%
  select(external_bias_model_id, external_data_label) %>%
  add_row(external_bias_model_id = "none", 
          external_data_label = "`No`~`external`~`data`~`in`~`model`",
          .before = 1)# %>%


external_data_models_index <- external_data_models_labels$external_bias_model_id  
external_data_labels_index <- external_data_models_labels$external_data_label

extra_knots_settings <- readRDS("extra_knots_settings.rds")
extra_knots_models <- readRDS("extra_knots_models.rds")
extra_knots_models <- extra_knots_models %>%
  slice(c(1,4,3,2))
extra_knots_models_index <- extra_knots_models$extra_knots_id
extra_knots_labels_index <- extra_knots_models$extra_knots_labels

estimand_labels <- readRDS("estimand_labels.rds")

rmst_estimand_vec <- estimand_labels %>%
  filter(estimand == "rmst", trt == 0) %>%
  pull(estimand_id)

for(rmst_estimand in rmst_estimand_vec){
  
  #rmst_estimand <- rmst_estimand_vec[1]
  # for(i in 1:2){
  #   
  # #i <- 1 
  #   #extra_knots_models
    
    scen_df <- scenarios %>%
      mutate(add_knots = if_else(is.na(add_knots), "default", add_knots)) %>%
      mutate(include_external_data = if_else(include_external_data, "yes", "no")) %>%
      mutate(bsmooth = if_else(bsmooth == T, "Smoothed", "Standard")) %>%
      # mutate(add_knots %in% c("default", "extra_knots5")) %>%
      filter(weibull_model_id == paste0("weibull_mod",1)) %>%
      filter(design_id == "single_arm",
             stan_fit_method == "mcmc",
             data_cut_model_id == "data_cut_mod2",
             external_bias_model_id %in% c("none", "external_bias_mod1"),
             backhaz == T) %>%
      left_join(performance_res, by = "scenario_fit_id")  %>%
      filter(estimand_id == rmst_estimand) %>%
      mutate("External data" = factor(external_bias_model_id,
                                      levels = external_data_models_index,
                                      labels = external_data_labels_index)) %>%
      mutate("Extra knots for external data" = factor(add_knots,
                                    levels = extra_knots_models_index,
                                    labels = extra_knots_labels_index)) %>%
      #mutate(`Extra_knots" = as.factor(`Extra knots`)) %>%
      mutate("GPM rates" = if_else(backhaz, "Included", "Not included")) %>%
      mutate("Stan fit method" = if_else(stan_fit_method == "mcmc", "MCMC", "Opt" )) %>%
      arrange(mspline_df) %>%
      mutate(`Trial data M-spline` = paste0("df = ", mspline_df)) %>%
      mutate(`Trial data M-spline` = fct_inorder(`Trial data M-spline`)) %>%
     # mutate(`External data` = as.character(`External data`)) %>%
      #mutate(`External data` = as.factor(`External data`)) %>%
      arrange(`Extra knots for external data`) %>%
      arrange(`Trial data M-spline`) %>%
      arrange(`External data`) %>%
      arrange(backhaz) %>%
      arrange(include_external_data) %>%
      group_by(stat) %>%
      mutate(y_height = -row_number()) %>%
      ungroup() 
    
    #View(scen_df)
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
      select(upper) %>%
      filter(!is.na(upper)) %>%
      pull() %>%
      max() %>%
      max(true_value) + 0.3
    
    y_min <- scen_df %>%
      select(y_height) %>%
      pull() %>%
      min() - 1
    
    y_max <- scen_df %>%
      select(y_height) %>%
      pull() %>%
      max() + 0.5
    
    x_axis_label <- scen_df %>%
      select(estimand, t) %>%
      distinct() %>%
      mutate(label = paste0(toupper(estimand), " at ", t, "-y")) %>%
      pull(label)
    
    #View(scen_df)
    
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
           legend.key.spacing.y = unit(1.5, "pt"),
           legend.text = element_text(size=8)) + 
      geom_segment(x = true_value, 
                   xend = true_value, 
                   y = y_min-0.5, 
                   yend = y_max, 
                   colour = "gray50",
                   alpha = 0.3)+
      geom_point(aes(#group = scenario_fit_id, 
        colour = `Trial data M-spline`,
        shape = `Extra knots for external data`,
        #  shape = `Extra knots`),
        fill = `External data`),
      na.rm = FALSE,
      alpha = 1,
      stroke = 1,
      size = 2) + 
      #  geom_segment(aes(x= true_value, xend = est, y = y_height, yend = y_height,
      #                   colour = `External data`),  alpha = 0.35)+
      #  geom_point(aes(x=lower, y = y_height, colour = `External data`), shape = 91, size = 3)+
      #  geom_point(aes(x=upper, y = y_height, colour = `External data`), shape = 93, size = 3)+
      geom_segment(aes(x= true_value, xend = est, y = y_height, yend = y_height,
                       colour = `Trial data M-spline`),  alpha = 0.35)+
      # geom_point(aes(x=lower, y = y_height, colour = `Trial data M-spline`),
      #            shape = 91, size = 3)+
      # geom_point(aes(x=upper, y = y_height, colour = `Trial data M-spline`),
      #            shape = 93, size = 3)+
      scale_colour_discrete()+
      scale_x_continuous(x_axis_label, limits = c(x_min, x_max))+
      scale_shape_manual(values = c(21,22,23,24))+
      # scale_shape_manual(name= "Basis", values = c(21,22))+
      scale_fill_manual(values = c("white","gray60"), labels = parse_format())+
      guides(fill = guide_legend(override.aes=list(shape=c(21,21),
                                              fill = c("white","gray70"),
                                               colour = "gray40"),
                            order = 1),
         shape = guide_legend(override.aes=list(colour = "gray40"),
                              order = 3),#)#,
         colour = guide_legend(override.aes=list(shape = 19,
                                                 size = 1.7),
                               order = 2))
    
    
    forest_plot  
    saveRDS(forest_plot, paste0("plots/single_arm/forest_",rmst_estimand,"_single_arm_plot_knots.rds"))
    
  tiff(file = paste0("plots/single_arm/forest_",rmst_estimand,"_single_arm_knots.tiff"),   
       width = 6.8, 
       height = 4.4,
       units = 'in',  
       res = 300, 
       compression = "lzw")
  print(forest_plot)
  dev.off()
  
  print(paste0("plots/single_arm/forest_",rmst_estimand,"_single_arm_knots.tiff"))
  
}



############################################################################
# Sensitivity on data cut.
############################################################################

scenarios <- readRDS("scenarios.rds")
performance_res <- readRDS("rmst_and_irmst_performance.rds")

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

new_model <-  
  scenarios %>%
  mutate(add_knots = if_else(is.na(add_knots), "default", add_knots)) %>%
  mutate(include_external_data = if_else(include_external_data, "yes", "no")) %>%
  mutate(bsmooth = if_else(bsmooth == T, "Smoothed", "Standard")) %>%
  # mutate(add_knots %in% c("default", "extra_knots5")) %>%
  filter(weibull_model_id == "weibull_mod1") %>%
  # filter based on analysis 1 .
  filter(design_id == "single_arm",
         stan_fit_method == "mcmc",
         mspline_df == 10,
         add_knots %in% c("none", "extra_knots1"),
         external_bias_model_id %in% c("none","external_bias_mod1"),
         backhaz == T) %>%
  mutate(new_model_id = paste0(external_bias_model_id, "_", add_knots)) %>%
  mutate(new_model_id_labels = 0)


for(i in 1:nrow(new_model)){
  # i <- 2
 # if(new_model$external_bias_model_id[i] == "none"){
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
    
  # } else{ 
  #   
  #   new_model$new_model_id_labels[i] <- external_data_models_labels$external_data_label[
  #     external_data_models_labels$external_bias_model_id == new_model$external_bias_model_id[i]]
  #   
  # }
}

# scenarios_test$new_model_id
new_model$new_model_id
new_model$new_model_id_labels
#new_model
reorder_vec <- c(2,3,1)
new_model_levels <- unique(new_model$new_model_id)[reorder_vec]
new_model_labels <- unique(new_model$new_model_id_labels)[reorder_vec]

# Find survival probabilities at datacut.

dgm_true <- readRDS("dgm_true.rds")

pars_dgm_surv_true_data_cut <- dgm_true %>%
  select(true_model_id, design_id) %>%
  filter(design_id == "single_arm") %>%
  mutate(big_data_file = paste0("../dgm_", true_model_id, "/big_data.rds")) %>%
  expand_grid(t = c(3,5,8)) %>%
  mutate(trt_group = 0) %>%
  mutate(save_file = paste0("../dgm_", true_model_id, "/survival_true_data_cut", row_number(),".rds")) %>%
  select(c(big_data_file, t, trt_group, save_file))

survival_big_data_cut <- function(big_data_file, t, trt_group, save_file) {
  
  big_data <- readRDS(big_data_file)
  
  res <- survival_big(big_data, t = t, trt_group = trt_group) 

  saveRDS(res, save_file)
  
}

package_attach <- readRDS("package_attach.rds")

sjob_data_cut_surv_true <- slurm_apply(survival_big_data_cut, 
                                  pars_dgm_surv_true_data_cut, 
                                  jobname = "surv_true_data_cut",
                                  nodes = 4, 
                                  cpus_per_node = 4, 
                                  submit = TRUE,
                                  pkgs = package_attach,
                                  global_objects = c("survival_big"),
                                  slurm_options = list(time='02:00:00',
                                                       partition='core',
                                                       "mem-per-cpu"= '16G'))

test <- readRDS("_rslurm_surv_true_data_cut/results_0.RDS")

data_cut1 <- readRDS("dgm_true_mod1/survival_true_data_cut1.rds")
data_cut2 <- readRDS("dgm_true_mod1/survival_true_data_cut2.rds")
data_cut3 <- readRDS("dgm_true_mod1/survival_true_data_cut3.rds")

data_cut1
data_cut2
data_cut3

trial_censoring_parameters <- readRDS("trial_censoring_parameters.rds")
trial_censoring_parameters <- trial_censoring_parameters %>%
  mutate(data_cut_labels = paste0(maxT, " years"))

#trial_censoring_parameters
trial_censoring_parameters <- trial_censoring_parameters %>%
  add_column(tibble(survival = round(c(data_cut1$value, data_cut2$value, data_cut3$value), 2)))

trial_censoring_parameters <- trial_censoring_parameters %>%
  mutate(survival = 100*survival) %>%
  mutate(survival = paste0(survival, "%"))

trial_censoring_parameters <- trial_censoring_parameters %>%
  mutate(data_cut_labels = paste0(data_cut_labels, " (", survival, " Survival)"))  
  
data_cut_models_index <- trial_censoring_parameters$data_cut_model_id
data_cut_labels_index <- trial_censoring_parameters$data_cut_labels 

estimand_labels <- readRDS("estimand_labels.rds")

rmst_estimand_vec <- estimand_labels %>%
  filter(estimand == "rmst", trt == 0) %>%
  pull(estimand_id)

for(rmst_estimand in rmst_estimand_vec){
  
 # rmst_estimand <- "rmst2_trt0"
#  
  # for(i in 1:2){
  #   

    #extra_knots_models
    
    scen_df <- new_model %>%
      mutate("Scenarios" = factor(new_model_id, 
                                  levels = new_model_levels,
                                  labels = new_model_labels)) %>%
      mutate(add_knots = if_else(is.na(add_knots), "default", add_knots)) %>%
      # mutate(add_knots %in% c("default", "extra_knots5")) %>%
      filter(weibull_model_id == paste0("weibull_mod",1)) %>%
      filter(design_id == "single_arm",
             stan_fit_method == "mcmc",
             mspline_df == 10,
             add_knots %in% c("none", "extra_knots1"),
             external_bias_model_id %in% c("none","external_bias_mod1"),
             backhaz == T) %>%
      left_join(performance_res, by = "scenario_fit_id")  %>%
      filter(estimand_id == rmst_estimand) %>%
      mutate("Data Cut-off" = factor(data_cut_model_id,
                                                      levels = data_cut_models_index ,
                                                      labels = data_cut_labels_index )) %>%
      #mutate(`Extra_knots" = as.factor(`Extra knots`)) %>%
      mutate("GPM rates" = if_else(backhaz, "Included", "Not included")) %>%
      mutate("Stan fit method" = if_else(stan_fit_method == "mcmc", "MCMC", "Opt" )) %>%
      arrange(Scenarios) %>%
      arrange(`Data Cut-off`) %>%
      group_by(stat) %>%
      mutate(y_height = -row_number()) %>%
      ungroup() 
    
    #View(scen_df)
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
      select(upper) %>%
      filter(!is.na(upper)) %>%
      pull() %>%
      max() %>%
      max(true_value) + 0.3
    
    y_min <- scen_df %>%
      select(y_height) %>%
      pull() %>%
      min() - 1
    
    y_max <- scen_df %>%
      select(y_height) %>%
      pull() %>%
      max() + 0.5
    
    x_axis_label <- scen_df %>%
      select(estimand, t) %>%
      distinct() %>%
      mutate(label = paste0(toupper(estimand), " at ", t, "-y")) %>%
      pull(label)
    
    #View(scen_df)
    
    col_values <- hue_pal()(3)
    fill_values <- hue_pal()(3)
    
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
            legend.key.spacing.y = unit(4, "pt"),
            legend.text = element_text(size=8)) + 
      geom_segment(x = true_value, 
                   xend = true_value, 
                   y = y_min-0.5, 
                   yend = y_max, 
                   colour = "gray50",
                   alpha = 0.3)+
      geom_point(aes(#group = scenario_fit_id, 
        colour = `Data Cut-off`,
        fill = `Data Cut-off`,
        shape = Scenarios),
        na.rm = FALSE,
        alpha = 1,
        stroke = 1,
        size = 2) + 
      geom_segment(aes(x= true_value, xend = est, y = y_height, yend = y_height,
                       colour = `Data Cut-off`),  alpha = 0.35)+
      # geom_point(aes(x=lower, y = y_height, colour = `Data Cut-off`),
      #            shape = 91, size = 3)+
      # geom_point(aes(x=upper, y = y_height, colour = `Data Cut-off`),
      #            shape = 93, size = 3)+
      scale_colour_manual(values = col_values)+
      scale_shape_manual("Settings", values = c(21,22,24))+
      scale_fill_manual(values = fill_values)+
      scale_x_continuous(x_axis_label, limits = c(x_min, x_max))+
      guides(                              
        shape = guide_legend(override.aes=list(colour = "gray60",
                                               fill = "gray60"),
                             order = 3),
        colour = guide_legend(override.aes=list(shape = 21,
                                                size = 1.7,
                                                stroke = 1,
                                                fill = fill_values),
                              order = 1),
        fill = guide_legend(override.aes=list(shape = 21,
                                              size = 1.7,
                                              stroke = 1,
                                              fill = c("white")),
                            order = 1))

    
    forest_plot  
    saveRDS(forest_plot, paste0("plots/single_arm/forest_",rmst_estimand,"_single_arm_plot_datacut.rds"))
    
    
  #   
  # }
  # 
  # 
  # 
  # plot1 <- readRDS(paste0("plots/single_arm/forest_",rmst_estimand,"_single_arm_plot1_datacut.rds"))
  # #plot1
  # 
  # plot2 <- readRDS(paste0("plots/single_arm/forest_",rmst_estimand,"_single_arm_plot2_datacut.rds"))
  # #plot2
  # 
  # plot_legend <- cowplot::get_legend(
  #   # create some space to the left of the legend
  #   plot1 + theme(legend.box.margin = margin(0, 12, 0, 12))
  # )
  # 
  # #plot_legend
  # plot_all <- plot_grid(plot_grid(plot1+
  #                                   theme(legend.position="none",
  #                                         plot.title = element_text(size=10, face="bold"))+
  #                                   labs(title = "(a) Weibull mixture"),
  #                                 plot2+
  #                                   theme(legend.position="none",
  #                                         plot.title = element_text(size=10, face="bold"))+
  #                                   labs(title = "(b) Weibull cure"),
  #                                 align = "h",
  #                                 #labels = c("(a) Mixture Weibull", 
  #                                 #           "(b) Weibull cure"),
  #                                 #label_y = 1.00,
  #                                 #label_size = 10,
  #                                 rel_widths=c(0.5,0.5)),
  #                       plot_legend,
  #                       rel_widths=c(1, 1),
  #                       axis = "l"
  # )
  # 
  
  tiff(file = paste0("plots/single_arm/forest_",rmst_estimand,"_single_arm_datacut.tiff"),   
       width = 5.6, 
       height = 3.3,
       units = 'in',  
       res = 300, 
       compression = "lzw")
  print(forest_plot)
  dev.off()
  
  print(paste0("plots/single_arm/forest_",rmst_estimand,"_single_arm_datacut.tiff"))
  
}




