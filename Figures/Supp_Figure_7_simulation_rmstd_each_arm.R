
#########################################################
# Supp Figure 7, Simulation study, rmst plot.
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
# RMST for each arm, in two arm trials (condensed)
################################################

scenarios <- readRDS("scenarios.rds")
scenarios_freq <- readRDS("scenarios_two_arm_freq.rds")

performance_res <- readRDS("rmst_and_irmst_performance.rds")
performance_freq_res <- readRDS("rmst_and_irmst_performance_two_arm_freq.rds")

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
         design_id == "two_arm",
         stan_fit_method == "mcmc",
         mspline_df == 10,
         data_cut_model_id == "data_cut_mod2",
         add_knots %in% c("none", "extra_knots1"),
         backhaz == T,
         waning_model_id == "waning_mod1") %>%
  mutate(new_model_id = paste0(external_bias_model_id, "_", add_knots)) %>%
  mutate(new_model_id_labels = 0) 


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
  filter(estimand == "rmst") %>%
  pull(estimand_id)

rmst_estimand_vec

for(rmst_estimand in rmst_estimand_vec){
  
  #rmst_estimand <- "rmst1_trt0"
  
  for(j in 1:3){
    #j <- 1
    scen_survextrap_all_df <- 
      new_model  %>%
      filter(external_bias_model_id %in% c("none", paste0("external_bias_mod", c(1,3,5)))) %>%
    #  filter(trt_effect_model_id == paste0("trt_effect_mod", j)) %>%
      left_join(performance_res, by = "scenario_fit_id")  %>%
      filter(estimand_id == rmst_estimand) %>%
      mutate("Scenarios" = factor(new_model_id, 
                                  levels = new_model_levels,
                                  labels = new_model_labels )) %>%
      mutate("External data" = factor(external_bias_model_id,
                                      levels = external_data_models_index,
                                      labels = external_data_labels_index)) %>%
      mutate("Extra knots" = factor(add_knots,
                                    levels = extra_knots_models_index,
                                    labels = extra_knots_labels_index)) %>%
      mutate("GPM rates" = if_else(backhaz, "Included", "Not included")) %>%
      mutate("All models" = paste0(external_bias_model_id, "_", add_knots)) %>%
      rename(Model = model) %>%
      arrange(Scenarios) %>%
      select(estimand,  t, trt, stat,est,lower, upper, scenario_fit_id, Scenarios, Model, trt_effect_model_id)
    
    #summary(as.factor(scen_survextrap_all_df$Model))
    
    scen_freq_all_df <- scenarios_freq %>%
    #  filter(trt_effect_model_id == paste0("trt_effect_mod", j)) %>%
      left_join(performance_freq_res %>% 
                  mutate(two_arm_freq_model_id = scenario_fit_id),
                by = "two_arm_freq_model_id")  %>%
      filter(estimand_id == rmst_estimand) %>%
      mutate(Scenarios = NA) %>%
      mutate(Model = parametric_model) %>%
      mutate(scenario_fit_id = two_arm_freq_model_id) %>%
      select(estimand, t, trt, stat, est, lower, upper, scenario_fit_id, Scenarios, Model, trt_effect_model_id)
    
    #summary(as.factor(scen_freq_all_df$Model))
    
    # Create new factors/labels..
    
   new_levels <- c("exp", "weibull", "spline", "ph", "nonph", "separate")
   new_labels <- c("Exponential, separate arms", "Weibull, separate arms", "Royston-Parmar (df = 3),\nseparate arms",
                   "survextrap, PH", "survextrap, Non-PH", "survextrap, separate arms") 
    
   scen_all_df <- bind_rows(scen_freq_all_df, 
                         scen_survextrap_all_df) %>%
     mutate(Model = factor(Model, levels = new_levels, labels = new_labels))
   
   #summary(scen_all_df$Model)
   
   scen_df <-  scen_all_df  %>%
      filter(trt_effect_model_id == paste0("trt_effect_mod", j)) %>%
      group_by(stat) %>%
      mutate(y_height = -row_number()) %>%
      ungroup() 
   # 
   # scen_df %>%
   #   filter(stat == "mean") %>%
   #   View()
    
   # summary(as.factor(scen_df$y_height)) 
   # summary(as.factor(scen_df$Scenarios))
   # summary(as.factor(scen_df$Model))
   # summary(as.factor(scen_df$est))

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
    
    x_min <- scen_all_df %>%
      filter(stat == "mean") %>% 
      select(lower) %>%
      filter(!is.na(lower)) %>%
      pull() %>%
      min() %>%
      min(true_value) - 0.3
    
    x_max <- scen_all_df %>%
      filter(stat == "mean") %>% 
      select(lower) %>%
      filter(!is.na(lower)) %>%
      pull() %>%
      max() %>%
      max(true_value) + 0.3
    
    trt_label <- scen_df %>%
      select(trt) %>%
      distinct() %>%
      mutate(label = if_else(trt == 0, "Control", "Active")) %>%
      pull() 
    
    x_axis_label <- scen_df %>%
      select(estimand, t) %>%
      distinct() %>%
      mutate(label = paste0("RMST at ", t, "-y \n ", trt_label, " arm" )) %>%
      pull(label)
    
    summary(scen_df$Scenarios)
    
    colour_values <- c(hue_pal()(6)[1], hue_pal()(6)[c(1,2,4,6)])
    fill_values <- c( "white", hue_pal()(6)[1],  hue_pal()(6)[c(2,4,6)])
    shape_values <- c( 3, 4, 8, 21,22,24)
    size_values <- c(1.2, 1.2, 1.2, 2, 2, 2)
    
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
            legend.box.background = element_rect(colour = "black")) + 
      geom_segment(x = true_value, 
                   xend = true_value, 
                   y = y_min-0.5, 
                   yend = y_max, 
                   colour = "gray50",
                   alpha = 0.3)+
      geom_point(aes(group = scenario_fit_id, 
                     colour = Scenarios,
                     fill = Scenarios,
                     shape = Model,
                     size = Model),
                 na.rm = FALSE,
                 alpha = 1,
                 stroke = 1) + 
      geom_segment(aes(x= true_value, xend = est, y = y_height, yend = y_height,
                       colour = Scenarios),  alpha = 0.35)+
      # geom_point(aes(x=lower, y = y_height, colour = Scenarios), 
      #            shape = 91, size = 3, alpha = 0.5)+
      # geom_point(aes(x=upper, y = y_height, colour = Scenarios), 
      #            shape = 93, size = 3, alpha = 0.5)+
      scale_colour_manual("Settings for survextrap\nmodels", values = colour_values)+
      scale_fill_manual("Settings for survextrap\nmodels", values = fill_values)+
      scale_shape_manual(values = shape_values)+
      scale_size_manual(values = size_values)+
      scale_x_continuous(x_axis_label, limits = c(x_min, x_max))+
      guides(                              
        shape = guide_legend(override.aes=list(colour = "gray60",
                                               fill = "gray60"),
                             order = 1),
        colour = guide_legend(override.aes=list(shape = 21,
                                                size = 1.7,
                                                stroke = 1,
                                                fill = c(fill_values, "white"),
                              order = 3)),
        size = "none")

    forest_plot
    
    saveRDS(forest_plot, paste0("plots/two_arm/forest_", rmst_estimand, "_model", j,"_two_arm_cond2.rds"))
    
  }
  
}


################################################
# Difference in RMST for each arm, in two arm trials (condensed)
################################################

estimand_labels <- readRDS("estimand_labels.rds")

irmst_estimand_vec <- estimand_labels %>%
  filter(estimand == "irmst") %>%
  pull(estimand_id)

for(irmst_estimand in irmst_estimand_vec){
  
  #irmst_estimand <- "irmst1"
  
  for(j in 1:3){
    #  j <- 1
    #names(scenarios)
    #View(scenarios)
    #summary(as.factor(scenarios$add_knots))
    #summary(as.factor(scenarios$backhaz))
    scen_survextrap_all_df <- 
      new_model  %>%
      filter(external_bias_model_id %in% c("none", paste0("external_bias_mod", c(1,3,5)))) %>%
      #  filter(trt_effect_model_id == paste0("trt_effect_mod", j)) %>%
      left_join(performance_res, by = "scenario_fit_id")  %>%
      filter(estimand_id == irmst_estimand) %>%
      mutate("Scenarios" = factor(new_model_id, 
                                  levels = new_model_levels,
                                  labels = new_model_labels )) %>%
      mutate("External data" = factor(external_bias_model_id,
                                      levels = external_data_models_index,
                                      labels = external_data_labels_index)) %>%
      mutate("Extra knots" = factor(add_knots,
                                    levels = extra_knots_models_index,
                                    labels = extra_knots_labels_index)) %>%
      mutate("GPM rates" = if_else(backhaz, "Included", "Not included")) %>%
      mutate("All models" = paste0(external_bias_model_id, "_", add_knots)) %>%
      rename(Model = model) %>%
      arrange(Scenarios) %>%
      select(estimand,  t, trt, stat,est,lower, upper, scenario_fit_id, Scenarios, Model, trt_effect_model_id)
    
    #summary(as.factor(scen_survextrap_all_df$Model))
    
    scen_freq_all_df <- scenarios_freq %>%
      #  filter(trt_effect_model_id == paste0("trt_effect_mod", j)) %>%
      left_join(performance_freq_res %>% 
                  mutate(two_arm_freq_model_id = scenario_fit_id),
                by = "two_arm_freq_model_id")  %>%
      filter(estimand_id == irmst_estimand) %>%
      mutate(Scenarios = NA) %>%
      mutate(Model = parametric_model) %>%
      mutate(scenario_fit_id = two_arm_freq_model_id) %>%
      select(estimand, t, trt, stat, est, lower, upper, scenario_fit_id, Scenarios, Model, trt_effect_model_id)
    
    #summary(as.factor(scen_freq_all_df$Model))
    
    # Create new factors/labels..
    
    new_levels <- c("exp", "weibull", "spline", "ph", "nonph", "separate")
    new_labels <- c("Exponential, separate arms", "Weibull, separate arms", "Royston-Parmar (df = 3),\nseparate arms",
                    "survextrap, PH", "survextrap, Non-PH", "survextrap, separate arms") 
    
    scen_all_df <- bind_rows(scen_freq_all_df, 
                             scen_survextrap_all_df) %>%
      mutate(Model = factor(Model, levels = new_levels, labels = new_labels))
    
    #summary(scen_all_df$Model)
    
    scen_df <-  scen_all_df  %>%
      filter(trt_effect_model_id == paste0("trt_effect_mod", j)) %>%
      group_by(stat) %>%
      mutate(y_height = -row_number()) %>%
      ungroup() 

    # summary(as.factor(scen_df$add_knots))
    # summary(as.factor(scen_df$))
    # View(scen_df)

    
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
    
    x_min <- scen_all_df %>%
      filter(stat == "mean") %>% 
      select(lower) %>%
      filter(!is.na(lower)) %>%
      pull() %>%
      min() %>%
      min(true_value) - 0.3
    
    x_max <- scen_all_df %>%
      filter(stat == "mean") %>% 
      select(lower) %>%
      filter(!is.na(lower)) %>%
      pull() %>%
      max() %>%
      max(true_value) + 0.3
    
    trt_label <- scen_df %>%
      select(trt) %>%
      distinct() %>%
      mutate(label = if_else(trt == 0, "Control", "Active")) %>%
      pull() 

    x_axis_label <- scen_df %>%
      select(estimand, t) %>%
      distinct() %>%
      mutate(label = paste0(" \u0394RMST", " at ", t, "-y")) %>%
      pull(label)
    
    
   
    
    colour_values <- c(hue_pal()(6)[1], hue_pal()(6)[c(1,2,4,6)])
    fill_values <- c( "white", hue_pal()(6)[1],  hue_pal()(6)[c(2,4,6)])
    shape_values <- c( 3, 4, 8, 21,22,24)
    size_values <- c(1.2, 1.2, 1.2, 2, 2,  2)
    
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
                     colour = Scenarios,
                     fill = Scenarios,
                     shape = Model,
                     size = Model),
                 na.rm = FALSE,
                 alpha = 1,
                 stroke = 1) + 
      geom_segment(aes(x= true_value, xend = est, y = y_height, yend = y_height,
                       colour = Scenarios),  alpha = 0.35)+
      # geom_point(aes(x=lower, y = y_height, colour = Scenarios), 
      #            shape = 91, size = 3, alpha = 0.5)+
      # geom_point(aes(x=upper, y = y_height, colour = Scenarios), 
      #            shape = 93, size = 3, alpha = 0.5)+
      scale_colour_manual("Settings for survextrap\nmodels", values = colour_values)+
      scale_fill_manual("Settings for survextrap\nmodels", values = fill_values)+
      scale_shape_manual(values = shape_values)+
      scale_size_manual(values = size_values)+
      scale_x_continuous(x_axis_label, limits = c(x_min, x_max), )+
      guides(                              
        shape = guide_legend(override.aes=list(colour = "gray60",
                                               fill = "gray60"),
                             order = 1),
        colour = guide_legend(override.aes=list(shape = 21,
                                                size = 1.7,
                                                stroke = 1,
                                                fill = c(fill_values, "white"),
                                                order = 3)),
        size = "none")
    
    forest_plot
    saveRDS(forest_plot, paste0("plots/two_arm/forest_", irmst_estimand, "_model", j, "_cond2.rds"))
    
    # create new plot with improved legend (ignore points in this plot)
    forest_plot_improved_legend <- 
      bind_rows(scen_df %>% 
                  filter(!is.na(Scenarios)),
                scen_df %>% 
                  filter(is.na(Scenarios)) %>%
                  mutate(Scenarios = "No external data,\nno extra knots"))  %>%
      mutate("Scenarios" = factor(Scenarios, 
                                  levels = new_model_labels,
                                  labels = new_model_labels )) %>%
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
                     colour = Scenarios,
                     fill = Scenarios,
                     shape = Model,
                     size = Model),
                 na.rm = FALSE,
                 alpha = 1,
                 stroke = 1) + 
      geom_segment(aes(x= true_value, xend = est, y = y_height, yend = y_height,
                       colour = Scenarios),  alpha = 0.35)+
      # geom_point(aes(x=lower, y = y_height, colour = Scenarios), 
      #            shape = 91, size = 3, alpha = 0.5)+
      # geom_point(aes(x=upper, y = y_height, colour = Scenarios), 
      #            shape = 93, size = 3, alpha = 0.5)+
      scale_colour_manual("Settings for survextrap\nmodels", 
                          values = colour_values,
                          na.translate = F)+
      scale_fill_manual("Settings for survextrap\nmodels", 
                        values = fill_values,
                        na.translate = F)+
      scale_shape_manual(values = shape_values)+
      scale_size_manual(values = size_values)+
      scale_x_continuous(x_axis_label, limits = c(x_min, x_max))+
      guides(                              
        shape = guide_legend(override.aes=list(colour = "gray60",
                                               fill = "gray60"),
                             order = 1),
        colour = guide_legend(override.aes=list(shape = 21,
                                                size = 1.7,
                                                stroke = 1,
                                                fill = fill_values,
                                                order = 3)),
        size = "none")
    
    
    forest_plot_improved_legend
    saveRDS(forest_plot_improved_legend, paste0("plots/two_arm/forest_", irmst_estimand, "_model", j, "_cond2_legend.rds")) 
  }
  
  
  
  plot1 <- readRDS(paste0("plots/two_arm/forest_", irmst_estimand, "_model", 1, "_cond2.rds"))
  plot2 <- readRDS(paste0("plots/two_arm/forest_", irmst_estimand, "_model", 2, "_cond2.rds"))
  plot3 <- readRDS(paste0("plots/two_arm/forest_", irmst_estimand, "_model", 3, "_cond2.rds"))
  #plot1
  # plot2
  #plot3
  plot_with_improved_legend <- readRDS(paste0("plots/two_arm/forest_", irmst_estimand, "_model", 1, "_cond2_legend.rds"))
  
  # import legend
  plot_legend <- cowplot::get_legend(
    # create some space to the left of the legend
    plot_with_improved_legend + theme(legend.box.margin = margin(0, 12, 0, 12))
  )
  
  plot_legend
  plot_all <- plot_grid(plot_grid(plot1+
                                    theme(legend.position="none",
                                          plot.title = element_text(size=10, face="bold"))+
                                    labs(title = "(b) Scenario 1: \n Constant effect"),
                                  plot2+
                                    theme(legend.position="none",
                                          plot.title = element_text(size=10, face="bold"))+
                                    labs(title = "(c) Scenario 2: \n Waning effect"),
                                  plot3+
                                    theme(legend.position="none",
                                          plot.title = element_text(size=10, face="bold"))+
                                    labs(title = "(d) Scenario 3: Delayed then \n waning effect"),
                                  align = "h",
                                  rel_widths=c(1,1,1),
                                  nrow = 1),
                        plot_legend,
                        rel_widths=c(1, 0.5),
                        axis = "l"
  )
  
  plot_all
  
  tiff(file = paste0("plots/two_arm/forest_", irmst_estimand, "_cond2.tiff"),   
       width = 6.5, 
       height = 4.8,
       units = 'in',  
       res = 300, 
       compression = "lzw")
  print(plot_all)
  dev.off()
  
  print(paste0("plots/two_arm/forest_", irmst_estimand, "_cond2.tiff"))
  
}



################################################
# Combine plots for RMST for each arm, and the
# difference in RMST, all for two arm trials.
################################################

estimand_labels <- readRDS("estimand_labels.rds")

# RMST for single arm trials.

rmst_estimand_number <- estimand_labels %>%
  filter(estimand == "rmst") %>%
  pull(estimand_id) %>%
  length()/2

for(rmst_estimand in paste0("rmst", 1:rmst_estimand_number)){
  
  #rmst_estimand <- "rmst1"
  plot1_control <- readRDS(paste0("plots/two_arm/forest_", rmst_estimand, "_trt0_model", 1,"_two_arm_cond2.rds"))
  plot1_active <- readRDS(paste0("plots/two_arm/forest_", rmst_estimand, "_trt1_model", 1,"_two_arm_cond2.rds"))
  plot1_diff <- readRDS(paste0("plots/two_arm/forest_i", rmst_estimand, "_model", 1,"_cond2.rds"))
  
  plot2_control <- readRDS(paste0("plots/two_arm/forest_", rmst_estimand, "_trt0_model", 2,"_two_arm_cond2.rds"))
  plot2_active <- readRDS(paste0("plots/two_arm/forest_", rmst_estimand, "_trt1_model", 2,"_two_arm_cond2.rds"))
  plot2_diff <- readRDS(paste0("plots/two_arm/forest_i", rmst_estimand, "_model", 2,"_cond2.rds"))
  
  plot3_control <- readRDS(paste0("plots/two_arm/forest_", rmst_estimand, "_trt0_model", 3,"_two_arm_cond2.rds"))
  plot3_active <- readRDS(paste0("plots/two_arm/forest_", rmst_estimand, "_trt1_model", 3,"_two_arm_cond2.rds"))
  plot3_diff <- readRDS(paste0("plots/two_arm/forest_i", rmst_estimand, "_model", 3,"_cond2.rds"))
  
  
  plot_with_improved_legend <- readRDS(paste0("plots/two_arm/forest_", irmst_estimand, "_model", 1, "_cond2_legend.rds"))
  
  
  plot_legend <- cowplot::get_legend(
    # create some space to the left of the legend
    plot_with_improved_legend + theme(legend.box.margin = margin(0, 12, 0, 12))
  )
  

  plot_scenario1 <- plot_grid(plot1_control+
                                theme(legend.position="none",
                                      plot.title = element_text(size=10, face="bold")),
                              plot1_active+
                                theme(legend.position="none",
                                      plot.title = element_text(size=10, face="bold")),
                              plot1_diff+
                                theme(legend.position="none",
                                      plot.title = element_text(size=10, face="bold")),
                              nrow = 1,
                              align = "h")+
    theme(legend.position="none",
          plot.title = element_text(hjust = 0.05, size=10, face="bold"))+
    labs(title = "(a) Scenario 1: Constant effect")
  #plot_scenario1
  
  plot_scenario2 <- plot_grid(plot2_control+
                                theme(legend.position="none",
                                      plot.title = element_text(size=10, face="bold")),
                              plot2_active+
                                theme(legend.position="none",
                                      plot.title = element_text(size=10, face="bold")),
                              plot2_diff+
                                theme(legend.position="none",
                                      plot.title = element_text(size=10, face="bold")),
                              nrow = 1,
                              align = "h")+
    theme(legend.position="none",
          plot.title = element_text(hjust = 0.05, size=10, face="bold"))+
    labs(title = "(b) Scenario 2: Waning effect")
  #plot_scenario2
  
  plot_scenario3 <- plot_grid(plot3_control+
                                theme(legend.position="none",
                                      plot.title = element_text(size=10, face="bold")),
                              plot3_active+
                                theme(legend.position="none",
                                      plot.title = element_text(size=10, face="bold")),
                              plot3_diff+
                                theme(legend.position="none",
                                      plot.title = element_text(size=10, face="bold")),
                              nrow = 1,
                              align = "h")+
    theme(legend.position="none",
          plot.title = element_text(hjust = 0.08, size=10, face="bold"))+
    labs(title = "(c) Scenario 3: Delayed then waning effect")
  
  plot_without_legend <- plot_grid(plot_scenario1,
                                   plot_scenario2,
                                   plot_scenario3,
                                   align = "h",
                                   rel_widths=c(1,1,1),
                                   nrow = 3)
  
  plot_all <- plot_grid(plot_without_legend,
                        plot_legend,
                        rel_widths=c(1, 0.55),
                        axis = "l"
  )
  
  plot_all
  
  tiff(file = paste0("plots/two_arm/forest_i", rmst_estimand, "_and_", 
                     rmst_estimand, "_all_two_arm_cond2.tiff"),   
       width = 7.4, 
       height = 7.8,
       units = 'in',  
       res = 300, 
       compression = "lzw")
  print(plot_all)
  dev.off()
  
  print(paste0("plots/two_arm/forest_i", rmst_estimand, "_and_", 
               rmst_estimand, "_all_two_arm_cond2.tiff"))
  
}

