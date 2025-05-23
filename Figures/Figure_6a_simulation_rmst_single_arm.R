#########################################################
# Figure 6a, Simulation study, rmst plot.
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

# Jobname where results are stored.
stores_res <- "directory/to/store/simulations"
setwd(store_res)

scenarios <- readRDS("scenarios.rds")
performance_res <- readRDS("rmst_and_irmst_performance.rds")

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
new_model

################################################
# RMST for single arm trials.
################################################

rmst_estimand_vec <- estimand_labels %>%
  filter(estimand == "rmst", trt == 0) %>%
  pull(estimand_id)

# Generate plot for each rmst time point.

for(rmst_estimand in rmst_estimand_vec){
  

    scen_df <-
      new_model  %>%
      left_join(performance_res, by = "scenario_fit_id")  %>%
      filter(estimand_id == rmst_estimand) %>%
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
      max() + 1
    
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
                      shape = Scenarios),
              na.rm = FALSE,
              alpha = 1,
              stroke = 1,
              size = 2)+ 
     geom_segment(aes(x= true_value, xend = est, y = y_height,
                         yend = y_height, colour = Scenarios),  alpha = 0.35)+
     # geom_point(aes(x=lower, y = y_height, colour = Scenarios), 
     #              shape = 91, size = 3, alpha = 0.5)+
     # geom_point(aes(x=upper, y = y_height, colour = Scenarios), 
     #                 shape = 93, size = 3, alpha = 0.5)+
     scale_colour_manual("Settings", values = colour_values)+
     scale_shape_manual("Settings", values = rep(21,7))+
     scale_x_continuous(x_axis_label, limits = c(x_min, x_max))+
     scale_fill_manual("Settings", values = fill_values)+         
   theme(plot.title = element_text(size=10, face="bold"))+
   labs(title = "(a) Weibull mixture")+
   guides(                              
     shape = guide_legend(override.aes=list(shape = 21,
                                            stroke = 1,
                                            size = 1.7,
                                            alpha = 1,
                                            colour = colour_values,
                                            fill = fill_values)))
      
 forest_plot 
 saveRDS(forest_plot, paste0("plots/single_arm/forest_",rmst_estimand,"_single_arm_plot.rds"))
    

  tiff(file = paste0("plots/single_arm/forest_",rmst_estimand,"_single_arm.tiff"),   
       width = 5.2, 
       height = 3.2,
       units = 'in',  
       res = 1200, 
       compression = "lzw")
  print(forest_plot)
  dev.off()
  
  print(paste0("plots/single_arm/forest_",rmst_estimand,"_single_arm.tiff"))
  
  
}
  
