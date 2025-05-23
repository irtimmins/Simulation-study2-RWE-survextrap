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

# Set directory to where simulation results are stored.
store_res <- "directory/to/store/simulations"
setwd(store_res)

scenarios <- readRDS("scenarios.rds")
performance_res <- readRDS("rmst_and_irmst_performance.rds")

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
  select(external_bias_model_id, external_data_label) %>%
  add_row(external_bias_model_id = "none", external_data_label = "No external data",
          .before = 1)

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
         Model %in% c("PH", "Non-PH")) %>%
  mutate(new_model_id = paste0(external_bias_model_id, "_", add_knots)) %>%
  mutate(new_model_id_labels = 0)

for(i in 1:nrow(new_model)){

  if(new_model$external_bias_model_id[i] == "none"){

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

two_arm_waning_settings <- readRDS("two_arm_waning_settings.rds")
two_arm_waning_settings$waning_model_labels <- 0

for(i in 1:nrow(two_arm_waning_settings)){
  if(is.na(two_arm_waning_settings$wane_period_start[i])){
    two_arm_waning_settings$waning_model_labels[i] <- "No treatment effect waning"
  } else {
    two_arm_waning_settings$waning_model_labels[i] <- paste0(
      "Waning from ",two_arm_waning_settings$wane_period_start[i], 
      " to ", two_arm_waning_settings$wane_period_stop[i], " years")
  }
}

waning_model_levels <- two_arm_waning_settings$waning_model_id
waning_model_labels <- two_arm_waning_settings$waning_model_labels

####################################################
# Difference in RMST across arms.
####################################################

estimand_labels <- readRDS("estimand_labels.rds")

irmst_estimand_vec <- estimand_labels %>%
  filter(estimand == "irmst") %>%
  pull(estimand_id)

for(irmst_estimand in irmst_estimand_vec){
  
  for(model in c("PH", "Non-PH")){
 
    for(j in 1:3){
    
    print(paste0("trt_effect_mod", j))

    scen_df <- 
      new_model %>%
      filter(external_bias_model_id %in% c("none", paste0("external_bias_mod", 1))) %>%
      filter(trt_effect_model_id == paste0("trt_effect_mod", j)) %>%
      filter(Model == model) %>%
      left_join(performance_res, by = "scenario_fit_id")  %>%
      filter(estimand_id == irmst_estimand) %>%
      mutate("Scenarios" = factor(new_model_id, 
                                  levels = new_model_levels,
                                  labels = new_model_labels)) %>%
      mutate("Waning assumption" = factor(waning_model_id, 
                                          levels = waning_model_levels,
                                          labels = waning_model_labels)) %>%
      mutate("GPM rates" = if_else(backhaz, "Included", "Not included")) %>%
      mutate("All models" = paste0(external_bias_model_id, "_", add_knots)) %>%
      arrange(Scenarios) %>%
      arrange(backhaz) %>%
      arrange(include_external_data) %>%
      group_by(stat) %>%
      mutate(y_height = -row_number()) %>%
      ungroup()

    scen_all_df <- 
      new_model %>%
      filter(external_bias_model_id %in% c("none", paste0("external_bias_mod", 1))) %>%
      left_join(performance_res, by = "scenario_fit_id")  %>%
      filter(estimand_id == irmst_estimand) %>%
      filter(stat == "mean")
    
    x_min <- scen_all_df %>%
      select(lower) %>%
      filter(!is.na(lower)) %>%
      pull() %>%
      min() - 0.3
    
    x_max <- scen_all_df %>%
      select(upper) %>%
      filter(!is.na(upper)) %>%
      pull() %>%
      max() + 0.3
    
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
    print(true_value)
    #true_value
    #names(scen_df)
    
    x_axis_label <- scen_df %>%
      select(estimand, t) %>%
      distinct() %>%
      mutate(label = paste0("RMSTD", " at ", t, "-y")) %>%
      pull(label)
    
    colour_values <- c(hue_pal()(6)[1], hue_pal()(6)[c(1,4)])
    fill_values <-c("white",hue_pal()(6)[1],  hue_pal()(6)[c(4)])

    forest_plot <- scen_df %>%
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
                     shape = `Waning assumption`),
                 na.rm = FALSE,
                 alpha = 1,
                 stroke = 1,
                 size = 2) + 
      geom_segment(aes(x= true_value, xend = est, y = y_height, yend = y_height,
                       colour = Scenarios),  alpha = 0.35)+
      # geom_point(aes(x=lower, y = y_height, colour = Scenarios), 
      #            shape = 91, size = 3, alpha = 0.5)+
      # geom_point(aes(x=upper, y = y_height, colour = Scenarios), 
      #            shape = 93, size = 3, alpha = 0.5)+
      scale_colour_manual("Settings", values = colour_values)+
      scale_fill_manual("Settings", values = fill_values)+
      scale_shape_manual(values = c(21,22,23,24))+
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
    
    model_str <- if_else(model == "PH", "ph", "nonph")
    saveRDS(forest_plot, paste0("plots/two_arm/forest_", irmst_estimand, "_model", j,"_", model_str, ".rds"))
    
    }
    
  plot1 <- readRDS(paste0("plots/two_arm/forest_", irmst_estimand,"_model", 1,"_", model_str, ".rds"))
  plot2 <- readRDS(paste0("plots/two_arm/forest_", irmst_estimand,"_model", 2,"_", model_str, ".rds"))
  plot3 <- readRDS(paste0("plots/two_arm/forest_", irmst_estimand,"_model", 3,"_", model_str, ".rds"))

  plot_legend <- cowplot::get_legend(
    # create some space to the left of the legend
    plot1 + theme(legend.box.margin = margin(0, 12, 0, 12))
  )
  
  # plot_legend
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
                        rel_widths=c(1, 0.7),
                        axis = "l"
  )

  saveRDS(plot_all, paste0("plots/two_arm/forest_", irmst_estimand, "_", model_str, ".rds"))
  
  tiff(file = paste0("plots/two_arm/forest_", irmst_estimand, "_", model_str, ".tiff"),   
       width = 7.5, 
       height = 6.8,
       units = 'in',  
       res = 300, 
       compression = "lzw")
  print(plot_all)
  dev.off()
  
  }
  
  plot1_ph <- readRDS(paste0("plots/two_arm/forest_", irmst_estimand,"_model", 1,"_", "ph", ".rds"))
  plot1_ph_no_legend <- plot1_ph +
    theme(legend.position="none",
          plot.title = element_text(size=10, face="bold"))+
    labs(title = "(a) Scenario 1: \nConstant effect,\nPH models")
  plot2_ph <- readRDS(paste0("plots/two_arm/forest_", irmst_estimand,"_model", 2,"_", "ph", ".rds"))
  plot2_ph_no_legend <- plot2_ph +
    theme(legend.position="none", 
          plot.title = element_text(size=10, face="bold"))+
    labs(title = "(b) Scenario 2: \nWaning effect,\nPH models")
  plot3_ph <- readRDS(paste0("plots/two_arm/forest_", irmst_estimand,"_model", 3,"_", "ph", ".rds"))
  plot3_ph_no_legend <- plot3_ph +
    theme(legend.position="none",
          plot.title = element_text(size=10, face="bold"))+
    labs(title = "(c) Scenario 3: \nDelayed then \nwaning effect,\nPH models")
  
  plot1_nonph <- readRDS(paste0("plots/two_arm/forest_", irmst_estimand,"_model", 1,"_", "nonph", ".rds"))
  plot1_nonph_no_legend <- plot1_nonph + 
    theme(legend.position="none",
          plot.title = element_text(size=10, face="bold"))+
    labs(title = "(d) Scenario 1: \nConstant effect,\nNon-PH models")

  plot2_nonph <- readRDS(paste0("plots/two_arm/forest_", irmst_estimand,"_model", 2,"_", "nonph", ".rds"))
  plot2_nonph_no_legend <- plot2_nonph + 
  theme(legend.position="none", 
        plot.title = element_text(size=10, face="bold"))+
    labs(title = "(e) Scenario 2: \nWaning effect,\nNon-PH models")
  plot3_nonph <- readRDS(paste0("plots/two_arm/forest_", irmst_estimand,"_model", 3,"_", "nonph", ".rds"))
  plot3_nonph_no_legend <- plot3_nonph + 
    theme(legend.position="none",
          plot.title = element_text(size=10, face="bold"))+
    labs(title = "(f) Scenario 3: \nDelayed then \nwaning effect,\nNon-PH models")
  
  
  plot_ph_legend <- cowplot::get_legend(
    # create some space to the left of the legend
    plot1_ph + theme(legend.box.margin = margin(0, 22, 0, 12))
  )
  

  plot_ph <- plot_grid(plot1_ph_no_legend,
                       plot2_ph_no_legend,
                       plot3_ph_no_legend,
                        rel_widths=c(1,1,1),
                        nrow = 1,
                        axis = "l",
                        align = "hv")
  
  
  plot_nonph <- plot_grid(plot1_nonph_no_legend,
                       plot2_nonph_no_legend,
                       plot3_nonph_no_legend,
                       rel_widths=c(1,1,1),
                       nrow = 1,
                       axis = "l",
                       align = "hv")
  
  
  plot_both <- plot_grid(plot_ph, 
                         plot_nonph,
                         ncol = 1,
                         axis = "l",
                         align = "hv")
  
  plot_both_with_legend <- plot_grid(plot_both,
                                     plot_ph_legend, 
                                     ncol = 2,
                                     rel_widths = c(1, 0.45),
                                     axis = "l",
                                     align = "hv")

  
  plot_both_with_legend  

  saveRDS(plot_both_with_legend, paste0("plots/two_arm/forest_", irmst_estimand, "_waning.rds"))
  
  tiff(file = paste0("plots/two_arm/forest_", irmst_estimand, "_waning.tiff"),   
       width = 7.6, 
       height = 5.4,
       units = 'in',  
       res = 300, 
       compression = "lzw")
  print(plot_both_with_legend)
  dev.off()
  
  
}

