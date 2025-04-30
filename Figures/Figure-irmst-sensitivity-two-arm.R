############################################
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
setwd(store_res)

scenarios <- readRDS("scenarios.rds")
performance_res <- readRDS("rmst_and_irmst_performance.rds")


############################################################################
# Sensitivity on MCMC vs Laplace "Opt" approximation
############################################################################


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
                                            "..",
                                           "...")

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
         mspline_df == 10,
         data_cut_model_id == "data_cut_mod2",
         add_knots %in% c("none", "extra_knots1"),
         external_bias_model_id %in% c("external_bias_mod1", "none"),
         backhaz == T,
         waning_model_id == "waning_mod1") %>%
  mutate(new_model_id = paste0(external_bias_model_id, "_", add_knots)) %>%
  mutate(new_model_id_labels = 0) #%>%
#slice(c(2,3,5,4,1,6,7))

summary(as.factor(new_model$external_bias_model_id))
#summary(as.factor(new_model$new_model_id))
#View(new_model)


for(i in 1:nrow(new_model)){

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
    
}
new_model_levels <- unique(new_model$new_model_id)[c(2,3,1)]
new_model_labels <- unique(new_model$new_model_id_labels)[c(2,3,1)]
new_model_levels
new_model_labels    

estimand_labels <- readRDS("estimand_labels.rds")

irmst_estimand_vec <- estimand_labels %>%
  filter(estimand == "irmst") %>%
  pull(estimand_id)

for(irmst_estimand in irmst_estimand_vec){
  
  #irmst_estimand <- "irmst1"

    for(j in 1:3){
      #names(scenarios)
      #View(scenarios)
      #summary(as.factor(scenarios$add_knots))
      #summary(as.factor(scenarios$backhaz))
      #j <- 1
      scen_df <- 
        new_model %>%
        filter(trt_effect_model_id == paste0("trt_effect_mod", j)) %>%
        left_join(performance_res, by = "scenario_fit_id")  %>%
        filter(estimand_id == irmst_estimand) %>%
        mutate("Scenarios" = factor(new_model_id, 
                                    levels = new_model_levels,
                                    labels = new_model_labels)) %>%
        mutate("Stan fit method" = if_else(stan_fit_method == "mcmc", "MCMC", "Opt" )) %>%
        mutate("GPM rates" = if_else(backhaz, "Included", "Not included")) %>%
        mutate("All models" = paste0(external_bias_model_id, "_", add_knots)) %>%
        #arrange(`Extra knots`) %>%
        arrange(backhaz) %>%
        arrange(include_external_data) %>%
        arrange(`Stan fit method`) %>%
        arrange(Model) %>%
        arrange(Scenarios) %>%
        group_by(stat) %>%
        mutate(y_height = -row_number()) %>%
        ungroup() #%>%
      #View()
      
      # summary(as.factor(scen_df$add_knots))
      # summary(as.factor(scen_df$))
      # View(scen_df)
      
      scen_all_df <- 
        new_model %>%
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
      #true_value
      #names(scen_df)
      
      x_axis_label <- scen_df %>%
        select(estimand, t) %>%
        distinct() %>%
        mutate(label = paste0("RMSTD", " at ", t, "-y")) %>%
        pull(label)
      
      colour_values <- hue_pal()(6)[1:3]
    #  fill_values <-c("white",hue_pal()(6)[1],  hue_pal()(6)[2:6])
      
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
              legend.text = element_text(size=8),
              legend.margin = margin(r = 15)) + 
        geom_segment(x = true_value, 
                     xend = true_value, 
                     y = y_min-0.5, 
                     yend = y_max, 
                     colour = "gray50",
                     alpha = 0.3)+
        geom_point(aes(group = scenario_fit_id, 
                       colour = Scenarios,
                       fill = `Stan fit method`,
                       shape = Model),
                   na.rm = FALSE,
                   alpha = 1,
                   stroke = 0.8,
                   size = 1.5) + 
        geom_segment(aes(x= true_value, xend = est, y = y_height, yend = y_height,
                         colour = Scenarios),  alpha = 0.35)+
        geom_point(aes(x=lower, y = y_height,  colour = Scenarios), shape = 91, size = 3,  alpha = 0.35)+
        geom_point(aes(x=upper, y = y_height,  colour = Scenarios), shape = 93, size = 3,  alpha = 0.35)+
        scale_colour_discrete("Settings")+
        scale_fill_manual(values = c("white","grey"))+
        scale_shape_manual(values = c(21,22,24))+
         scale_x_continuous(x_axis_label, limit = c(x_min, x_max))+
         guides(    
           colour = guide_legend(order = 1),
           shape =  guide_legend(order = 2),
           fill = guide_legend(override.aes=list(shape=c(21,21),
                                                 fill = c("white","gray70"),
                                                 colour = "gray40"),
                               order = 3))
      
      forest_plot
      saveRDS(forest_plot, paste0("plots/two_arm/forest_", irmst_estimand, "_model", j, "_stan.rds"))
      
      
    }
    
    
    plot1 <- readRDS(paste0("plots/two_arm/forest_", irmst_estimand,"_model", 1, "_stan.rds"))
    plot2 <- readRDS(paste0("plots/two_arm/forest_", irmst_estimand,"_model", 2, "_stan.rds"))
    plot3 <- readRDS(paste0("plots/two_arm/forest_", irmst_estimand,"_model", 3, "_stan.rds"))

    #plot1
    #plot2
    #plot3
    
    
    plot_legend <- cowplot::get_legend(
      # create some space to the left of the legend
      plot1 + theme(legend.box.margin = margin(12, 12, 12, 12))
    )
    
    plot_legend
    plot_all <- plot_grid(plot_grid(plot1+
                                      theme(legend.position="none",
                                            plot.title = element_text(size=10, 
                                                                      face="bold",
                                                                      lineheight = 1.2))+
                                      labs(title = "(c) Scenario 1: \n Constant effect"),
                                    plot2+
                                      theme(legend.position="none",
                                            plot.title = element_text(size=10, 
                                                                      face="bold",
                                                                      lineheight = 1.2))+
                                      labs(title = "(d) Scenario 2: \n Waning effect"),
                                    plot3+
                                      theme(legend.position="none",
                                            plot.title = element_text(size=10, 
                                                                      face="bold",
                                                                      lineheight = 1.2))+
                                      labs(title = "(e) Scenario 3: Delayed then \n waning effect"),
                                    align = "h",
                                    rel_widths=c(1,1,1),
                                    nrow = 1),
                          plot_legend,
                          rel_widths=c(1.2, 0.8),
                          axis = "l"
    )

    
    tiff(file = paste0("plots/two_arm/forest_", irmst_estimand, "_stan.tiff"),   
         width = 7.8, 
         height = 5.2,
         units = 'in',  
         res = 300, 
         compression = "lzw")
    print(plot_all)
    dev.off()
    
    print(paste0("plots/two_arm/forest_", irmst_estimand, "_stan.tiff"))

    
    
    plot_single <- readRDS(paste0("plots/single_arm/forest_",substr(irmst_estimand, start = 2, stop = 6),"_trt0_single_arm_plot_stan.rds"))
    
    plot_single
    
    plot_single_legend <-  cowplot::get_legend(
      # create some space to the left of the legend
      plot_single + theme(legend.title = element_text(margin = margin(3, 0, 3, 0)),
                          legend.box.margin = margin(0, 20, 0, 12),
                          #    legend.text = element_text(size=7),
                          legend.key.spacing.y = unit(1, "pt"))
    )
    
    plot_single_without_legend <- plot_single+
      theme(legend.position="none",
            plot.title = element_text(size=10, face="bold",
                                      lineheight = 1.2))+
      labs(title = "Single-arm trial \n  (a) Weibull mixture")
    
    plot_single_mod <- plot_grid(NULL, plot_single_without_legend, NULL, plot_single_legend, NULL, rel_widths = c(0.3, 0.9,-0.02,1.1, 0.1), nrow = 1)
    plot_single_mod  
    
    
    #  plot_single
    
    plot_both <- plot_grid(plot_single_mod,plot_all, ncol = 1, rel_heights = c(1.2, 1.4))
    
    tiff(file = paste0("plots/single_and_two_arm/forest_", irmst_estimand, "_all_stan.tiff"),   
         width = 7.2, 
         height = 8.6,
         units = 'in',  
         res = 300, 
         compression = "lzw")
    print(plot_both)
    dev.off()
    
    
    
    
}



############################################################################
# Sensitivity on knot locations.
############################################################################


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

reorder_vec <- c(1,4,3,2)
extra_knots_models_index <- extra_knots_models$extra_knots_id[reorder_vec]
extra_knots_labels_index <- extra_knots_models$extra_knots_labels[reorder_vec]


estimand_labels <- readRDS("estimand_labels.rds")

irmst_estimand_vec <- estimand_labels %>%
  filter(estimand == "irmst") %>%
  pull(estimand_id)


for(irmst_estimand in irmst_estimand_vec){
  
  #irmst_estimand <- "irmst1"
  
  for(j in 1:3){
   # j <- 1
    
    scen_df <- 
      scenarios %>%
      mutate(add_knots = if_else(is.na(add_knots), "default", add_knots)) %>%
      mutate(include_external_data = if_else(include_external_data, "yes", "no")) %>%
      mutate("External data" = factor(external_bias_model_id,
                                      levels = external_data_models_index,
                                      labels = external_data_labels_index)) %>%
      mutate("Extra knots for\n external data" = factor(add_knots,
                                                      levels = extra_knots_models_index,
                                                      labels = extra_knots_labels_index)) %>%
      #summary(as.factor(scen_df$`Extra knots for external data`))  
      mutate("GPM rates" = if_else(backhaz, "Included", "Not included")) %>%
      rename(Model = model) %>%
      mutate(Model = case_when(Model == "ph" ~ 1,
                               Model == "nonph" ~ 2,
                               Model == "separate" ~ 3)) %>%
      mutate(Model = factor(Model, levels = c(1,2,3), 
                            labels = c("PH", 
                                       "Non-PH", 
                                       "Separate arms" ))) %>%
      mutate("Stan fit method" = if_else(stan_fit_method == "mcmc", "MCMC", "Opt" )) %>%
      filter(trt_effect_model_id == paste0("trt_effect_mod", j)) %>%
  filter(design_id == "two_arm",
         weibull_model_id == "weibull_mod1",
         stan_fit_method == "mcmc",
         data_cut_model_id == "data_cut_mod2",
         external_bias_model_id %in% c("external_bias_mod1"),
         backhaz == T,
         waning_model_id == "waning_mod1",
         mspline_df == 10) %>%
      arrange(mspline_df) %>%
      mutate(`Trial data M-spline` = paste0("df = ", mspline_df)) %>%
      mutate(`Trial data M-spline` = fct_inorder(`Trial data M-spline`)) %>%
      #arrange(`Trial data M-spline`) %>%
      arrange(`Extra knots for\n external data`) %>%
      left_join(performance_res, by = "scenario_fit_id")  %>%
      filter(estimand_id == irmst_estimand) %>%
   #   arrange(`External data`) %>%
   #   arrange(backhaz) %>%
  #    arrange(include_external_data) %>%
 #     arrange(mspline_df) %>%
      arrange(Model) %>%
      arrange(`Extra knots for\n external data`) %>%
      group_by(stat) %>%
      mutate(y_height = -row_number()) %>%
      ungroup() #%>%
    #View()
    
    # summary(as.factor(scen_df$add_knots))
    # summary(as.factor(scen_df$))
    # View(scen_df)
    
    scen_all_df <- 
      scenarios %>%
      mutate(add_knots = if_else(is.na(add_knots), "default", add_knots)) %>%
      mutate(include_external_data = if_else(include_external_data, "yes", "no")) %>%
      filter(design_id == "two_arm",
             weibull_model_id == "weibull_mod1",
             stan_fit_method == "mcmc",
             data_cut_model_id == "data_cut_mod2",
             external_bias_model_id %in% c("external_bias_mod1"),
             backhaz == T,
             waning_model_id == "waning_mod1") %>%
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
    #true_value
    #names(scen_df)
    
    col_values <- hue_pal()(3)
    fill_values <- hue_pal()(3)
    
    
    x_axis_label <- scen_df %>%
      select(estimand, t) %>%
      distinct() %>%
      mutate(label = paste0("RMSTD", " at ", t, "-y")) %>%
      pull(label)
    
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
            legend.margin = margin(t = 3, r = 15, b = 3, l = 3),
            legend.text = element_text(size=8)) + 
      geom_segment(x = true_value, 
                   xend = true_value, 
                   y = y_min-0.5, 
                   yend = y_max, 
                   colour = "gray50",
                   alpha = 0.3)+
      geom_point(aes(#group = scenario_fit_id, 
        colour = `Extra knots for\n external data`,
        shape = Model,
        fill = `Extra knots for\n external data`) ,
        #  shape = `Extra knots`),
        # fill = `External data`),
        na.rm = FALSE,
        alpha = 1,
        stroke = 1,
        size = 2) + 
      geom_segment(aes(x= true_value, xend = est, y = y_height, yend = y_height,
                       colour = `Extra knots for\n external data`),  alpha = 0.35)+
      # geom_point(aes(x=lower, y = y_height, colour = `Extra knots for\n external data`), 
      #            shape = 91, size = 3,  alpha = 0.35)+
      # geom_point(aes(x=upper, y = y_height, colour = `Extra knots for\n external data`), 
      #            shape = 93, size = 3,  alpha = 0.35)+
      scale_colour_manual( values = col_values)+
      scale_fill_manual(values = fill_values)+
      scale_shape_manual(values = c(21,22,24))+
          #  scale_fill_manual(values = c("white","gray60"))+
      scale_x_continuous(x_axis_label, limit = c(x_min, x_max))+
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
    saveRDS(forest_plot, paste0("plots/two_arm/forest_", irmst_estimand, "_model", j, "_knots.rds"))
    
    
  }

  plot1 <- readRDS(paste0("plots/two_arm/forest_", irmst_estimand, "_model", 1, "_knots.rds"))
  plot2 <- readRDS(paste0("plots/two_arm/forest_", irmst_estimand, "_model", 2, "_knots.rds"))
  plot3 <- readRDS(paste0("plots/two_arm/forest_", irmst_estimand, "_model", 3, "_knots.rds"))

  plot_legend <- cowplot::get_legend(
    # create some space to the left of the legend
    plot1 + theme(legend.box.margin = margin(0, 6, 0, 6))
  )
  
  plot_legend
  plot_all <- plot_grid(plot_grid(plot1+
                                    theme(legend.position="none",
                                          plot.title = element_text(size=10, 
                                                                    face="bold",
                                                                    lineheight = 1.2))+
                                    labs(title = "Two-arm trials \n  (b) Scenario 1: \n Constant effect"),
                                  plot2+
                                    theme(legend.position="none",
                                          plot.title = element_text(size=10, 
                                                                    face="bold",
                                                                    lineheight = 1.2))+
                                    labs(title = "\n (c) Scenario 2: \n Waning effect"),
                                  plot3+
                                    theme(legend.position="none",
                                          plot.title = element_text(size=10, 
                                                                    face="bold",
                                                                    lineheight = 1.2))+
                                    labs(title = "\n (d) Scenario 3: Delayed then \n waning effect"),
                                  align = "h",
                                  rel_widths=c(1,1,1),
                                  nrow = 1),
                        plot_legend,
                        rel_widths=c(1, 0.75),
                        axis = "l"
  )
  
  plot_all
  
  tiff(file = paste0("plots/two_arm/forest_", irmst_estimand, "_knots.tiff"),   
       width = 7.2, 
       height = 4.2,
       units = 'in',  
       res = 300, 
       compression = "lzw")
  print(plot_all)
  dev.off()
  
  
  print(paste0("plots/two_arm/forest_", irmst_estimand, "_knots.tiff"))

  
  plot_single <- readRDS(paste0("plots/single_arm/forest_",substr(irmst_estimand, start = 2, stop = 6),"_trt0_single_arm_plot_knots.rds"))
 
  plot_single
  
  plot_single_legend <-  cowplot::get_legend(
    # create some space to the left of the legend
    plot_single + theme(legend.title = element_text(margin = margin(3, 0, 3, 0)),
                            legend.box.margin = margin(0, 20, 0, 12),
                        #    legend.text = element_text(size=7),
                            legend.key.spacing.y = unit(1, "pt"))
  )
  
  plot_single_without_legend <- plot_single+
    theme(legend.position="none",
          plot.title = element_text(size=10, face="bold",
                                    lineheight = 1.2))+
    labs(title = "Single-arm trial \n  (a) Weibull mixture")
  
  plot_single_mod <- plot_grid(NULL, plot_single_without_legend, NULL, plot_single_legend, NULL, rel_widths = c(0.3, 0.9,-0.02,1.1, 0.1), nrow = 1)
  plot_single_mod  
  
  
#  plot_single
 
  plot_both <- plot_grid(plot_single_mod,plot_all, ncol = 1, rel_heights = c(1.7, 1.4))
  
  tiff(file = paste0("plots/single_and_two_arm/forest_", irmst_estimand, "_all_knots.tiff"),   
       width = 7.2, 
       height = 7.3,
       units = 'in',  
       res = 300, 
       compression = "lzw")
  print(plot_both)
  dev.off()
  
  
  
}



