#########################################################
# Sup Figure: Hazard ratio plots.
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
library(ggh4x)

# Jobname where results are stored.
stores_res <- "directory/to/store/simulations"
setwd("/scratch/klvq491/simsurvextrap_slurm_mix_weib_full1")
setwd(store_res)

# Load in data on scenarios and performance metrics.
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
scenarios <- readRDS("scenarios.rds")

new_model <-  scenarios %>%
  mutate(add_knots = if_else(is.na(add_knots), "default", add_knots)) %>%
  mutate(include_external_data = if_else(include_external_data, "yes", "no")) %>%
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



###########################################################
# Combine true and estimated hazards.
###########################################################


for(j in 1:3){
  
  scen_df <- scenarios %>%
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
    # mutate(add_knots %in% c("default", "extra_knots5")) %>%
    filter(weibull_model_id == "weibull_mod1") %>%
    # filter based on analysis 1 .
    filter(trt_effect_model_id == paste0("trt_effect_mod",j),
           design_id == "two_arm", 
           stan_fit_method == "mcmc",
           mspline_df == 10,
           data_cut_model_id == "data_cut_mod2",
           add_knots %in% c("none", "extra_knots1"),
           external_bias_model_id %in%  c("none", 
                                          "external_bias_mod1", "external_bias_mod3",
                                          "external_bias_mod5"),
           waning_model_id == "waning_mod1",
           backhaz == T) 
  
  dgm_true <- readRDS("dgm_true.rds")
  
  for(i in 1:nrow(scen_df)){
    #i <- 200
    temp <- readRDS(paste0(scen_df$scenario_fit_id[i], "/all_res.rds"))
    #summary(as.factor(temp$estimand))
    temp <- temp %>%
      filter(estimand == "hr",
             isim <= 50)
    
    if(i == 1){
      res <- temp
    }  else  {
      res <- rbindlist(list(res, temp))
    }
    print(paste0("file ", i, "/", nrow(scen_df)))
    res <- as_tibble(res)
  }
  
  
  #i <- 1
  # Extract scenarios to plot.
  scen_df2 <- scen_df  %>%
    mutate(new_model_id = paste0(external_bias_model_id, "_", add_knots)) %>%
    mutate("Scenarios" = factor(new_model_id, 
                                levels = new_model_levels,
                                labels = new_model_labels)) %>%
    arrange(Model) %>%
    arrange(Scenarios)
  
  hr_plot <- 
    scen_df2 %>%  
    left_join(res, by = "scenario_fit_id")  %>%
    
    arrange(add_knots) %>%
    arrange(external_bias_model_id) %>%
    arrange(backhaz) %>%
    arrange(include_external_data) %>%
    ungroup() %>%
    filter(isim <= 25) %>%
    mutate(line_alpha = if_else(isim  == 0, 1, 0)) %>%
    mutate(line_width = if_else(isim  == 0, 1, 0)) %>%
    mutate(line_colour = if_else(isim  == 0, 1, 0),
           line_colour = as.factor(line_colour)) %>%
    #summary(as.factor(haz_plot$line_alpha))  
    #summary(as.factor(haz_plot$isim))  
    filter(t > 1e-2) %>%
    ggplot(aes(x = t, y = value, alpha = line_alpha, colour = line_colour, 
               group = isim, linewidth = line_width))+
    theme_bw()+
    theme(axis.text.y=element_text(size = 6),
          axis.title.y=element_text(size = 8),
          axis.text.x = element_text(size = 6),
          axis.title.x = element_text(size = 8),
          strip.text.x = element_text(size = 6, 
                                      margin = margin(0.085,0,0.085,0, "cm")),
          strip.text.y = element_text(size = 6, 
                                      margin = margin(0.085,0.2,0.085,0.2, "cm"),
                                      angle  = 0),
          strip.background = element_rect(fill = "white", 
                                          colour = "black", linewidth = rel(2))) +
    geom_line()+
    geom_vline(xintercept = 5, colour = "gray30", linetype = "dashed", alpha = 0.6) +
    geom_hline(yintercept = 1, linewidth = 1, 
                  alpha = 0.55, colour = "gray40",linetype = "solid")+
    scale_linewidth(range = c(0.4,1))+
    scale_alpha(range = c(0.15,1))+
    scale_color_manual(values = c("#830051", "black")) + 
    scale_x_continuous("Time (years)")+
    #scale_y_continuous("Survival",limits = c(0,1), labels = scales::percent)+
    scale_y_continuous("Hazard ratio", breaks = c(0.25, 0.5, 1, 2, 4)) +
    coord_trans(y = "log10", ylim=c(0.2, 4))+
   #facet_wrap(Model~Scenarios, nrow = 5)+
    facet_grid(Scenarios~Model)+
    guides(linewidth = "none", alpha = "none", colour = "none") 
  
  
  hr_plot
  saveRDS(hr_plot, paste0("plots/two_arm/hr_trt_effect_mod", j, ".rds")) 
  

}

plot1 <- readRDS(paste0("plots/two_arm/hr_trt_effect_mod", 1, ".rds"))
plot2 <- readRDS(paste0("plots/two_arm/hr_trt_effect_mod", 2, ".rds"))
plot3 <- readRDS(paste0("plots/two_arm/hr_trt_effect_mod", 3, ".rds"))

# 
for(i in 1:3){
tiff(file = paste0("plots/two_arm/hr_trt_effect_mod", i, ".tiff"),   
     width = 7.5, 
     height = 5.5,
     units = 'in',  
     res = 300, 
     compression = "lzw")
print(get(paste0("plot",i)))
dev.off()
}

# Attempt to combine all.

plot_all <- plot_grid(plot1+
                        theme(legend.position="none",
                              plot.title = element_text(size=10, face="bold"))+
                        labs(title = "(a) Scenario 1: \n Constant effect"),
                      plot2+
                        theme(legend.position="none",
                              plot.title = element_text(size=10, face="bold"))+
                        labs(title = "(b) Scenario 2: \n Waning effect"),
                      plot3+
                        theme(legend.position="none",
                              plot.title = element_text(size=10, face="bold"))+
                        labs(title = "(c) Scenario 3: Delayed then \n waning effect"),
                      align = "v",
                      rel_heights=c(1,1,1),
                      nrow = 3)

plot_all

tiff(file = "plots/two_arm/hr_all.tiff",   
     width = 6.5, 
     height = 8,
     units = 'in',  
     res = 300, 
     compression = "lzw")
print(plot_all)
dev.off()

