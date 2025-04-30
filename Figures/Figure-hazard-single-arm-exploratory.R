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
library(data.table)

jobname <- "mix_weib_test19"
user <- Sys.info()["user"]
project_directory <- paste0("/projects/aa/", user, "/")
store_res <- paste0(project_directory, "simsurvextrap_slurm_", jobname, "/")
setwd(store_res)

scenarios <- readRDS("scenarios.rds")
scenarios <- scenarios %>%
  filter(backhaz == T) %>%
  filter(design_id == "single_arm") %>%
  filter(mspline_df == 10) %>%
  filter(stan_fit_method == "mcmc",
         data_cut_model_id == "data_cut_mod2") %>%
  filter(external_id == "none" | add_knots %in% paste0("extra_knots",1:7)) %>%
  filter(external_id != "none" | add_knots == "none") 


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
  mutate(haz_bias = if_else(haz_bias_temp > 0, paste0("External~data~with~`+`*`",haz_bias_temp, "`*`%`~`bias`"), haz_bias))  %>%
  mutate(haz_bias = if_else(haz_bias_temp < 0, paste0("External~data~with~`\u2212`*`",abs(haz_bias_temp), "`*`%`~`bias`"), haz_bias))  %>%
  mutate(haz_bias = if_else(haz_bias == "0", "Unbiased~external~data", haz_bias)) %>%
  select(-haz_bias_temp)  %>%
  rename(external_data_label = haz_bias) %>%
  # arrange(abs(loghaz_bias)) %>%
  select(external_bias_model_id, external_data_label) %>%
  add_row(external_bias_model_id = "none", external_data_label = "`No`~`external`~`data`",
          .before = 1)# %>%

external_data_models_index <- external_data_models_labels$external_bias_model_id  
external_data_labels_index <- external_data_models_labels$external_data_label


################################################
# Extra knots labels.
################################################


extra_knots_settings <- readRDS("extra_knots_settings.rds")
extra_knots_models <- readRDS("extra_knots_models.rds")
extra_knots_models_index <- extra_knots_models$extra_knots_id
extra_knots_labels_index <- extra_knots_models$extra_knots_labels

estimand_labels <- readRDS("estimand_labels.rds")

##################################################
# Prep hazard times.
##################################################

dgm_true <- readRDS("dgm_true.rds")
#temp <- readRDS("dgm_true_mod3/survival_true.rds")
#head(temp)

###########################################################
# Combine true and estimated hazards.
###########################################################

for(i in 1:nrow(scenarios)){
  temp <- readRDS(paste0(scenarios$scenario_fit_id[i], "/all_res.rds"))
  #summary(as.factor(temp$isim))
  #sum(temp$isim == 0)
  temp <- temp %>%
    filter(estimand == "hazard",
           isim <= 50)

  if(i == 1){
    res <- temp
  }  else  {
    res <- rbindlist(list(res, temp))
  }
  print(paste0("file ", i, "/", nrow(scenarios)))
  res <- as_tibble(res)
}




for(i in 1:2){
  #i <- 1
  scen_df <- scenarios %>%
    mutate(add_knots = if_else(add_knots == "none", "default", add_knots)) %>%
    mutate(include_external_data = if_else(include_external_data, "yes", "no")) %>%
    mutate(bsmooth = if_else(bsmooth == T, "Smoothed basis", "Standard basis")) %>%
    mutate("External data" = factor(external_bias_model_id,
                                    levels = external_data_models_index,
                                    labels = external_data_labels_index)) %>%
    mutate("Extra knots" = factor(add_knots,
                                  levels = extra_knots_models_index,
                                  labels = extra_knots_labels_index)) %>%
    mutate("GPM rates" = if_else(backhaz, "GPM included", "GPM not included")) %>%
    arrange(`GPM rates`) %>%
    mutate(external_data = 
             if_else(external_bias_model_id == "none", 
                     "No external data", "With external data"))
  
  haz_plot <- scen_df %>%  
    filter(weibull_model_id == paste0("weibull_mod",i)) %>%
    left_join(res, by = "scenario_fit_id") %>%
    arrange(add_knots) %>%
    arrange(external_bias_model_id) %>%
    arrange(backhaz) %>%
    arrange(include_external_data) %>%
    ungroup() %>%
    filter(isim <= 20) %>%
    mutate(line_alpha = if_else(isim  == 0, 1, 0)) %>%
    mutate(line_width = if_else(isim  == 0, 1, 0)) %>%
    mutate(line_colour = if_else(isim  == 0, 1, 0),
           line_colour = as.factor(line_colour)) %>%
    #summary(as.factor(haz_plot$line_alpha))  
    #summary(as.factor(haz_plot$isim))  
    filter(t > 1e-2) %>%
    ggplot(aes(x = t, y = value, alpha = line_alpha, colour = line_colour, 
               group = isim, linewidth = line_width))+
    theme_classic()+
    theme(axis.text.y=element_text(size = 6),
          axis.title.y=element_text(size = 8),
          axis.text.x = element_text(size = 6),
          axis.title.x = element_text(size = 8),
          strip.text.x = element_text(size = 6)) +
    #  geom_vline(xintercept = c(5), colour = c("gray60"))+
    #  geom_vline(xintercept = c(25), colour = c("gray70"))+
    geom_line()+
    scale_linewidth(range = c(0.4,1))+
    scale_alpha(range = c(0.15,1))+
    scale_color_manual(values = c("#830051", "black")) + 
    scale_x_continuous("Time (years)")+
    scale_y_continuous("Hazard",limits = c(0,0.45))+
    facet_wrap(~`External data`, ncol = 3, labeller = label_parsed)+
    guides(linewidth = "none", alpha = "none", colour = "none") 
  
  haz_plot
  saveRDS(haz_plot, paste0("plots/single_arm/hazard_plot",i,".rds"))
  
}

####################################################
# Combine both DGMs.
####################################################

scen_df <- scenarios %>%
  mutate(DGM = if_else(weibull_model_id == "weibull_mod1", "Weibull~mixture", "Weibull~cure")) %>%
  mutate(DGM = fct_inorder(DGM)) %>%
  mutate(add_knots = if_else(add_knots == "none", "default", add_knots)) %>%
  mutate(include_external_data = if_else(include_external_data, "yes", "no")) %>%
  mutate(bsmooth = if_else(bsmooth == T, "Smoothed basis", "Standard basis")) %>%
  mutate("External data" = factor(external_bias_model_id,
                                  levels = external_data_models_index,
                                  labels = external_data_labels_index)) %>%
  mutate("Extra knots" = factor(add_knots,
                                levels = extra_knots_models_index,
                                labels = extra_knots_labels_index)) %>%
  mutate("GPM rates" = if_else(backhaz, "GPM included", "GPM not included")) %>%
  arrange(`GPM rates`) %>%
  mutate(external_data = 
           if_else(external_bias_model_id == "none", 
                   "No external data", "With external data"))

haz_plot <- scen_df %>%  
  left_join(res, by = "scenario_fit_id") %>%
  arrange(add_knots) %>%
  arrange(external_bias_model_id) %>%
  arrange(backhaz) %>%
  arrange(include_external_data) %>%
  ungroup() %>%
  filter(isim <= 50) %>%
  mutate(line_alpha = if_else(isim  == 0, 1, 0)) %>%
  mutate(line_width = if_else(isim  == 0, 1, 0)) %>%
  mutate(line_colour = if_else(weibull_model_id == "weibull_mod1", 1, 2)) %>%
  mutate(line_colour = line_colour*(isim  != 0)) %>%
  mutate(line_colour = as.factor(line_colour)) %>%
  #summary(as.factor(haz_plot$line_alpha))  
  #summary(as.factor(haz_plot$isim))  
  filter(t > 1e-2) %>%
  ggplot(aes(x = t, y = value, alpha = line_alpha, colour = line_colour, 
             group = isim, linewidth = line_width))+
  theme_classic()+
  theme(axis.text.y=element_text(size = 6),
        axis.title.y=element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.title.x = element_text(size = 8),
        strip.text.x = element_text(size = 5, 
                                    margin = margin(0.1,0,0.1,0, "cm"))) +
  #  geom_vline(xintercept = c(5), colour = c("gray60"))+
  #  geom_vline(xintercept = c(25), colour = c("gray70"))+
  geom_line()+
  scale_linewidth(range = c(0.4,1))+
  scale_alpha(range = c(0.075,1))+
  scale_color_manual(values = c("black", "#830051", "#00441b")) + 
  scale_x_continuous("Time (years)")+
  scale_y_continuous("Hazard",limits = c(0,0.45))+
  facet_wrap(DGM~`External data`, ncol = 3, labeller = label_parsed)+
  guides(linewidth = "none", alpha = "none", colour = "none") 

haz_plot
saveRDS(haz_plot, paste0("plots/single_arm/hazard_plot_all.rds"))



