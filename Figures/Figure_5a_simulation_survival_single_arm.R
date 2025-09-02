#########################################################
# Figure 5a, Simulation study,
# plot survival functions.
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
library(data.table)
library(forcats)
library(ggh4x)

# Jobname where results are stored.
stores_res <- "directory/to/store/simulations"
setwd(store_res)
#setwd("/scratch/klvq491/simsurvextrap_slurm_mix_weib_full1")

# Read in scenarios data
scenarios <- readRDS("scenarios.rds")
# Filter scenarios to plot.
scenarios <- scenarios %>%
  filter(design_id == "single_arm",
         stan_fit_method == "mcmc",
         mspline_df == 10,
         data_cut_model_id == "data_cut_mod2",
         add_knots %in% c("none", "extra_knots1"),
         backhaz == T)

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
# Model labels.
################################################

new_model <-  scenarios %>%
  mutate(add_knots = if_else(is.na(add_knots), "default", add_knots)) %>%
  mutate(include_external_data = if_else(include_external_data, "yes", "no")) %>%
  mutate(bsmooth = if_else(bsmooth == T, "Smoothed", "Standard")) %>%
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

new_model_levels <- unique(new_model$new_model_id)
new_model_labels <- unique(new_model$new_model_id_labels)

###########################################################
# Combine true and estimated survival.
###########################################################

dgm_true <- readRDS("dgm_true.rds")

for(i in 1:nrow(scenarios)){
  
  temp <- readRDS(paste0(scenarios$scenario_fit_id[i], "/all_res.rds"))
  
  temp <- temp %>%
    filter(estimand == "survival",
           isim <= 50)
  
  if(i == 1){
    res <- temp
  }  else  {
    res <- rbindlist(list(res, temp))
  }
  print(paste0("file ", i, "/", nrow(scenarios)))
  res <- as_tibble(res)
}


# Extract scenarios to plot.
scen_df <- scenarios %>%
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

surv_plot <-  scen_df %>%  
  left_join(res, by = "scenario_fit_id")  %>%
  mutate("Scenarios" = factor(new_model_id, 
                              levels = new_model_levels,
                              labels = new_model_labels)) %>%
  mutate("External data" = factor(external_bias_model_id,
                                  levels = external_data_models_index,
                                  labels = external_data_labels_index)) %>%
  mutate("Extra knots" = factor(add_knots,
                                levels = extra_knots_models_index,
                                labels = extra_knots_labels_index)) %>%
  filter(external_bias_model_id %in%  c("none", 
                                        "external_bias_mod1", "external_bias_mod3",
                                        "external_bias_mod5")) %>%
  arrange(add_knots) %>%
  arrange(external_bias_model_id) %>%
  arrange(backhaz) %>%
  arrange(include_external_data) %>%
  ungroup() %>%
  filter(isim <= 25) %>%
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
        strip.text.x = element_text(size = 6, 
                                    margin = margin(0.085,0,0.085,0, "cm"))) +
  geom_line()+
  scale_linewidth(range = c(0.4,1))+
  scale_alpha(range = c(0.15,1))+
  scale_color_manual(values = c("#830051", "black")) + 
  scale_x_continuous("Time (years)")+
  scale_y_continuous("Survival",limits = c(0,1), labels = scales::percent)+
  facet_manual(~Scenarios, design = c(
    "#AABB#
        CCDDEE"
  ))+
  guides(linewidth = "none", alpha = "none", colour = "none") 

saveRDS(surv_plot, "plots/Figure_5a.rds")
