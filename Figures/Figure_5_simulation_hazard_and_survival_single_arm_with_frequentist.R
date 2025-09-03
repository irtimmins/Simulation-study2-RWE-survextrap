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
setwd("/scratch/klvq491/simsurvextrap_slurm_mix_weib_full1")

###########################################################
# Frequentist models: combine true and estimated survival.
###########################################################

# Read in scenarios data for frequentist (exponential, Weibull etc models)
scenarios_freq <- readRDS("scenarios_single_arm_freq.rds")
dgm_true <- readRDS("dgm_true.rds")

for(i in 1:nrow(scenarios_freq)){

  temp <- readRDS(paste0(scenarios_freq$single_arm_freq_model_id[i], "/all_res.rds"))
  
  temp <- temp %>%
    filter(estimand %in% c("survival", "hazard"),
           isim <= 50)
  
  if(i == 1){
    res_freq <- temp
  }  else  {
    res_freq <- rbindlist(list(res_freq, temp))
  }
  print(paste0("file ", i, "/", nrow(scenarios_freq)))
  res_freq <- as_tibble(res_freq)
}


###########################################################
# survextrap models: combine true and estimated survival.
###########################################################

# Read in scenarios data for survextrap models
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
# Survextrap models: combine true and estimated survival.
###########################################################

dgm_true <- readRDS("dgm_true.rds")

for(i in 1:nrow(scenarios)){
  
  temp <- readRDS(paste0(scenarios$scenario_fit_id[i], "/all_res.rds"))
  
  temp <- temp %>%
    filter(estimand %in% c("survival", "hazard"),
           isim <= 50)
  
  if(i == 1){
    res_survextrap <- temp
  }  else  {
    res_survextrap <- rbindlist(list(res_survextrap, temp))
  }
  print(paste0("file ", i, "/", nrow(scenarios)))
  res_survextrap <- as_tibble(res_survextrap)
}


# Extract survextrap scenarios to plot.
scen_survextrap_df <- scenarios %>%
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
  slice(c(2,3,5,4,1,6,7)) %>%
  filter(external_bias_model_id %in%  c("none", 
                                        "external_bias_mod1", "external_bias_mod3",
                                        "external_bias_mod5")) %>%
  select(scenario_fit_id,new_model_id) %>%
  left_join(res_survextrap, by = "scenario_fit_id") %>%
  rename(model_id = new_model_id)
  
#scen_survextrap_dfnew_model_idscen_survextrap_df$new_model_id_labels

model_levels <- c("exp",
                  "weibull", "spline", 
                  "none_none", 
                  "none_extra_knots1", 
                  "external_bias_mod3_extra_knots1", 
                  "external_bias_mod1_extra_knots1", 
                  "external_bias_mod5_extra_knots1")

model_labels <- c( "exponential model", 
                   "Weibull model",
                   "Royston-Parmar (df = 3) model",
                   "survextrap model, no external data,\nno extra knots", 
                   "survextrap model, no external data,\nwith extra knots at t=5,10,25",
                   "survextrap model, \nwith external data with +20% bias", 
                   "survextrap model, \nwith unbiased external data" , 
                   "survextrap model, \nwith external data with −20% bias" )

scen_freq_df <- scenarios_freq %>%
  select(single_arm_freq_model_id, parametric_model) %>%
  rename(scenario_fit_id = single_arm_freq_model_id) %>%
  left_join(res_freq, by = "scenario_fit_id") %>%
  rename(model_id = parametric_model)

scen_all_df <- bind_rows(scen_freq_df, 
                         scen_survextrap_df) 
  
surv_plot <-  
  scen_all_df %>%
  filter(estimand == "survival") %>%
  mutate("Scenarios" = factor(model_id, 
                              levels = model_levels ,
                              labels = model_labels)) %>%
  ungroup() %>%
  filter(isim <= 25) %>%
  mutate(line_alpha = if_else(isim  == 0, 1, 0)) %>%
  mutate(line_width = if_else(isim  == 0, 1, 0)) %>%
  mutate(line_colour = if_else(isim  == 0, 1, 0),
         line_colour = as.factor(line_colour)) %>%
  filter(t > 1e-2) %>%
  ggplot(aes(x = t, y = value, alpha = line_alpha, colour = line_colour, 
             group = isim, linewidth = line_width))+
  theme_classic()+
  theme(axis.text.y=element_text(size = 6),
        axis.title.y=element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.title.x = element_text(size = 8),
        strip.text.x = element_text(size = 5, 
                                    margin = margin(0.065,0,0.065,0, "cm"))) +
  geom_vline(xintercept = c(5), alpha = 0.6,
             colour = c("gray30"), linewidth = 0.8, linetype = "dashed")+
  geom_line()+
  scale_linewidth(range = c(0.4,1))+
  scale_alpha(range = c(0.15,1))+
  scale_color_manual(values = c("#830051", "black")) + 
  scale_x_continuous("Time (years)")+
  scale_y_continuous("Survival", 
                     limits = c(0,1),
                     labels = scales::percent)+
  facet_manual(~Scenarios, design = c("AABBCC 
                                      #DDEE# 
                                      FFGGHH"))+
  guides(linewidth = "none", alpha = "none", colour = "none") 

surv_plot

saveRDS(surv_plot, "plots/Figure_5a_with_frequentist.rds")





# Plot hazard functions.
haz_plot <- scen_all_df %>%
  filter(estimand == "hazard") %>%
  mutate("Scenarios" = factor(model_id, 
                              levels = model_levels ,
                              labels = model_labels)) %>%
  ungroup() %>%
  filter(isim <= 25) %>%
  mutate(line_alpha = if_else(isim  == 0, 1, 0)) %>%
  mutate(line_width = if_else(isim  == 0, 1, 0)) %>%
  mutate(line_colour = if_else(isim  == 0, 1, 0),
         line_colour = as.factor(line_colour)) %>%
  filter(t > 1e-2) %>%
  ggplot(aes(x = t, y = value, alpha = line_alpha, colour = line_colour, 
             group = isim, linewidth = line_width))+
  theme_classic()+
  theme(axis.text.y=element_text(size = 6),
        axis.title.y=element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.title.x = element_text(size = 8),
        strip.text.x = element_text(size = 5, 
                                    margin = margin(0.065,0,0.065,0, "cm"))) +
  geom_vline(xintercept = c(5), alpha = 0.6,
             colour = c("gray30"), linewidth = 0.8, linetype = "dashed")+
  geom_line()+
  scale_linewidth(range = c(0.4,1))+
  scale_alpha(range = c(0.15,1))+
  scale_color_manual(values = c("#830051", "black")) + 
  scale_x_continuous("Time (years)")+
  scale_y_continuous("Hazard",limits = c(0,0.5))+
  facet_manual(~Scenarios, design = c("AABBCC 
                                      #DDEE# 
                                      FFGGHH"))+
  guides(linewidth = "none", alpha = "none", colour = "none") 

haz_plot
saveRDS(haz_plot, "plots/Figure_5b_with_frequentist.rds")

# Combine survival and hazard plots using plot_grid.
plot_all <- plot_grid(surv_plot +
                        theme(  plot.title = element_text(hjust = -0.1,
                                                          size=8, face="bold"))+
                        labs(title = "(a) Survival curves"),
                      haz_plot+
                        theme(  plot.title = element_text(hjust = -0.1,
                                                          size=8, face="bold"))+
                        labs(title = "(b) Hazard curves"),
                      align = "v",
                      rel_heights=c(0.5,0.5),
                      ncol = 1)

tiff(file = "plots/Figure_5_with_frequentist.tiff",   
     width = 5.0, 
     height = 6.5,
     units = 'in',  
     res = 300, 
     compression = "lzw")
print(plot_all)
dev.off()

