################
setwd("~/survival_extrapolation/simsurvextrap")

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
library(survminer)
library(RColorBrewer)

jobname <- "mix_weib_full1"
user <- Sys.info()["user"]
project_directory <- paste0("/projects/aa/", user, "/")
store_res <- paste0(project_directory, "simsurvextrap_slurm_", jobname, "/")

source("R/simulate_dgm_mixture_weibull.R")

setwd(store_res)

#####################################################
# Plot options
##################################################### 

margins <- unit(c(0,0,0,0), "cm")
legend_text <- element_text(size=8)
legend_key <- unit(0.3, "lines")
legend_space_y <- unit(1.5, "pt")


#####################################################
# Single arm DGMs. 
##################################################### 

dgm_true <- readRDS("dgm_true.rds")

scenarios <- readRDS("scenarios.rds")

scenarios <- scenarios %>%
  filter(design_id == "single_arm",
         stan_fit_method == "mcmc",
         mspline_df == 10,
         data_cut_model_id == "data_cut_mod2",
         add_knots %in% c("none"),
         backhaz == T)


for(i in 1:nrow(scenarios)){
  
  temp <- readRDS(paste0(scenarios$scenario_fit_id[i], "/all_res.rds"))
  
  temp <- temp %>%
    filter(isim == 0)
  
  if(i == 1){
    res <- temp
  }  else  {
    res <- rbindlist(list(res, temp))
  }
  print(paste0("file ", i, "/", nrow(scenarios)))
  res <- as_tibble(res)
}


##################################################
# Prep cetuximab control Kaplan-Meier.
##################################################

trial_data <- cetux %>%
  as_tibble() %>%
  mutate(trt = if_else(treat == "Cetuximab", 1, 0)) %>%
  filter(trt == 0) %>%
  select(-months) %>%
  rename(event = d,
         time = years)

km_fit <- survfit(Surv(time, event) ~ 1, data=trial_data)
#plot(km_fit)

km_plot <- ggsurvplot(km_fit, data=trial_data)
km_surv_df <- km_plot$data.survplot 
max_trial_data <- max(km_surv_df$time)


##################################################
# Survival plot.
##################################################

trial_data_km <- km_surv_df %>% 
  mutate(Dataset = "Bonner trial Kaplan-Meier") %>%
  select(time, surv, Dataset) %>%
  rename(t = time, value = surv)

margins1 <- unit(c(0.1,0.1,0,0), "cm")
margins2 <- unit(c(0,0,0,0), "cm")

surv_plot <- scenarios %>%  
  filter(weibull_model_id == paste0("weibull_mod",1)) %>%
  left_join(res, by = "scenario_fit_id") %>%
  mutate(Dataset = if_else(weibull_model_id == "weibull_mod1", 
                           "Weibull mixture DGM",
                           "Weibull cure DGM")) %>%
  filter(estimand == "survival") %>%
  select(t, value, Dataset) %>%
  bind_rows(trial_data_km) %>%
  ungroup() %>%
  filter(t > 1e-2) %>%
  ggplot(aes(x = t, y = value,
             colour = Dataset, alpha = Dataset, linewidth = Dataset))+
  theme_classic()+
  theme(axis.text.y=element_text(size = 6),
        axis.title.y=element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.title.x = element_text(size = 8),
        strip.text.x = element_text(size = 4),
        legend.position = "inside",
        legend.position.inside  = c(.6,0.8),
        legend.text=legend_text,
        legend.key.size=legend_key,
        legend.key.spacing.y = legend_space_y,
        plot.margin = get(paste0("margins",i)))+
  geom_vline(xintercept = c(5), colour = c("gray30"), linewidth = 0.8, linetype = "dashed")+
  #  geom_vline(xintercept = c(25), colour = c("gray70"))+
  scale_colour_manual("", values = c("gray20", "firebrick")) + 
  scale_alpha_manual("",values = c(0.9, 0.6))+
  scale_linewidth_manual("",values = c(0.9, 0.5))+
  geom_line(linewidth = 0.9)+
  scale_x_continuous("Time (years)")+
  scale_y_continuous("Survival", limits = c(0,1), labels = scales::percent)#+

# 
surv_plot
# 
all_haz_plot <- 
  scenarios %>%  
  left_join(res, by = "scenario_fit_id") %>%
  filter(estimand %in% c("hazard","excess_hazard", "gpm_hazard")) %>%
  filter(weibull_model_id == paste0("weibull_mod",1)) %>%
  mutate(estimand = factor(estimand, labels = c("Hazard", "Excess hazard", "GPM rate"),
                           levels = c("hazard","excess_hazard", "gpm_hazard"))) %>%
  filter(estimand != "Excess hazard") %>%
  filter(t > 1e-2) %>%
  ggplot(aes(x = t, y = value, colour = estimand, linetype = estimand))+
  theme_classic()+
  theme(axis.text.y=element_text(size = 6),
        axis.title.y=element_text(size = 8, vjust = -2),
        axis.text.x = element_text(size = 6),
        axis.title.x = element_text(size = 8),
        strip.text.x = element_text(size = 4),
        legend.position = "inside",
        legend.position.inside  = c(.65,.82),
        legend.text=legend_text,
        legend.key.size=legend_key,
        legend.key.spacing.y = legend_space_y,
        plot.margin = margins) +
  geom_vline(xintercept = c(5), colour = c("gray30"), linewidth = 0.8, linetype = "dashed")+
  #  geom_vline(xintercept = c(25), colour = c("gray70"))+
  geom_line(linewidth = 0.9)+
  scale_colour_manual("", values = c("firebrick", "green3")) + 
  scale_linetype_manual("",values = c("solid", "dashed")) +
  scale_x_continuous("Time (years)")+
  scale_y_continuous("Hazard", limits = c(0,0.4))

all_haz_plot

title_vec <- c("Single-arm trials\n(a) Weibull mixture")

plot_all <- plot_grid(
  surv_plot,
  all_haz_plot,
  rel_widths=c(0.5,0.5),
  ncol = 2,
  align = "hv")+
  theme(plot.title = element_text(hjust = 0.03,
                                  size=8, face="bold",
                                  lineheight = 1.2))+
  labs(title = title_vec[1])

single_arm <- plot_all


#####################################################
# Two arm DGMs. 
##################################################### 

dgm_true <- readRDS("dgm_true.rds")
scenarios2 <- readRDS("scenarios.rds")


# Select three DGM scenarios.
scenarios2 <- scenarios2 %>%
  filter(design_id == "two_arm",
         stan_fit_method == "mcmc",
         mspline_df == 10,
         data_cut_model_id == "data_cut_mod2",
         add_knots %in% c("none"),
         model == "separate",
         backhaz == T)

for(i in 1:nrow(scenarios2)){
  
  temp <- readRDS(paste0(scenarios2$scenario_fit_id[i], "/all_res.rds"))
  
  temp <- temp %>%
    filter(isim == 0)
  
  if(i == 1){
    res2 <- temp
  }  else  {
    res2 <- rbindlist(list(res2, temp))
  }
  print(paste0("file ", i, "/", nrow(scenarios2)))
  res2 <- as_tibble(res2)
}

##############################################
margins1 <- unit(c(0.1,0.1,0,0), "cm")
margins2 <- unit(c(0,0.1,0,0), "cm")
margins3 <- unit(c(0,0.1,0,0), "cm")

for(i in 1:3){
  #i <- 1
  surv_plot_two_arms <- scenarios2 %>%  
    filter(trt_effect_model_id == paste0("trt_effect_mod",i)) %>%
    left_join(res2, by = "scenario_fit_id") %>%
    filter(estimand == "survival") %>%
    mutate(trt = factor(trt, levels = c(1,0), labels = c("Active","Control"))) %>%
    ungroup() %>%
    filter(t > 1e-2) %>%
    ggplot(aes(x = t, y = value, colour = trt))+
    theme_classic()+
    theme(axis.text.y=element_text(size = 6),
          axis.title.y=element_text(size = 8),
          axis.text.x = element_text(size = 6),
          axis.title.x = element_text(size = 8),
          strip.text.x = element_text(size = 6),
          legend.position = "inside",
          legend.position.inside  = c(.6,0.6),
          legend.title=element_blank(),
          legend.text=legend_text,
          legend.key.size=legend_key,
          legend.key.spacing.y = legend_space_y,
          plot.margin = get(paste0("margins",i)))+
    geom_vline(xintercept = c(5), colour = c("gray30"), linewidth = 0.8, linetype = "dashed")+
    # geom_vline(xintercept = c(25), colour = c("gray70"))+
    geom_line(linewidth =  1)+
    scale_color_manual(values = c("firebrick","steelblue")) + 
    scale_x_continuous("Time (years)")+
    scale_y_continuous("Survival", limits = c(0,1), labels = scales::percent)#+

  
 hr_plot_two_arms <- 
    scenarios2 %>%  
    filter(trt_effect_model_id == paste0("trt_effect_mod",i)) %>%
    left_join(res2, by = "scenario_fit_id") %>%
    filter(estimand %in% c("hr", "excess_hazard")) %>%
    pivot_wider(names_from = "estimand", values_from = value) %>%
    select(c(t, trt, isim, excess_hazard, hr)) %>%
    pivot_wider(names_from = trt, values_from = c("hr", "excess_hazard")) %>%
    mutate(excess_hr = excess_hazard_1/excess_hazard_0) %>%
    rename(hr = hr_NA) %>%
    pivot_longer(cols = c(hr, excess_hr), names_to = "estimand")  %>%
    mutate(estimand = factor(estimand, levels = c("hr", "excess_hr"),
                             labels = c("Hazard ratio", "Excess hazard ratio"))) %>%
    filter(t > 1e-2) %>%
    ggplot(aes(x = t, y = value, colour = estimand, linetype = estimand))+
    theme_classic()+
    theme(axis.text.y=element_text(size = 6),
          axis.title.y=element_text(size = 8, vjust = -2),
          axis.text.x = element_text(size = 6),
          axis.title.x = element_text(size = 8),
          strip.text.x = element_text(size = 6),
          legend.position = "inside",
          legend.position.inside  = c(.65,0.83),
          legend.title=element_blank(),
          legend.text=legend_text,
          legend.key.size=legend_key,
          legend.key.spacing.y = legend_space_y,
          plot.margin = margins)+
    geom_vline(xintercept = c(5), colour = c("gray30"), linewidth = 0.8, linetype = "dashed")+
    # geom_vline(xintercept = c(25), colour = c("gray70"))+
    geom_hline(yintercept = 1, linewidth = 1, 
               alpha = 0.55, colour = "gray40",linetype = "solid")+
    geom_line(linewidth = 1.0)+
    scale_color_manual(values = brewer.pal(n = 12, name = "Paired")[c(10,8)]) + 
    scale_linetype_manual(values = c("solid", "dashed"))+
    scale_x_continuous("Time (years)")+
    scale_y_continuous("Hazard ratio", limits = c(0.4, 1.45), 
                       breaks = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4))

  
  hr_plot_two_arms
  assign(paste0("surv_fig",i), surv_plot_two_arms)
  assign(paste0("hr_fig",i), hr_plot_two_arms)
  
  title_vec <- c("Two-arm trials\n(b) Scenario 1: Constant effect on excess hazard",
                 "(c) Scenario 2: Waning effect",
                 "(d) Scenario 3: Delayed then waning effect")
  
  hjust_vec <- c(0.03, 0.03, 0.04)
  
  plot_all_two_arms <- plot_grid(
    surv_plot_two_arms ,
    hr_plot_two_arms,
    align = "hv",
    rel_widths=c(0.5,0.5),
    ncol = 2)+ 
    theme(plot.title = element_text(hjust = hjust_vec[i],
                                    size=8, face="bold",
                                    lineheight = 1.2))+
    labs(title = title_vec[i])
  
  
  assign(paste0("all_fig",i), plot_all_two_arms)
  
}


plot_all_dgm <- plot_grid(
  single_arm,
  all_fig1,
  all_fig2,
  all_fig3,
  align = "hv",
  rel_heights = c(1.08,1.08,1,1),
  ncol = 1)

plot_all_dgm

tiff(file = "plots/single_and_two_arm/figure_dgm_survival_hazard_hr.tiff",   
     width = 6.2, 
     height = 7.1,
     units = 'in',  
     res = 300, 
     compression = "lzw")
print(plot_all_dgm)
dev.off()


