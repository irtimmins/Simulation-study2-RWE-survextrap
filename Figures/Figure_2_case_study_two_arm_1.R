#########################################################
# Two arms, manuscript plot with both survival and HRs.
#########################################################

library(ggplot2)
library(dplyr)
library(viridis) # for colour palettes
library(stringr)
library(gridExtra)
library(pracma)
library(tidyr)
library(flextable)
library(abind)
library(cowplot)
library(magrittr)
library(purrr)
library(wrapr)
library(patchwork)
library(readr)

# Load survival and hazard results.

# Load from single arms results.
load("surv_plots.Rdata") # load survival stats
load("haz_plots.Rdata") # load survival stats
# Load two arm results.
load("surv_plots_two_arms.Rdata") # load survival stats
load("surv_plots_cetux.Rdata") # load survival stats
load("haz_plots_two_arms.Rdata") # load hazard stats
load("haz_plots_cetux.Rdata") # load hazard stats
load("haz_ratio_plots.Rdata") # load hazard stats



two_arm_models <- c("Proportional \nhazards", 
                    "Non-proportional \nhazards",
                    "Separate arms")

two_arm_dataset <- c("Trial only", 
                     "Trial + \nPopulation rates",
                     "Trial + \nPopulation rates + \nRegistry")

two_arm_title_vec <- c("(a) Proportional excess hazards",
                       "(b) Flexible non-proportional hazards           ",
                       "(c) Separate arms modelling")

two_arm_colour_fill <- "darkgreen"

for(i in 1:3){
  
  # Specify legend parameters
  legend_text <- element_text(size=8)
  legend_key <- unit(0.3, "lines")
  legend_space_y <- unit(1.5, "pt")
  
  # Specify figure margins
  two_arm_margins1 <- unit(c(0.0,0.0,0.0,0), "cm")
  two_arm_margins2 <- unit(c(0.0,0.4,0.0,0), "cm")
  two_arm_hjust_vec <- c(-1.9/.pt, 1.5/.pt,-1.1/.pt)
  
  surv_plot_two_arm <- surv_plots_combined %>% 
    filter(df == 10, prior_rate==1, model == two_arm_models[i])  %>%
    filter(`Extra knots` == "10,15,25") %>%
    filter(`External data` == "Trial + \nPopulation rates + \nRegistry") %>%
    mutate(treat = factor(treat, levels = c("Cetuximab", "Control"),
                          labels = c("Cetuximab+\nRadiotherapy", "Radiotherapy"))) %>%
    mutate(treat_dummy = treat) %>%
    ggplot(aes(colour = treat)) +
    theme_classic()+
    theme(axis.text.y=element_text(size = 6),
          axis.title.y=element_text(size = 8),
          axis.text.x = element_text(size = 6),
          axis.title.x = element_text(size = 8),
          strip.text.x = element_text(size = 4),
          plot.margin = two_arm_margins1,
          legend.position = "inside" ,
          legend.position.inside = c(0.72, 0.75),
          legend.text=legend_text,
          legend.key.size=legend_key,
          legend.key.spacing.y = legend_space_y)+
    geom_ribbon(aes(x = t, ymin = lower, ymax = upper, group = treat, fill = treat), alpha = 0.1, colour = NA)+
    geom_line(aes(x=t, y=median), lwd=0.85) +
    geom_line(aes(x=t, y=lower, group = treat, alpha = treat), linetype = "solid") +
    geom_line(aes(x=t, y=upper, group = treat, alpha = treat), linetype = "solid") +
    geom_vline(xintercept = max(control$years), colour = "gray30", linetype = "dashed") +
    scale_colour_manual(NULL, values = c("chocolate3","royalblue1")) +
    scale_fill_manual(NULL, values =  c("chocolate3","royalblue1")) +
    scale_alpha_manual(NULL, values = c(0.25, 0.5))+
    xlab("Time (Years)") + 
    scale_y_continuous("Overall Survival", limits = c(0,1), labels = scales::percent) 
  
  haz_ratio_plot_two_arm <- haz_ratio_plots %>%
    mutate(model = factor(model, levels = c("PH", "NPH", "Separate fits"),
                          labels = c("Proportional \nhazards", 
                                     "Non-proportional \nhazards",
                                     "Separate arms"))) %>%
    mutate("External data" = factor(
      ed, 
      levels = c("Trial only", "Trial + Population rates", "Trial + Population rates + Registry"),
      labels = c("Trial only", "Trial + \nPopulation rates", "Trial + \nPopulation rates + \nRegistry"))) %>%
    filter(df == 10, prior_rate==1, model == two_arm_models[i])  %>%
    filter(`Extra knots` == "10,15,25") %>%
    filter(`External data` == "Trial + \nPopulation rates + \nRegistry") %>%
    filter(t>0) %>%
    ggplot() +
    theme_classic()+
    theme(axis.text.y=element_text(size = 6),
          axis.title.y=element_text(size = 8),
          axis.text.x = element_text(size = 6),
          axis.title.x = element_text(size = 8),
          strip.text.x = element_text(size = 4),
          plot.margin = two_arm_margins2)+
    geom_ribbon(aes(x = t, ymin = lower, ymax = upper), fill= two_arm_colour_fill, alpha = 0.12)+
    geom_line(aes(x=t, y=median), lwd=1, colour = two_arm_colour_fill) +
    geom_line(aes(x=t, y=lower), linetype = "solid", alpha = 0.2, colour = two_arm_colour_fill) +
    geom_line(aes(x=t, y=upper), linetype = "solid", alpha = 0.2, colour = two_arm_colour_fill) +
    geom_vline(xintercept = max(control$years), colour = "gray30", linetype = "dashed") +
    geom_hline(yintercept = 1, linewidth = 1, 
               alpha = 0.55, colour = "gray40",linetype = "solid")+
    xlab("Time (Years)") + 
    scale_y_continuous("Hazard ratio", breaks = c(0.25, 0.5, 1, 2, 4)) +
    coord_trans(y = "log10", ylim=c(0.2, 4))
  
  assign(paste0("surv_plot_two_arm", i), surv_plot_two_arm)
  assign(paste0(" haz_ratio_plot_two_arm", i),  haz_ratio_plot_two_arm )
  
  plot_surv_hr_two_arm <- surv_plot_two_arm + 
    theme(plot.title = element_text(hjust = two_arm_hjust_vec[i],
                                    size=8, face="bold",
                                    lineheight = 1.2),
          plot.tag = element_text(size = 8))+
    labs(title = two_arm_title_vec[i])+
    haz_ratio_plot_two_arm +
    plot_layout(nrow = 1, ncol = 2)
  
  plot_surv_hr_two_arm
  assign(paste0("plot_surv_hr_two_arm", i), plot_surv_hr_two_arm)
  
}

plot_all_two_arm <- plot_grid(
  plot_surv_hr_two_arm1,
  NULL,
  plot_surv_hr_two_arm2,
  NULL,
  plot_surv_hr_two_arm3,
  align = "hv",
  rel_heights = c(1,-0.05,1,-0.05,1),
  ncol = 1)

tiff(file = "Plots/figure2.tiff",   
     width = 5.8, 
     height = 5.3,
     units = 'in',  
     res = 300, 
     compression = "lzw")
print(plot_all_two_arm)
dev.off()
