#########################################################
# Figure 1, Illustration of structural uncertainty.
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

#############################################################
# specify path to case study data
#############################################################

store_res <- "directory/to/store/case_study"
setwd(store_res)

# Load datasets.
load("surv_plots.Rdata") # load survival stats
load("haz_plots.Rdata") # load hazard stats

# Extract scenarios without external data to plot.

surv_df <- surv_plots %>%
  filter(`External data` == "Trial only",
         df == 10,
         `Prior for sigma` == 1,
         `Extra knots` %in% c("None", "10,15,25"))

haz_df <- haz_plots %>%
  filter(`External data` == "Trial only",
         df == 10,
         `Prior for sigma` == 1,
         `Extra knots` %in% c("None", "10,15,25"))

# Run case_study_script.R to get models_control.
knots1 <- models_control[[1]]$mspline$knots
knots12 <- models_control[[12]]$mspline$knots

colour_fill <- "deepskyblue3"
colour_KM <- "gray20"

margins1 <- unit(c(0.0,0.0,0.0,0), "cm")
margins2 <- unit(c(0.0,0.4,0.0,0), "cm")

hjust1 <- c(0/.pt)
hjust12 <- c(0/.pt)

title1  <- "(a) Model fitted to trial data, without extra knots"
title12  <-  "(b) Model fitted to trial data, with extra knots"

# Model numbers 1 and 12 from case_study_script.1.
for(i in c(1,12)){
  # i <- 12
  knots <- get(paste0("knots", i))
  title <- get(paste0("title", i)) 
  
  surv_plot_single_arm <-
  surv_df %>% 
    filter(model_number == i)  %>%
    ggplot() +
    theme_classic()+
    theme(axis.text.y=element_text(size = 6),
          axis.title.y=element_text(size = 8),
          axis.text.x = element_text(size = 6),
          axis.title.x = element_text(size = 8),
          strip.text.x = element_text(size = 4),
          plot.margin = margins1)+
    geom_vline(xintercept = knots, colour = "gray30", alpha = 0.3)+
    geom_ribbon(aes(x = t, ymin = lower, ymax = upper), fill= colour_fill, alpha = 0.12)+
    geom_line(aes(x=t, y=median), lwd=1, colour = colour_fill) +
    geom_line(aes(x=t, y=lower), linetype = "solid", alpha = 0.2, colour = colour_fill) +
    geom_line(aes(x=t, y=upper), linetype = "solid", alpha = 0.2, colour = colour_fill) +
    geom_step(data=km_control, aes(x=time, y=surv), lwd=1, colour = colour_KM, 
              inherit.aes = FALSE, alpha = 0.75) +
    geom_vline(xintercept = max(control$years), colour = "gray30", linetype = "dashed") +
    xlab("Time (Years)") + 
    scale_y_continuous("Overall Survival", limits = c(0,1), labels = scales::percent) 
  
  haz_plot_single_arm <- haz_df %>% 
    filter(model_number == i)  %>%
    ggplot() +
    theme_classic()+
    theme(axis.text.y=element_text(size = 6),
          axis.title.y=element_text(size = 8),
          axis.text.x = element_text(size = 6),
          axis.title.x = element_text(size = 8),
          strip.text.x = element_text(size = 4),
          plot.margin = margins2)+
    geom_vline(xintercept = knots, colour = "gray30", alpha = 0.3)+
    geom_ribbon(aes(x = t, ymin = lower, ymax = upper), fill= colour_fill, alpha = 0.12)+
    geom_line(aes(x=t, y=median), lwd=1, colour = colour_fill) +
    geom_line(aes(x=t, y=lower), linetype = "solid", alpha = 0.2, colour = colour_fill) +
    geom_line(aes(x=t, y=upper), linetype = "solid", alpha = 0.2, colour = colour_fill) +
    geom_vline(xintercept = max(control$years), colour = "gray30", linetype = "dashed")  +
    xlab("Time (Years)") + 
    scale_y_continuous("Hazard", limits = c(0, 0.75))
  
  assign(paste0("surv_plot_single_arm", i), surv_plot_single_arm)
  assign(paste0("haz_plot_single_arm", i), haz_plot_single_arm)
  
  plot_surv_haz_single_arm <- surv_plot_single_arm + 
    theme(plot.title = element_text(hjust = get(paste0("hjust",i)),
                                    size=8, face="bold",
                                    lineheight = 1.2),
          plot.tag = element_text(size = 8))+
    labs(title = title)+
    haz_plot_single_arm +
    plot_layout(nrow = 1, ncol = 2)
  
  plot_surv_haz_single_arm 
  assign(paste0("plot_surv_haz_single_arm", i), plot_surv_haz_single_arm)
  
}

plot_all_single_arm <- plot_grid(
  plot_surv_haz_single_arm1,
  NULL,
  plot_surv_haz_single_arm12,
  align = "hv",
  rel_heights = c(1,-0.02, 1),
  ncol = 1)

pdf(file = "Plots/Figure_1.pdf",   
     width = 5.8, 
     height = 3.8)  
print(plot_all_single_arm)
dev.off()

ggsave("Plots/Figure_1.svg",
       width = 5.8, 
       height = 3.8)  
