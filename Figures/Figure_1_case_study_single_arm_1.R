#########################################################
# Single arm, manuscript plot with both survival and hazard.
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

# Load datasets.
load("surv_plots.Rdata") # load survival stats
load("haz_plots.Rdata") # load hazard stats

colour_fill <- "deepskyblue3"
colour_KM <- "gray20"
margins1 <- unit(c(0.0,0.0,0.0,0), "cm")
margins2 <- unit(c(0.0,0.4,0.0,0), "cm")
hjust_vec <- c(-0.8/.pt, -4.5/.pt, 1.5/.pt)

dataset <- c("Trial only", 
             "Trial + Population rates",
             "Trial + Population rates + Registry")

title_vec <- c("(a) Trial data only",
               "(b) Trial data and population rates",
               "(c) Trial, population rates and SEER registry data")

population_data_plot <- function(plot, dataset){
  
  bh_label <- tibble(x=37, y=0.078, lab = "Population \nmortality")
  seer_label <- tibble(x=15, y=0.3, lab = "SEER registry \ndata")
  
  if(dataset == "Trial only"){
    
    plot 
    
  } else if (dataset == "Trial + Population rates") {
    
    plot+
      geom_step(data=cetux_bh %>% filter(time < 40), 
                aes(x=time, y=hazard), inherit.aes = FALSE, colour="gray30") +
      geom_text(data = bh_label, aes(x=x, y=y, label = lab), 
                colour = "gray30", 
                size = 8/.pt, lineheight = 0.85) 
    
  } else if (dataset == "Trial + Population rates + Registry") {
    
    plot+
      geom_step(data=cetux_seer, 
                aes(x=start, y=haz), inherit.aes = FALSE, colour= "brown4") +
      geom_step(data=cetux_bh %>% filter(time < 40), 
                aes(x=time, y=hazard), inherit.aes = FALSE, colour="gray30") +
      geom_text(data = bh_label, aes(x=x, y=y, label = lab), 
                size = 8/.pt, colour="gray30", lineheight = 0.85) +
      geom_text(data = seer_label, aes(x=x, y=y, label = lab), 
                size = 8/.pt, colour= "brown4", lineheight = 0.85 ) 
    
    
  }
  
}



for(i in 1:3){
  
 surv_plot_single_arm <- surv_plots %>% 
    filter(df == 10, prior_rate==1)  %>%
    filter(`Extra knots` == "10,25") %>%
    filter(`External data` == dataset[i]) %>%
    ggplot() +
    theme_classic()+
    theme(axis.text.y=element_text(size = 6),
          axis.title.y=element_text(size = 8),
          axis.text.x = element_text(size = 6),
          axis.title.x = element_text(size = 8),
          strip.text.x = element_text(size = 4),
          plot.margin = margins1)+
    geom_ribbon(aes(x = t, ymin = lower, ymax = upper), fill= colour_fill, alpha = 0.12)+
    geom_line(aes(x=t, y=median), lwd=1, colour = colour_fill) +
    geom_line(aes(x=t, y=lower), linetype = "solid", alpha = 0.2, colour = colour_fill) +
    geom_line(aes(x=t, y=upper), linetype = "solid", alpha = 0.2, colour = colour_fill) +
    geom_step(data=km_control, aes(x=time, y=surv), lwd=1, colour = colour_KM, 
              inherit.aes = FALSE, alpha = 0.75) +
    geom_vline(xintercept = max(control$years), colour = "gray30", linetype = "dashed") +
    xlab("Time (Years)") + 
    scale_y_continuous("Overall Survival", limits = c(0,1), labels = scales::percent) 
  
  haz_plot_single_arm <- haz_plots %>% 
    filter(t > 0, t <= 40) %>%
    filter(df == 10, prior_rate==1)  %>%
    filter(`Extra knots` == "10,25") %>%
    filter(`External data` == dataset[i]) %>%
    ggplot() +
    theme_classic()+
    theme(axis.text.y=element_text(size = 6),
          axis.title.y=element_text(size = 8),
          axis.text.x = element_text(size = 6),
          axis.title.x = element_text(size = 8),
          strip.text.x = element_text(size = 4),
          plot.margin = margins2)+
    geom_ribbon(aes(x = t, ymin = lower, ymax = upper), fill= colour_fill, alpha = 0.12)+
    geom_line(aes(x=t, y=median), lwd=1, colour = colour_fill) +
    geom_line(aes(x=t, y=lower), linetype = "solid", alpha = 0.2, colour = colour_fill) +
    geom_line(aes(x=t, y=upper), linetype = "solid", alpha = 0.2, colour = colour_fill) +
    geom_vline(xintercept = max(control$years), colour = "gray30", linetype = "dashed")#%.>%
  
  haz_plot_single_arm <- haz_plot_single_arm %.>%
    population_data_plot(plot = ., dataset = dataset[i]) +
    xlab("Time (Years)") + 
    scale_y_continuous("Hazard", limits = c(0, 0.75))
  
  
  assign(paste0("surv_plot_single_arm", i), surv_plot_single_arm)
  assign(paste0("haz_plot_single_arm", i), haz_plot_single_arm)
  
  plot_surv_haz_single_arm <- surv_plot_single_arm + 
    theme(plot.title = element_text(hjust = hjust_vec[i],
                                    size=8, face="bold",
                                    lineheight = 1.2),
          plot.tag = element_text(size = 8))+
    labs(title = title_vec[i])+
    haz_plot_single_arm +
    plot_layout(nrow = 1, ncol = 2)
  
  #plot_surv_haz_single_arm 
  assign(paste0("plot_surv_haz_single_arm", i), plot_surv_haz_single_arm)
  
}

plot_all_single_arm <- plot_grid(
  plot_surv_haz_single_arm1,
  NULL,
  plot_surv_haz_single_arm2,
  NULL,
  plot_surv_haz_single_arm3,
  align = "hv",
  rel_heights = c(1,-0.05,1,-0.05,1),
  ncol = 1)

plot_all_single_arm

tiff(file = "figure1.tiff",   
     width = 5.8, 
     height = 5.3,
     units = 'in',  
     res = 300, 
     compression = "lzw")
print(plot_all_single_arm)
dev.off()

