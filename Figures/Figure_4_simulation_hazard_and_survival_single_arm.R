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

jobname <- "mix_weib_full1"
user <- Sys.info()["user"]
project_directory <- paste0("/projects/aa/", user, "/")
store_res <- paste0(project_directory, "simsurvextrap_slurm_", jobname, "/")
setwd(store_res)

 haz_plot <- readRDS(paste0("plots/single_arm/hazard_plot.rds"))
 surv_plot <- readRDS(paste0("plots/single_arm/survival_plot.rds"))
 plot_grid(surv_plot +
             theme(  plot.title = element_text(hjust = -0.2,
                                               size=10, face="bold"))+
             labs(title = "(a) Survival curves"),
           haz_plot +
             theme(plot.title = element_text(hjust = -0.2,
                                             size=10, face="bold"))+
             labs(title = "(b) Hazard curves"),
           align = "v",
           rel_heights=c(0.5,0.5),
           ncol = 1)
 
 plot_all <- plot_grid(surv_plot +
                         theme(  plot.title = element_text(hjust = -0.1,
                                                           size=8, face="bold"))+
                         labs(title = "(a) Survival curves"),
                       haz_plot+
                         theme(  plot.title = element_text(hjust = -0.1,
                                                           size=8, face="bold"))+
                         labs(title = "(b) Hazard curves"),
                       #  rel_widths=c(1, 0.5), 
                       align = "v",
                       rel_heights=c(0.5,0.5),
                       ncol = 1)
 
 #plot_all
 tiff(file = "plots/single_arm/hazard_and_survival_plot.tiff",   
      width = 5.0, 
      height = 6.2,
      units = 'in',  
      res = 300, 
      compression = "lzw")
 print(plot_all)
 dev.off()
 
 
