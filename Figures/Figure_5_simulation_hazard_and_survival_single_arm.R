#########################################################
# Figure 5, Simulation study,
# survival and hazard functions.
# Combine 5a and 5b together.
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

# Jobname where results are stored.
stores_res <- "directory/to/store/simulations"
setwd(store_res)

# Run Figure_5a and Figure_5b scripts first and import.
surv_plot <- readRDS("plots/Figure_5a.rds")
haz_plot <- readRDS("plots/Figure_5b.rds")

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

tiff(file = "plots/Figure_5.tiff",   
     width = 5.0, 
     height = 6.2,
     units = 'in',  
     res = 300, 
     compression = "lzw")
print(plot_all)
dev.off()

