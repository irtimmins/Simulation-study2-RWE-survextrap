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
library(pracma)
library(posterior)
library(stats)
library(data.table)
library(scales)
library(ggpubr)
library(grid)
library(gridExtra) 


jobname <- "mix_weib_full1"
user <- Sys.info()["user"]
project_directory <- paste0("/projects/aa/", user, "/")
store_res <- paste0(project_directory, "simsurvextrap_slurm_", jobname, "/")
setwd(store_res)

scenarios <- readRDS("scenarios.rds")
performance_res <- readRDS("rmst_and_irmst_performance.rds")

estimand_labels <- readRDS("estimand_labels.rds") 

irmst_estimand_vec <- estimand_labels %>%
  filter(estimand == "irmst") %>%
  pull(estimand_id)

for(irmst_estimand in irmst_estimand_vec){
#  irmst_estimand <- "irmst1"
  estimand_number <- substr(irmst_estimand, 6, nchar(irmst_estimand))
  
  single_arm_plot <- readRDS(paste0("plots/single_arm/forest_rmst",estimand_number,"_trt0_single_arm_plot.rds"))

  single_arm_plot_legend <-  cowplot::get_legend(
    # create some space to the left of the legend
    single_arm_plot + theme(legend.title = element_text(margin = margin(3, 0, 3, 0)),
                            legend.box.margin = margin(0, 12, 0, 12),
                            legend.text = element_text(size=7),
                            legend.key.spacing.y = unit(1, "pt"))
  )
  
 single_arm_without_legend <- single_arm_plot+
     theme(legend.position="none",
           plot.title = element_text(size=10, face="bold",
                                     lineheight = 1.2))+
   labs(title = "Single-arm trial \n  (a) Weibull mixture")

 
 single_arm_plot_all <- plot_grid(NULL, 
                                  single_arm_without_legend, 
                                  NULL, 
                                  single_arm_plot_legend, 
                                  NULL, 
                                  rel_widths = c(0.5, 0.8,-0.02,1, 0.3), 
                                  nrow = 1)
 
 single_arm_plot_all 
  
  two_arm_plot1 <- readRDS(paste0("plots/two_arm/forest_", irmst_estimand,"_model", 1, ".rds"))
  two_arm_plot2 <- readRDS(paste0("plots/two_arm/forest_", irmst_estimand, "_model", 2, ".rds"))
  two_arm_plot3 <- readRDS(paste0("plots/two_arm/forest_", irmst_estimand, "_model", 3, ".rds"))
  #plot1
    # plot2
    #plot3
    
    
    plot_two_arm_legend <- cowplot::get_legend(
      # create some space to the left of the legend
      two_arm_plot1 + theme(legend.title = element_text(margin = margin(3, 0, 3, 0)),
                            legend.box.margin = margin(0, 12, 0, 12),
                            legend.text = element_text(size=7),
                            legend.key.spacing.y = unit(1, "pt"),
                            legend.spacing.y = unit(-10, "pt"))
    )
    
    plot_two_arm_without_legend <- plot_grid(
      two_arm_plot1+
        theme(legend.position="none",
              plot.title = element_text(size=10, face="bold",
                                        lineheight = 1.2))+
        labs(title = "Two-arm trials \n  (b) Scenario 1: \n Constant effect"),
      two_arm_plot2+
        theme(legend.position="none",
              plot.title = element_text(size=10, face="bold",
                                        lineheight = 1.2))+
        labs(title = "\n (c) Scenario 2: \n Waning effect"),
      two_arm_plot3+
        theme(legend.position="none",
              plot.title = element_text(size=10, face="bold",
                                        lineheight = 1.2))+
        labs(title = "\n (d) Scenario 3: Delayed then \n waning effect"),
      align = "h",
      rel_widths=c(1,1,1),
      nrow = 1)
   
    
   two_arm_plot_all <- plot_grid(plot_two_arm_without_legend,
      plot_two_arm_legend,
      rel_widths=c(1, 0.6),
      axis = "l"
    )
    
  two_arm_plot_all


  plot_all <- plot_grid(single_arm_plot_all,
            two_arm_plot_all,
          #  rel_widths=c(1, 0.5),
            axis = "l",
            nrow = 2,
            rel_heights=c(0.62, 1))
  
  plot_all
  
  tiff(file = paste0("plots/single_and_two_arm/forest_rmst_and_irmst", 
                     estimand_number, "_all.tiff"),   
       width = 6.6, 
       height = 6.8,
       units = 'in',  
       res = 300, 
       compression = "lzw")
  print(plot_all)
  dev.off()
  
  paste0("plots/single_and_two_arm/forest_rmst_and_irmst", 
         estimand_number, "_all.tiff")

}


