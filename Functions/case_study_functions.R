
library(survextrap)
library(ggplot2)
library(dplyr)
library(viridis) # for colour palettes
library(memoise)
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

########################################################################################
# Case study helper functions.
########################################################################################

# Memoise of survextrap functions.

survextrap_mem <- memoise(survextrap, cache = cachem::cache_disk(dir = ".cache",
                                                                 max_age = 31557600, # keep for 1-year
                                                                 max_size = 2e10)) # 20Gb of space 
rmst_mem <- memoise(rmst_fast, cache = cachem::cache_disk(dir = ".cache",
                                                          max_age = 31557600,
                                                          max_size = 2e10))
irmst_mem <- memoise(irmst_fast, cache = cachem::cache_disk(dir = ".cache",
                                                            max_age = 31557600,
                                                            max_size = 2e10))
survival_mem <- memoise(survival, cache = cachem::cache_disk(dir = ".cache",
                                                             max_age = 31557600,
                                                             max_size = 2e10))
hazard_mem <- memoise(hazard, cache = cachem::cache_disk(dir = ".cache",
                                                         max_age = 31557600,
                                                         max_size = 2e10))
hazard_ratio_mem <- memoise(hazard_ratio, cache = cachem::cache_disk(dir = ".cache",
                                                                     max_age = 31557600,
                                                                     max_size = 2e10))


# Functions for figure style.

theme_paper <- function(){ 
  
  theme_classic() + 
    theme(
      plot.title = element_text(hjust = 0.5, size=20),
      legend.position = "right",
      #legend.position.inside = c(0.6, 0.95),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 14),
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 14),
      strip.text = element_text(size = 20, color = "dark green")
      
    )
}

theme_paper2 <- function(){
  theme(axis.text.y=element_text(size = 6),
        axis.title.y=element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.title.x = element_text(size = 8),
        strip.text.x = element_text(size = 4))
}


