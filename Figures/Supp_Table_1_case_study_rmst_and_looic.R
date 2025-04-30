#####################################################
# LOOIC table, calculated within trial follow-up
#####################################################

library(dplyr)
library(stringr)
library(tidyr)
library(flextable)
library(abind)
library(magrittr)
library(purrr)
library(wrapr)
library(patchwork)
library(readr)


load("looic.Rdata") # load looic stats 

looic_table <- looic %>% 
  mutate(looic = round(looic, 1), 
         rmst = sprintf("%.2f (%.2f, %.2f)", round(median,2), 
                        round(`lower`, 2), round(`upper`, 2)),
         prior_rate = paste0("Gamma(2, ", prior_rate, ")")
  ) %>%
  select(df, ek, prior_rate , looic, rmst) 

write_csv(looic_table, "single_arm_looic_table.csv")