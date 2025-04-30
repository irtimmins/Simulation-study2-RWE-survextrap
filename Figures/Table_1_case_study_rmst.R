#########################################################
# Two arms, manuscript table with RMST and iRMST.
#########################################################

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

load("rmst_two_arms.Rdata") # load rmst stats
load("irmst_two_arms.Rdata") # load rmst stats

summary(as.factor(irmst_two_arms$ek))

control_waning <- rmst_two_arms %>% 
  filter(ek == "10,15,25") %>%
  filter(model == "PH") %>%
  filter(treat == "Control") %>%
  select(-wane) %>%
  expand_grid(wane = c(6, 10, 20))

summary(as.factor(irmst_two_arms$ed))

two_arm_table <- irmst_two_arms %>% 
  mutate(treat = "Contrast") %>%
  bind_rows(rmst_two_arms) %>%
  bind_rows(rmst_cont %>% mutate(model = "Separate fits", treat = "Control")) %>%
  bind_rows(rmst_cetux %>% mutate(model = "Separate fits", treat = "Cetuximab")) %>%
  bind_rows(control_waning) %>%
  mutate(treat = factor(treat, levels = c("Cetuximab", "Control", "Contrast"))) %>%
  filter(ek == "10,15,25") %>%
  mutate(estimand = sprintf("%.2f (%.2f, %.2f)", round(median,2), 
                            round(lower,2), round(upper,2))) %>%
  mutate(ed = factor(ed, 
                     levels = c("Trial only", "Trial + Population rates",  "Trial + Population rates + Registry"),
                     labels = c("trial_only", "trial_gpms", "trial_gpm_registry")))  %>%
  select(-c(df, ek, prior_rate, variable, t, median, lower, upper))%>%
  filter(treat %in% c("Control", "Contrast")) %>%
  arrange(ed) %>%
  arrange(treat) %>%
  #pivot_wider(names_from = treat, values_from = c("trial_only", "trial_gpms", "trial_gpm_registry"))# %>%
  pivot_wider(names_from = c(ed, treat), values_from = estimand) %>%
  mutate(model = factor(model, levels = c("PH", "NPH", "Separate fits"))) %>%
  arrange(wane, model) %>%
  mutate(wane = ifelse(wane==0, NA, wane)) %>%
  rename(Model = model) %>%
  rename(`Waning time (years)` = wane)

two_arm_table
write_csv(two_arm_table, "two_arm_irmst_table.csv")

