######################
library(survival)
library(flexsurv)
library(survextrap)
library(tidyr)
library(dplyr)
library(readr)
library(pracma)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(scales)
library(survminer)
library(ggsurvfit)
require(HMDHFDplus)
require(relsurv) 
require(dplyr)

sessionInfo()

# MS: Shouldn't need to set working directory if in an R project
# setwd("/home/klvq491/survival_extrapolation/simsurvextrap")

# jobname <- "mix_weib_test14"
# user <- Sys.info()["user"]
# #scratch_directory <- paste0("/scratch/", user, "/")
# project_directory <- paste0("/projects/aa/statistical_innovation/itimmins/simsurvextrap/aim2_simulations/slurm")
# store_res <- paste0(project_directory, "/simsurvextrap_slurm_", jobname, "/")

#############################################################
# call scripts to create functions.
#############################################################

source("R/estimands_ratetable.R")

##############################################################
country <- c("USA")  ## US rates
for(ctry in country){
  print(ctry)
  # Read in popmort file from HMD if the file does not exist
  if(!file.exists(paste0("popmort_files/",ctry,".f.txt"))){
    # Females
    HMD.f <- readHMDweb(CNTRY = ctry, item = "fltper_1x1", username = myHMDusername,
                        password = myHMDpassword)
    write.table(HMD.f, file=paste0("popmort_files/",ctry,".f.txt"))
    # Males
    HMD.m <- readHMDweb(CNTRY = ctry, item = "mltper_1x1", username = myHMDusername,
                        password = myHMDpassword)
    write.table(HMD.m, file=paste0("popmort_files/",ctry,".m.txt"))
  }
  ratetable <- relsurv::transrate.hmd(male = paste0("popmort_files/",ctry,".m.txt"),
                                      female = paste0("popmort_files/",ctry,".f.txt"))
  
  ## reshape lifetable to be a tidy data.frame, and convert rates to per person-year 
  # (needed if time scale of survival model is in years)
  
  ratetable.df <- as.data.frame.table(ratetable, responseName = "exprate") %>%
    mutate(exprate = 365.25 * exprate) %>%
    mutate(age = as.numeric(as.character(age)),
           year = as.numeric(as.character(year)))
  
  assign(paste0("ratetable.",ctry), ratetable)
  assign(paste0("ratetable.",ctry,".df"), ratetable.df)
}

####################################################

# setwd(store_res)
set.seed(130513)

# Using the cetuximab trial data, we artificially create ages, calendar year and sex
# These are not correlated, nor are they correlated with survival time
# All are assumed male (MS: Why?)
trial_data <- cetux %>%
  as_tibble() %>%
  mutate(trt = if_else(treat == "Cetuximab", 1, 0)) %>%
  mutate(trt = as.factor(trt)) %>%
  select(-months) %>%
  rename(event = d,
         time = years) %>%
  mutate(sex = factor("male"), 
         age = rnorm(n = n(), mean = 60, sd = 9),
         year = sample(seq(as.Date('1999/04/01'), # randomisation year
                           as.Date('2002/03/31'), 
                           by="day"), 
                       n())) %>%
  mutate(attained_age = floor(age + time),
         attained_year = lubridate::year(year + time*365.25)) %>%
  left_join(
    ratetable.USA.df,
    by = c("attained_age"="age", 
           "attained_year"="year", 
           "sex"="sex")
  ) %>%
  rename(backhaz = exprate)


control_data <- trial_data %>%
  filter(trt == 0)

active_data <- trial_data %>%
  filter(trt == 1)

#View(trial_data)


##############################################################
# Calculate expected survival rates for registry data.
##############################################################

set.seed(130513)
# MS: The object below is never used. Why??
sim_ipd_SEER_data <-
  expand_grid(age_years = rnorm(n = 1e3, mean = 60, sd = 9),
              time_years = seq(from = 5, to = 30),
              event = 1) %>%
  mutate(
    age = 365.25*age_years,
    time = 365.25*time_years,
    sex = factor("male"),
    year = sample(seq( # randomisation year
      as.Date('1980/04/01'),
      as.Date('2002/03/31'),
      by = "day"
    ),
    n(),
    replace = TRUE)
  )  


expected_survival <- survexp(time ~ 1, 
                             data= sim_trial_data,
                             rmap = list(year=year, 
                                         age=age,
                                         sex=sex),
                             ratetable=ratetable.USA,
                             method='individual.s')


sim_trial_standarised <- sim_trial_data %>%
  mutate(expsurv = expected_survival) %>%
  group_by(time_years) %>%
  summarise(expsurv = mean(expsurv),
            count = n()) %>%
  rename(time = time_years) %>%
  select(-count) %>%
  mutate(start = time,
         stop = time,
         backsurv_start = expsurv,
         backsurv_stop = expsurv)

sim_trial_standarised_start <- sim_trial_standarised %>%
  select(start, backsurv_start)

sim_trial_standarised_stop <- sim_trial_standarised %>%
  select(stop, backsurv_stop)

# add background rates to cetux_seer.

cetux_seer_backhaz <- cetux_seer %>%
  left_join(sim_trial_standarised_start, join_by("start")) %>%
  left_join(sim_trial_standarised_stop, join_by("stop")) %>%
  mutate(trt = if_else(treat == "Cetuximab", 1, 0)) %>%
  mutate(trt = as.factor(trt))

#############################################################################
# Models without external data.
#############################################################################

set.seed(130513)
no_ext_data_ph_mod <- survextrap(
  Surv(time, event) ~ trt,
  data = trial_data,
 # fit_method = "opt",
  smooth_model = "random_walk",
  backhaz = "backhaz"
)

plot(no_ext_data_ph_mod)

set.seed(130513)
no_ext_data_non_ph_mod <- survextrap(
  Surv(time, event) ~ trt,
  data = trial_data,
  nonprop = T,
#  fit_method = "opt",
  smooth_model = "random_walk",
  backhaz = "backhaz"
)

plot(no_ext_data_non_ph_mod)

set.seed(130513)
no_ext_data_sep_arms_control <- survextrap(
  Surv(time, event) ~ 1,
  data = trial_data %>%
    filter(trt == 0),
 # fit_method = "opt",
  smooth_model = "random_walk",
  backhaz = "backhaz"
)

plot(no_ext_data_sep_arms_control)

set.seed(130513)
no_ext_data_sep_arms_active <- survextrap(
  Surv(time, event) ~ 1,
  data = trial_data %>%
    filter(trt == 1),
 # fit_method = "opt",
  smooth_model = "random_walk",
  backhaz = "backhaz"
)

plot(no_ext_data_sep_arms_active)


#############################################################################
# Models with external data.
#############################################################################

set.seed(130513)
ext_data_ph_mod <- survextrap(
  Surv(time, event) ~ trt,
  data = trial_data,
  external = cetux_seer_backhaz,
  mspline = list(add_knots = c(10,15,25)),
 # fit_method = "opt",
  smooth_model = "random_walk",
  backhaz = "backhaz"
)

plot(ext_data_ph_mod)

set.seed(130513)
ext_data_non_ph_mod <- survextrap(
  Surv(time, event) ~ trt,
  data = trial_data,
  external = cetux_seer_backhaz,
  mspline = list(add_knots = c(10,15,25)),
  nonprop = T,
#  fit_method = "opt",
  smooth_model = "random_walk",   
  backhaz=cetux_bh
)

plot(ext_data_non_ph_mod)

set.seed(130513)
ext_data_sep_arms_control <- survextrap(
  Surv(time, event) ~ 1,
  data = trial_data %>%
    filter(trt == 0),
  external = cetux_seer_backhaz,
  mspline = list(add_knots = c(10,15,25)),
 # fit_method = "opt",
  smooth_model = "random_walk",
  backhaz = "backhaz"
)

plot(ext_data_sep_arms_control)

set.seed(130513)
ext_data_sep_arms_active <- survextrap(
  Surv(time, event) ~ 1,
  data = trial_data %>%
    filter(trt == 1),
  fit_method = "opt",
  smooth_model = "random_walk",
  backhaz = "backhaz"
)

plot(ext_data_sep_arms_active)


#############################################################################
# Estimates for hazard, survival and rmst for each model.
#############################################################################


t_vec <- c(seq(from = 0, to = 5, length.out = 1e2),seq(from = 5+0.001, to = 40, length.out = 1e2))
rmst_vec <- c(5,10,20,30,40)

newdata_trt <- tibble(trt = as.factor(c(0,1))) 

#############################################################################
# Estimates for models without external data.
#############################################################################


## PH models.


est_no_ext_data_ph_mod_control <- estimates_ratetable(fitted_mod = no_ext_data_ph_mod,
                            t = t_vec, 
                            rmst_t = rmst_vec, 
                            newdata = newdata_trt %>%
                              filter(trt == 0), 
                            excess_hazard = TRUE,
                            ratetable = ratetable.USA,
                            scale_ratetable = 365.25,
                            trial_data = trial_data %>%
                              filter(trt == 0))

est_no_ext_data_ph_mod_control <- est_no_ext_data_ph_mod_control %>%
  mutate(model = "ph", external_data = "none", backhaz = T)

est_no_ext_data_ph_mod_control  %>%
  filter(estimand == "rmst")

est_no_ext_data_ph_mod_active <- estimates_ratetable(fitted_mod = no_ext_data_ph_mod,
                                                      t = t_vec, 
                                                      rmst_t = rmst_vec, 
                                                      newdata = newdata_trt %>%
                                                        filter(trt == 1), 
                                                      excess_hazard = TRUE,
                                                      ratetable = ratetable.USA,
                                                      scale_ratetable = 365.25,
                                                      trial_data = trial_data %>%
                                                        filter(trt == 1))

est_no_ext_data_ph_mod_active <- est_no_ext_data_ph_mod_active %>%
  mutate(model = "ph", external_data = "none", backhaz = T)

est_no_ext_data_ph_mod_active %>%
  filter(estimand == "rmst")
est_no_ext_data_ph_mod_control %>%
  filter(estimand == "rmst")

est_irmst_no_ext_data_ph_mod <- estimates_irmst_ratetable(fitted_mod = no_ext_data_ph_mod,
                                                          rmst_t = rmst_vec, 
                                                          newdata = newdata_trt, 
                                                          excess_hazard = TRUE,
                                                          ratetable = ratetable.USA,
                                                          scale_ratetable = 365.25,
                                                          trial_data = trial_data)

est_irmst_no_ext_data_ph_mod <- est_irmst_no_ext_data_ph_mod %>%
  mutate(model = "ph", external_data = "none", backhaz = T)


## Non-PH models.

est_no_ext_data_non_ph_mod_control <- estimates_ratetable(fitted_mod = no_ext_data_non_ph_mod,
                                                      t = t_vec, 
                                                      rmst_t = rmst_vec, 
                                                      newdata = newdata_trt %>%
                                                        filter(trt == 0), 
                                                      excess_hazard = TRUE,
                                                      ratetable = ratetable.USA,
                                                      scale_ratetable = 365.25,
                                                      trial_data = trial_data %>%
                                                        filter(trt == 0))

est_no_ext_data_non_ph_mod_control <- est_no_ext_data_non_ph_mod_control %>%
  mutate(model = "nonph", external_data = "none", backhaz = T)

est_no_ext_data_non_ph_mod_active <- estimates_ratetable(fitted_mod = no_ext_data_non_ph_mod,
                                                     t = t_vec, 
                                                     rmst_t = rmst_vec, 
                                                     newdata = newdata_trt %>%
                                                       filter(trt == 1), 
                                                     excess_hazard = TRUE,
                                                     ratetable = ratetable.USA,
                                                     scale_ratetable = 365.25,
                                                     trial_data = trial_data %>%
                                                       filter(trt == 1))

est_no_ext_data_non_ph_mod_active <- est_no_ext_data_non_ph_mod_active %>%
  mutate(model = "nonph", external_data = "none", backhaz = T)


est_irmst_no_ext_data_non_ph_mod <- estimates_irmst_ratetable(fitted_mod = no_ext_data_non_ph_mod,
                                                          rmst_t = rmst_vec, 
                                                          newdata = newdata_trt, 
                                                          excess_hazard = TRUE,
                                                          ratetable = ratetable.USA,
                                                          scale_ratetable = 365.25,
                                                          trial_data = trial_data)

est_irmst_no_ext_data_non_ph_mod <- est_irmst_no_ext_data_non_ph_mod %>%
  mutate(model = "nonph", external_data = "none", backhaz = T)


## Separate arm models.

est_no_ext_data_sep_arms_mod_control <- estimates_ratetable(fitted_mod = no_ext_data_sep_arms_control,
                                                          t = t_vec, 
                                                          rmst_t = rmst_vec, 
                                                          newdata = newdata_trt %>%
                                                            filter(trt == 0), 
                                                          excess_hazard = TRUE,
                                                          ratetable = ratetable.USA,
                                                          scale_ratetable = 365.25,
                                                          trial_data = trial_data %>%
                                                            filter(trt == 0))

est_no_ext_data_sep_arms_mod_control <- est_no_ext_data_sep_arms_mod_control %>%
  mutate(model = "sep_arms", external_data = "none", backhaz = T)

est_no_ext_data_sep_arms_mod_active <- estimates_ratetable(fitted_mod = no_ext_data_sep_arms_active,
                                                         t = t_vec, 
                                                         rmst_t = rmst_vec, 
                                                         newdata = newdata_trt %>%
                                                           filter(trt == 1), 
                                                         excess_hazard = TRUE,
                                                         ratetable = ratetable.USA,
                                                         scale_ratetable = 365.25,
                                                         trial_data = trial_data %>%
                                                           filter(trt == 1))

est_no_ext_data_sep_arms_mod_active <- est_no_ext_data_sep_arms_mod_active %>%
  mutate(model = "sep_arms", external_data = "none", backhaz = T)


est_irmst_no_ext_data_sep_arms_mod <- estimates_irmst_separate_ratetable(
  control_model = no_ext_data_sep_arms_control,
  active_model = no_ext_data_sep_arms_active,
  rmst_t = rmst_vec, 
  newdata = tibble(trt = as.factor(c(0,1))), 
  excess_hazard = TRUE,
  ratetable = ratetable.USA,
  scale_ratetable = 365.25,
  trial_data = trial_data)

est_irmst_no_ext_data_sep_arms_mod <- est_irmst_no_ext_data_sep_arms_mod %>%
  mutate(model = "sep_arms", external_data = "none", backhaz = T)


#############################################################################
# Estimates for models with external data.
#############################################################################

est_ext_data_ph_mod_control <- estimates_ratetable(fitted_mod = ext_data_ph_mod,
                                                      t = t_vec, 
                                                      rmst_t = rmst_vec, 
                                                      newdata = newdata_trt %>%
                                                        filter(trt == 0), 
                                                      excess_hazard = TRUE,
                                                      ratetable = ratetable.USA,
                                                      scale_ratetable = 365.25,
                                                      trial_data = trial_data %>%
                                                        filter(trt == 0))

est_ext_data_ph_mod_control <- est_ext_data_ph_mod_control %>%
  mutate(model = "ph", external_data = "SEER", backhaz = T)


est_ext_data_ph_mod_active <- estimates_ratetable(fitted_mod = ext_data_ph_mod,
                                                     t = t_vec, 
                                                     rmst_t = rmst_vec, 
                                                     newdata = newdata_trt %>%
                                                       filter(trt == 1), 
                                                     excess_hazard = TRUE,
                                                     ratetable = ratetable.USA,
                                                     scale_ratetable = 365.25,
                                                     trial_data = trial_data %>%
                                                       filter(trt == 1))

est_ext_data_ph_mod_active <- est_ext_data_ph_mod_active %>%
  mutate(model = "ph", external_data = "SEER", backhaz = T)

est_irmst_ext_data_ph_mod <- estimates_irmst_ratetable(fitted_mod = ext_data_ph_mod,
                                                          rmst_t = rmst_vec, 
                                                          newdata = newdata_trt, 
                                                          excess_hazard = TRUE,
                                                          ratetable = ratetable.USA,
                                                          scale_ratetable = 365.25,
                                                          trial_data = trial_data)

est_irmst_ext_data_ph_mod <- est_irmst_ext_data_ph_mod %>%
  mutate(model = "ph", external_data = "SEER", backhaz = T)


est_ext_data_non_ph_mod_control <- estimates_ratetable(fitted_mod = ext_data_non_ph_mod,
                                                          t = t_vec, 
                                                          rmst_t = rmst_vec, 
                                                          newdata = newdata_trt %>%
                                                            filter(trt == 0), 
                                                          excess_hazard = TRUE,
                                                          ratetable = ratetable.USA,
                                                          scale_ratetable = 365.25,
                                                          trial_data = trial_data %>%
                                                            filter(trt == 0))

est_ext_data_non_ph_mod_control <- est_ext_data_non_ph_mod_control %>%
  mutate(model = "nonph", external_data = "SEER", backhaz = T)

est_ext_data_non_ph_mod_active <- estimates_ratetable(fitted_mod = ext_data_non_ph_mod,
                                                         t = t_vec, 
                                                         rmst_t = rmst_vec, 
                                                         newdata = newdata_trt %>%
                                                           filter(trt == 1), 
                                                         excess_hazard = TRUE,
                                                         ratetable = ratetable.USA,
                                                         scale_ratetable = 365.25,
                                                         trial_data = trial_data %>%
                                                           filter(trt == 1))

est_ext_data_non_ph_mod_active <- est_ext_data_non_ph_mod_active %>%
  mutate(model = "nonph", external_data = "SEER", backhaz = T)

est_irmst_ext_data_non_ph_mod <- estimates_irmst_ratetable(fitted_mod = ext_data_non_ph_mod,
                                                       rmst_t = rmst_vec, 
                                                       newdata = newdata_trt, 
                                                       excess_hazard = TRUE,
                                                       ratetable = ratetable.USA,
                                                       scale_ratetable = 365.25,
                                                       trial_data = trial_data)

est_irmst_ext_data_non_ph_mod <- est_irmst_ext_data_non_ph_mod %>%
  mutate(model = "nonph", external_data = "SEER", backhaz = T)


est_ext_data_sep_arms_mod_control <- estimates_ratetable(fitted_mod = ext_data_sep_arms_control,
                                                            t = t_vec, 
                                                            rmst_t = rmst_vec, 
                                                            newdata = newdata_trt %>%
                                                              filter(trt == 0), 
                                                            excess_hazard = TRUE,
                                                            ratetable = ratetable.USA,
                                                            scale_ratetable = 365.25,
                                                            trial_data = trial_data %>%
                                                              filter(trt == 0))

est_ext_data_sep_arms_mod_control <- est_ext_data_sep_arms_mod_control %>%
  mutate(model = "sep_arms", external_data = "SEER", backhaz = T)

est_ext_data_sep_arms_mod_active <- estimates_ratetable(fitted_mod = ext_data_sep_arms_active,
                                                           t = t_vec, 
                                                           rmst_t = rmst_vec, 
                                                           newdata = newdata_trt %>%
                                                             filter(trt == 1), 
                                                           excess_hazard = TRUE,
                                                           ratetable = ratetable.USA,
                                                           scale_ratetable = 365.25,
                                                           trial_data = trial_data %>%
                                                             filter(trt == 1))

est_ext_data_sep_arms_mod_active <- est_ext_data_sep_arms_mod_active %>%
  mutate(model = "sep_arms", external_data = "SEER", backhaz = T)

est_irmst_ext_data_sep_arms_mod <- estimates_irmst_separate_ratetable(
  control_model = ext_data_sep_arms_control,
  active_model = ext_data_sep_arms_active,
  rmst_t = rmst_vec, 
  newdata = tibble(trt = as.factor(c(0,1))), 
  excess_hazard = TRUE,
  ratetable = ratetable.USA,
  scale_ratetable = 365.25,
  trial_data = trial_data)

est_irmst_ext_data_sep_arms_mod <- est_irmst_ext_data_sep_arms_mod %>%
  mutate(model = "sep_arms", external_data = "SEER", backhaz = T)



#############################################################################
# Combine data for figures.
#############################################################################

model_levels <-  c("ph", "nonph", "sep_arms")

model_labels <- c("Proportional hazards",
                  "Non-proportional hazards",
                  "Separate arms")

est_all <- bind_rows(est_no_ext_data_ph_mod_control,
                     est_no_ext_data_ph_mod_active,
                     est_no_ext_data_non_ph_mod_control,
                     est_no_ext_data_non_ph_mod_active,
                     est_no_ext_data_sep_arms_mod_control,
                     est_no_ext_data_sep_arms_mod_active,
                     est_ext_data_ph_mod_control,
                     est_ext_data_ph_mod_active,
                     est_ext_data_non_ph_mod_control,
                     est_ext_data_non_ph_mod_active,
                     est_ext_data_sep_arms_mod_control,
                     est_ext_data_sep_arms_mod_active) %>%
  mutate(external_data = factor(external_data,
                                labels = c("No external data", "Incorporating SEER data"),
                                levels = c("none", "SEER"))) %>%
  mutate(Model = factor(model,
                                labels = model_labels,
                                levels = model_levels)) %>%
  mutate(Treatment = factor(trt,
                        labels = c("Radiotherapy", "Cetuximab"),
                        levels = c(0, 1)))


km_fit <- survfit(Surv(time, event) ~ trt, data=trial_data)
#plot(km_fit)

km_plot <- ggsurvplot(km_fit, data=trial_data)
#km_plot <- ggsurvplot(km_fit, data=big_df %>% slice(1:1000))

km_surv_df <- km_plot$data.survplot %>%
  mutate(Treatment = factor(trt,
                            labels = c("Radiotherapy", "Cetuximab"),
                            levels = c(0, 1))) %>%
  mutate(Model = NA)

max_trial_data <- max(km_surv_df$time)

# Survival plots.

surv_plot <- est_all %>%
  filter(estimand == "survival") %>%
  ggplot(aes(x = t, y = value, 
             colour = Model, 
             linetype = Treatment))+
  theme_classic()+
  theme(axis.text.y=element_text(size = 6),
        axis.title.y=element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.title.x = element_text(size = 8),
        strip.text.x = element_text(size = 6),
        legend.title=element_text(size=8), 
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.4, "cm")) +
    geom_vline(xintercept = c(max_trial_data), colour = c("gray60"))+
    geom_vline(xintercept = c(25), colour = c("gray70"))+
  geom_line(alpha = 1)+
  geom_step(data = km_surv_df %>% filter(trt == 0), 
            aes(x = time, y = surv),
            colour = "gray30",
            linewidth = 0.6)+
  geom_step(data = km_surv_df %>% filter(trt == 1), 
            aes(x = time, y = surv),
            colour = "gray20",
            alpha = 0.8,
            linetype = "dashed",
            linewidth = 0.4)+
  geom_step(data = km_surv_df %>% filter(trt == 1), 
            aes(x = time, y = surv),
            colour = "gray20",
            alpha = 0.6,
            linetype = "dotted",
            linewidth = 0.4)+
  #scale_linewidth(range = c(0.4,1))+
  #scale_alpha(range = c(0.15,1))+
 # scale_color_manual(values = c(hue_pal()(3), "black"),
#                       breaks = model_labels) + 
  scale_x_continuous("Time (years)")+
  scale_y_continuous("Survival",limits = c(0,1), labels = scales::percent)+
  facet_wrap(~external_data, ncol = 1)
 # guides(linewidth = "none", alpha = "none", colour = "none") 

surv_plot


haz_plot <- est_all %>%
  filter(estimand == "hazard") %>%
  ggplot(aes(x = t, y = value, 
             colour = Model, 
             linetype = Treatment))+
  theme_classic()+
  theme(axis.text.y=element_text(size = 6),
        axis.title.y=element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.title.x = element_text(size = 8),
        strip.text.x = element_text(size = 6),
        legend.title=element_text(size=8), 
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.4, "cm")) +
  geom_vline(xintercept = c(max_trial_data), colour = c("gray60"))+
  geom_vline(xintercept = c(25), colour = c("gray70"))+
  geom_line(alpha = 1)+
  geom_step(data=cetux_seer %>% mutate(external_data = as.factor("Incorporating SEER data")),
            aes(x=start, y=haz), inherit.aes = FALSE,
            colour = "coral4",
            alpha = 1,
            linewidth = 0.5) +
  #scale_linewidth(range = c(0.4,1))+
  #scale_alpha(range = c(0.15,1))+
  #scale_color_manual(values = c("#830051", "black")) + 
  scale_x_continuous("Time (years)")+
  scale_y_continuous("Hazard")+
  facet_wrap(~external_data, ncol = 1)

haz_plot
surv_plot
# guides(linewidth = "none", alpha = "none", colour = "none") 

  #head()

plot_all <- plot_grid(surv_plot +
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

#?plot_grid()
plot_all

saveRDS(plot_all, "plots/figure_case_study_two_arms.rds")


tiff(file = "plots/figure_case_study_two_arms.tiff",   
     width = 5.5, 
     height = 6.2,
     units = 'in',  
     res = 600, 
     compression = "lzw")
print(plot_all)
dev.off()




###############################################
# Create table for RMST results.
###############################################


est_rmst_all <- bind_rows(est_no_ext_data_ph_mod_control,
                     est_no_ext_data_ph_mod_active,
                     est_no_ext_data_non_ph_mod_control,
                     est_no_ext_data_non_ph_mod_active,
                     est_no_ext_data_sep_arms_mod_control,
                     est_no_ext_data_sep_arms_mod_active,
                     est_ext_data_ph_mod_control,
                     est_ext_data_ph_mod_active,
                     est_ext_data_non_ph_mod_control,
                     est_ext_data_non_ph_mod_active,
                     est_ext_data_sep_arms_mod_control,
                     est_ext_data_sep_arms_mod_active) %>%
  filter(estimand == "rmst")

est_irmst_all <- bind_rows(est_irmst_no_ext_data_ph_mod,
                           est_irmst_no_ext_data_non_ph_mod,
                           est_irmst_no_ext_data_sep_arms_mod,
                           est_irmst_ext_data_ph_mod,
                           est_irmst_ext_data_non_ph_mod,
                           est_irmst_ext_data_sep_arms_mod)

est_table <- bind_rows(est_rmst_all,
                       est_irmst_all) %>%
  filter(t == 30)

################################################

est_table_wide <- est_table %>%
  mutate(value_chr = sprintf("%.2f", round(value, 2)),
         lower_chr = sprintf("%.2f", round(lower, 2)),
         upper_chr = sprintf("%.2f", round(upper, 2))) %>%
  mutate(value_str = paste0(value_chr, " (", lower_chr, ", ", upper_chr, ")")) %>%
  mutate(category = 
           case_when(
             trt == 0 ~ "Control",
             trt == 1 ~ "Active",
             is.na(trt) == T ~ "Difference",
             .default = as.character(trt)
           )
  ) %>%
  mutate(Model = 
           case_when(
             model == "ph" ~ "Proportional hazards",
             model == "nonph" ~ "Non-proportional hazards",
             model == "sep_arms" ~ "Separate arms",
             .default = as.character(model)
           )
  ) %>%
  mutate(External_data = 
           case_when(
             external_data == "none" ~ "Trial data alone",
             external_data == "SEER" ~ "Trial and SEER registry data",
             .default = as.character(external_data)
           )
  ) %>%
  select(category, External_data,Model, value_str) %>%
  pivot_wider(names_from = category,
              values_from = value_str)

write_csv(est_table_wide, "tables/case_study_rmst.csv")

################################################







