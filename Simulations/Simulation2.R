#############################################################
# Load packages
#############################################################

library(tidyr)
library(dplyr)
library(ggplot2)
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
library(emg)
library(survminer)
library(purrr)

#############################################################
# specify path where results are to be stored.
#############################################################

store_res <- "directory/to/store/simulations"
setwd("/home/klvq491/simulation_survival/Simulation-study2-RWE-survextrap/")
setwd("/scratch/klvq491/simsurvextrap_slurm_mix_weib_full2")

#############################################################
# call scripts to create functions.
#############################################################

source("Functions/simulate_dgm_mixture_weibull.R")
source("Functions/simulation_functions.R")
source("Functions/estimands.R")

if(!dir.exists(store_res)){
  dir.create(store_res)
}

setwd(store_res)

#########################################################
# Specify data generating mechanisms (DGMs).
#########################################################

# write DGMs as tibble

set.seed(2035235)

# Number of simuation iterations

nsim <- 1000

# Mixture Weibull models for cause-specific hazard.
weibull_models <- tibble(lambda1 = c(0.52), 
                         lambda2 = c(0.13),
                         gamma1 = c(1.53),
                         gamma2 = c(0.82),
                         pmix = c(0.41)) %>%
  mutate(weibull_model_id = paste0("weibull_mod", row_number()))

# GPM rates.
backhaz_models <- tibble(lambda_gpm = 0.000028, 
                         gamma_gpm = 0.0936) %>%
  mutate(backhaz_model_id = paste0("backhaz_mod", row_number())) 

# Age distribution
# For cetuximab Bonner trial.

age_dist <- tibble(age_mean = 60,
                   age_sd = 9) 

# Frailty effect sizes
# Simulate without frailty effects.
frailty_models <- tibble(alpha_frailty = c(0)) %>%
  mutate(frailty_model_id = paste0("frailty_mod", row_number()))

# Define treatment effect models 
single_arm_trt_effect_models <- tibble(arms = "control",  
                            hr_function = NA,  
                            beta = NA,  
                            design_id = "single_arm",
                            trt_effect_model_id = NA)
# Define coefficients for the treatment effect functions
# (see Supplementary Table 2)
beta1 <- log(0.7)
beta2 <- c(log(0.5)/(1-tanh(-1.2)), 0.8, -1.2)
beta3 <- c(0, -2.8, 0.8, 0.4, 0.35)
saveRDS(beta1, "beta1.rds")
saveRDS(beta2, "beta2.rds")
saveRDS(beta3, "beta3.rds")

two_arm_trt_effect_models <- tibble(arms = "both",
                         hr_function = c("ph","wane", "delay_wane"),
                         beta = c("beta1", "beta2", "beta3"),
                         design_id = "two_arm") %>%
  mutate(trt_effect_model_id = paste0("trt_effect_mod", row_number()))

trt_effect_models <- bind_rows(single_arm_trt_effect_models, 
                               two_arm_trt_effect_models)
 
# Sample size scenarios.

single_arm_sample_models <-  tibble(N = c(200), 
                                    design_id = c("single_arm")) %>%
  mutate(sample_size_model_id = paste0("single_arm_sample_mod", row_number()))

two_arm_sample_models <-  tibble(N = c(400), 
                                 design_id = c("two_arm")) %>%
  mutate(sample_size_model_id = paste0("two_arm_sample_mod", row_number()))

sample_size_models <-  bind_rows(single_arm_sample_models, 
                             two_arm_sample_models) 

# further trial survival parameters
# including censoring distribution and administrative censoring time.

trial_censoring_parameters <- tibble(cens_min = c(1,3,6),
                                    cens_max = c(3,5,8),
                                    left_trunc = c(0,0,0),
                                    maxT = c(3,5,8),
                                    count_format = FALSE) %>%
  mutate(data_cut_model_id = paste0("data_cut_mod", row_number()))
saveRDS(trial_censoring_parameters, "trial_censoring_parameters.rds")

# external data parameters
# including censoring distribution and administrative censoring time.

external_data_models <- tibble(N = 600, # external data sample size.
                               loghaz_bias =  c(0, log(1.1), log(1.2), log(0.9), log(0.8)),
                               arms = "control",
                               hr_function = NA, 
                               beta = NA)  %>%
  mutate(external_bias_model_id = paste0("external_bias_mod", row_number()))

saveRDS(external_data_models, "external_data_models.rds")
external_data_models <- readRDS("external_data_models.rds")

# Labeling needed for figures / tables.

external_data_models_labels <- 
  external_data_models %>%
  arrange(-loghaz_bias) %>%
  mutate(haz_bias = -100+100*exp(loghaz_bias),
         haz_bias_temp = round(haz_bias)) %>%
  mutate(haz_bias = as.character(haz_bias)) %>%
  mutate(haz_bias = if_else(haz_bias_temp > 0, paste0("`+`*`",haz_bias_temp, "`*`%`~`bias`"), haz_bias))  %>%
  mutate(haz_bias = if_else(haz_bias_temp < 0, paste0("`\u2212`*`",abs(haz_bias_temp), "`*`%`~`bias`"), haz_bias))  %>%
  mutate(haz_bias = if_else(haz_bias == "0", "`Unbiased`", haz_bias)) %>%
  select(-haz_bias_temp)  %>%
  rename(external_data_label = haz_bias) %>%
  select(external_bias_model_id, external_data_label) %>%
  add_row(external_bias_model_id = "none", external_data_label = "`No`~`external`~`data`",
          .before = 1)# %>%

saveRDS(external_data_models_labels, "external_data_models_labels.rds")
external_data_models_labels <- readRDS("external_data_models_labels.rds")

external_censoring_parameters <- tibble(cens_min = NA, # no censoring for external data.
                                       cens_max = NA,
                                       left_trunc = 5, # starts at 5-years.
                                       maxT = 25, # ends at 25 years
                                       count_format = TRUE) 

dgm_true <- tibble(nsim = nsim) %>%
  cross_join(weibull_models)  %>%
  cross_join(backhaz_models) %>%
  cross_join(age_dist) %>%
  cross_join(frailty_models) %>%
  mutate(true_control_model_id = paste0("true_control_mod", row_number())) %>% 
  cross_join(trt_effect_models) %>%
  mutate(true_model_id = paste0("true_mod", row_number())) %>%
  relocate(ends_with("_id"), .after = nsim)

saveRDS(dgm_true, "dgm_true.rds")
dgm_true <- readRDS("dgm_true.rds")

dgm_trial <- dgm_true %>%
  full_join(sample_size_models, 
            by = join_by(design_id),
            relationship = "many-to-many") %>%
  cross_join(trial_censoring_parameters) %>%
  mutate(trial_id = paste0("dgm_trial", row_number())) %>%
  relocate(sample_size_model_id, .after = trt_effect_model_id) %>%
  relocate(trial_id)

saveRDS(dgm_trial, "dgm_trial.rds")
dgm_trial <- readRDS("dgm_trial.rds")

dgm_external <- dgm_true %>%
  select(-c(arms, hr_function, beta, design_id, trt_effect_model_id, true_model_id)) %>%
  distinct() %>%
  cross_join(external_data_models) %>%
  cross_join(external_censoring_parameters) %>%
  mutate(external_id = paste0("dgm_external", row_number() )) %>%
  relocate(ends_with("_id"), .after = nsim)

saveRDS(dgm_external, "dgm_external.rds")
dgm_trial <- readRDS("dgm_trial.rds")

#################################################
# General population hazard rates
# needed for excess hazard models.
#################################################

dgm_backhaz_models <- dgm_true %>%
  select(backhaz_model_id, age_mean, age_sd, lambda_gpm, gamma_gpm) %>%
  distinct() 

for(i in 1:nrow(dgm_backhaz_models)){
  if(!dir.exists(dgm_backhaz_models$backhaz_model_id[i])){
    dir.create(dgm_backhaz_models$backhaz_model_id[i])
  }  
}

for(i in 1:nrow(dgm_backhaz_models)){
  dgm_backhaz_models %>%
    filter(backhaz_model_id == dgm_backhaz_models$backhaz_model_id[i]) %>%
    select(lambda_gpm, gamma_gpm) %>%
    saveRDS(paste0(dgm_backhaz_models$backhaz_model_id[i], "/backhaz.rds"))
}

###########################################
# Create pars for slurm script.
###########################################

pars_trial <- dgm_trial %>%
  uncount(nsim) %>% # create individual row for each dataset.
  relocate(N) %>% 
  mutate(nsim = 1, .after = N) %>%
  mutate(loghaz_bias = 0) %>%
  mutate(lower = 1e-08, 
         upper = 10000, 
         nodes = 100) %>%
  mutate(seed = sample(10^6:10^8, n(), replace = FALSE),
         return_data = FALSE) %>%
  relocate(c(hr_function, beta, alpha_frailty, loghaz_bias), .after = pmix) %>%
  relocate(c(cens_min, cens_max, maxT, arms, left_trunc, lower, upper, nodes, 
             seed, count_format, return_data), .after = age_sd) %>%
  group_by(trial_id) %>%
  mutate(save_file = paste0("../", trial_id, "/trial_data", row_number(), ".rds"))  %>%  
  ungroup() %>%
  slice(sample(1:n())) # randomise row positions for even slurm runtimes

saveRDS(pars_trial, "pars_trial.rds")
pars_trial <- readRDS("pars_trial.rds")

pars_single_arm_trial <- pars_trial %>%
  filter(design_id == "single_arm") %>%
  select(-ends_with("_id")) 

saveRDS(pars_single_arm_trial, "pars_single_arm_trial.rds")
pars_single_arm_trial <- readRDS("pars_single_arm_trial.rds")

pars_two_arm_trial <- pars_trial %>%
  filter(design_id == "two_arm") %>%
  select(-ends_with("_id")) 

saveRDS(pars_two_arm_trial, "pars_two_arm_trial.rds")
pars_two_arm_trial <- readRDS("pars_two_arm_trial.rds")

pars_external <- dgm_external %>%
  uncount(nsim) %>% # create individual row for each dataset.
  relocate(N) %>% 
  mutate(nsim = 1, .after = N) %>%
  mutate(lower = 1e-08, 
         upper = 10000, 
         nodes = 100) %>%
  mutate(seed = sample(10^6:10^8, n(), replace = FALSE),
         return_data = FALSE) %>%
  relocate(c(hr_function, beta, alpha_frailty, loghaz_bias), .after = pmix) %>%
  relocate(c(cens_min, cens_max, maxT, arms, left_trunc, lower, upper, nodes, 
             seed, count_format, return_data), .after = age_sd) %>%
  group_by(external_id) %>%
  mutate(save_file = paste0("../", external_id, "/trial_data", row_number(), ".rds"))  %>%  
  ungroup() %>%
  select(-ends_with("_id")) %>%
  slice(sample(1:n())) # randomise row positions for even slurm runtimes

saveRDS(pars_external, "pars_external.rds")
pars_external <- readRDS("pars_external.rds")

########################################################
# Create directories for trial data and external data.
########################################################

for(i in 1:nrow(dgm_trial)){
  if(!dir.exists(dgm_trial$trial_id[i])){
    dir.create(dgm_trial$trial_id[i])
  }  
}

for(i in 1:nrow(dgm_external)){
  if(!dir.exists(dgm_external$external_id[i])){
    dir.create(dgm_external$external_id[i])
  }  
}

#################################################################
## Generate data from each dgm using slurm script.
#################################################################

objects_attach <- c(paste0("beta",1:3),
                 "vecquad_gl", "f_ch", "haz")

package_attach <- c("dplyr", "tidyr", "ggplot2", "readr",
  "flexsurv", "survextrap", "rstan", "emg", "survival",
  "pracma","rstpm2", "posterior", "data.table" , "rsimsum")
saveRDS(package_attach, "package_attach.rds")

sjob_dgm_trial_single_arm <- slurm_apply(sim_dgm_trt_mix, 
                        pars_single_arm_trial, 
                        jobname = "dgm_single",
                        nodes = 4, 
                        cpus_per_node = 4, 
                        submit = TRUE,
                        global_objects = objects_attach,
                        pkgs = package_attach,
                        slurm_options = list(time='02:00:00',
                                             partition='core',
                                             "mem-per-cpu"= '16G'))

saveRDS(sjob_dgm_trial_single_arm, "sjob_dgm_trial_single_arm.rds")

sjob_dgm_trial_two_arm <- slurm_apply(sim_dgm_trt_mix, 
                                      pars_two_arm_trial, 
                                      jobname = "dgm_two",
                                      nodes = 20, 
                                      cpus_per_node = 4, 
                                      submit = TRUE,
                                      global_objects = objects_attach,
                                      pkgs = package_attach,
                                      slurm_options = list(time='02:00:00',
                                                           partition='core',
                                                           "mem-per-cpu"= '16G'))

saveRDS(sjob_dgm_trial_two_arm, "sjob_dgm_trial_two_arm.rds")

sjob_dgm_external <- slurm_apply(sim_dgm_trt_mix_external, 
                                      pars_external, 
                                      jobname = "dgm_ext",
                                      nodes = 40, 
                                      cpus_per_node = 4, 
                                      submit = TRUE,
                                      global_objects = c(objects_attach, 
                                                         "sim_dgm_trt_mix",
                                                         "expected_survival"),
                                      pkgs = package_attach,
                                      slurm_options = list(time='02:00:00',
                                                           partition='core',
                                                           "mem-per-cpu"= '16G'))

saveRDS(sjob_dgm_external, "sjob_dgm_external.rds")


# Check trial/external data has been generated.
pars_single_arm_trial_rerun <- pars_single_arm_trial %>%
  mutate(output_file = substr(save_file, start = 4, stop = nchar(save_file))) %>%
  mutate(job_run = if_else(file.exists(output_file), TRUE, FALSE))
summary(as.factor(pars_single_arm_trial_rerun$job_run))

pars_two_arm_trial_rerun <- pars_two_arm_trial %>%
  mutate(output_file = substr(save_file, start = 4, stop = nchar(save_file))) %>%
  mutate(job_run = if_else(file.exists(output_file), TRUE, FALSE))

pars_external_rerun <- pars_external %>%
  mutate(output_file = substr(save_file, start = 4, stop = nchar(save_file))) %>%
  mutate(job_run = if_else(file.exists(output_file), TRUE, FALSE))

# Rerun any failed jobs.
sjob_dgm_trial_two_arm_rerun <- slurm_apply(sim_dgm_trt_mix, 
                                      pars_two_arm_trial_rerun %>%
                                        filter(job_run == FALSE) %>%
                                        select(-c(output_file, job_run)),# %>%
                                        #mutate(seed = sample(10^4:10^5, n(), replace = FALSE)), # reseed  if necessary. 
                                      jobname = "dgm_two_rerun",
                                      nodes = 12, 
                                      cpus_per_node = 4, 
                                      submit = TRUE,
                                      global_objects = objects_attach,
                                      pkgs = package_attach,
                                      slurm_options = list(time='02:00:00',
                                                           partition='core',
                                                           "mem-per-cpu"= '16G'))

saveRDS(sjob_dgm_trial_two_arm_rerun, "sjob_dgm_trial_two_arm_rerun.rds")


##################################################
# generate very large datasets to obtain true values.
##################################################

for(i in 1:nrow(dgm_true)){
  if(!dir.exists(paste0("dgm_",dgm_true$true_model_id[i]))){
    dir.create(paste0("dgm_",dgm_true$true_model_id[i]))
  }  
}

pars_true_big <- dgm_true %>% 
  mutate(cens_min = NA,
         cens_max = NA,
         left_trunc = 0,
         maxT = 100,
         loghaz_bias = 0,
         count_format = FALSE,
         lower = 1e-08, 
         upper = 10000, 
         nodes = 100) %>%
  mutate(N = 1e5, nsim = 1e3) %>%
  uncount(nsim) %>% # create individual row for each dataset.
  relocate(N) %>% 
  mutate(nsim = 1, .after = N) %>%
  mutate(seed = sample(10^6:10^8, n(), replace = FALSE),
         return_data = FALSE) %>%
  relocate(c(hr_function, beta, alpha_frailty, loghaz_bias), .after = pmix) %>%
  relocate(c(cens_min, cens_max, maxT, arms, left_trunc, lower, upper, nodes, 
             seed, count_format, return_data), .after = age_sd) %>%
  group_by(true_model_id) %>%
  mutate(save_file = paste0("../dgm_", true_model_id, "/data", row_number(), ".rds"))  %>%  
  ungroup() %>%
  slice(sample(1:n())) %>% # randomise row positions for even slurm runtimes
  select(-ends_with("_id")) 
saveRDS(pars_true_big, "pars_true_big.rds")
pars_true_big <- readRDS("pars_true_big.rds")

sjob_dgm_big <- slurm_apply(sim_dgm_trt_mix, 
                                pars_true_big, 
                                 jobname = "dgm_big",
                                 nodes = 20, 
                                 cpus_per_node = 4, 
                                 submit = TRUE,
                                 global_objects = objects_attach,
                                 pkgs = package_attach,
                                 slurm_options = list(time='36:00:00',
                                                      partition='core',
                                                      "mem-per-cpu"= '16G'))

saveRDS(sjob_dgm_big , "sjob_dgm_big.rds")

pars_true_big_rerun <- pars_true_big %>%
  mutate(output_file = substr(save_file, start = 4, stop = nchar(save_file))) %>%
  mutate(job_run = if_else(file.exists(output_file), TRUE, FALSE))

sjob_dgm_big_rerun <- slurm_apply(sim_dgm_trt_mix, 
                            pars_true_big_rerun %>%
                              filter(job_run == FALSE) %>%
                              select(-c(output_file, job_run)), 
                            jobname = "dgm_big_rerun",
                            nodes = 10, 
                            cpus_per_node = 4, 
                            submit = TRUE,
                            global_objects = objects_attach,
                            pkgs = package_attach,
                            slurm_options = list(time='36:00:00',
                                                 partition='core',
                                                 "mem-per-cpu"= '16G'))


############################################
# survextrap model scenarios
############################################

extra_knots_settings <- list(NA,
                             c(5,10,25),
                             c(10,25),
                             c(25)
                             )

count <- 1

for(i in 1:length(extra_knots_settings)){

    if(is.na(extra_knots_settings[[i]][1])) {
      
      names(extra_knots_settings)[i] <- "none"
  
    } else {
     
      names(extra_knots_settings)[i] <- paste0("extra_knots", count)
      
      count <-  count + 1
      
    }
}

#extra_knots_settings
saveRDS(extra_knots_settings, "extra_knots_settings.rds")
extra_knots_settings <- readRDS("extra_knots_settings.rds")

extra_knots_models <- tibble(extra_knots_id = names(extra_knots_settings)) %>%
  mutate(extra_knots_labels = 0)

for(i in 1:nrow(extra_knots_models)){
  #length(extra_knots_settings[[i]])
  if(any(is.na(extra_knots_settings[[i]]))){
    extra_knots_models$extra_knots_labels[i] <- "No extra knots"
  } else if(length(extra_knots_settings[[i]]) ==1 & !any(is.na(extra_knots_settings[[i]]))){
    extra_knots_models$extra_knots_labels[i] <- paste0("Extra knot at t=",paste(extra_knots_settings[[i]], collapse = ","))
  } else if(length(extra_knots_settings[[i]]) > 1){
    extra_knots_models$extra_knots_labels[i] <- paste0("Extra knots at t=",paste(extra_knots_settings[[i]], collapse = ","))
  }
}

saveRDS(extra_knots_models, "extra_knots_models.rds")
extra_knots_models <- readRDS("extra_knots_models.rds")
knot_location_settings <- tibble(knot_settings = c("default"))

single_arm_settings <- tibble(model = "single_arm") %>%
    mutate(design_id = "single_arm")

two_arm_settings <- tibble(model = c("ph","nonph", "separate")) %>%
  mutate(design_id = "two_arm")

two_arm_waning_settings <- tibble(wane_period_start = c(NA, 5, 5, 5),
                                  wane_period_stop = c(NA, 6, 10, 20),
                                  wane_nt = 20)  %>%
  mutate(waning_model_id = paste0("waning_mod", row_number()))
saveRDS(two_arm_waning_settings, "two_arm_waning_settings.rds")

survextrap_settings <- expand_grid(mspline_df = c(3,6,10),
                                    smooth_model = c("random_walk"),
                                    bsmooth = c(T)) %>%
  cross_join(rbind(single_arm_settings,
                   two_arm_settings)) %>%
  cross_join(tibble(
    stan_fit_method = c("mcmc","opt"),
    chains = c(4,NA),
    iter = c(2000, NA))
  ) %>%
  cross_join(two_arm_waning_settings) %>%
  filter(is.na(wane_period_start) | model %in% c("ph", "nonph")) %>%
  filter(stan_fit_method == "mcmc" | mspline_df == 10) %>%  
  cross_join(knot_location_settings)

saveRDS(survextrap_settings, "survextrap_settings.rds")
survextrap_settings <- readRDS("survextrap_settings.rds")

external_data_settings <- expand_grid(
  backhaz = c(F,T),
  include_external_data = c(F,T),
  add_knots = names(extra_knots_settings)) %>%
  filter(add_knots != "none" | include_external_data == F) %>%
  filter(add_knots == "none" | backhaz == T) %>%
  filter(!(backhaz == F & include_external_data == T))
external_data_settings
saveRDS(external_data_settings, "external_data_settings.rds")
external_data_settings <- readRDS("external_data_settings.rds")

# external_data_settings
trial_data_map <- dgm_trial %>%
  select(ends_with("_id"))
saveRDS(trial_data_map, "trial_data_map.rds")
trial_data_map <- readRDS("trial_data_map.rds")

# Create files to map together various ids.
external_data_map <- dgm_external %>%
              select(true_control_model_id,
                     external_bias_model_id,
                     external_id) %>%
  expand_grid(include_external_data = c(T,F)) %>%
  mutate(external_bias_model_id = if_else(include_external_data == F, "none", external_bias_model_id),
         external_id = if_else(include_external_data == F, "none", external_id)) %>%
  distinct() 
saveRDS(external_data_map, "external_data_map.rds")
external_data_map <- readRDS("external_data_map.rds")

# Combine trial, external data and survextrap models
# to encompass all factorial scenarios.

scenarios_all <-  
  trial_data_map %>%
    full_join(external_data_map,  
              join_by(true_control_model_id),
              relationship = "many-to-many") %>%
    distinct() %>%
    full_join(survextrap_settings,
              by = join_by(design_id),
              relationship = "many-to-many") %>%
    full_join(external_data_settings,
              by = join_by(include_external_data),
              relationship = "many-to-many") %>%
    filter(!is.na(backhaz))  %>%
    ungroup() %>%
    mutate(scenario_all_id = row_number())  
  
## Define scenarios of interest

# Single arm sensitivity analysis on data cut.
scenarios_analysis1 <- scenarios_all %>%
  filter(design_id == "single_arm",
         stan_fit_method == "mcmc",
         mspline_df == 10,
         add_knots %in% c("none", "extra_knots1"),
         backhaz == T)

scenarios_analysis_all <- scenarios_analysis1 %>%
  distinct() %>%
  arrange(scenario_all_id)

# Only include scenarios of interest
scenarios <- scenarios_analysis_all %>%
  select(-scenario_all_id) %>%
  mutate(scenario_fit_id = paste0("scen_fit", row_number()))  

# scenarios <- scenarios %>%
#   filter(design_id == "single_arm")
saveRDS(scenarios, "scenarios.rds")
scenarios <- readRDS("scenarios.rds")

scenarios_nsim <- scenarios %>%
  uncount(nsim) %>% # create individual row for each simulation iteration
  group_by(scenario_fit_id) %>% 
  mutate(isim = row_number()) %>%
  ungroup() %>%
  mutate(trial_data = paste0("../", trial_id, "/trial_data", isim, ".rds")) %>%
  mutate(external_data = paste0("../", external_id, "/trial_data", isim, ".rds") ) %>%
  mutate(external_data = if_else(external_id != "none", external_data, NA)) %>%
  mutate(backhaz_data = if_else(backhaz == T, paste0("../", backhaz_model_id, "/backhaz.rds"), NA)) %>%
  mutate(save_file = paste0("../", scenario_fit_id, "/fit_est", isim, ".rds")) %>%
  mutate(save_mspline = NA) %>%
  mutate(seed = sample(10^6:10^8, n(), replace = FALSE)) %>%
  slice(sample(1:n())) 
saveRDS(scenarios_nsim, "scenarios_nsim.rds")
scenarios_nsim <- readRDS("scenarios_nsim.rds")

pars_single_arm_fit <- scenarios_nsim %>%
  filter(design_id == "single_arm") %>%
  select(-c(include_external_data,backhaz, isim)) %>%
  select(-ends_with("_id")) %>%
  relocate(starts_with("wane"), .after = "add_knots")
saveRDS(pars_single_arm_fit, "pars_single_arm_fit.rds")
pars_single_arm_fit <- readRDS("pars_single_arm_fit.rds")

pars_two_arm_fit <- scenarios_nsim %>%
  filter(design_id == "two_arm") %>%
  select(-c(include_external_data,backhaz, isim)) %>%
  select(-ends_with("_id")) %>%
  relocate(starts_with("wane"), .after = "add_knots")
saveRDS(pars_two_arm_fit, "pars_two_arm_fit.rds")
pars_two_arm_fit <- readRDS("pars_two_arm_fit.rds")

############################################
# create directories to store results
############################################

for(i in 1:nrow(scenarios)){
  if(!dir.exists(scenarios$scenario_fit_id[i])){
    dir.create(scenarios$scenario_fit_id[i])
  }
}

############################################
# Run slurm jobs for model fit and estimation.
############################################

t_vec <- c(seq(from = 0, to = 5, length.out = 100), seq(from = 5.01, to = 40, length.out = 100))
rmst_vec <- c(30, 40)
saveRDS(t_vec, "t_vec.rds")
saveRDS(rmst_vec, "rmst_vec.rds")


fit_objects_attach <- c("extra_knots_settings",
                        "rmst_vec",
                        "t_vec",
                        "estimates",
                        "estimates_rmst",
                        "estimates_irmst",
                        "estimates_irmst_separate",
                        "estimates_hr",
                        "estimates_hr_separate",
                        "rmst_samples_GL",
                        "set_knots")

sjob_single_arm_fit <- slurm_apply(fit_est_slurm , 
                            pars_single_arm_fit, 
                            jobname = "fit_single",
                            nodes = 80, 
                            cpus_per_node = 4, 
                            submit = TRUE,
                            global_objects = fit_objects_attach,
                            pkgs = package_attach,
                            slurm_options = list(time='12:00:00',
                                                 partition='core',
                                                 "mem-per-cpu"= '16G'))

saveRDS(sjob_single_arm_fit , "sjob_single_arm_fit.rds")

check_status <- paste0("sacct -S ", as.character(Sys.Date()-4),
                       " -u ",  Sys.info()["user"] ," --format=JobID,Jobname,partition,state,elapsed,ncpus -X")
system(check_status)


saveRDS(sjob_two_arm_fit , "sjob_two_arm_fit.rds")
test <- readRDS("_rslurm_fit_single/results_0.RDS")
test

pars_single_arm_fit %>%
  slice(1:20) %>%
  pull(save_file)

####################################################
## Identify and re-run failed jobs.
#####################################################

pars_single_arm_fit_rerun <- pars_single_arm_fit %>%
  mutate(output_file = substr(save_file, start = 4, stop = nchar(save_file))) %>%
  mutate(job_run = if_else(file.exists(output_file), TRUE, FALSE))
summary(as.factor(pars_single_arm_fit_rerun$job_run))

pars_two_arm_fit_rerun <- pars_two_arm_fit %>%
  mutate(output_file = substr(save_file, start = 4, stop = nchar(save_file))) %>%
  mutate(job_run = if_else(file.exists(output_file), TRUE, FALSE))  %>%
  filter(job_run == F) %>%
  select(-c(output_file, job_run)) 
nrow(pars_two_arm_fit_rerun)

sjob_two_arm_fit_rerun <- slurm_apply(fit_est_slurm, 
                                      pars_two_arm_fit_rerun,
                                      jobname = "fit_two_rerun",
                                      nodes = 10, 
                                      cpus_per_node = 1, 
                                      submit = TRUE,
                                      global_objects = fit_objects_attach,
                                      pkgs = package_attach,
                                      slurm_options = list(time='02:00:00',
                                                           partition='core',
                                                           "mem-per-cpu"= '16G'))

if(!dir.exists("test_reruns")){
  dir.create("test_reruns")
}

# Explore performance of failed jobs.

pars_rerun_test <- pars_two_arm_fit_rerun %>% 
  filter(job_run == FALSE | row_number() <= 100 ) %>%
  select(-c(output_file, job_run)) %>%
  mutate(save_file = paste0("../test_reruns/test", row_number(), ".rds"),
         save_mspline = paste0("../test_reruns/test_mspline", row_number(), ".rds")) %>%
  slice(sample(1:n())) 

sjob_two_arm_fit_rerun_test <- slurm_apply(fit_est_slurm, 
                                pars_rerun_test,
                                jobname = "fit_two_rerun_test",
                                nodes = 10, 
                                cpus_per_node = 1, 
                                submit = TRUE,
                                global_objects = fit_objects_attach,
                                pkgs = package_attach,
                                slurm_options = list(time='00:30:00',
                                                     partition='core',
                                                     "mem-per-cpu"= '16G'))

saveRDS(sjob_two_arm_fit_rerun_test, "sjob_two_arm_fit_rerun_test.rds")
system(check_status)

#######################################################
## Label the estimands
#######################################################

hazard_estimands <- tibble(estimand = "hazard", 
                           t = t_vec) %>%
  expand_grid(trt = c(0,1)) %>%
  group_by(trt) %>%
  mutate(estimand_id = paste0(estimand, row_number(),"_trt", trt ), .after = estimand)

excess_hazard_estimands <- tibble(estimand = "excess_hazard", 
                                 t = t_vec) %>%
  expand_grid(trt = c(0,1)) %>%
  group_by(trt) %>%
  mutate(estimand_id = paste0(estimand, row_number(),"_trt", trt ), .after = estimand)

gpm_hazard_estimands <- tibble(estimand = "gpm_hazard", 
                                  t = t_vec) %>%
  expand_grid(trt = c(0,1)) %>%
  group_by(trt) %>%
  mutate(estimand_id = paste0(estimand, row_number(),"_trt", trt ), .after = estimand)

survival_estimands <- tibble(estimand = "survival", 
                             t = t_vec) %>%
  expand_grid(trt = c(0,1)) %>%
  group_by(trt) %>%
  mutate(estimand_id = paste0(estimand, row_number(),"_trt", trt ), .after = estimand)

hr_estimands <- tibble(estimand = "hr", 
                          t = t_vec) %>%
  mutate(estimand_id = paste0(estimand, row_number()), .after = estimand) %>%
  mutate(trt= NA, .after= estimand_id)

rmst_estimands <- tibble(estimand = "rmst", 
                         t = rmst_vec) %>%
  expand_grid(trt = c(0,1)) %>%
  group_by(trt) %>%
  mutate(estimand_id = paste0(estimand, row_number(),"_trt", trt ), .after = estimand)

irmst_estimands <- tibble(estimand = "irmst", 
                          t = rmst_vec) %>%
  mutate(estimand_id = paste0(estimand, row_number()), .after = estimand) %>%
  mutate(trt= NA, .after= estimand_id)

estimand_labels <- bind_rows(hazard_estimands,
                             excess_hazard_estimands,
                             gpm_hazard_estimands,
                             survival_estimands,
                             hr_estimands,
                             rmst_estimands,
                             irmst_estimands) %>%
  ungroup()

saveRDS(estimand_labels, "estimand_labels.rds")
estimand_labels <- readRDS("estimand_labels.rds")


################################################
# Get true hazard/HR.
################################################

beta1 <- readRDS("beta1.rds")
beta2 <- readRDS("beta2.rds")
beta3 <- readRDS("beta3.rds")
dgm_true <- readRDS("dgm_true.rds")
estimand_labels <- readRDS("estimand_labels.rds")
package_attach <- readRDS("package_attach.rds")

pars_dgm_haz_true <- dgm_true %>%
  select(true_model_id, design_id, lambda1, lambda2, gamma1, gamma2, pmix, lambda_gpm, gamma_gpm, 
         age_mean, age_sd, alpha_frailty, hr_function, beta) %>%
  mutate(N_large = 1e7,
         save_file = paste0("../dgm_",true_model_id, "/hazard_true.rds")) %>%
  select(-true_model_id)

sjob_dgm_haz_true <- slurm_apply(dgm_haz_true, 
                                pars_dgm_haz_true, 
                                jobname = "dgm_haz_true",
                                nodes = 2, 
                                cpus_per_node = 4, 
                                submit = TRUE,
                                pkgs = package_attach,
                                global_objects = c("estimand_labels",
                                                   paste0("beta",1:3), 
                                                   "haz"),
                                slurm_options = list(time='04:00:00',
                                                     partition='core',
                                                     "mem-per-cpu"= '128G'))

################################################
# Get true survival/RMST/iRMST
################################################

# Get big datasets for each true model.

dgm_true <- readRDS("dgm_true.rds")
pars_true_big <- readRDS("pars_true_big.rds")
package_attach <- readRDS("package_attach.rds")

pars_dgm_combine <- dgm_true %>%
  select(nsim, true_model_id) %>%
  mutate(save_file = paste0("../dgm_",true_model_id, "/big_data.rds"))

sjob_dgm_combine <- slurm_apply(dgm_combine, 
                                pars_dgm_combine, 
                                jobname = "dgm_combine",
                                nodes = 10, 
                                cpus_per_node = 4, 
                                submit = TRUE,
                                pkgs = package_attach,
                                slurm_options = list(time='01:00:00',
                                                     partition='core',
                                                     "mem-per-cpu"= '16G'))

pars_dgm_surv_true <- dgm_true %>%
  select(true_model_id, design_id) %>%
  mutate(big_data_file = paste0("../dgm_", true_model_id, "/big_data.rds")) %>%
  mutate(save_file = paste0("../dgm_", true_model_id, "/survival_true.rds"))

sjob_dgm_surv_true <- slurm_apply(dgm_surv_true, 
                                pars_dgm_surv_true, 
                                jobname = "dgm_surv_true",
                                nodes = 10, 
                                cpus_per_node = 4, 
                                submit = TRUE,
                                pkgs = package_attach,
                                global_objects = c("estimand_labels",
                                                   "survival_big"),
                                slurm_options = list(time='02:00:00',
                                                     partition='core',
                                                     "mem-per-cpu"= '16G'))

pars_dgm_rmst_true <- dgm_true %>%
  select(true_model_id, design_id) %>%
  mutate(big_data_file = paste0("../dgm_", true_model_id, "/big_data.rds")) %>%
  mutate(save_file = paste0("../dgm_", true_model_id, "/rmst_true.rds"))

sjob_dgm_rmst_true <- slurm_apply(dgm_rmst_true, 
                                  pars_dgm_rmst_true, 
                                  jobname = "dgm_rmst_true",
                                  nodes = 10, 
                                  cpus_per_node = 4, 
                                  submit = TRUE,
                                  pkgs = package_attach,
                                  global_objects = c("estimand_labels",
                                                     "rmst_big"),
                                  slurm_options = list(time='02:00:00',
                                                       partition='core',
                                                       "mem-per-cpu"= '16G'))


pars_dgm_irmst_true <- dgm_true %>%
  filter(design_id == "two_arm") %>%
  select(true_model_id, design_id) %>%
  mutate(big_data_file = paste0("../dgm_", true_model_id, "/big_data.rds")) %>%
  mutate(save_file = paste0("../dgm_", true_model_id, "/irmst_true.rds"))

sjob_dgm_irmst_true <- slurm_apply(dgm_irmst_true, 
                                  pars_dgm_irmst_true, 
                                  jobname = "dgm_irmst_true",
                                  nodes = 10, 
                                  cpus_per_node = 4, 
                                  submit = TRUE,
                                  pkgs = package_attach,
                                  global_objects = c("estimand_labels",
                                                     "irmst_big"),
                                  slurm_options = list(time='02:00:00',
                                                       partition='core',
                                                       "mem-per-cpu"= '16G'))


################################################
# Combine true for hazard, survival and rmst.
################################################

estimand_labels <- readRDS("estimand_labels.rds")

for(i in 1:nrow(dgm_true)){

  haz_true <- readRDS(paste0("dgm_true_mod", i, "/hazard_true.rds"))
  surv_true <- readRDS(paste0("dgm_true_mod", i, "/survival_true.rds"))
  rmst_true <- readRDS(paste0("dgm_true_mod", i, "/rmst_true.rds"))
  
  all_true <- haz_true %>%
    bind_rows(surv_true) %>%
    bind_rows(rmst_true)
  
  #  View(all_true)
  if(dgm_true$design_id[i] == "two_arm"){
    
    irmst_true <- readRDS(paste0("dgm_true_mod", i, "/irmst_true.rds")) 
    
    all_true <- bind_rows(all_true, irmst_true)
  } 
  
  all_true <- all_true %>%
    select(-t) %>%
    left_join(estimand_labels %>% 
                select(estimand_id, t), 
              join_by(estimand_id) ) %>%
    relocate(t, .after = trt)

  saveRDS(all_true, paste0("dgm_true_mod", i, "/all_true.rds"))   
  
}

####################################################
# combine simulated and true results for estimands. 
####################################################

scenarios <- readRDS("scenarios.rds")
dgm_true <- readRDS("dgm_true.rds")
estimand_labels <- readRDS("estimand_labels.rds")

pars_est_combine <- scenarios %>%
# filter(design_id == "single_arm") %>%
  select(scenario_fit_id, true_model_id) %>%
  left_join(dgm_true, join_by(true_model_id)) %>%
  select(nsim, scenario_fit_id, true_model_id) %>%
  mutate(true_file = paste0("../dgm_",true_model_id, "/all_true.rds" ),
         save_estimates_file = paste0("../", scenario_fit_id, "/fit_est"), 
         save_results_file = paste0("../", scenario_fit_id, "/all_res.rds")) %>%
  select(-true_model_id)

sjob_all_est_combine <- slurm_apply(combine_sim_and_true, 
                             pars_est_combine,
                             jobname = "all_est_combine",
                             nodes = 4, 
                             cpus_per_node = 4, 
                             submit = TRUE,
                             pkgs = package_attach,
                             global_objects = c("estimand_labels"),
                             slurm_options = list(time='00:10:00',
                                                  partition='core',
                                                  "mem-per-cpu"= '4G'))

system(check_status)


####################################################
# Use rsimsum to compute RMST/iRMST performance measures.
####################################################

scenarios <- readRDS("scenarios.rds")
estimand_labels <- readRDS("estimand_labels.rds")

pars_simsum_rmst <- scenarios %>%
#  filter(design_id == "single_arm") %>%
  select(design_id, scenario_fit_id) %>%
  mutate(results_file = paste0("../", scenario_fit_id,"/all_res.rds")) %>%
  mutate(save_file = paste0("../", scenario_fit_id,"/rmst_and_irmst_performance.rds"))

sjob_simsum_rmst <- slurm_apply(simsum_rmst, 
                                pars_simsum_rmst,
                                jobname = "simsum_rmst",
                                nodes = 4, 
                                cpus_per_node = 4, 
                                submit = TRUE,
                                pkgs = package_attach,
                                global_objects = c("estimand_labels"),
                                slurm_options = list(time='02:00:00',
                                                     partition='core',
                                                     "mem-per-cpu"= '16G'))



for(scen_fit in scenarios$scenario_fit_id){
  temp <- readRDS(paste0(scen_fit, "/rmst_and_irmst_performance.rds"))
 if(scen_fit == scenarios$scenario_fit_id[1]){
   comb <- temp
 } else{
   comb <- bind_rows(comb, temp)
 }
  print(scen_fit)
}

performance_all <- comb
performance_all <- as_tibble(performance_all) %>%
  mutate(scenario_fit_id = as.character(scenario_fit_id))
saveRDS(performance_all, "rmst_and_irmst_performance.rds")

performance_all %>%
  select(estimand_id) %>%
  distinct() %>%
  summarise(count = n())

performance_all %>%
  filter(stat == "nsim") %>%
  select(est) %>%
  distinct()


########################################
# Combine mspline knots locations.
########################################

comb <- NULL
for(i in 1:nrow(scenarios)){
  
  for(j in 1:nsim){
 
  mspline <- readRDS(paste0("scen_fit",i, "/mspline", j, ".rds"))  
  mspline$knots
  
  temp <- tibble(knots = 1:length(mspline$knots),
                 value = mspline$knots,
                 isim = j,
                 scenario_fit_id = paste0("scen_fit",i))
  comb <- rbind(comb, temp)
  
  }

 print(i)
}


saveRDS(comb, "knots_all.rds")

########################################
# Create plots and log directory.
########################################

if(!dir.exists("plots")){
  dir.create("plots")
  dir.create("plots/single_arm")
  dir.create("plots/two_arm")
  dir.create("plots/single_and_two_arm")
}  

