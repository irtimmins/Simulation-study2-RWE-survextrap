##################################################
library(dplyr)
library(tibble)
library(tidyr)
library(flexsurv)
library(relsurv)
library(rslurm)


source("Functions/simulation_frequentist_functions.R")

setwd("/home/klvq491/simulation_survival/Simulation-study2-RWE-survextrap/")
setwd("/scratch/klvq491/simsurvextrap_slurm_mix_weib_full1")

trial_data_map <- readRDS("trial_data_map.rds")

scenarios_single_arm_freq <- trial_data_map %>%
  filter(design_id == "single_arm",
         data_cut_model_id == "data_cut_mod2") %>%
  expand_grid(parametric_model = c("exp", "weibull", "spline"),
              excess_hazard = T) %>%
  mutate(single_arm_freq_model_id = paste0("single_arm_freq_mod", row_number()))

for(i in 1:nrow(scenarios_single_arm_freq)){
  if(!dir.exists(scenarios_single_arm_freq$single_arm_freq_model_id[i])){
    dir.create(scenarios_single_arm_freq$single_arm_freq_model_id[i])
  }  
}

scenarios_two_arm_freq <- trial_data_map %>%
  filter(design_id == "two_arm",
         data_cut_model_id == "data_cut_mod2") %>%
  expand_grid(parametric_model = c("exp", "weibull", "spline"),
              excess_hazard = T) %>%
  mutate(two_arm_freq_model_id = paste0("two_arm_freq_mod", row_number()))

for(i in 1:nrow(scenarios_two_arm_freq)){
  if(!dir.exists(scenarios_two_arm_freq$two_arm_freq_model_id[i])){
    dir.create(scenarios_two_arm_freq$two_arm_freq_model_id[i])
  }  
}


nsim <- 10

pars_single_arm_freq <- scenarios_single_arm_freq %>%
  uncount(nsim) %>% 
  group_by(single_arm_freq_model_id) %>% 
  mutate(isim = row_number()) %>%
  ungroup() %>%
  mutate(trial_data = paste0("../", trial_id, "/trial_data", isim, ".rds")) %>%
  mutate(save_file = paste0("../", single_arm_freq_model_id, "/fit_est", isim, ".rds"))  %>%
  mutate(backhaz_data = paste0("../", backhaz_model_id, "/backhaz.rds")) %>%
  mutate(seed = sample(10^6:10^8, n(), replace = FALSE)) %>%
  select(parametric_model, trial_data, excess_hazard, backhaz_data, save_file, seed) %>%
  slice(sample(1:n())) 


pars_two_arm_freq <- scenarios_two_arm_freq %>%
  uncount(nsim) %>% 
  group_by(two_arm_freq_model_id) %>% 
  mutate(isim = row_number()) %>%
  ungroup() %>%
  mutate(trial_data = paste0("../", trial_id, "/trial_data", isim, ".rds")) %>%
  mutate(save_file = paste0("../", two_arm_freq_model_id, "/fit_est", isim, ".rds"))  %>%
  mutate(backhaz_data = paste0("../", backhaz_model_id, "/backhaz.rds")) %>%
  mutate(seed = sample(10^6:10^8, n(), replace = FALSE)) %>%
  select(parametric_model, trial_data, excess_hazard, backhaz_data, save_file, seed) %>%
  slice(sample(1:n())) 




t_vec <- readRDS("t_vec.rds")
rmst_vec <- readRDS("rmst_vec.rds")

fit_objects_attach <- c("rmst_vec",
                        "t_vec")

package_attach <- c( "dplyr", "tidyr", "readr" ,
                     "flexsurv","data.table" , "relsurv")


sjob_single_arm_freq <- slurm_apply(flexsurv_single_arm_fit, 
                                    pars_single_arm_freq,
                                    jobname = "fit_single_freq",
                                    nodes = 20, 
                                    cpus_per_node = 1, 
                                    submit = TRUE,
                                    global_objects = fit_objects_attach,
                                    pkgs = package_attach,
                                    slurm_options = list(time='01:00:00',
                                                         partition='core',
                                                         "mem-per-cpu"= '16G'))


sjob_two_arm_freq <- slurm_apply(flexsurv_two_arm_fit, 
                                    pars_two_arm_freq,
                                    jobname = "fit_two_freq",
                                    nodes = 10, 
                                    cpus_per_node = 1, 
                                    submit = TRUE,
                                    global_objects = fit_objects_attach,
                                    pkgs = package_attach,
                                    slurm_options = list(time='01:00:00',
                                                         partition='core',
                                                         "mem-per-cpu"= '16G'))

check_status <- paste0("sacct -S ", as.character(Sys.Date()-4),
                       " -u ",  Sys.info()["user"] ," --format=JobID,Jobname,partition,state,elapsed,ncpus -X")
system(check_status)

# test <- readRDS("_rslurm_fit_single_freq/results_5.RDS")
# test <- readRDS("_rslurm_fit_single_freq/results_5.RDS")
test <- readRDS("two_arm_freq_mod1/fit_est5.rds")
summary(as.factor(test$estimand))


View(test)

