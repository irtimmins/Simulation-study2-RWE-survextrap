
library(dplyr)
library(tibble)
library(tidyr)
library(flexsurv)
library(relsurv)
library(rslurm)
library(purrr)
library(rsimsum)

source("Functions/simulation_frequentist_functions.R")
source("Functions/estimands.R")


setwd("/home/klvq491/simulation_survival/Simulation-study2-RWE-survextrap/")
setwd("/scratch/klvq491/simsurvextrap_slurm_mix_weib_full1")

trial_data_map <- readRDS("trial_data_map.rds")

scenarios_single_arm_freq <- trial_data_map %>%
  filter(design_id == "single_arm",
         data_cut_model_id == "data_cut_mod2") %>%
  expand_grid(parametric_model = c("exp", "weibull", "spline"),
              excess_hazard = T) %>%
  mutate(single_arm_freq_model_id = paste0("single_arm_freq_mod", row_number()))

saveRDS(scenarios_single_arm_freq, "scenarios_single_arm_freq.rds")

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

saveRDS(scenarios_two_arm_freq, "scenarios_two_arm_freq.rds")


for(i in 1:nrow(scenarios_two_arm_freq)){
  if(!dir.exists(scenarios_two_arm_freq$two_arm_freq_model_id[i])){
    dir.create(scenarios_two_arm_freq$two_arm_freq_model_id[i])
  }  
}


nsim <- 1000

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
                                    nodes = 30, 
                                    cpus_per_node = 4, 
                                    submit = TRUE,
                                    global_objects = fit_objects_attach,
                                    pkgs = package_attach,
                                    slurm_options = list(time='24:00:00',
                                                         partition='core',
                                                         "mem-per-cpu"= '16G'))


sjob_two_arm_freq <- slurm_apply(flexsurv_two_arm_fit, 
                                    pars_two_arm_freq,
                                    jobname = "fit_two_freq",
                                    nodes = 100, 
                                    cpus_per_node = 4, 
                                    submit = TRUE,
                                    global_objects = fit_objects_attach,
                                    pkgs = package_attach,
                                    slurm_options = list(time='48:00:00',
                                                         partition='core',
                                                         "mem-per-cpu"= '16G'))

check_status <- paste0("sacct -S ", as.character(Sys.Date()-4),
                       " -u ",  Sys.info()["user"] ," --format=JobID,Jobname,partition,state,elapsed,ncpus -X")
system(check_status)

#test <- readRDS("_rslurm_fit_single_freq/results_5.RDS")
#test <- readRDS("_rslurm_fit_single_freq/results_5.RDS")
#test <- readRDS("two_arm_freq_mod1/fit_est5.rds")
#summary(as.factor(test$estimand))


##################################################
# Summarise results.
##################################################


estimand_labels <- readRDS("estimand_labels.rds")

pars_est_combine_single_arm_freq <- scenarios_single_arm_freq %>%
  mutate(true_file = paste0("../dgm_",true_model_id, "/all_true.rds" ),
         save_estimates_file = paste0("../", single_arm_freq_model_id, "/fit_est"), 
         save_results_file = paste0("../", single_arm_freq_model_id, "/all_res.rds")) %>%
  mutate(nsim = nsim) %>%
  relocate(nsim) %>%
  select(nsim, single_arm_freq_model_id,  true_file, save_estimates_file, save_results_file ) %>%
  rename(scenario_fit_id = single_arm_freq_model_id)

sjob_all_est_combine_single_arm_freq <- slurm_apply(
  combine_sim_and_true, 
  pars_est_combine_single_arm_freq,
  jobname = "all_est_single_freq",
  nodes = 5, 
  cpus_per_node = 4, 
  submit = TRUE,
  pkgs = package_attach,
  global_objects = c("estimand_labels"),
  slurm_options = list(time='00:10:00',
                       partition='core',
                       "mem-per-cpu"= '4G'))

pars_est_combine_two_arm_freq <- scenarios_two_arm_freq %>%
  mutate(true_file = paste0("../dgm_",true_model_id, "/all_true.rds" ),
         save_estimates_file = paste0("../", two_arm_freq_model_id, "/fit_est"), 
         save_results_file = paste0("../", two_arm_freq_model_id, "/all_res.rds")) %>%
  mutate(nsim = nsim) %>%
  relocate(nsim) %>%
  select(nsim, two_arm_freq_model_id,  true_file, save_estimates_file, save_results_file ) %>%
  rename(scenario_fit_id = two_arm_freq_model_id)

sjob_all_est_combine_two_arm_freq <- slurm_apply(
  combine_sim_and_true, 
  pars_est_combine_two_arm_freq,
  jobname = "all_est_two_freq",
  nodes = 5, 
  cpus_per_node = 4, 
  submit = TRUE,
  pkgs = package_attach,
  global_objects = c("estimand_labels"),
  slurm_options = list(time='00:10:00',
                       partition='core',
                       "mem-per-cpu"= '4G'))

system(check_status)

#View(pars_est_combine_single_arm_freq)

# test <- readRDS("two_arm_freq_mod1/all_res.rds")
# test
# summary(as.factor(test$model_type))
# test2 <- readRDS("_rslurm_all_est_two_freq/results_0.RDS")
# test2
# View(test)
# test3 <- readRDS("dgm_true_mod2/all_true.rds")
# test3
# combine_sim_and_true


####################################################
# Use rsimsum to compute RMST/iRMST performance measures.
####################################################

estimand_labels <- readRDS("estimand_labels.rds")

pars_simsum_rmst_single_arm_freq  <- scenarios_single_arm_freq %>%
  select(design_id, single_arm_freq_model_id) %>%
  mutate(results_file = paste0( single_arm_freq_model_id, "/all_res.rds")) %>%
  mutate(save_file = paste0(single_arm_freq_model_id, "/rmst_and_irmst_performance.rds")) %>%
  rename(scenario_fit_id = single_arm_freq_model_id)

pars_simsum_rmst_two_arm_freq  <- scenarios_two_arm_freq %>%
  select(design_id, two_arm_freq_model_id) %>%
  mutate(results_file = paste0(two_arm_freq_model_id, "/all_res.rds")) %>%
  mutate(save_file = paste0(two_arm_freq_model_id, "/rmst_and_irmst_performance.rds")) %>%
  rename(scenario_fit_id = two_arm_freq_model_id)

pmap(pars_simsum_rmst_single_arm_freq, simsum_rmst)
pmap(pars_simsum_rmst_two_arm_freq, simsum_rmst)

# Combine simsum results across single arm scenarios.

for(scen_fit in scenarios_single_arm_freq$single_arm_freq_model_id){
  temp <- readRDS(paste0(scen_fit, "/rmst_and_irmst_performance.rds"))
  if(scen_fit == scenarios_single_arm_freq$single_arm_freq_model_id[1]){
    comb <- temp
  } else{
    comb <- bind_rows(comb, temp)
  }
  print(scen_fit)
}

performance_all_single_arm_freq <- comb %>%
  as_tibble() %>%
  mutate(scenario_fit_id = as.character(scenario_fit_id))
saveRDS(performance_all_single_arm_freq, "rmst_and_irmst_performance_single_arm_freq.rds")

performance_all_single_arm_freq %>%
  select(estimand_id) %>%
  distinct() %>%
  summarise(count = n())

# Combine simsum results across single arm scenarios.

for(scen_fit in scenarios_two_arm_freq$two_arm_freq_model_id){
  temp <- readRDS(paste0(scen_fit, "/rmst_and_irmst_performance.rds"))
  if(scen_fit == scenarios_two_arm_freq$two_arm_freq_model_id[1]){
    comb <- temp
  } else{
    comb <- bind_rows(comb, temp)
  }
  print(scen_fit)
}

performance_all_two_arm_freq <- comb %>%
  as_tibble() %>%
  mutate(scenario_fit_id = as.character(scenario_fit_id))
saveRDS(performance_all_two_arm_freq, "rmst_and_irmst_performance_two_arm_freq.rds")

performance_all_two_arm_freq %>%
  select(estimand_id) %>%
  distinct() %>%
  summarise(count = n())

performance_all_two_arm_freq %>%
  filter(t == 40) %>%
  filter(estimand == "irmst") %>%
  filter(stat == "bias")

performance_all_single_arm_freq %>%
  # filter(t == 40) %>%
  #  filter(estimand == "irmst") %>%
  filter(stat == "nsim") %>%
  select(est) %>%
  distinct() %>%
  pull()

performance_all_two_arm_freq %>%
 # filter(t == 40) %>%
#  filter(estimand == "irmst") %>%
  filter(stat == "nsim") %>%
  select(est) %>%
  distinct() %>%
  pull()

##########################################





