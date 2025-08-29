##################################################
library(dplyr)
library(tibble)
library(flexsurv)
library(relsurv)


setwd("/home/klvq491/simulation_survival/Simulation-study2-RWE-survextrap/")
setwd("/scratch/klvq491/simsurvextrap_slurm_mix_weib_full1")

dgm_true <- readRDS("dgm_true.rds")
dgm_true

?transrate


##################################################

test <- readRDS("dgm_trial1/trial_data1.rds")
test_res <- readRDS("scen_fit1/fit_est1.rds")
test_res

summary(as.factor(test_res$estimand))

test_ratetable <- survexp.us
dim(test_ratetable)
test_ratetable[,,]
t_vec <- readRDS("t_vec.rds")
rmst_vec <- readRDS("rmst_vec.rds")

##################################

dgm_backhaz <- readRDS("/projects/aa/klvq491/simsurvextrap_slurm_mix_weib_full1/backhaz_mod1/backhaz.rds")
#test_backhaz
saveRDS(test_backhaz,"backhaz_mod1/backhaz.rds")

survexp.us[,1,1]
############ Create life table.

men <- cbind(exp(-365.241*exp(-14.5+.08*(0:100))),exp(-365*exp(-14.7+.085*(0:100))))
women <- cbind(exp(-365.241*exp(-15.5+.085*(0:100))),exp(-365*exp(-15.7+.09*(0:100))))
table <- transrate(men,women,yearlim=c(1980,1990),int.length=10)

#men <- tibble(age = )

dgm_backhaz <- readRDS("/projects/aa/klvq491/simsurvextrap_slurm_mix_weib_full1/backhaz_mod1/backhaz.rds")
gamma_gpm <- dgm_backhaz$gamma_gpm
lambda_gpm <- dgm_backhaz$lambda_gpm

men <- women <- tibble(age = 0:120) %>%
  mutate(survival = (1-pgompertz(q = age+1,  shape = gamma_gpm, rate = lambda_gpm))/
           (1-pgompertz(q = age,  shape = gamma_gpm, rate = lambda_gpm)))  %>%
  select(survival) %>%
  as.matrix()


rate_table_gompertz <- transrate(men,women,yearlim=c(2000, 2000),int.length=1)

#rate_table
#survexp.us
#?transrate()

test_data <- readRDS("dgm_trial1/trial_data1.rds")

test_mod <- flexsurvspline(Surv(time, event)~1, 
               data = test_data, 
           #   bhazard = backhaz,
               k = 3,
               spline = "splines2ns")

plot(test_mod)


test_est <- standsurv(test_mod, type = "rmst", t = 20, ci = TRUE, se = TRUE, boot = F)

test_est


test_mod_bh <- flexsurvspline(Surv(time, event)~1, 
                           data = test_data, 
                           bhazard = backhaz,
                           k = 3,
                           spline = "splines2ns")

test_data_bh <- test_data %>%
  mutate(sex = factor("male")) %>%
  mutate(diag = as.Date("01/01/2000", "%d/%m/%Y")) %>%
  mutate(agedays = floor(age * 365.25))

test_bh_est <- standsurv(test_mod_bh, 
                         type = "rmst",
                         t = 20, 
                         ci = TRUE, 
                         se = TRUE, 
                         boot = F,
                         ratetable = rate_table_gompertz,
                         rmap = list(sex = sex,
                                   year = diag,
                                   age = agedays),
                         scale.ratetable = 365.25,
                         newdata = test_data  %>%
                           mutate(sex = factor("male")) %>%
                           mutate(diag = as.Date("01/01/2000", "%d/%m/%Y")) %>%
                           mutate(agedays = floor(age * 365.25))%>%
                           select(agedays, sex, diag))
#?standsurv()
test_bh_est

flexsurv_single_arm_fit <- function(parametric_model,
                               trial_data,
                               excess_hazard = F,
                               backhaz_data = NA,
                               save_file,
                               seed = NULL){
  
  t_vec <- c(seq(from = 0, to = 5, length.out = 100), seq(from = 5.01, to = 40, length.out = 100))
  rmst_vec <- c(30, 40)
  parametric_model <- "spline"
  excess_hazard <- T
  trial_data <- "dgm_trial1/trial_data1.rds"
  backhaz_data <- "backhaz_mod1/backhaz.rds"
  seed <- 1042
  
  set.seed(seed)
  
  trial_data <- readRDS(trial_data) %>%
    mutate(trt = as.factor(trt))

  
  if(excess_hazard == T){
    
    backhaz <- readRDS(backhaz_data)

    # Gompertz distribution GPM parameters
    lambda_gpm <- backhaz$lambda_gpm
    gamma_gpm <- backhaz$gamma_gpm  
    
    trial_data <- trial_data %>%
      mutate(sex = factor("male")) %>%
      mutate(diag = as.Date("01/01/2000", "%d/%m/%Y")) %>%
      mutate(agedays = floor(age * 365.25))
    
    # Create lifetable based on Gompertz parameters.
    
    men <- women <- tibble(age = 0:120) %>%
      mutate(survival = (1-pgompertz(q = age+1,  shape = gamma_gpm, rate = lambda_gpm))/
               (1-pgompertz(q = age,  shape = gamma_gpm, rate = lambda_gpm)))  %>%
      select(survival) %>%
      as.matrix()
    
    rate_table_gompertz <- transrate(men,women,yearlim=c(2000, 2000),int.length=1)
    
    
  }
  
  
  if(parametric_model %in% c("exp", "weibull") ){
    
    if(excess_hazard == F){
      
      flexmod <- try(flexsurvreg(Surv(time, event)~1, 
                                  data = trial_data, 
                                  dist = parametric_model))  
      
    } else {
      
      flexmod <- try(flexsurvreg(Surv(time, event)~1, 
                                  data = trial_data, 
                                  bhazard = backhaz,
                                  dist = parametric_model))  
      
    }
    
    
    
  } else if(parametric_model == "spline"){
    
    if(excess_hazard == F){
      
      flexmod <- try(flexsurvspline(Surv(time, event)~1, 
                                     data = trial_data, 
                                     k = 3,
                                     spline = "splines2ns"))
      
    } else {
      
      flexmod <- try(flexsurvspline(Surv(time, event)~1, 
                                     data = trial_data, 
                                     bhazard = backhaz,
                                     k = 3,
                                     spline = "splines2ns"))
      
    }
    
  }
  
  mod_run_status <- inherits(flexmod, "try-error")  
#######
  
  if(excess_hazard == F){
    
  est_rmst <- standsurv(flexmod, 
                        type = "rmst", 
                        t = rmst_vec, 
                        ci = TRUE, 
                        se = TRUE, 
                        boot = F)
  
  
  est_survival <- standsurv(flexmod, 
                        type = "survival", 
                        t = t_vec, 
                        ci = TRUE, 
                        se = TRUE, 
                        boot = F)
  
  
  est_hazard <- standsurv(flexmod, 
                            type = "hazard", 
                            t = t_vec, 
                            ci = TRUE, 
                            se = TRUE, 
                            boot = F)
  
  
  } else {
    
    est_rmst <- 
      standsurv(flexmod, 
              type = "rmst",
              t = rmst_vec, 
              ci = TRUE, 
              se = TRUE, 
              boot = F,
              ratetable = rate_table_gompertz,
              rmap = list(sex = sex,
                          year = diag,
                          age = agedays),
              scale.ratetable = 365.25,
              newdata = trial_data  %>%
                mutate(sex = factor("male")) %>%
                mutate(diag = as.Date("01/01/2000", "%d/%m/%Y")) %>%
                mutate(agedays = floor(age * 365.25))%>%
                select(agedays, sex, diag))
   
    est_survival <- 
      standsurv(flexmod, 
                type = "survival",
                t = t_vec, 
                ci = TRUE, 
                se = TRUE, 
                boot = F,
                ratetable = rate_table_gompertz,
                rmap = list(sex = sex,
                            year = diag,
                            age = agedays),
                scale.ratetable = 365.25,
                newdata = trial_data  %>%
                  mutate(sex = factor("male")) %>%
                  mutate(diag = as.Date("01/01/2000", "%d/%m/%Y")) %>%
                  mutate(agedays = floor(age * 365.25))%>%
                  select(agedays, sex, diag))
    
     
    est_hazard <- 
      standsurv(flexmod, 
                type = "hazard",
                t = t_vec, 
                ci = TRUE, 
                se = TRUE, 
                boot = F,
                ratetable = rate_table_gompertz,
                rmap = list(sex = sex,
                            year = diag,
                            age = agedays),
                scale.ratetable = 365.25,
                newdata = trial_data  %>%
                  mutate(sex = factor("male")) %>%
                  mutate(diag = as.Date("01/01/2000", "%d/%m/%Y")) %>%
                  mutate(agedays = floor(age * 365.25)) %>%
                  select(agedays, sex, diag))
    
    

  }
  
  est <- 
    est_rmst %>%
    mutate(estimand = "rmst") %>%
    bind_rows(est_survival %>%
                mutate(estimand = "survival"),
              est_hazard %>%
                mutate(estimand = "hazard")) %>%
    mutate(trt = factor(0, levels = c(0,1))) %>%
    relocate(estimand, trt) %>%
    rename(value = at1, 
           lower = at1_lci,
           upper = at1_uci, 
           sd = at1_se)
  
  

  saveRDS(est, save_file)
  
  
  
  
  
}


?flexsurv









