
# Function for model fit and estimation of exponential, weibull and spline
# single arm models from flexsurv.

flexsurv_single_arm_fit <- function(parametric_model,
                                    trial_data,
                                    excess_hazard = F,
                                    backhaz_data = NA,
                                    save_file,
                                    seed = NULL){
  
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
  
  if(!mod_run_status){
    
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
      rename(t = time,
             value = at1, 
             lower = at1_lci,
             upper = at1_uci, 
             sd = at1_se) %>%
      relocate(-sd)
    
    
    
    saveRDS(est, save_file)
    
    
    
  }
  
  
} 


# Function for model fit and estimation of exponential, weibull and spline
# two arm models from flexsurv.


flexsurv_two_arm_fit <- function(parametric_model,
                                 trial_data,
                                 excess_hazard = F,
                                 backhaz_data = NA,
                                 save_file,
                                 seed = NULL){
  
  # t_vec <- c(seq(from = 0, to = 5, length.out = 100), seq(from = 5.01, to = 40, length.out = 100))
  # rmst_vec <- c(30, 40)
  # 
  # excess_hazard <- T
  # trial_data <- "dgm_trial4/trial_data1.rds"
  # backhaz_data <- "backhaz_mod1/backhaz.rds"
  # seed <- 1042
  # 
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
  
  
  # parametric_model <- "weibull"
  # 
  # test0 <- try(flexsurvreg(Surv(time, event)~1,
  #                          data = trial_data %>%
  #                            filter(trt == 0), 
  #                          #   bhazard = backhaz,
  #                          dist = parametric_model))
  # 
  # test1 <- try(flexsurvreg(Surv(time, event)~1,
  #                          data = trial_data %>%
  #                            filter(trt == 1), 
  #                          #  bhazard = backhaz,
  #                          dist = parametric_model))
  # 
  # test0$loglik+test1$loglik
  # 
  # flexmod$loglik
  # 
  # 
  # test0 <- try(flexsurvspline(Surv(time, event)~1,
  #                             data = trial_data %>%
  #                               filter(trt == 0), 
  #                             bknots =  c(-3.3, 1),
  #                             knots = c(log(1), log(1.5), log(2)),
  #                             k = 3,
  #                             spline = "splines2ns"))
  # plot(test0)
  # 
  # test1 <- try(flexsurvspline(Surv(time, event)~1,
  #                             data = trial_data %>%
  #                               filter(trt == 1), 
  #                             bknots =  c(-3.3, 1),
  #                             knots = c(log(1), log(1.5), log(2)),
  #                             #    bhazard = backhaz,
  #                             k = 3,
  #                             spline = "splines2ns"))
  # 
  # plot(test1)
  # test0$loglik+test1$loglik
  # #test1$knots
  # 
  # # flexmod$loglik
  # # 
  # standsurv(test0, 
  #           type = "rmst", 
  #           t = rmst_vec, 
  #           ci = TRUE, 
  #           se = TRUE, 
  #           boot = F)
  # standsurv(test1, 
  #           type = "rmst", 
  #           t = rmst_vec, 
  #           ci = TRUE, 
  #           se = TRUE, 
  #           boot = F)
  
  
  if(parametric_model == "exp"){
    
    if(excess_hazard == F){
      
      flexmod <- try(flexsurvreg(Surv(time, event)~trt,
                                 data = trial_data, 
                                 dist = parametric_model))
      flexmod
      
    } else {
      
      flexmod <- try(flexsurvreg(Surv(time, event)~trt, 
                                 data = trial_data, 
                                 bhazard = backhaz,
                                 dist = parametric_model)) 
      
    } 
    
  }
  
  else if(parametric_model == "weibull"){
    
    if(excess_hazard == F){
      
      flexmod <- try(flexsurvreg(Surv(time, event)~trt,
                                 data = trial_data, 
                                 dist = parametric_model,
                                 anc = list(shape=~trt)))
      flexmod
      
    } else {
      
      flexmod <- try(flexsurvreg(Surv(time, event)~trt, 
                                 data = trial_data, 
                                 bhazard = backhaz,
                                 dist = parametric_model,
                                 anc = list(shape=~trt)))  
      
    } 
    
  }
  
  else if(parametric_model == "spline"){
    
    if(excess_hazard == F){
      
      flexmod <- try(flexsurvspline(Surv(time, event)~trt, 
                                    data = trial_data, 
                                    k = 3,
                                    bknots =  c(-3.3, 1),
                                    knots = c(log(1), log(1.5), log(2)),
                                    spline = "splines2ns",
                                    anc = list(gamma1=~trt,
                                               gamma2=~trt,
                                               gamma3=~trt,
                                               gamma4=~trt)))
      
    } else {
      
      flexmod <- try(flexsurvspline(Surv(time, event)~trt, 
                                    data = trial_data, 
                                    bhazard = backhaz,
                                    k = 3,
                                    spline = "splines2ns",
                                    anc = list(gamma1=~trt,
                                               gamma2=~trt,
                                               gamma3=~trt,
                                               gamma4=~trt)))
      
      
    }
    
  }
  
  mod_run_status <- inherits(flexmod, "try-error")  
  
  if(!mod_run_status){
    
    if(excess_hazard == F){
      
      est_rmst <- standsurv(flexmod, 
                            type = "rmst", 
                            at = list(list(trt=as.factor(0)), list(trt=as.factor(1))),
                            t = rmst_vec, 
                            ci = TRUE, 
                            se = TRUE, 
                            boot = F)
      
      est_rmst_control <- est_rmst %>%
        rename(value = at1, 
               lower = at1_lci, 
               upper = at1_uci,
               sd = at1_se,
               t = time) %>%
        mutate(estimand = "rmst", 
               trt = factor(0, levels = c(0,1))) %>%
        select(estimand, trt, t, value, lower, upper, sd)
      
      est_rmst_active <- est_rmst %>%
        rename(value = at2, 
               lower = at2_lci, 
               upper = at2_uci,
               sd = at2_se,
               t = time) %>%
        mutate(estimand = "rmst", 
               trt = factor(1, levels = c(0,1))) %>%
        select(estimand, trt, t, value, lower, upper, sd)
      
      est_irmst <- standsurv(flexmod, 
                             type = "rmst", 
                             at = list(list(trt=as.factor(0)), list(trt=as.factor(1))),
                             t = rmst_vec, 
                             contrast = "difference",
                             ci = TRUE, 
                             se = TRUE, 
                             boot = F) %>%
        rename(value = contrast2_1, 
               lower = contrast2_1_lci, 
               upper = contrast2_1_uci,
               sd = contrast2_1_se,
               t = time) %>%
        mutate(estimand = "irmst", trt = NA) %>%
        select(estimand, trt, t, value, lower, upper, sd)
      
      est_survival <- standsurv(flexmod, 
                                type = "survival", 
                                at = list(list(trt=as.factor(0)), list(trt=as.factor(1))),
                                t = t_vec, 
                                ci = TRUE, 
                                se = TRUE, 
                                boot = F)
      
      est_survival_control <- est_survival %>%
        rename(value = at1, 
               lower = at1_lci, 
               upper = at1_uci,
               sd = at1_se,
               t = time) %>%
        mutate(estimand = "survival", 
               trt = factor(0, levels = c(0,1))) %>%
        select(estimand, trt, t, value, lower, upper, sd)
      
      
      est_survival_active <- est_survival %>%
        rename(value = at2, 
               lower = at2_lci, 
               upper = at2_uci,
               sd = at2_se,
               t = time) %>%
        mutate(estimand = "survival", 
               trt = factor(1, levels = c(0,1))) %>%
        select(estimand, trt, t, value, lower, upper, sd)
      
      
      est_hazard <- standsurv(flexmod, 
                              type = "hazard", 
                              at = list(list(trt=as.factor(0)), list(trt=as.factor(1))),
                              t = t_vec, 
                              ci = TRUE, 
                              se = TRUE, 
                              boot = F)
      
      est_hazard_control <- est_hazard %>%
        rename(value = at1, 
               lower = at1_lci, 
               upper = at1_uci,
               sd = at1_se,
               t = time) %>%
        mutate(estimand = "hazard", 
               trt = factor(0, levels = c(0,1))) %>%
        select(estimand, trt, t, value, lower, upper, sd)
      
      est_hazard_active <- est_hazard %>%
        rename(value = at2, 
               lower = at2_lci, 
               upper = at2_uci,
               sd = at2_se,
               t = time) %>%
        mutate(estimand = "hazard", 
               trt = factor(1, levels = c(0,1))) %>%
        select(estimand, trt, t, value, lower, upper, sd)
      
      est_hr <- standsurv(flexmod, 
                          type = "hazard", 
                          contrast = "ratio",
                          at = list(list(trt=as.factor(0)), list(trt=as.factor(1))),
                          t = t_vec, 
                          ci = TRUE, 
                          se = TRUE, 
                          boot = F) %>%
        rename(value = contrast2_1, 
               lower = contrast2_1_lci, 
               upper = contrast2_1_uci,
               sd = contrast2_1_se,
               t = time) %>%
        mutate(estimand = "hr", trt = NA) %>%
        select(estimand, trt, t, value, lower, upper, sd) 
      
      
    } else {
      
      est_rmst_control <- 
        standsurv(flexmod, 
                  type = "rmst",
                  at = list(list(trt=as.factor(0))),
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
                    filter(trt == 0) %>%
                    mutate(sex = factor("male")) %>%
                    mutate(diag = as.Date("01/01/2000", "%d/%m/%Y")) %>%
                    mutate(agedays = floor(age * 365.25))%>%
                    select(trt, agedays, sex, diag)) %>%
        rename(value = at1, 
               lower = at1_lci, 
               upper = at1_uci,
               sd = at1_se,
               t = time) %>%
        mutate(estimand = "rmst", 
               trt = factor(0, levels = c(0,1))) %>%
        select(estimand, trt, t, value, lower, upper, sd)
      
      
      est_rmst_active <- 
        standsurv(flexmod, 
                  type = "rmst",
                  at = list(list(trt=as.factor(1))),
                  t = rmst_vec, 
                  ci = TRUE, 
                  se = TRUE, 
                  boot = F,
                  ratetable = rate_table_gompertz,
                  rmap = list(sex = sex,
                              year = diag,
                              age = agedays),
                  scale.ratetable = 365.25,
                  newdata = trial_data %>%
                    filter(trt == 1) %>%
                    mutate(sex = factor("male")) %>%
                    mutate(diag = as.Date("01/01/2000", "%d/%m/%Y")) %>%
                    mutate(agedays = floor(age * 365.25))%>%
                    select(trt, agedays, sex, diag)) %>%
        rename(value = at1, 
               lower = at1_lci, 
               upper = at1_uci,
               sd = at1_se,
               t = time) %>%
        mutate(estimand = "rmst", 
               trt = factor(1, levels = c(0,1))) %>%
        select(estimand, trt, t, value, lower, upper, sd)
      
      est_irmst <- standsurv(flexmod, 
                             type = "rmst", 
                             at = list(list(trt=as.factor(0)), list(trt=as.factor(1))),
                             t = rmst_vec, 
                             contrast = "difference",
                             ci = TRUE, 
                             se = TRUE, 
                             boot = F,
                             ratetable = rate_table_gompertz,
                             rmap = list(sex = sex,
                                         year = diag,
                                         age = agedays),
                             scale.ratetable = 365.25,
                             newdata = trial_data %>%
                               mutate(sex = factor("male")) %>%
                               mutate(diag = as.Date("01/01/2000", "%d/%m/%Y")) %>%
                               mutate(agedays = floor(age * 365.25)) %>%
                               select(trt, agedays, sex, diag)) %>%
        rename(value = contrast2_1, 
               lower = contrast2_1_lci, 
               upper = contrast2_1_uci,
               sd = contrast2_1_se,
               t = time) %>%
        mutate(estimand = "irmst", trt = NA) %>%
        select(estimand, trt, t, value, lower, upper, sd)
      
      
      est_survival_control <- 
        standsurv(flexmod, 
                  type = "survival",
                  at = list(list(trt=as.factor(0))),
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
                    filter(trt == 0) %>%
                    mutate(sex = factor("male")) %>%
                    mutate(diag = as.Date("01/01/2000", "%d/%m/%Y")) %>%
                    mutate(agedays = floor(age * 365.25))%>%
                    select(trt, agedays, sex, diag)) %>%
        rename(value = at1, 
               lower = at1_lci, 
               upper = at1_uci,
               sd = at1_se,
               t = time) %>%
        mutate(estimand = "survival", 
               trt = factor(0, levels = c(0,1))) %>%
        select(estimand, trt, t, value, lower, upper, sd)
      
      est_survival_active <- 
        standsurv(flexmod, 
                  type = "survival",
                  at = list(list(trt=as.factor(1))),
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
                    filter(trt == 1) %>%
                    mutate(sex = factor("male")) %>%
                    mutate(diag = as.Date("01/01/2000", "%d/%m/%Y")) %>%
                    mutate(agedays = floor(age * 365.25))%>%
                    select(trt, agedays, sex, diag)) %>%
        rename(value = at1, 
               lower = at1_lci, 
               upper = at1_uci,
               sd = at1_se,
               t = time) %>%
        mutate(estimand = "survival", 
               trt = factor(1, levels = c(0,1))) %>%
        select(estimand, trt, t, value, lower, upper, sd)
      
      
      est_hazard_control <- 
        standsurv(flexmod, 
                  type = "hazard",
                  at = list(list(trt=as.factor(0))),
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
                    filter(trt == 0) %>%
                    mutate(sex = factor("male")) %>%
                    mutate(diag = as.Date("01/01/2000", "%d/%m/%Y")) %>%
                    mutate(agedays = floor(age * 365.25))%>%
                    select(trt, agedays, sex, diag)) %>%
        rename(value = at1, 
               lower = at1_lci, 
               upper = at1_uci,
               sd = at1_se,
               t = time) %>%
        mutate(estimand = "hazard", 
               trt = factor(0, levels = c(0,1))) %>%
        select(estimand, trt, t, value, lower, upper, sd)
      
      est_hazard_active <- 
        standsurv(flexmod, 
                  type = "hazard",
                  at = list(list(trt=as.factor(1))),
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
                    filter(trt == 1) %>%
                    mutate(sex = factor("male")) %>%
                    mutate(diag = as.Date("01/01/2000", "%d/%m/%Y")) %>%
                    mutate(agedays = floor(age * 365.25))%>%
                    select(trt, agedays, sex, diag)) %>%
        rename(value = at1, 
               lower = at1_lci, 
               upper = at1_uci,
               sd = at1_se,
               t = time) %>%
        mutate(estimand = "hazard", 
               trt = factor(1, levels = c(0,1))) %>%
        select(estimand, trt, t, value, lower, upper, sd)
      
      est_hr <- standsurv(flexmod, 
                          type = "hazard", 
                          contrast = "ratio",
                          at = list(list(trt=as.factor(0)), list(trt=as.factor(1))),
                          t = t_vec, 
                          ci = TRUE, 
                          se = TRUE, 
                          boot = F,
                          ratetable = rate_table_gompertz,
                          rmap = list(sex = sex,
                                      year = diag,
                                      age = agedays),
                          scale.ratetable = 365.25,
                          newdata = trial_data %>%
                            mutate(sex = factor("male")) %>%
                            mutate(diag = as.Date("01/01/2000", "%d/%m/%Y")) %>%
                            mutate(agedays = floor(age * 365.25)) %>%
                            select(trt, agedays, sex, diag)) %>%
        rename(value = contrast2_1, 
               lower = contrast2_1_lci, 
               upper = contrast2_1_uci,
               sd = contrast2_1_se,
               t = time) %>%
        mutate(estimand = "hr", trt = NA) %>%
        select(estimand, trt, t, value, lower, upper, sd) 
      
      
    }
    
    est <- 
      bind_rows(est_rmst_control,
                est_rmst_active,
                est_irmst,
                est_survival_control,
                est_survival_active,
                est_hazard_control,
                est_hazard_active,
                est_hr)
    
    saveRDS(est, save_file)
    
  }
  
} 


