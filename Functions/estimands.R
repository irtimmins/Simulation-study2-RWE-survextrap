
# Derive estimates of hazard, survival and RMST for a single arm from a survextrap model.

estimates <- function(fitted_model,
                      t,
                      rmst_t, 
                      newdata, 
                      excess_hazard = FALSE,
                      lambda_gpm = NULL,
                      gamma_gpm = NULL,
                      trial_data = NULL,
                      wane_period = NULL,
                      wane_nt = NULL,
                      newdata0 = NULL){
  
  summ_fns <- list(median = ~median(.x, na.rm = T),
                   ~quantile(.x, probs=c(0.025, 0.975), na.rm = T),
                   sd = ~sd(.x, na.rm = T))
  
  if(excess_hazard == FALSE){
    
    haz_est <- hazard(fitted_model, t=t,  summ_fns=summ_fns,
                      wane_period = wane_period, wane_nt = wane_nt,
                      newdata = newdata, newdata0= newdata0) %>%
      mutate(estimand = "hazard") %>%
      rename(value=median, lower =`2.5%`, upper=`97.5%`) %>%
      relocate(estimand) 
    
    surv_est <- survival(fitted_model, t=t, summ_fns=summ_fns,
                         wane_period = wane_period, wane_nt = wane_nt,
                         newdata = newdata, newdata0= newdata0) %>%
      mutate(estimand = "survival") %>%
      rename(value=median, lower =`2.5%`, upper=`97.5%`) %>%
      relocate(estimand) 
    
  } else if(excess_hazard == TRUE) {
    
    trial_data <- trial_data %>%
      filter(trt %in% newdata[["trt"]])
    
    # Compute marginal background hazard and expected survival rates.
    backhaz_rates <- expand_grid(t = t,
                                 age = trial_data[["age"]]) %>%
      mutate(backhaz = hgompertz(x = age+t,  shape = gamma_gpm, rate = lambda_gpm)) %>%
      mutate(expsurv = (1-pgompertz(q = age+t,  shape = gamma_gpm, rate = lambda_gpm))/
               (1-pgompertz(q = age,  shape = gamma_gpm, rate = lambda_gpm))) %>%
      group_by(t) %>% 
      summarise(backhaz = sum(backhaz*expsurv)/sum(expsurv),
                expsurv = mean(expsurv)) %>%
      ungroup()
    
    haz_est <- hazard(fitted_model, t=t, summ_fns=summ_fns,
                      wane_period = wane_period, wane_nt = wane_nt,
                      newdata = newdata, newdata0= newdata0) %>%
      mutate(estimand = "hazard") %>%
      relocate(estimand) %>%
      rename(value=median, lower =`2.5%`, upper=`97.5%`) %>%
      left_join(backhaz_rates, join_by("t")) %>%
      mutate(across(c(value, lower, upper), ~ .x + backhaz)) %>%
      select(-c(backhaz, expsurv)) 
    
    surv_est <- survival(fitted_model, t=t, summ_fns=summ_fns,
                         wane_period = wane_period, wane_nt = wane_nt,
                         newdata = newdata, newdata0= newdata0) %>%
      mutate(estimand = "survival") %>%
      relocate(estimand) %>%
      rename(value=median, lower =`2.5%`, upper=`97.5%`) %>%
      left_join(backhaz_rates, join_by("t")) %>%
      mutate(across(c(value, lower, upper, sd), ~ .x*expsurv)) %>%
      select(-c(backhaz, expsurv)) 
    
  } 
  
  rmst_est <- estimates_rmst(fitted_model = fitted_model,
                                  rmst_t = rmst_t, 
                                  newdata = newdata,
                                  excess_hazard = excess_hazard,
                                  lambda_gpm = lambda_gpm,
                                  gamma_gpm = gamma_gpm,
                                  trial_data = trial_data,
                                  wane_period = wane_period,
                                  wane_nt = wane_nt,
                                  newdata0 = newdata0) 
  
  est_df <- haz_est %>% 
    bind_rows(surv_est, 
              rmst_est) 
  
  return(est_df)
  
}



# Derive rmst estimands using Gauss-Legendre approximation for speed.

estimates_rmst <- function(fitted_model,
                                rmst_t, 
                                newdata,
                                excess_hazard = FALSE,
                                lambda_gpm = NULL,
                                gamma_gpm = NULL,
                                trial_data,
                                wane_period = NULL,
                                wane_nt = NULL,
                                newdata0 = NULL){
  
  summ_fns <- list(median = ~median(.x, na.rm = T),
                   ~quantile(.x, probs=c(0.025, 0.975), na.rm = T),
                   sd = ~sd(.x, na.rm = T))
  
  # previous version used survextrap rmst function.
  #  rmst_est <- rmst(fitted_model, t=rmst_t, newdata = newdata, summ_fns=summ_fns) %>%
  #    rename(value=median, lower =`2.5%`, upper=`97.5%`) %>%
  #    select(trt, t, value, lower, upper, sd) %>%
  #    mutate(estimand = "rmst") %>%
  #    relocate(estimand) 
  
  # RMST loops over timepoints.
  
  rmst_est <- NULL
  
  for(time_rmst in rmst_t){
    
    # Get the RMST MCMC iterations to extract quantiles from.
    
    rmst_samples <- rmst_samples_GL(fitted_model, 
                                         time_rmst = time_rmst,
                                         newdata = newdata,
                                         excess_hazard = excess_hazard,
                                         lambda_gpm = lambda_gpm,
                                         gamma_gpm = gamma_gpm,
                                         trial_data = trial_data,
                                         wane_period = wane_period, 
                                         wane_nt = wane_nt,
                                         newdata0 = newdata0,
                                         nodes = 100)   
    
    rmst_quantiles <- quantile(rmst_samples, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
    
    rmst_est_temp <- newdata %>%
      mutate(estimand = "rmst") %>%
      mutate(t = time_rmst, 
             value = rmst_quantiles[["50%"]],
             lower = rmst_quantiles[["2.5%"]],
             upper = rmst_quantiles[["97.5%"]],
             sd = sd(rmst_samples, na.rm= TRUE)) %>%
      relocate(estimand)
    
    if(is.null(rmst_est)){
      
      rmst_est <- rmst_est_temp
      
    } else{
      
      rmst_est <- bind_rows(rmst_est,
                            rmst_est_temp)
    }
    
  }
  
  #  }
  
  
  return(rmst_est)
  
}


# rmst MCMC samples using Gauss-Legendre integration on the survival function.

rmst_samples_GL <- function(fitted_model, 
                            time_rmst,
                            newdata,
                            excess_hazard = FALSE,
                            lambda_gpm = NULL,
                            gamma_gpm = NULL,
                            trial_data,
                            wane_period = NULL, 
                            wane_nt = NULL,
                            newdata0= NULL,
                            nodes = 100){
  
  GL <- gaussLegendre(n = nodes, 0, time_rmst)
  GL_times <- GL$x
  GL_weights <- GL$w
  
  relative_survival_chains <- survival(fitted_model, t = GL_times, 
                                       newdata = newdata,                                    
                                       wane_period = wane_period, 
                                       wane_nt = wane_nt,
                                       newdata0= newdata0,
                                       sample = TRUE)[,,]
  N <- nrow(trial_data)
  e_vec <- matrix(1, ncol = N, nrow = 1) # vector of ones
  W <- diag(GL_weights)
  
  if(excess_hazard == FALSE){
    
    expected_survival_grid <- expand_grid(time = GL_times,
                                          age = trial_data[["age"]]) %>%
      mutate(expected_survival = 1) %>%
      pull(expected_survival) %>%
      matrix(nrow = N)
    
  }else if(excess_hazard == TRUE){
    
    expected_survival_grid <- expand_grid(time = GL_times,
                                          age = trial_data[["age"]]) %>%
      mutate(expected_survival = (1-pgompertz(q = age+time,  shape = gamma_gpm, rate = lambda_gpm))/
               (1-pgompertz(q = age,  shape = gamma_gpm, rate = lambda_gpm))) %>%
      pull(expected_survival) %>%
      matrix(nrow = N)
    
  }
  
  rmst_samples <- (1/N) * e_vec %*% expected_survival_grid %*% W %*% relative_survival_chains %>% 
    as.vector()
  
  return(rmst_samples)
  
}

# irmst for PH/non-PH models with treatment covariate.

estimates_irmst <- function(fitted_model,
                            rmst_t, 
                            newdata,
                            excess_hazard = FALSE,
                            lambda_gpm = NULL,
                            gamma_gpm = NULL,
                            trial_data,
                            wane_period = NULL,
                            wane_nt = NULL,
                            newdata0 = NULL){
  
  summ_fns <- list(median = ~median(.x, na.rm = T),
                   ~quantile(.x, probs=c(0.025, 0.975), na.rm = T),
                   sd = ~sd(.x, na.rm = T))
  
  #  irmst_est <-  irmst(fitted_model, t=rmst_t, newdata = newdata, summ_fns=summ_fns) %>%
  #     mutate(trt = NA) %>%
  #      rename(value=median, lower =`2.5%`, upper=`97.5%`) %>%
  #    select(trt, t, value, lower, upper, sd) %>%
  #    mutate(estimand = "irmst") %>%
  #    relocate(estimand)  
  # 
  irmst_est <- NULL
  
  for(time_rmst in rmst_t){
    
    rmst_samples_control <- rmst_samples_GL(fitted_model, 
                                                 time_rmst = time_rmst,
                                                 newdata = newdata %>%
                                                   filter(trt == 0),
                                                 excess_hazard = excess_hazard,
                                                 lambda_gpm = lambda_gpm,
                                                 gamma_gpm = gamma_gpm,
                                                 trial_data = trial_data %>%
                                                   filter(trt == 0),
                                                 nodes = 100)    
    
    rmst_samples_active <- rmst_samples_GL(fitted_model, 
                                                time_rmst = time_rmst,
                                                newdata = newdata %>%
                                                  filter(trt == 1),
                                                excess_hazard = excess_hazard,
                                                lambda_gpm = lambda_gpm,
                                                gamma_gpm = gamma_gpm,
                                                trial_data = trial_data %>%
                                                  filter(trt == 1),
                                                wane_period = wane_period,
                                                wane_nt = wane_nt,
                                                newdata0 = newdata0,
                                                nodes = 100)    
    
    irmst_samples <- rmst_samples_active - rmst_samples_control

    irmst_quantiles <- quantile(irmst_samples, probs = c(0.025, 0.5, 0.975), na.rm= TRUE)
    
    irmst_est_temp <- tibble(trt = NA) %>%
      mutate(estimand = "irmst") %>%
      mutate(t = time_rmst, 
             value = irmst_quantiles[["50%"]],
             lower = irmst_quantiles[["2.5%"]],
             upper = irmst_quantiles[["97.5%"]],
             sd = sd(irmst_samples, na.rm= TRUE)) %>%
      relocate(estimand)
    
    if(is.null(irmst_est)){
      
      irmst_est <- irmst_est_temp
      
    } else{
      
      irmst_est <- bind_rows(irmst_est,
                             irmst_est_temp)
    }
    
  }
  
  
  
  return(irmst_est)
  
}

# irmst estimates for separately modeled arms.

estimates_irmst_separate <- function(control_model,
                                     active_model,
                                     rmst_t, 
                                     newdata = tibble(trt = as.factor(c(0,1))), 
                                     excess_hazard = FALSE,
                                     lambda_gpm = NULL,
                                     gamma_gpm = NULL,
                                     trial_data = NULL){
  

  summ_fns <- list(median = ~median(.x, na.rm = T),
                   ~quantile(.x, probs=c(0.025, 0.975), na.rm = T),
                   sd = ~sd(.x, na.rm = T))
  
  # if(excess_hazard == FALSE){
  #   
  #   irmst_samples <- rmst(active_model, t=rmst_t, sample=TRUE) - 
  #     rmst(control_model, t=rmst_t, sample=TRUE)
  #   
  #   irmst_est <- irmst_samples %>%
  #     aperm(c(2,3,1)) %>% # reorder as iterations, chains, variables
  #     as_draws_array() %>% # convert to format for summarise_draws
  #     posterior::summarise_draws(median,
  #                                ~quantile(.x, probs=c(0.025, 0.975)),
  #                                sd) %>%
  #     rename(estimand = variable, value = median, lower="2.5%", upper="97.5%") %>% 
  #     mutate(estimand = "irmst", trt = NA) %>%
  #     mutate(t = rmst_t) %>%
  #     relocate(estimand, trt, t)
  #   
  
  #} else {
  
  irmst_est <- NULL
  
  for(time_rmst in rmst_t){
    
    rmst_samples_control <- rmst_samples_GL(control_model, 
                                            time_rmst = time_rmst,
                                            newdata = newdata %>%
                                              filter(trt == 0),
                                            excess_hazard = excess_hazard,
                                            lambda_gpm = lambda_gpm,
                                            gamma_gpm = gamma_gpm,
                                            trial_data = trial_data %>%
                                              filter(trt == 0),
                                            nodes = 100)    
    
    rmst_samples_active <- rmst_samples_GL(active_model, 
                                           time_rmst = time_rmst,
                                           newdata = newdata %>%
                                             filter(trt == 1),
                                           excess_hazard = excess_hazard,
                                           lambda_gpm = lambda_gpm,
                                           gamma_gpm = gamma_gpm,
                                           trial_data = trial_data %>%
                                             filter(trt == 1),
                                           nodes = 100)    
    
    # Make sure separate models have same number of iterations.
    # otherwise return NAs.
    
    if(length(rmst_samples_active) == length(rmst_samples_control)){
      
      irmst_samples <- rmst_samples_active - rmst_samples_control
      
      irmst_quantiles <- quantile(irmst_samples, probs = c(0.025, 0.5, 0.975),  na.rm= TRUE)
      
      irmst_est_temp <- tibble(trt = NA) %>%
        mutate(estimand = "irmst") %>%
        mutate(t = time_rmst, 
               value = irmst_quantiles[["50%"]],
               lower = irmst_quantiles[["2.5%"]],
               upper = irmst_quantiles[["97.5%"]],
               sd = sd(irmst_samples, na.rm= TRUE)) %>%
        relocate(estimand)
      
    } else{
      
      irmst_est_temp <- tibble(trt = NA) %>%
        mutate(estimand = "irmst") %>%
        mutate(t = time_rmst, 
               value = NA,
               lower = NA,
               upper = NA,
               sd = NA) %>%
        relocate(estimand)
      
    }
    
    if(is.null(irmst_est)){
      
      irmst_est <- irmst_est_temp
      
    } else{
      
      irmst_est <- bind_rows(irmst_est,
                             irmst_est_temp)
    }
    
    #}
    
    
    
  }
  
  
  
  return(irmst_est)
  
}

# hazard ratio estimates for PH/non-PH models with treatment covariate.

estimates_hr <- function(fitted_model,
                         t, 
                         newdata =  tibble(trt = as.factor(c(0,1))),
                         excess_hazard = FALSE,
                         lambda_gpm = NULL,
                         gamma_gpm = NULL,
                         trial_data = NULL,
                         wane_period = NULL,
                         wane_nt = NULL,
                         newdata0 = NULL){
  
  summ_fns <- list(median = ~median(.x, na.rm = T),
                   ~quantile(.x, probs=c(0.025, 0.975), na.rm = T),
                   sd = ~sd(.x, na.rm = T))
  
  if(excess_hazard == FALSE){
    
    # If no standardisation needed, then use survextrap hazard_ratio function.
    hr_est <- hazard_ratio(fitted_model, 
                           t=t, 
                           newdata = newdata, 
                           summ_fns=summ_fns,  
                           wane_period = wane_period,
                           wane_nt = wane_nt,
                           newdata0 = newdata0) %>%
      mutate(trt = NA) %>%
      rename(value=median, lower =`2.5%`, upper=`97.5%`) %>%
      select(trt, t, value, lower, upper, sd) %>%
      mutate(estimand = "irmst") %>%
      relocate(estimand) 
    
  } else if(excess_hazard == TRUE){
    
    hazard_active_samples <-  hazard(fitted_model, 
                                     t=t, 
                                     newdata = newdata %>%
                                       filter(trt == 1),
                                     sample = TRUE,
                                     wane_period = wane_period,
                                     wane_nt = wane_nt,
                                     newdata0 = newdata0) 
    
    hazard_control_samples <-  hazard(fitted_model, 
                                      t=t, 
                                      newdata = newdata %>%
                                        filter(trt == 0),
                                      sample = TRUE)
    
    hazard_active_draws <- hazard_active_samples %>%
      aperm(c(2,3,1)) %>% # reorder as iterations, chains, variables
      as_draws_array()
    
    hazard_control_draws <- hazard_control_samples %>%
      aperm(c(2,3,1)) %>% # reorder as iterations, chains, variables
      as_draws_array()
    
    active_trial_data <- trial_data %>%
      filter(trt == 1)
    
    control_trial_data <- trial_data %>%
      filter(trt == 0)
    
    backhaz_rates_active <- expand_grid("t" = t,
                                        age = active_trial_data[["age"]]) %>%
      mutate(backhaz = hgompertz(x = age+t,  shape = gamma_gpm, rate = lambda_gpm)) %>%
      mutate(expsurv = (1-pgompertz(q = age+t,  shape = gamma_gpm, rate = lambda_gpm))/
               (1-pgompertz(q = age,  shape = gamma_gpm, rate = lambda_gpm))) %>%
      group_by(t) %>%
      summarise(backhaz = sum(backhaz*expsurv)/sum(expsurv),
                expsurv = mean(expsurv)) %>%
      ungroup()
    
    backhaz_rates_control <- expand_grid("t" = t,
                                         age = control_trial_data[["age"]]) %>%
      mutate(backhaz = hgompertz(x = age+t,  shape = gamma_gpm, rate = lambda_gpm)) %>%
      mutate(expsurv = (1-pgompertz(q = age+t,  shape = gamma_gpm, rate = lambda_gpm))/
               (1-pgompertz(q = age,  shape = gamma_gpm, rate = lambda_gpm))) %>%
      group_by(t) %>%
      # marginal excess hazards are weighted by expected survival
      summarise(backhaz = sum(backhaz*expsurv)/sum(expsurv), 
                expsurv = mean(expsurv)) %>%
      ungroup()
    
    # Add excess hazard rates to get all-cause hazard.
    for(i in 1:length(t)) {
      hazard_active_draws[,,i] <- hazard_active_draws[,,i] + backhaz_rates_active$backhaz[i] 
      hazard_control_draws[,,i] <- hazard_control_draws[,,i] + backhaz_rates_control$backhaz[i] 
    }
    
    hr_samples <- hazard_active_draws / hazard_control_draws 
    
    hr_est <- hr_samples %>% # reorder as iterations, chains, variables
      as_draws_array() %>% # convert to format for summarise_draws
      posterior::summarise_draws(median = ~median(.x, na.rm = T),
                                 ~quantile(.x, probs=c(0.025, 0.975), na.rm = T),
                                 sd = ~sd(.x, na.rm = T)) %>%
      rename(estimand = variable, value = median, lower="2.5%", upper="97.5%") %>% 
      mutate(estimand = "hr", trt = NA) %>%
      mutate("t" = t) %>%
      relocate(estimand, trt, t)
    
  }
  
  return(hr_est)
  
}

# hazard ratio estimates for separately modeled arms.

estimates_hr_separate <- function(control_model,
                                  active_model,
                                  t,
                                  newdata = tibble(trt = as.factor(c(0,1))), 
                                  excess_hazard = FALSE,
                                  lambda_gpm = NULL,
                                  gamma_gpm = NULL,
                                  trial_data = NULL){
  
  summ_fns <- list(median = ~median(.x, na.rm = T),
                   ~quantile(.x, probs=c(0.025, 0.975), na.rm = T),
                   sd = ~sd(.x, na.rm = T))
  
  hazard_active_samples <-  hazard(active_model, 
                                   t=t, 
                                   sample = TRUE) 
  
  hazard_control_samples <-  hazard(control_model, 
                                    t=t, 
                                    sample = TRUE)
  
  hazard_active_draws <- hazard_active_samples %>%
    aperm(c(2,3,1)) %>% # reorder as iterations, chains, variables
    as_draws_array()
  
  hazard_control_draws <- hazard_control_samples %>%
    aperm(c(2,3,1)) %>% # reorder as iterations, chains, variables
    as_draws_array()
  
 if (excess_hazard == TRUE){
    
    # For excess hazard models we have extra step
    # of needing to standardise over trial arms by patient age.
    
    # trial data on active arm and control arm
    active_trial_data <- trial_data %>%
      filter(trt == 1)
    
    control_trial_data <- trial_data %>%
      filter(trt == 0)
    
    # calculate other-cause hazard at each time point for each trial patient, 
    # and standardised over patient age.
    backhaz_rates_active <- expand_grid("t" = t,
                                        age = active_trial_data[["age"]]) %>%
      mutate(backhaz = hgompertz(x = age+t,  shape = gamma_gpm, rate = lambda_gpm)) %>%
      mutate(expsurv = (1-pgompertz(q = age+t,  shape = gamma_gpm, rate = lambda_gpm))/
               (1-pgompertz(q = age,  shape = gamma_gpm, rate = lambda_gpm))) %>%
      group_by(t) %>%
      summarise(backhaz = sum(backhaz*expsurv)/sum(expsurv),
                expsurv = mean(expsurv)) %>%
      ungroup()
    
    backhaz_rates_control <- expand_grid("t" = t,
                                         age = control_trial_data[["age"]]) %>%
      mutate(backhaz = hgompertz(x = age+t,  shape = gamma_gpm, rate = lambda_gpm)) %>%
      mutate(expsurv = (1-pgompertz(q = age+t,  shape = gamma_gpm, rate = lambda_gpm))/
               (1-pgompertz(q = age,  shape = gamma_gpm, rate = lambda_gpm))) %>%
      group_by(t) %>%
      summarise(backhaz = sum(backhaz*expsurv)/sum(expsurv),
                expsurv = mean(expsurv)) %>%
      ungroup()
    
    # Add background hazard to each MCMC draw.
    # Loop over timepoint vector.
    for(i in 1:length(t)) {
      hazard_active_draws[,,i] <- hazard_active_draws[,,i] + backhaz_rates_active$backhaz[i] 
      hazard_control_draws[,,i] <- hazard_control_draws[,,i] + backhaz_rates_control$backhaz[i] 
    }
    

  }
  
  if(all(dim(hazard_active_draws) == dim(hazard_control_draws))){
    
    hr_samples <- hazard_active_draws / hazard_control_draws 
    
  } else{
    # If dimensions misaligned then give NAs.   
    hr_samples <- hazard_active_draws
    hr_samples[,,] <- NA 
  }
  
  hr_est <- hr_samples %>% # reorder as iterations, chains, variables
    as_draws_array() %>% # convert to format for summarise_draws
    posterior::summarise_draws(median = ~median(.x, na.rm = T),
                               ~quantile(.x, probs=c(0.025, 0.975), na.rm = T),
                               sd = ~sd(.x, na.rm = T)) %>%
    rename(estimand = variable, value = median, lower="2.5%", upper="97.5%") %>% 
    mutate(estimand = "hr", trt = NA) %>%
    mutate("t" = t) %>%
    relocate(estimand, trt, t)
  
  return(hr_est)
  
}





# Functions to extract true values.

survival_big <- function(big_data, t, trt_group){
  
  res <- big_data %>%
    filter(trt == trt_group) %>%
    mutate(survival = 1*(time > t)) %>%
    summarise(value = mean(survival),
              value_se = sqrt(value*(1-value)/n())) #%>%
   # relocate(trt)
  
  return(res)
  
}


rmst_big <- function(big_data, t, trt_group){
  
  res <- big_data %>%
    filter(trt == trt_group) %>%
    mutate(truncT = pmin(time, t)) %>%
    summarise(value = mean(truncT),
              value_se = sd(truncT)/sqrt(n())) #%>%
  #  relocate(trt)
  
  return(res)
  
}


irmst_big <- function(big_data, t){
  
  res <- big_data %>%
    mutate(truncT = pmin(time, t)) %>%
    group_by(trt) %>%
    summarise(rmst_mean = mean(truncT),
              rmst_var = sd(truncT)^2/n()) %>%
    mutate(trt_sense = if_else(trt == 1, 1, -1)) %>%
    summarise(value = sum(rmst_mean*trt_sense),
              value_se = sqrt(sum(rmst_var))) 
  
  return(res)
  
}


# Function to combine trial datasets into big data file (e.g. N = 10^7).

dgm_combine <- function(nsim, true_model_id, save_file){
  
  for(j in 1:nsim){
    
    file_temp <- paste0("../dgm_", true_model_id, "/data", j, ".rds")
    
    if(file.exists(file_temp)){
      if(j == 1){ 
        comb <- readRDS(file_temp)
      } else {
        temp <- readRDS(file_temp)
        comb <- rbindlist(list(comb, temp))
      }
      if((j %% 10) == 0) print(paste0(j, "/", nsim, " ", true_model_id))
    }
  }
  saveRDS(comb, save_file)
}



# calculates true survival times from a large trial dataset.

dgm_surv_true <- function(true_model_id, design_id, big_data_file, save_file){
  
  big_data <- readRDS(big_data_file)
  
  true_surv_df <- estimand_labels %>%
    filter(estimand == "survival")
  
  if(design_id == "single_arm"){
    
    true_surv_df <- true_surv_df  %>%
      filter(trt == 0)
  }
  
  surv_comb <- NULL
  
  for(j in 1:nrow(true_surv_df)){
    
    id <- true_surv_df %>%
      select(estimand_id) %>%
      slice(j) %>%
      pull()
    
    t <- true_surv_df %>%
      select(t) %>%
      slice(j) %>%
      pull()
    
    trt_j <- true_surv_df %>%
      select(trt) %>%
      slice(j) %>%
      pull()
    
    surv_temp <- survival_big(big_data, t = t, trt_group = trt_j)  %>%
      mutate(estimand_id = id,
             estimand = "survival",
             trt = trt_j) %>%
      relocate(estimand, estimand_id, trt)
    
    surv_comb <- surv_comb %>%
      bind_rows(surv_temp)
    
    if((j %% 10) == 0) print(paste0(j, "/", nrow(true_surv_df), ", ",true_model_id))
  }
  
  surv_comb <- as_tibble(surv_comb)
  saveRDS(surv_comb, save_file)
  
}


# calculates true rmst times from a large trial dataset.

dgm_rmst_true <- function(true_model_id, design_id, big_data_file, save_file){
  
  big_data <- readRDS(big_data_file)
  
  true_rmst_df <- estimand_labels %>%
    filter(estimand == "rmst")
  
  if(design_id == "single_arm"){
    
    true_rmst_df <- true_rmst_df  %>%
      filter(trt == 0)
  }
  
  rmst_comb <- NULL
  
  for(j in 1:nrow(true_rmst_df)){
    
    id <- true_rmst_df %>%
      select(estimand_id) %>%
      slice(j) %>%
      pull()
    
    t <- true_rmst_df %>%
      select(t) %>%
      slice(j) %>%
      pull()
    
    trt_j <- true_rmst_df %>%
      select(trt) %>%
      slice(j) %>%
      pull()
    
    rmst_temp <- rmst_big(big_data, t = t, trt_group = trt_j)  %>%
      mutate(estimand_id = id,
             estimand = "rmst",
             trt = trt_j) %>%
      relocate(estimand, estimand_id, trt)
    
    rmst_comb <- rmst_comb %>%
      bind_rows(rmst_temp)
    
    print(paste0(j, "/", nrow(true_rmst_df), ", ", true_model_id))
    
  }
  
  rmst_comb <- as_tibble(rmst_comb)
  saveRDS(rmst_comb, save_file)

  
}


# calculates true irmst times from a large trial dataset.

dgm_irmst_true <- function(true_model_id, design_id, big_data_file, save_file){
  
  big_data <- readRDS(big_data_file)
  
  true_irmst_df <- estimand_labels %>%
    filter(estimand == "irmst")
  
  irmst_comb <- NULL
  
  for(j in 1:nrow(true_irmst_df)){
    #  j <- 1
    id <- true_irmst_df %>%
      select(estimand_id) %>%
      slice(j) %>%
      pull()
    
    t <- true_irmst_df %>%
      select(t) %>%
      slice(j) %>%
      pull()
    
    irmst_temp <- irmst_big(big_data, t)  %>%
      mutate(estimand_id = id,
             estimand = "irmst",
             trt = NA) %>%
      relocate(estimand, estimand_id, trt)
    #head(irmst_temp)
    irmst_comb <- irmst_comb %>%
      bind_rows(irmst_temp)
    
    print(paste0(j, "/", nrow(true_irmst_df), ", ", true_model_id))
    
  }
  
  saveRDS(irmst_comb, save_file)
  
}

# Function to evaluate true all-cause hazard and excess/gpm hazards using
# analytical formulae where possible.

dgm_haz_true <- function(design_id, lambda1, lambda2, gamma1, gamma2, pmix, lambda_gpm, gamma_gpm, 
                         age_mean, age_sd, alpha_frailty, hr_function, beta, N_large, save_file){
  
  left_trunc <- 0
  
  if(is.character(beta) & !is.na(beta)) beta <- get(beta)
  
  if(design_id == "single_arm"){
    
    beta <- 0
    hr_function <- "ph"
    
    estimand_labels <- estimand_labels %>%
      filter(trt == 0)
  }
  
  # Identify true (marginal) excess hazards.
  
  if(alpha_frailty == 0){
    
    true_excess_haz_df <- estimand_labels %>%
      filter(estimand == "excess_hazard") %>%
      mutate(frailty = 0) %>%
      mutate(value = haz(t, lambda1, lambda2, gamma1, gamma2, pmix, hr_function, beta, trt, alpha_frailty, frailty, left_trunc)) %>%
      select(-c(frailty)) %>%
      mutate(value_se = NA)
    
  } else{
    
    true_excess_haz_df <- estimand_labels %>%
      filter(estimand == "excess_hazard") %>%
      expand_grid(frailty = rnorm(N_large,0,1)) %>%
      mutate(value = haz(t, lambda1, lambda2, gamma1, gamma2, pmix, hr_function, beta, trt, alpha_frailty, frailty, left_trunc)) %>%
      group_by(estimand_id) %>%
      summarise(value = mean(value),
                count= n()) %>%
      left_join(estimand_labels, join_by(estimand_id)) %>%
      select(estimand, estimand_id, t, trt, value) %>%
      arrange(t) %>%
      mutate(value_se = NA)
    
  }
  
  # Identify true (marginal) background GPM hazards.
  # Do each time point separately to reduce memory.
  for(i in 1:nrow(estimand_labels)){
  
  true_gpm_haz_df_temp <- estimand_labels %>%
    filter(estimand == "gpm_hazard") %>%
    slice(i) %>%
    expand_grid(age = rnorm(N_large,age_mean,age_sd)) %>%
    mutate(backhaz = hgompertz(x = age+t,  shape = gamma_gpm, rate = lambda_gpm)) %>%
    mutate(expsurv = (1-pgompertz(q = age+t,  shape = gamma_gpm, rate = lambda_gpm))/
             (1-pgompertz(q = age,  shape = gamma_gpm, rate = lambda_gpm))) %>%
    group_by(estimand_id) %>%
    summarise(value = sum(backhaz*expsurv)/sum(expsurv),
              count = n()) %>%
    left_join(estimand_labels, join_by(estimand_id)) %>%
    select(estimand, estimand_id, t, trt, value) %>%
    arrange(t) %>%
    mutate(value_se = NA)
  
  if(i == 1){
    true_gpm_haz_df <- true_gpm_haz_df_temp 
  } else{
    true_gpm_haz_df <- true_gpm_haz_df %>%
      bind_rows(true_gpm_haz_df_temp)
  }
  
  if((i %% 100) == 0) print(paste0(i, "/", nrow(estimand_labels), " completed" ))
  
  }
  
  # Derive true all-cause hazards.
  
  true_haz_df <- estimand_labels %>%
    filter(estimand == "hazard") 
  
  if(!(all(true_haz_df$t == true_gpm_haz_df$t) & 
       all(true_gpm_haz_df$t == true_excess_haz_df$t) & 
       all(true_haz_df$trt == true_gpm_haz_df$trt) &
       all(true_gpm_haz_df$trt == true_excess_haz_df$trt))){
    stop("Check consistency of estimand timepoints")
  } 
  
  true_haz_df <- true_haz_df %>%
    mutate(value = true_excess_haz_df[["value"]]+
             true_gpm_haz_df[["value"]]) %>%
    mutate(value_se = NA)
  
  if(design_id == "two_arm"){
    
    true_hr_df <- 
      true_haz_df %>%
      separate_wider_delim(
        estimand_id,
        delim = "_",
        names = c("estimand_id_num", "trt_character")
      ) %>%
      separate_wider_regex(
        estimand_id_num,
        patterns = c(
          estimand_check = "[A-Za-z]+", 
          id_num = "[0-9]+"
        )
      ) %>%
      select(-c(estimand_check,trt))  %>%
      pivot_wider(names_from = trt_character, values_from = value) %>%
      mutate(value = trt1/trt0) %>%
      select(c(-trt0,trt1)) %>%
      mutate(estimand = "hr",
             estimand_id = paste0(estimand, id_num),
             trt = NA) %>%
      select(estimand, estimand_id, t, trt, value, value_se)
    
    true_hr_df_test <- estimand_labels %>%
      filter(estimand == "hr")  
    
    if(!(all(true_hr_df$estimand_id == true_hr_df_test$estimand_id) &
         all(true_hr_df$t == true_hr_df_test$t))){
      stop("Check consistency of HR estimand timepoints")
    }
    
  }
  
  est_df <-  true_haz_df %>%
    bind_rows(true_excess_haz_df) %>%
    bind_rows(true_gpm_haz_df)
  
  if(design_id == "two_arm"){
    
    est_df <- est_df %>%
      bind_rows(true_hr_df)
    
  }
  
  saveRDS(est_df, save_file)
  print(save_file)
}



# Function to approximate true (all-cause) hazard using splines.
 
dgm_haz_true_model <- function(big_data_file, design_id, N, df, return_model) {
  
  big_data <- readRDS(big_data_file)
  
  if(design_id == "single_arm"){
    trt_group_vec <- c(0)
  } else if(design_id == "two_arm"){
    trt_group_vec <- c(0,1)
  }
  
  #trt_group_vec
  
  for(trt_group in trt_group_vec){
    
    big_data_temp <- as_tibble(big_data) %>%
      filter(trt == trt_group) %>%
      slice(1:N) 
    
    rstpm_mod <- stpm2(Surv(time, event) ~ 1, 
                       data = big_data_temp, 
                       df = df)
    
    saveRDS(rstpm_mod, paste0(return_model, "_trt", trt_group, ".rds"))
    
  }
  
}


# Evaluate true (all-cause) hazard using splines.

dgm_haz_true_rstpm2 <- function(true_model_id, design_id, model_file){
  
  
  if(design_id == "single_arm"){
    trt_group_vec <- c(0)
  } else if(design_id == "two_arm"){
    trt_group_vec <- c(0,1)
  }
  
  
  for(trt_group in trt_group_vec){
    
    #trt_group <- 0
    
    rstpm_mod <- readRDS(paste0(model_file, "_trt", trt_group, ".rds"))
    
    est_haz_mod <- predict(rstpm_mod, 
                           type = "hazard", 
                           newdata = tibble(time = time_vec_mod),
                           se.fit = TRUE ) %>%
      as_tibble() %>%
      mutate(estimand = "hazard",
             trt = trt_group,
             t = time_vec_mod,
             value_se = (upper-lower)/(2*1.96)) %>%
      relocate(estimand, trt, t) %>%
      rename(value = Estimate) %>%
      relocate(value_se, .after= value)  %>%
      select(-c(lower, upper)) %>%
      filter(t > 0)
    
    saveRDS(est_haz_mod, paste0(model_file, "_trt", trt_group, "_hazard.rds"))
    
    est_surv_mod <- predict(rstpm_mod, 
                            type = "surv", 
                            newdata = tibble(time = time_vec_mod),
                            se.fit = TRUE ) %>%
      as_tibble() %>%
      mutate(estimand = "survival",
             trt = trt_group,
             t = time_vec_mod,
             value_se = (upper-lower)/(2*1.96)) %>%
      relocate(estimand, trt, t) %>%
      rename(value = Estimate) %>%
      relocate(value_se, .after= value)  %>%
      select(-c(lower, upper)) %>%
      filter(t > 0)
    
    saveRDS(est_surv_mod, paste0(model_file, "_trt", trt_group, "_survival.rds"))
    print(paste0(model_file, "_trt", trt_group ))
    
  }
  
  
}


# Function to combine true and simulated estimates for each scenario i.

combine_sim_and_true <- function(nsim, scenario_fit_id, 
                                 true_file, save_estimates_file, save_results_file){
  
  est <- NULL
  
  for(j in 1:nsim){
    #j <- 1
    # print progress track every 10th iteration
    if((j %% 10) == 0) print(paste0(j, "/", nsim, ", ", scenario_fit_id))
    
    est_file <- paste0(save_estimates_file,  j, ".rds")
    
    if(file.exists(est_file)){
      
      temp <- readRDS(est_file)
      #temp
      temp <- temp %>%
        mutate(trt = if_else(trt == 0, 0, 1)) %>% # convert from factor
        mutate("scenario_fit_id" = scenario_fit_id, 
               isim = j,
               model_type = "simulation")  %>%
        left_join(estimand_labels, by = c("estimand", "t", "trt")) %>%
        relocate(estimand_id, .after = estimand)
      
      est <-  rbindlist(list(est, temp))
      
    }
  }
  
  est <- as_tibble(est)
  saveRDS(est, paste0(save_estimates_file, "_all.rds"))
  
  # read in true values of estimands.
  true <- readRDS(true_file)
  #head(true)
  res <- true %>%
    mutate("scenario_fit_id" = scenario_fit_id, 
           isim = 0, 
           model_type = "true") %>%
    select(-value_se) %>%
    full_join(est, by = c("estimand", "estimand_id","trt", "t" ,"value", 
                          "scenario_fit_id", "isim", "model_type")) %>%
    relocate(c( "isim", "model_type", "scenario_fit_id"), .after = last_col())

  saveRDS(res, save_results_file)
  
}

# Function to combine files (used across scenarios in the blocks folder).
# read_file must be correctly formatted 

combine_scenarios <- function(start, stop, read_file, save_file){
  
  for(i in start:stop){
    
    # add i.rds to read_file prefix.
    temp <- readRDS(paste0(read_file,i,".rds"))
    
    if(i == start){
      res <- temp
    }  else  {
      res <- rbindlist(list(res, temp))
    }
    print(paste0("file ", i, "/", stop-start+1))
    res <- as_tibble(res)
  }
  
  saveRDS(res, save_file)
# saveRDS(res, "../all_res.rds")
  
}


# Run rsimsum for rmst and irmst for each scenario and save results.

simsum_rmst <- function(design_id, scenario_fit_id, results_file, save_file){
  
  res <- readRDS(results_file)
  
  performance_combine <- NULL
  
  estimand_temp <- estimand_labels %>%
    filter(estimand %in% c("rmst","irmst"))
  
  if(design_id == "single_arm"){
    
    estimand_temp <- estimand_temp %>%
      filter(trt == 0 & !is.na(trt))
  }
  
  
  for(j in 1:nrow(estimand_temp)){
    #j <- 1
    estimand_id_temp <- estimand_temp$estimand_id[j]
    
    true_value <- res %>%
      filter(scenario_fit_id == scenario_fit_id) %>%
      filter(estimand_id == estimand_id_temp) %>%
      filter(isim == 0) %>%
      select(value) %>%
      pull()

    performance_temp <- res %>%
      filter(scenario_fit_id == scenario_fit_id) %>%
      filter(estimand_id == estimand_id_temp) %>%
      filter(isim != 0) %>%
      rename(high = upper, low = lower) %>%
      simsum(estvarname = "value", true = true_value, se =  "sd", 
             ci.limits = c("low", "high"),  by = "scenario_fit_id") %>%
      summary() %>%
      tidy() %>%
      mutate(scenario_fit_id = as.character(scenario_fit_id)) %>%
      add_row(stat = "true", est = true_value, scenario_fit_id = scenario_fit_id) # %>%
    

    performance_temp_row <- performance_temp %>%
      filter(stat == "bias") %>%
      mutate(across(c(est,lower,upper), ~ .x + true_value)) %>%
      mutate(stat = "mean")
    
    performance_temp <- performance_temp %>%
      add_row(performance_temp_row) %>%
      mutate(estimand_id = estimand_id_temp) %>%
      left_join(estimand_labels, 
                join_by(estimand_id)) %>%
      relocate(c(scenario_fit_id, estimand, estimand_id, trt, t), .after=last_col())

    performance_combine <- performance_combine %>%
      bind_rows(performance_temp) 
    
  }
 
 print(save_file)
 
 saveRDS(performance_combine, save_file)
  
  
}





