
# Function for model fit and estimation.

fit_est_slurm <- function(mspline_df,
                          smooth_model,
                          bsmooth,
                          model,
                          stan_fit_method,
                          chains,
                          iter,
                          knot_settings,
                          add_knots,
                          wane_period_start,
                          wane_period_stop,
                          wane_nt,
                          trial_data,
                          external_data,
                          backhaz_data,
                          save_file,
                          save_mspline = NA,
                          seed = NULL){
  
  start_all <- Sys.time()
  
  set.seed(seed)
  
  trial_data <- readRDS(trial_data) %>%
    mutate(trt = as.factor(trt))
  
  knots <- set_knots(trial_data, setting = knot_settings, mspline_df = mspline_df)
  
  add_knots <- extra_knots_settings[[add_knots]]
  
  mspline <- list(df = mspline_df,
                  knots = knots,
                  add_knots = add_knots,
                  bsmooth = bsmooth)
  
  if(all(is.na(add_knots))) mspline <- within(mspline, rm("add_knots")) 
  
  # Remove knots from mspline object to ensure default knots are used.
  
  if(is.na(knots)) mspline <- within(mspline, rm("knots")) 
  
  if(!is.na(save_mspline)) saveRDS(mspline, save_mspline)
  
  
  if(is.na(external_data)) {
    
    external_data <- NULL 
    
  } else {
    
    external_data <- readRDS(external_data) %>%
      mutate(trt = factor(trt, levels = c(0,1)))
    
  }
  
  if(is.na(backhaz_data)) {
    
    backhaz_col_name <- NULL
    lambda_gpm <- NULL
    gamma_gpm <- NULL
    
  } else {
    
    backhaz <- readRDS(backhaz_data)
    backhaz_col_name <- "backhaz"
    
    lambda_gpm <- backhaz$lambda_gpm
    gamma_gpm <- backhaz$gamma_gpm
    
  }
  
  start_fit <- Sys.time()
  
  if(model == "single_arm"){
    
    fitted_model <- survextrap(Surv(time, event) ~ 1, 
                               mspline = mspline, 
                               data = trial_data,
                               fit_method= stan_fit_method,
                               chains = chains,
                               iter = iter,
                               smooth_model = smooth_model,
                               external = external_data,
                               backhaz = backhaz_col_name) 
    
  } else if (model == "ph" | model == "nonph"){
    
    nonprop <- ifelse(model == "ph", F, T)
    
    fitted_model <- survextrap(Surv(time, event) ~ trt, 
                               nonprop = nonprop,
                               mspline = mspline, 
                               data = trial_data,
                               fit_method= stan_fit_method,
                               chains = chains,
                               iter = iter,
                               smooth_model = smooth_model,
                               external = external_data,
                               backhaz = backhaz_col_name)    
    
  } else if (model == "separate") {
    
    control_model <- survextrap(Surv(time, event) ~ 1, 
                                mspline = mspline, 
                                data = trial_data %>%
                                  filter(trt == 0),
                                fit_method= stan_fit_method,
                                chains = chains,
                                iter = iter,
                                smooth_model = smooth_model,
                                external = external_data,
                                backhaz = backhaz_col_name)   
    
    # Ensure that active arm separate arm model
    # has no extra knots
    if("add_knots" %in% names(mspline)){
      mspline_active <- within(mspline, rm("add_knots")) 
    } else{
      mspline_active <- mspline  
    }
    
    active_model <- survextrap(Surv(time, event) ~ 1, 
                               mspline = mspline_active,
                               data = trial_data %>%
                                 filter(trt == 1),
                               fit_method = stan_fit_method,
                               chains = chains,
                               iter = iter,
                               smooth_model = smooth_model,
                               backhaz = backhaz_col_name)   
    
  }
  
  end_fit <- Sys.time()
  
  start_est <- Sys.time()
  
  newdata_trt <- tibble(trt = as.factor(c(0,1)))
  
  excess_hazard <- if_else(is.na(backhaz_data), FALSE, TRUE) 
  
  if(model == "single_arm") { 
    
    est_df <- estimates(fitted_model, 
                        t = t_vec, 
                        rmst_t = rmst_vec, 
                        newdata = newdata_trt %>%
                          filter(trt == 0), 
                        excess_hazard = excess_hazard,
                        lambda_gpm = lambda_gpm,
                        gamma_gpm = gamma_gpm,
                        trial_data = trial_data)
    
  } else if( model == "ph" | model == "nonph") {
    

    if(!is.na(wane_period_start) & !is.na(wane_period_stop)){
      wane_period <- c(wane_period_start, wane_period_stop)
      newdata0 <- newdata_trt %>%
        filter(trt == 0)
    } else{
      wane_period <- NULL
      newdata0 <- NULL
      wane_nt <- NULL
    }
    
    # hazard, survival and rmst in each arm.
    
    est_control <- estimates(fitted_model, 
                             t = t_vec, 
                             rmst_t = rmst_vec, 
                             newdata = newdata_trt %>%
                               filter(trt == 0), 
                             excess_hazard = excess_hazard,
                             lambda_gpm = lambda_gpm,
                             gamma_gpm = gamma_gpm,
                             trial_data = trial_data %>%
                               filter(trt == 0))
    
    est_active <- estimates(fitted_model, 
                            t = t_vec, 
                            rmst_t = rmst_vec, 
                            newdata = newdata_trt %>%
                              filter(trt == 1), 
                            excess_hazard = excess_hazard,
                            lambda_gpm = lambda_gpm,
                            gamma_gpm = gamma_gpm,
                            trial_data = trial_data %>%
                              filter(trt == 1),
                            wane_period = wane_period,
                            wane_nt = wane_nt,
                            newdata0 = newdata0)
    
    #### irmst and hazard ratio, use separate functions.
    
    est_irmst <- estimates_irmst(fitted_model,
                                 rmst_t = rmst_vec, 
                                 newdata = newdata_trt,
                                 excess_hazard = excess_hazard,
                                 lambda_gpm = lambda_gpm,
                                 gamma_gpm = gamma_gpm,
                                 trial_data = trial_data,
                                 wane_period = wane_period,
                                 wane_nt = wane_nt,
                                 newdata0 = newdata0)
    
    est_hr <- estimates_hr(fitted_model,
                           t = t_vec, 
                           newdata = newdata_trt,
                           excess_hazard = excess_hazard,
                           lambda_gpm = lambda_gpm,
                           gamma_gpm = gamma_gpm,
                           trial_data = trial_data,
                           wane_period = wane_period,
                           wane_nt = wane_nt,
                           newdata0 = newdata0)
    
    est_df <- bind_rows(est_control,
                        est_active,
                        est_irmst,
                        est_hr) 
    
  }
  
  else if (model == "separate") {
    
    # hazard, survival and rmst in each arm.
    
    est_control <- estimates(control_model,
                             t = t_vec, 
                             rmst_t = rmst_vec, 
                             newdata = newdata_trt %>%
                               filter(trt == 0), 
                             excess_hazard = excess_hazard,
                             lambda_gpm = lambda_gpm,
                             gamma_gpm = gamma_gpm,
                             trial_data = trial_data %>%
                               filter(trt == 0))
    
    est_active <- estimates(active_model,
                            t = t_vec, 
                            rmst_t = rmst_vec, 
                            newdata = newdata_trt %>%
                              filter(trt == 1), 
                            excess_hazard = excess_hazard,
                            lambda_gpm = lambda_gpm,
                            gamma_gpm = gamma_gpm,
                            trial_data = trial_data %>%
                              filter(trt == 1))
    
    #### irmst and hazard ratio estimates.
    
    est_irmst <- estimates_irmst_separate(control_model = control_model,
                                          active_model = active_model,
                                          rmst_t = rmst_vec, 
                                          newdata = newdata_trt, 
                                          excess_hazard = excess_hazard,
                                          lambda_gpm = lambda_gpm,
                                          gamma_gpm = gamma_gpm,
                                          trial_data = trial_data)
    
    est_hr <- estimates_hr_separate(control_model = control_model,
                                    active_model = active_model,
                                    t = t_vec, 
                                    newdata = newdata_trt, 
                                    excess_hazard = excess_hazard,
                                    lambda_gpm = lambda_gpm,
                                    gamma_gpm = gamma_gpm,
                                    trial_data = trial_data)
    
    est_df <- bind_rows(est_control,
                        est_active,
                        est_irmst,
                        est_hr) 
    
  }
  
  end_est <- Sys.time()
  
  end_all <- Sys.time()
  
  timings <- tibble("estimand" = c("run_time", "fit_time", "est_time"),
                    "value" = c(as.numeric(difftime(end_all, start_all, units='mins')),
                                as.numeric(difftime(end_fit, start_fit, units='mins')),
                                as.numeric(difftime(end_est, start_est, units='mins'))
                    )  )
  
  diags_df <- NULL
  
  res <- est_df %>% 
    #  full_join(diags_df, by=c("estimand","value")) %>%
    full_join(timings, by=c("estimand","value") )
  
  if (!is.na(save_file)) saveRDS(res, file=save_file)
  
  return(res)
}



# function for defining knot location setting.


set_knots <- function(trial_data, setting = "default" , mspline_df = 10){
  
  if(setting == "default"){
    
  knots <- NA
    
  } else if(setting == "knot_setting1"){
   
    probs_setting1 <- c(seq(from = 0, to = 0.9, length.out = (mspline_df-2))[-1], 0.95, 1)

    knots <- quantile(trial_data %>% 
                              filter(event == 1) %>% 
                              pull(time),
                            probs = probs_setting1)

  } else if(setting == "knot_setting2"){
    
    probs_setting2 <- seq(from = 0, to = 1, length.out = mspline_df)[-1]
    
    knots <- quantile(trial_data %>% 
                        pull(time),
                      probs = probs_setting2)
    
  } else if(setting == "knot_setting3"){
    
    probs_setting3 <- c(seq(from = 0, to = 0.9, length.out = (mspline_df-2))[-1], 0.95, 1)
    
    knots <- quantile(trial_data %>% 
                        pull(time),
                      probs = probs_setting3)
    
  } else if(setting == "knot_setting4"){
    
    probs_setting4 <- c(seq(from = 0, to = 0.95, length.out = (mspline_df-2))[-1], 0.975, 1)
    
    knots <- quantile(trial_data %>% 
                        pull(time),
                      probs = probs_setting4)
    
  } else if(setting == "knot_setting5"){
    
    knots <- seq(from = 0, to = 5, length.out = mspline_df)[-1]
    
  } else if(setting == "knot_setting6"){
  
    knots <- c(seq(from = 0, to = 4.5, length.out = mspline_df-2)[-1], 4.75, 5)
    
  } else if(setting == "knot_setting7"){
    
    knots <- c(seq(from = 0, to = 4.75, length.out = mspline_df-2)[-1], 4.875, 5)
    
  }
  
  return(knots)
  
}




# rmst with Gauss-Legendre approximation.
rmst_fast <- function(x,t,newdata=NULL,
                      newdata0=NULL, wane_period=NULL, wane_nt=10, disc_rate=0,
                      niter=NULL,summ_fns=NULL,sample=FALSE, 
                      nodes = 100){
  
  newdata <- survextrap:::default_newdata(x, newdata)
  res <- NULL
  
  # Perform for each row of newdata at a time. 
  # If newdata is NULL, perform loop once
  N <- max(nrow(newdata), 1)
  for(i in 1:N){
    
    if(is.null(newdata)){
      newdata_temp <- NULL 
    } else {
      newdata_temp <- newdata %>% slice(i)
    }    
    
    # Perform each timepoint separately
    res_j <- NULL
    
    for(j in 1:length(t)){
      
      time_rmst <- t[j]
      
      GL <- gaussLegendre(n = nodes, 0, time_rmst)
      GL_times <- GL$x
      GL_weights <- GL$w 
      
      W <- diag(GL_weights)
      
      survival_chains <- survival(x, t=GL_times, 
                                  newdata=newdata_temp, newdata0=newdata0,
                                  wane_period = wane_period, wane_nt=wane_nt, 
                                  niter=niter, sample=TRUE, 
                                  disc_rate)[,,]
      
      iterations <- dim(survival_chains)[2]
      
      # evaluate Gauss-Legendre integrals across chains using matrix multiplication
      rmst_samples <- rep(1,nodes) %*% W %*% survival_chains
      rmst_samples_array <- as.array(rmst_samples)
      dim(rmst_samples_array) <- c(1,iterations,1)
      
      temp_j <- survextrap:::summarise_output(rmst_samples_array, t = time_rmst, summ_fns, 
                                              newdata=newdata_temp, sample=sample, summ_name = "rmst")
      
      if(is.null(res_j)){
        res_j <- temp_j
      } else{
        
        if(!sample) 
          res_j <- bind_rows(res_j, temp_j) 
        else { 
          # The array is built up along
          # first dimension for timepoints.
          res_j <- abind(res_j, temp_j, along = 1) 
        }
      }
    }  
    
    
    if(is.null(res)){
      res <- res_j
    } else{
      if(!sample) 
        res <- bind_rows(res, res_j) 
      else { 
        # The array is built up along
        # third dimension for covariates.
        res <- abind(res, res_j, along = 3) 
      }
    }
  }
  
  
  return(res)
  
}

irmst_fast <- function (x, t, newdata = NULL, newdata0 = NULL, wane_period = NULL, 
                        wane_nt = 10, niter = NULL, summ_fns = NULL, sample = FALSE, 
                        disc_rate = 0) 
{
  newdata <- survextrap:::default_newdata_comparison(x, newdata)
  
  rmst1 <- rmst_fast(x, t = t, newdata = newdata[1, , drop = FALSE], 
                     newdata0 = newdata0[1, , drop = FALSE], wane_period = wane_period, 
                     wane_nt = wane_nt, niter = niter, sample = TRUE, disc_rate = disc_rate)
  
  rmst2 <- rmst_fast(x, t = t, newdata = newdata[2, , drop = FALSE], 
                     newdata0 = newdata0[2, , drop = FALSE], wane_period = wane_period, 
                     wane_nt = wane_nt, niter = niter, sample = TRUE, disc_rate = disc_rate)
  
  irmst_sam <- rmst2 - rmst1
  
  res <- survextrap:::summarise_output(irmst_sam, summ_fns, t, newdata = NULL, 
                                       summ_name = "irmst", sample = sample)
  res
}


