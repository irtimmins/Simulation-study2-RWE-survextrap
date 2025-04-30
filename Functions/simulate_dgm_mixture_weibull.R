
# Simulate survival times from mixture Weibull model with Gompertz GPM rates.

sim_dgm_trt_mix <- function(N,
                            nsim = 1,
                            lambda1,
                            lambda2,
                            gamma1,
                            gamma2,
                            pmix,
                            hr_function = "ph", # proportional hazards default. 
                            beta, # treatment effect
                            alpha_frailty = 0, # frailty on log-hazard scale
                            loghaz_bias = 0, # bias for external data
                            lambda_gpm,
                            gamma_gpm,
                            age_mean,
                            age_sd,
                            cens_min = NA, #  use NA for no uniform censoring
                            cens_max = NA, #  use NA for no uniform censoring
                            maxT = Inf,
                            arms = "both", # specify both arms or "active" or "control".
                            left_trunc = 0, # time of left censoring.
                            lower = 1e-08, 
                            upper = 10000, 
                            nodes = 100,
                            seed,
                            count_format = FALSE,
                            return_data = TRUE,
                            save_file = NA
                            ){
  
  
  set.seed(seed)
  
  if(is.character(beta) & !is.na(beta)) beta <- get(beta)
  
  if(arms == "both"){
    
    sim_data <-  expand_grid(id=1:floor(N/2),
                             i=1:nsim,
                             trt = c(0,1)) %>%
      group_by(i) %>%
      mutate(id = row_number()) %>%
      ungroup()
      
    
  } else if(arms == "control"){
    
    sim_data <- expand_grid(id=1:N,
                            i=1:nsim,
                            trt = 0)
    
    beta <- 0
    hr_function <- "ph"
    
  } else if(arms == "active"){
    
    sim_data <- expand_grid(id=1:N,
                            i=1:nsim,
                            trt = 1)
    
  }
  
  # simulate age distribution and frailty.
  
  sim_data <- sim_data %>%
    mutate(age = rnorm(n(), mean = age_mean, sd = age_sd)) %>%
    mutate(time = 0,
           event = 0) %>%
    mutate(frailty = rnorm(n()))
  
  # Censoring times.
  
  if(is.na(maxT)) maxT <- Inf
  
  sim_data$censtime <- maxT
  
  if(!is.na(cens_min) & !is.na(cens_max)) {
    sim_data <- sim_data %>%
      mutate(censtime = runif(n(), min = cens_min, max = cens_max)) %>%
      mutate(censtime = pmin(censtime, maxT))
  }
  
  
  # Simulate log-survival times from uniform[0,1] distribution
  
  lnU  <- log(runif(nrow(sim_data)))
  
  # Root findings interval for vuniroot.
  
  root_interval <- cbind(rep(lower,nrow(sim_data)),rep(upper,nrow(sim_data)))
  
  # survival times
  
  t <- vuniroot(f_ch,
                interval=root_interval, 
                lnU=lnU, 
                nodes=nodes, 
                lambda1 = lambda1, 
                lambda2 = lambda2, 
                gamma1 = gamma1, 
                gamma2 = gamma2, 
                pmix = pmix,
                hr_function = hr_function,
                beta = beta, 
                trt = sim_data[["trt"]],
                alpha_frailty = alpha_frailty,
                frailty = sim_data[["frailty"]],
                loghaz_bias = loghaz_bias,
                lambda_gpm = lambda_gpm,
                gamma_gpm = gamma_gpm,
                age = sim_data[["age"]],
                left_trunc = left_trunc)$root
  
  
  sim_data <- sim_data %>%
    mutate(time = t,
    # ensure no survival times are left of truncation point (due to rounding errors etc.)
                      time = pmax(t, left_trunc)) %>% 
    mutate(time = pmin(time, censtime),
           event = as.numeric(time<censtime)) %>%
    relocate(time, event, censtime, 
             .after = frailty) %>%
    mutate(backhaz = hgompertz(x = age+time,  shape = gamma_gpm, rate = lambda_gpm)) 
  
  
  if(count_format == TRUE){

    # Convert to aggregate form. 
    
    # First use survSplit
    min_t <- floor(min(sim_data$time))
    max_t <- ceiling(max(sim_data$time))

    sim_data_split <- survSplit(Surv(time, event) ~.,
                                sim_data,
                                cut=min_t:max_t, 
                                episode ="timegroup") %>%
      as_tibble() %>%
      rename(start = tstart)
    
    sim_data <-  sim_data_split %>%
      mutate(stop = start + 1)  %>%
      # Add expected survival rates based on Gompertz model.
      mutate(across(c("start", "stop"), 
                    ~expected_survival(., age, lambda_gpm, gamma_gpm), 
                    .names = "backsurv_{.col}")) %>%
      group_by(i,trt,start) %>%
      summarise(n = n(),
                r = n-sum(event),
                backsurv_start = mean(backsurv_start),
                backsurv_stop = mean(backsurv_stop),
                .groups = "keep") %>%
      mutate(stop = start + 1) %>%
      relocate(c(start, stop), .after = r) %>%
      ungroup() %>%
      filter(start >= left_trunc) %>%
      select(i, trt, start, stop, n, r, backsurv_start, backsurv_stop) %>%
      ungroup()

  }
  
  if (!is.na(save_file)) saveRDS(sim_data, file=save_file)
  
  if (return_data) return(sim_data)
}



# Alternate approach to left-truncation for simulating count data.
# This preserves the age distribution.

sim_dgm_trt_mix_external <- function(N,
                                          nsim = 1,
                                          lambda1,
                                          lambda2,
                                          gamma1,
                                          gamma2,
                                          pmix,
                                          hr_function = "ph", # proportional hazards default. 
                                          beta, # treatment effect
                                          alpha_frailty = 0, # frailty on log-hazard scale
                                          loghaz_bias = 0, # bias for external data
                                          lambda_gpm,
                                          gamma_gpm,
                                          age_mean,
                                          age_sd,
                                          cens_min = NA, #  use NA for no uniform censoring
                                          cens_max = NA, #  use NA for no uniform censoring
                                          maxT = Inf,
                                          arms = "both", # specify both arms or "active" or "control".
                                          left_trunc = 0, # time of left censoring.
                                          lower = 1e-08, 
                                          upper = 10000, 
                                          nodes = 100,
                                          seed,
                                          count_format = FALSE,
                                          return_data = TRUE,
                                          save_file = NA
){
  
  sim_data <- sim_dgm_trt_mix(N*10,
                                   nsim,
                                   lambda1,
                                   lambda2,
                                   gamma1,
                                   gamma2,
                                   pmix,
                                   hr_function, # proportional hazards default. 
                                   beta, # treatment effect
                                   alpha_frailty, # frailty on log-hazard scale
                                   loghaz_bias, # bias for external data
                                   lambda_gpm,
                                   gamma_gpm,
                                   age_mean,
                                   age_sd,
                                   cens_min, #  use NA for no uniform censoring
                                   cens_max, #  use NA for no uniform censoring
                                   maxT,
                                   arms, # specify both arms or "active" or "control".
                                   left_trunc = 0, # time of left censoring.
                                   lower, 
                                   upper, 
                                   nodes,
                                   seed,
                                   count_format = FALSE,
                                   return_data = TRUE,
                                   save_file = NA) 
  
  # Left truncate and reduce to specified number of patients.
  sim_data <- sim_data %>%
    filter(time >= left_trunc) %>%
    slice(1:N)
  
  # Convert to aggregate form. 
  
  # First use survSplit
  min_t <- floor(min(sim_data$time))
  max_t <- ceiling(max(sim_data$time))
  
  sim_data_split <- survSplit(Surv(time, event) ~.,
                              sim_data,
                              cut=min_t:max_t, 
                              episode ="timegroup") %>%
    as_tibble() %>%
    rename(start = tstart)
  
  sim_data <-  sim_data_split %>%
    mutate(stop = start + 1)  %>%
    # Add expected survival rates based on Gompertz model.
    mutate(across(c("start", "stop"), 
                  ~expected_survival(., age, lambda_gpm, gamma_gpm), 
                  .names = "backsurv_{.col}")) %>%
    group_by(i,trt,start) %>%
    summarise(n = n(),
              r = n-sum(event),
              backsurv_start = mean(backsurv_start),
              backsurv_stop = mean(backsurv_stop),
              .groups = "keep") %>%
    mutate(stop = start + 1) %>%
    relocate(c(start, stop), .after = r) %>%
    ungroup() %>%
    filter(start >= left_trunc) %>%
    select(i, trt, start, stop, n, r, backsurv_start, backsurv_stop) %>%
    ungroup()
  
  if (!is.na(save_file)) saveRDS(sim_data, file=save_file)
  if (return_data) return(sim_data)

}


# Expected background survival from conditional Weibull model.
expected_survival <- function(time, age, lambda_gpm, gamma_gpm){
  exp_surv <- (1-pgompertz(q = age+time,  shape = gamma_gpm, rate = lambda_gpm))/
    (1-pgompertz(q = age,  shape = gamma_gpm, rate = lambda_gpm))
  return(exp_surv)
}

# Hazard function for mixture Weibull with treatment effects.

haz <- function(t, lambda1, lambda2, gamma1, gamma2, pmix, hr_function, beta, trt, alpha_frailty, frailty, left_trunc) {
  
  pdf1 <- lambda1*gamma1*(t^(gamma1-1))*exp(-lambda1*t^gamma1)
  pdf2 <- lambda2*gamma2*(t^(gamma2-1))*exp(-lambda2*t^gamma2)
  S1 <- exp(-lambda1*t^gamma1)
  S2 <- exp(-lambda2*t^gamma2)
  
  # Mixture Weibull baseline log hazard.
  
  loghaz <- log(pmix*pdf1 + (1-pmix)*pdf2)-log(pmix*S1 + (1-pmix)*S2)
  
  # Add frailty and bias terms
  
  loghaz <- loghaz + alpha_frailty*frailty
  
  # Add treatment effect 
  
  if(hr_function == "ph"){
    
    loghaz <- loghaz + beta[1]*trt
    
  } else if (hr_function == "wane") {
    
    loghaz <- loghaz +beta[1]*(1-tanh(x = (beta[2]*t+beta[3])))*trt
    
  } else if (hr_function == "delay_wane") {
    
    loghaz <- loghaz +
      (beta[1]+beta[2]*demg(x = t, mu = beta[3], sigma = beta[4], lambda = beta[5]))*trt
    
  }
  
  haz <- exp(loghaz)
  
  haz[t <= left_trunc] <- 0 
  
  return(haz)
  
}

## Cumulative hazard function combining cause-specific and other-cause.

f_ch <- function(t,lnU,nodes, lambda1, lambda2, gamma1, gamma2, pmix, hr_function, 
                 beta, trt, alpha_frailty, frailty, loghaz_bias, 
                 lambda_gpm, gamma_gpm, age, left_trunc) {
  
  
  # use analytic result if proportional hazards.
  if(hr_function == "ph"){ 
  
  S1 <- exp(-lambda1*t^gamma1)
  S2 <- exp(-lambda2*t^gamma2)
  S1_left_trunc <- exp(-lambda1*left_trunc^gamma1)
  S2_left_trunc <- exp(-lambda2*left_trunc^gamma2)
  
  ch_cause <- -log(pmix*S1+(1-pmix)*S2)*exp(alpha_frailty*frailty+beta[1]*trt)
  # cumulative hazard at left truncation point.
  ch_cause_left_trunc <- -log(pmix*S1_left_trunc+(1-pmix)*S2_left_trunc)*exp(alpha_frailty*frailty+beta[1]*trt)
  # take the difference
  ch_cause <- ch_cause - ch_cause_left_trunc
  
  } else {   # else use numerical integration
  
  ch_cause <- vecquad_gl(haz,0,t,nodes, lambda1, lambda2, gamma1, gamma2, pmix, hr_function,
                         beta, trt, alpha_frailty, frailty, left_trunc) 
  }
  
  # note the cancellation of terms in the formula for other cause with left truncation:
  # ch(t) denotes the cumulative hazard at time t.
  #
  # ch_other = ch(time+age | time>left_trunc)-ch(left_trunc+age | time>left_trunc)
  #          = ch(time+age)-ch(left_trunc)-(ch(left_trunc+age)-ch(left_trunc)) # terms cancel
  #          = ch(time+age)-ch(age)
  
  ch_other <-  lambda_gpm*(gamma_gpm^(-1))*(exp(gamma_gpm*(t+age))-1)-
    lambda_gpm*(gamma_gpm^(-1))*(exp(gamma_gpm*(left_trunc+age))-1)
  
  ch <- (ch_cause+ch_other)*exp(loghaz_bias)
  
  ch[t <= left_trunc] <- 0
  
  return(-ch-lnU)
}

# Gaussian quadrature function for integration

vecquad_gl <- function(fn,a,b,nodes,...) {
  nw  <- gaussLegendre(n = nodes,-1,1)
  z1 <- fn(t = 0.5*as.matrix(b-a)%*%t(nw$x)+as.vector(0.5*(a+b)),...)
  
  return((0.5*(b-a))*rowSums(sweep(z1,2, nw$w,"*")))
}  



## Cause-specific cumulative hazard function (for Weibull mixture) using Gaussian quadrature.

f1 <- function(t,lnU,nodes, lambda1, lambda2, gamma1, gamma2, pmix, hr_function, beta, trt, alpha_frailty, frailty, left_trunc) {
  
  ch <- vecquad_gl(haz,0,t,nodes, lambda1, lambda2, gamma1, gamma2, pmix, hr_function, beta, trt, alpha_frailty, frailty, left_trunc) 
  
  return(-ch-lnU)
}



## Cumulative hazard function for Gompertz GPM using closed-form.
## Conditioned on baseline age, with left truncation

f1.gpm <- function(t,lnU,lambda_gpm, gamma_gpm, age, left_trunc) {
  
  cumh.gpm <- lambda_gpm*(gamma_gpm^(-1))*(exp(gamma_gpm*(t+age))-1)-
    lambda_gpm*(gamma_gpm^(-1))*(exp(gamma_gpm*age)-1)
  
  cumh.gpm.left_trunc <- lambda_gpm*(gamma_gpm^(-1))*(exp(gamma_gpm*(left_trunc+age))-1)-
    lambda_gpm*(gamma_gpm^(-1))*(exp(gamma_gpm*age)-1)
  
  ch <- cumh.gpm - cumh.gpm.left_trunc
  
  ch[t <= left_trunc] <- 0
  
  return(-ch-lnU)
}

