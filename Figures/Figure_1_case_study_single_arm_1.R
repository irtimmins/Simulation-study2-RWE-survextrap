######################
library(survival)
library(flexsurv)
library(survextrap)
library(tidyr)
library(dplyr)
require(HMDHFDplus)
require(relsurv) 
require(dplyr)


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

set.seed(130513)

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

#############################################################################


control_mod1 <- survextrap(
  Surv(time, event) ~ 1,
  data = control_data,
  fit_method = "opt",
  smooth_model = "random_walk"
)

plot(control_mod1)

control_mod2 <- survextrap(
  Surv(time, event) ~ 1,
  data = control_data,
  fit_method = "opt",
  smooth_model = "random_walk",
  backhaz = "backhaz"
)


plot(control_mod2)



#############################################################################


sim_trial_data <-
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

summary(ratetable.USA)
summary(survexp.us)

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

cetux_seer_backhaz <- cetux_seer %>%
  left_join(sim_trial_standarised_start, join_by("start")) %>%
  left_join(sim_trial_standarised_stop, join_by("stop"))

#cetux
control_mod3 <- survextrap(
  Surv(time, event) ~ 1,
  data = control_data,
  mspline = list(add_knots = c(10,15,25)),
  external=cetux_seer_backhaz,
  fit_method = "opt",
  smooth_model = "random_walk",
  backhaz = "backhaz"
)

plot(control_mod3)

##############################################################
##############################################################

newdata_trt <- tibble(trt = as.factor(c(0,1))) %>%
  filter(trt == 0)

trial_data <- trial_data %>%
  filter(trt %in% newdata_trt[["trt"]])

t_vec <- seq(from = 0 , to = 40, length.out = 10)
rmst_vec <- c(5,10,20,30,40)

View(backhaz_rates[is.na(backhaz_rates$backhaz),])

backhaz_rates <- trial_data %>%
  select(age, year, sex) %>%
  expand_grid(t = t_vec) %>%
  mutate(expsurv = expsurv_ratetable(t = t, 
                                     age = age, 
                                     year = year, 
                                     sex = sex, 
                                     ratetable = ratetable.USA,
                                     scale_ratetable = 365.25)) %>%
  mutate(backhaz = backhaz_ratetable(t = t, 
                                     age = age, 
                                     year = year, 
                                     sex = sex,
                                     ratetable = ratetable.USA,
                                     scale_ratetable = 365.25)) %>%
  group_by(t) %>%
  summarise(backhaz = sum(backhaz*expsurv)/sum(expsurv),
            expsurv = mean(expsurv)) %>%
  ungroup()

#View(backhaz_rates)

summ_fns <- list(median=median, ~quantile(.x, probs=c(0.025, 0.975)), sd = sd)

test <- hazard(control_mod3, t=t_vec, newdata = newdata_trt, summ_fns=summ_fns)
View(test)

haz_est <- hazard(control_mod3, t=t_vec, newdata = newdata_trt, summ_fns=summ_fns) %>%
  mutate(estimand = "hazard") %>%
  relocate(estimand) %>%
  rename(value=median, lower =`2.5%`, upper=`97.5%`) %>%
  left_join(backhaz_rates, join_by("t")) %>%
  mutate(across(c(value, lower, upper), ~ .x + backhaz)) %>%
  select(-c(backhaz, expsurv)) 

View(haz_est)

surv_est <- survival(control_mod3, t=t_vec, newdata = newdata_trt, summ_fns=summ_fns) %>%
  mutate(estimand = "survival") %>%
  relocate(estimand) %>%
  rename(value=median, lower =`2.5%`, upper=`97.5%`) %>%
  left_join(backhaz_rates, join_by("t")) %>%
  mutate(across(c(value, lower, upper, sd), ~ .x*expsurv)) %>%
  select(-c(backhaz, expsurv)) 

head(surv_est)
View(surv_est)


est3 <- estimates_ratetable(fitted_mod = control_mod3,
                            t = t_vec, 
                            rmst_t = rmst_vec, 
                            newdata = newdata_trt %>%
                              filter(trt == 0), 
                            excess_hazard = TRUE,
                            ratetable = ratetable.USA,
                            scale_ratetable = 365.25,
                            trial_data = trial_data)
View(est3)




