# If running interactively, requires 64Gb of memory
# (or remove models once their summary stats have been calculated)

library(survextrap)
library(ggplot2)
library(dplyr)
library(viridis) # for colour palettes
library(memoise)
library(stringr)
library(gridExtra)
library(pracma)
library(tidyr)
library(flextable)
library(abind)
library(cowplot)
library(magrittr)
library(purrr)
library(wrapr)
library(patchwork)
library(readr)

user <- Sys.info()["user"]

# Load rmst_fast and irmst_fast.

source("Functions/functions.R")

# Load datasets from Bonner trial. 
# Control and Cetuximab arms.
control <- readRDS("Data/Bonner_control.RDS")
cetuximab <- readRDS("Data/Bonner_active.RDS")

# Both arms data.
cetux <- readRDS("Data/Bonner_all.RDS")

########################################################################################
# Single arm models (Control arm only)
########################################################################################

# Load datasets from Bonner trial. 
# Control and Cetuximab arms.
control <- readRDS("Data/Bonner_control.RDS")
cetuximab <- readRDS("Data/Bonner_active.RDS")

# Both arms data.
cetux <- readRDS("Data/Bonner_all.RDS")

# Derive Kaplan-Meier curves of Bonner trial.
# Both arms
survminer::ggsurvplot(survfit(Surv(years, d) ~ treat, data=cetux)) + xlab("Years")
# Control
km_control <- survminer::ggsurvplot(survfit(Surv(years, d) ~ 1, data=control))$data.survplot
# Cetuximab
km_cetux <- survminer::ggsurvplot(survfit(Surv(years, d) ~ 1, data=cetuximab))$data.survplot

tmax <- max(cetux$years)

# Derive priors.
mspline_base <- mspline_spec(Surv(years, d) ~ 1, data=cetux, df=10)
prior_hscale <- p_meansurv(median=25, upper=100, mspline=mspline_base) # prior for the baseline log hazard scale parameter. We'll leave this untouched once calculated (from 10df mspline)
prior_haz_const(mspline_base, prior_hscale = prior_hscale)


### Running models
options(mc.cores = 1) 
# options(mc.cores = parallel::detectCores())
chains <- 4; iter <- 2000
smooth_model <- "exchangeable"

# Whether to run single arm models or load in predictions from disk
run_single_arm_models <- FALSE

if(run_single_arm_models){
  models_control <- list()
  # We define the default model and investigate 1-way sensitivity analyses to this set-up
  df_base <- 10 # degrees of freedom
  ek_base <- NULL # extra knots after end of trial follow-up
  l_base <- 1 # rate for sigma prior 
  prior_hsd_base <- p_gamma(2, l_base) # changing prior for sigma parameter
  mspline_base <- mspline_spec(Surv(years, d) ~ 1, data=cetux, df=df_base, add_knots=ek_base)
  external_data_base <- "Trial only" # no background rates and no SEER data for default
  external_base <- NULL
  backhaz_base <- NULL
  models_control[[paste0("df=",df_base,"; extra_knots=", ek_base,"; prior_rate=", l_base, 
                         "; external_data=", external_data_base)]] <- 
    survextrap_mem(Surv(years, d) ~ 1, data=control, mspline=mspline_base,
                   chains=chains, iter=iter,
                   prior_hscale=prior_hscale, prior_hsd = prior_hsd_base,
                   external=external_base, 
                   backhaz=backhaz_base,
                   smooth_model = smooth_model)
  
  # loop over no. df for trial data
  for(df in c(3,6,10)){
    print(paste("df=", df))
    mspline <- mspline_spec(Surv(years, d) ~ 1, data=cetux, df=df, add_knots=ek_base)
    models_control[[paste0("df=",df,"; extra_knots=", ek_base,"; prior_rate=", l_base,
                           "; external_data=", external_data_base)]] <- 
      survextrap_mem(Surv(years, d) ~ 1, data=control, mspline=mspline,
                     chains=chains, iter=iter,
                     prior_hscale=prior_hscale, prior_hsd = prior_hsd_base,
                     external=external_base, 
                     backhaz=backhaz_base,
                     smooth_model = smooth_model
      )
  }
  
  
  # extra knots, for trial data alone, trial data+BH, trial data+BH+seer
  extra_knots <- list(NULL, 25, c(10, 25), c(10, 15, 25))
  external_data <- c("Trial only", "Trial + Population rates", "Trial + Population rates + Registry")
  for(ek in extra_knots){
    for(ed in external_data){
      print(paste("extra knots=",paste(ek, collapse=",")))
      # Spline specification
      mspline <- mspline_spec(Surv(years, d) ~ 1, data=cetux, df=df_base, add_knots=ek)
      
      # external data
      if(ed == "Trial + Population rates + Registry"){
        external <- cetux_seer
        backhaz <- cetux_bh
      } else if(ed == "Trial + Population rates"){
        external <- NULL   
        backhaz <- cetux_bh
      } else if(ed == "Trial only"){
        external <- NULL
        backhaz <- NULL
      } else stop("Options of external data do not match")
      models_control[[paste0("df=",df_base,"; extra_knots=", paste(ek, collapse=","),"; prior_rate=", l_base,
                             "; external_data=", ed)]] <-  
        survextrap_mem(Surv(years, d) ~ 1, data=control, mspline=mspline,
                       chains=chains, iter=iter,
                       prior_hscale=prior_hscale, prior_hsd = prior_hsd_base,
                       external=external, 
                       backhaz=backhaz,
                       smooth_model = smooth_model
        )
    }
  }
  
  # priors (trial data only)
  for(l in c(1, 5, 20)){ # changing prior for sigma parameter
    print(paste("prior rate=", l))
    prior_hsd <- p_gamma(2, l) 
    
    models_control[[paste0("df=",df_base,"; extra_knots=", ek_base,"; prior_rate=", l,
                           "; external_data=", external_data_base)]] <- 
      survextrap_mem(Surv(years, d) ~ 1, data=control, mspline=mspline_base,
                     chains=chains, iter=iter,
                     prior_hscale=prior_hscale, prior_hsd = prior_hsd,
                     external=external_base, 
                     backhaz=backhaz_base,
                     smooth_model = smooth_model
      )
  }
  
  
  
  ############################################
  # LOOIC stats from trial data only
  looic <- data.frame()
  for(k in seq_along(models_control)){
    print(k)
    name <- names(models_control)[k]
    # Use a regular expression to extract the number following "df="
    df <- as.numeric(str_extract(name, "(?<=df=)[^;]+"))
    ek <- str_extract(name, "(?<=extra_knots=)[^;]+")
    prior_rate <- as.numeric(str_extract(name, "(?<=prior_rate=)[^;]+"))
    ed <- str_extract(name, "(?<=external_data=)[^;]+")
    
    if(ed == "Trial only"){
      looic <- rbind(looic,
                     data.frame(df = df, ek = ek, prior_rate = prior_rate, ed = ed,
                                looic = models_control[[k]]$loo$estimates["looic","Estimate"],
                                rmst_mem(models_control[[k]], t=5))
      )
    }                                  
  }
  save(looic, file =  "looic.Rdata") # save looic stats to disk
  
  # rmst metrics at 40 years
  rmst_control <- data.frame()
  for(k in seq_along(models_control)){
    print(k)
    name <- names(models_control)[k]
    # Use a regular expression to extract the number following "df="
    df <- as.numeric(str_extract(name, "(?<=df=)[^;]+"))
    ek <- str_extract(name, "(?<=extra_knots=)[^;]+")
    prior_rate <- as.numeric(str_extract(name, "(?<=prior_rate=)[^;]+"))
    ed <- str_extract(name, "(?<=external_data=)[^;]+")
    
    rmst_control <- rbind(rmst_control,
                          data.frame(df = df, ek = ek, prior_rate = prior_rate, ed = ed,
                                     rmst_mem(models_control[[k]], t=40))
    )
  }
  save(rmst_control, file =  "rmst_control.Rdata") # save rmst stats to disk
  # saveRDS(haz_plot, paste0("plots/single_arm/hazard_plot",i,".rds"))
  
  
  
  # survival metrics
  surv_plots <- data.frame()
  for(k in seq_along(models_control)){
    print(k)
    name <- names(models_control)[k]
    # Use a regular expression to extract the number following "df="
    df <- as.numeric(str_extract(name, "(?<=df=)[^;]+"))
    ek <- str_extract(name, "(?<=extra_knots=)[^;]+")
    prior_rate <- as.numeric(str_extract(name, "(?<=prior_rate=)[^;]+"))
    ed <- str_extract(name, "(?<=external_data=)[^;]+")
    
    surv_plots <- rbind(surv_plots,
                        data.frame(df = df, ek = ek, prior_rate = prior_rate, ed = ed,
                                   survival_mem(models_control[[k]], tmax=40), 
                                   model_number = k)
    )
    
  }
  surv_plots$df <- factor(surv_plots$df)
  surv_plots$`Extra knots` <- factor(ifelse(is.na(surv_plots$ek), "None", surv_plots$ek))
  surv_plots$`Prior for sigma` <- factor(surv_plots$prior_rate)
  surv_plots$`External data` <- factor(surv_plots$ed, levels = c("Trial only", "Trial + Population rates", "Trial + Population rates + Registry"))
  save(surv_plots, file =  "surv_plots.Rdata") # save survival stats to disk
  
  # hazard metrics
  haz_plots <- data.frame()
  for(k in seq_along(models_control)){
    print(k)
    name <- names(models_control)[k]
    # Use a regular expression to extract the number following "df="
    df <- as.numeric(str_extract(name, "(?<=df=)[^;]+"))
    ek <- str_extract(name, "(?<=extra_knots=)[^;]+")
    prior_rate <- as.numeric(str_extract(name, "(?<=prior_rate=)[^;]+"))
    ed <- str_extract(name, "(?<=external_data=)[^;]+")
    
    haz_plots <- rbind(haz_plots,
                       data.frame(df = df, ek = ek, prior_rate = prior_rate, ed=ed,
                                  hazard_mem(models_control[[k]], t=c(seq(0,6,length=100),7:40)), 
                                  model_number = k)
    )
    
  }
  haz_plots$df <- factor(haz_plots$df)
  haz_plots$`Extra knots` <- factor(ifelse(is.na(haz_plots$ek), "None", haz_plots$ek))
  haz_plots$`Prior for sigma` <- factor(haz_plots$prior_rate)
  haz_plots$`External data` <- factor(haz_plots$ed, levels = c("Trial only", "Trial + Population rates", "Trial + Population rates + Registry"))
  save(haz_plots, file =  "haz_plots.Rdata") # save hazard stats to disk
  
} else {
  load("looic.Rdata") # load looic stats 
  load("surv_plots.Rdata") # load survival stats
  load("haz_plots.Rdata") # load hazard stats
  load("rmst_control.Rdata") # load rmst stats
  
}




############################################
# Survival Plots
############################################

# Varying number of degrees of freedom (Trial only)
s1 <- ggplot(surv_plots %>% filter(`Extra knots` == "None" & prior_rate==1 & `External data` == "Trial only")) +
  geom_line(aes(x=t, y=median, 
                col=df), lwd=1.3) +
  geom_line(aes(x=t, y=lower, col=df), linetype = "dashed") +
  geom_line(aes(x=t, y=upper, col=df), linetype = "dashed") +
  geom_step(data=km_control, aes(x=time, y=surv), lwd=1.3,
            inherit.aes = FALSE) +
  xlab("Time (Years)") + 
  scale_y_continuous("Overall Survival", limits = c(0,1), labels = scales::percent) +
  geom_vline(xintercept = max(control$years)) +
  ggtitle("(a) Varying number of degrees of freedom") +
  theme_classic()+
  theme_paper2()+
  theme(plot.title = element_text(hjust = 0,
                                  size=10,
                                  face="bold",
                                  lineheight = 1.2))
s1
# Varying number of extra knots (Trial only)
# s2 <- ggplot(surv_plots %>% filter(df == 10 & prior_rate==1 & `External data` == "Trial only")) +
#   geom_line(aes(x=t, y=median, 
#                 col=`Extra knots`), lwd=1.3) +
#   geom_line(aes(x=t, y=lower, col=`Extra knots`), linetype = "dashed") +
#   geom_line(aes(x=t, y=upper, col=`Extra knots`), linetype = "dashed") +
#   geom_step(data=km_control, aes(x=time, y=surv), lwd=1.3,
#             inherit.aes = FALSE) +
#   xlab("Time (Years)") + 
#   scale_y_continuous("Overall Survival", limits = c(0,1), labels = scales::percent) +
#   geom_vline(xintercept = max(control$years)) +
#   ggtitle("(b) Varying extra knot locations after trial follow-up")+
#   theme_classic()+
#   theme_paper2()+
#   theme(plot.title = element_text(hjust = 0))
# s2
# Varying prior for sigma (Trial only)
s2 <- ggplot(surv_plots %>% filter(df == 10 & `Extra knots` == "None" & `External data` == "Trial only")) +
  geom_line(aes(x=t, y=median, 
                col=`Prior for sigma`), lwd=1.3) +
  geom_line(aes(x=t, y=lower, col=`Prior for sigma`), linetype = "dashed") +
  geom_line(aes(x=t, y=upper, col=`Prior for sigma`), linetype = "dashed") +
  geom_step(data=km_control, aes(x=time, y=surv), lwd=1.3,
            inherit.aes = FALSE) +
  xlab("Time (Years)") + 
  scale_y_continuous("Overall Survival", limits = c(0,1), labels = scales::percent) +
  geom_vline(xintercept = max(control$years)) +
  ggtitle("(b) Varying prior for sigma")+
  theme_classic()+
  theme_paper2()+
  theme(plot.title = element_text(hjust = 0,
                                  size=10,
                                  face="bold",
                                  lineheight = 1.2))


plot_single_sensitivity <- plot_grid(s1, s2, ncol=1, align = "hv")

tiff(file = "figure_survival_cetux_single_arm_sensitivity.tiff",   
     width = 5.8, 
     height = 4.5,
     units = 'in',  
     res = 300, 
     compression = "lzw")
print(plot_single_sensitivity)
dev.off()



############## 
# Survival varying extra knots across trial only, trial + BH, trial + BH + seer
s4 <- ggplot(surv_plots %>% filter(df == 10 & prior_rate==1)) +
  geom_line(aes(x=t, y=median, 
                col=`Extra knots`), lwd=1.3) +
  geom_line(aes(x=t, y=lower, col=`Extra knots`), linetype = "dashed") +
  geom_line(aes(x=t, y=upper, col=`Extra knots`), linetype = "dashed") +
  facet_wrap(~`External data`, nrow = 3) +
  geom_step(data=km_control, aes(x=time, y=surv), lwd=1.3,
            inherit.aes = FALSE) +
  xlab("Time (Years)") + 
  scale_y_continuous("Overall Survival", limits = c(0,1), labels = scales::percent) +
  geom_vline(xintercept = max(control$years)) +
 # ggtitle("Varying extra knot locations after trial follow-up")+
  theme_classic()+
  theme_paper2()+
  theme(plot.title = element_text(hjust = 0))+
  theme(axis.title.x = element_text(size = 8),
        strip.text.x = element_text(size = 10))+  
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "grey20"))
s4

tiff(file = "figure_survival_cetux_single_arm_sensitivity_knots.tiff",   
     width = 5.8, 
     height = 5.3,
     units = 'in',  
     res = 300, 
     compression = "lzw")
print(s4)
dev.off()


############################################
# Hazard Plots

# Varying number of degrees of freedom
h1 <- ggplot(haz_plots %>% filter(`Extra knots` == "None" & prior_rate==1 & `External data` == "Trial only")) +
  geom_line(aes(x=t, y=median, 
                col=df), lwd=1.3) +
  geom_line(aes(x=t, y=lower, col=df), linetype = "dashed") +
  geom_line(aes(x=t, y=upper, col=df), linetype = "dashed") +
  xlab("Time (Years)") + 
  ylab("Hazard") +
  geom_vline(xintercept = max(control$years)) +
  ggtitle("(a) Varying number of degrees of freedom") +
  theme_classic()+
  theme_paper2()+
  theme(plot.title = element_text(hjust = 0))
h1
# Varying number of extra knots
h2 <- ggplot(haz_plots %>% filter(df == 10 & prior_rate==1 & `External data` == "Trial only")) +
  geom_line(aes(x=t, y=median, 
                col=`Extra knots`), lwd=1.3) +
  geom_line(aes(x=t, y=lower, col=`Extra knots`), linetype = "dashed") +
  geom_line(aes(x=t, y=upper, col=`Extra knots`), linetype = "dashed") +
  xlab("Time (Years)") + 
  ylab("Hazard") +
  geom_vline(xintercept = max(control$years)) +
  ggtitle("(b) Varying extra knot locations after trial follow-up")+
  theme_classic()+
  theme_paper2()+
  theme(plot.title = element_text(hjust = 0))
h2
# Varying prior for sigma
h3 <- ggplot(haz_plots %>% filter(df == 10 & `Extra knots` == "None" & `External data` == "Trial only")) +
  geom_line(aes(x=t, y=median, 
                col=`Prior for sigma`), lwd=1.3) +
  geom_line(aes(x=t, y=lower, col=`Prior for sigma`), linetype = "dashed") +
  geom_line(aes(x=t, y=upper, col=`Prior for sigma`), linetype = "dashed") +
  xlab("Time (Years)") + 
  ylab("Hazard") +
  geom_vline(xintercept = max(control$years)) +
  ggtitle("(c) Varying prior for sigma")+
  theme_classic()+
  theme_paper2()+
  theme(plot.title = element_text(hjust = 0))
h3

plot_single_sensitivity_hazard <- plot_grid(h1, h2, h3, ncol=1, align = "hv")


tiff(file = "figure_hazard_cetux_single_arm_sensitivity.tiff",   
     width = 5.8, 
     height = 5.3,
     units = 'in',  
     res = 300, 
     compression = "lzw")
print(plot_single_sensitivity_hazard)
dev.off()



################################ 
# Hazard varying extra knots across trial only, trial + BH, trial + BH + seer
cetux_seer_plot <- cbind(cetux_seer, `External data` =  factor("Trial + Population rates + Registry", levels = c("Trial only", "Trial + Population rates", "Trial + Population rates + Registry")))
cetux_bh_plot <- rbind(
  cbind(cetux_bh, `External data` =  factor("Trial + Population rates", levels = c("Trial only", "Trial + Population rates", "Trial + Population rates + Registry"))),
  cbind(cetux_bh, `External data` =  factor("Trial + Population rates + Registry", levels = c("Trial only", "Trial + Population rates", "Trial + Population rates + Registry")))
)
bh_text <- tibble(x=35, y=0.1, lab = "Population mortality",
                      `External data` =  factor(c("Trial + Population rates + Registry", "Trial + Population rates"), levels = c("Trial only", "Trial + Population rates", "Trial + Population rates + Registry")))
seer_text <- tibble(x=20, y=0.3, lab = "Registry data",
                  `External data` =  factor(c("Trial + Population rates + Registry"), levels = c("Trial only", "Trial + Population rates", "Trial + Population rates + Registry")))

h4 <- ggplot(haz_plots %>% filter(df == 10 & prior_rate==1)) +
  geom_line(aes(x=t, y=median, 
                col=`Extra knots`), lwd=1.3) +
  geom_line(aes(x=t, y=lower, col=`Extra knots`), linetype = "dashed") +
  geom_line(aes(x=t, y=upper, col=`Extra knots`), linetype = "dashed") +
  geom_step(data=cetux_seer_plot, 
            aes(x=start, y=haz), inherit.aes = FALSE) +
  geom_step(data=cetux_bh_plot %>% filter(time < 40), 
            aes(x=time, y=hazard), inherit.aes = FALSE, colour="darkgrey") +
  geom_text(data = bh_text, aes(x=x, y=y, label = lab), colour="darkgrey") +
  geom_text(data = seer_text, aes(x=x, y=y, label = lab)) +
  facet_wrap(~`External data`, ncol = 1) +
  xlab("Time (Years)") + 
  ylab("Hazard") +
  geom_vline(xintercept = max(control$years), colour = "gray30", linetype = "dashed")+
 # ggtitle("Varying extra knot locations after trial follow-up")+
  theme_classic()+
  theme_paper2()+
  theme(plot.title = element_text(hjust = 0))+
  theme(axis.title.x = element_text(size = 8),
        strip.text.x = element_text(size = 10))+  
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "grey20"))
h4

tiff(file = "figure_hazard_cetux_single_arm_sensitivity_knots.tiff",   
     width = 5.8, 
     height = 5.3,
     units = 'in',  
     res = 300, 
     compression = "lzw")
print(h4)
dev.off()






