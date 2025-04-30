# If running interactively, you cannot get all models into memory you need 64Gb of memory
# (or remove models once their summary stats have been calculated)

rm(list=ls())

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

source("R/functions_aim2.R")

# Choose working directory.
# For testing scripts.
# if(user == "klvq491"){
#   jobname <- "mix_weib_full1"
#   scratch_directory <- paste0("/scratch/", user, "/")
#   store_res <- paste0(scratch_directory, "simsurvextrap_slurm_", jobname, "/")
#   if(!dir.exists("case_study")){
#     dir.create("case_study")
#   }
#   setwd(paste0(store_res, "/case_study"))
#   
# } else {
  
wd <- "/projects/aa/statistical_innovation/survextrap/cetux_case_study"

if(!dir.exists(wd)){
  dir.create(wd)
 }
setwd(wd)

# Data 
control <- cetux[cetux$treat=="Control",]
cetuximab <- cetux[cetux$treat=="Cetuximab",]

theme_paper <- function(){ 
  
  theme_classic() + 
    theme(
      plot.title = element_text(hjust = 0.5, size=20),
      legend.position = "right",
      #legend.position.inside = c(0.6, 0.95),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 14),
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 14),
      strip.text = element_text(size = 20, color = "dark green")
      
    )
}

theme_paper2 <- function(){
  theme(axis.text.y=element_text(size = 6),
        axis.title.y=element_text(size = 8),
        axis.text.x = element_text(size = 6),
        axis.title.x = element_text(size = 8),
        strip.text.x = element_text(size = 4))
}




# KM
# Both arms
survminer::ggsurvplot(survfit(Surv(years, d) ~ treat, data=cetux)) + xlab("Years")
# Control
km_control <- survminer::ggsurvplot(survfit(Surv(years, d) ~ 1, data=control))$data.survplot
# Cetuximab
km_cetux <- survminer::ggsurvplot(survfit(Surv(years, d) ~ 1, data=cetuximab))$data.survplot

tmax <- max(cetux$years)

# Whether to run single arm models or load in predictions from disk
run_single_arm_models <- FALSE

# Whether to run two arm models or load in predictions from disk
run_two_arm_models <- FALSE


# Prior specification
mspline_base <- mspline_spec(Surv(years, d) ~ 1, data=cetux, df=10)
prior_hscale <- p_meansurv(median=25, upper=100, mspline=mspline_base) # prior for the baseline log hazard scale parameter. We'll leave this untouched once calculated (from 10df mspline)
prior_haz_const(mspline_base, prior_hscale = prior_hscale)

# memoise
survextrap_mem <- memoise(survextrap, cache = cachem::cache_disk(dir = ".cache",
                                                                 max_age = 31557600, # keep for 1-year
                                                                 max_size = 2e10)) # 20Gb of space 
rmst_mem <- memoise(rmst_fast, cache = cachem::cache_disk(dir = ".cache",
                                                     max_age = 31557600,
                                                     max_size = 2e10))
irmst_mem <- memoise(irmst_fast, cache = cachem::cache_disk(dir = ".cache",
                                                     max_age = 31557600,
                                                     max_size = 2e10))
survival_mem <- memoise(survival, cache = cachem::cache_disk(dir = ".cache",
                                                             max_age = 31557600,
                                                             max_size = 2e10))
hazard_mem <- memoise(hazard, cache = cachem::cache_disk(dir = ".cache",
                                                         max_age = 31557600,
                                                         max_size = 2e10))
hazard_ratio_mem <- memoise(hazard_ratio, cache = cachem::cache_disk(dir = ".cache",
                                                         max_age = 31557600,
                                                         max_size = 2e10))

### Running models
options(mc.cores = 1) 
# options(mc.cores = parallel::detectCores())
chains <- 4; iter <- 2000
smooth_model <- "exchangeable"

########################################################################################
# Single arm models (Control arm only)
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
                                   survival_mem(models_control[[k]], tmax=40))
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
                                  hazard_mem(models_control[[k]], t=c(seq(0,6,length=100),7:40)))
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

#####################################################
# LOOIC table, calculated within trial follow-up
#####################################################

looic_table <- looic %>% 
  mutate(looic = round(looic, 1), 
         rmst = sprintf("%.2f (%.2f, %.2f)", round(median,2), 
                        round(`lower`, 2), round(`upper`, 2)),
         prior_rate = paste0("Gamma(2, ", prior_rate, ")")
  ) %>%
  select(df, ek, prior_rate , looic, rmst) 

write_csv(looic_table, "single_arm_looic_table.csv")


looic_table %>%
  knitr::kable(col.names = c("df", "Extra knots", "Prior for sigma", "LOOIC", "Restricted mean survival (5 years)"),
               format = "rst")



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



#########################################################
# Single arm, manuscript plot with both survival and hazard.
#########################################################

colour_fill <- "deepskyblue3"
colour_KM <- "gray20"
margins1 <- unit(c(0.0,0.0,0.0,0), "cm")
margins2 <- unit(c(0.0,0.4,0.0,0), "cm")
hjust_vec <- c(-0.8/.pt, -4.5/.pt, 1.5/.pt)

dataset <- c("Trial only", 
             "Trial + Population rates",
             "Trial + Population rates + Registry")

title_vec <- c("(a) Trial data only",
               "(b) Trial data and population rates",
               "(c) Trial, population rates and SEER registry data")

population_data_plot <- function(plot, dataset){
  
  bh_label <- tibble(x=37, y=0.078, lab = "Population \nmortality")
  seer_label <- tibble(x=15, y=0.3, lab = "SEER registry \ndata")
  
  if(dataset == "Trial only"){
    
    plot 
    
  } else if (dataset == "Trial + Population rates") {
    
    plot+
      geom_step(data=cetux_bh %>% filter(time < 40), 
                aes(x=time, y=hazard), inherit.aes = FALSE, colour="gray30") +
      geom_text(data = bh_label, aes(x=x, y=y, label = lab), 
                colour = "gray30", 
                size = 8/.pt, lineheight = 0.85) 
    
  } else if (dataset == "Trial + Population rates + Registry") {
    
    plot+
      geom_step(data=cetux_seer, 
                aes(x=start, y=haz), inherit.aes = FALSE, colour= "brown4") +
      geom_step(data=cetux_bh %>% filter(time < 40), 
                aes(x=time, y=hazard), inherit.aes = FALSE, colour="gray30") +
      geom_text(data = bh_label, aes(x=x, y=y, label = lab), 
                size = 8/.pt, colour="gray30", lineheight = 0.85) +
      geom_text(data = seer_label, aes(x=x, y=y, label = lab), 
                size = 8/.pt, colour= "brown4", lineheight = 0.85 ) 
    
    
  }
  
}

  
  
for(i in 1:3){
  
  
  surv_plot_single_arm <- surv_plots %>% 
    filter(df == 10, prior_rate==1)  %>%
    filter(`Extra knots` == "10,25") %>%
    filter(`External data` == dataset[i]) %>%
    ggplot() +
    theme_classic()+
    theme(axis.text.y=element_text(size = 6),
          axis.title.y=element_text(size = 8),
          axis.text.x = element_text(size = 6),
          axis.title.x = element_text(size = 8),
          strip.text.x = element_text(size = 4),
          plot.margin = margins1)+
    geom_ribbon(aes(x = t, ymin = lower, ymax = upper), fill= colour_fill, alpha = 0.12)+
    geom_line(aes(x=t, y=median), lwd=1, colour = colour_fill) +
    geom_line(aes(x=t, y=lower), linetype = "solid", alpha = 0.2, colour = colour_fill) +
    geom_line(aes(x=t, y=upper), linetype = "solid", alpha = 0.2, colour = colour_fill) +
    geom_step(data=km_control, aes(x=time, y=surv), lwd=1, colour = colour_KM, 
              inherit.aes = FALSE, alpha = 0.75) +
    geom_vline(xintercept = max(control$years), colour = "gray30", linetype = "dashed") +
    xlab("Time (Years)") + 
    scale_y_continuous("Overall Survival", limits = c(0,1), labels = scales::percent) 
  
  
  haz_plot_single_arm <- haz_plots %>% 
    filter(t > 0, t <= 40) %>%
    filter(df == 10, prior_rate==1)  %>%
    filter(`Extra knots` == "10,25") %>%
    filter(`External data` == dataset[i]) %>%
    ggplot() +
    theme_classic()+
    theme(axis.text.y=element_text(size = 6),
          axis.title.y=element_text(size = 8),
          axis.text.x = element_text(size = 6),
          axis.title.x = element_text(size = 8),
          strip.text.x = element_text(size = 4),
          plot.margin = margins2)+
    geom_ribbon(aes(x = t, ymin = lower, ymax = upper), fill= colour_fill, alpha = 0.12)+
    geom_line(aes(x=t, y=median), lwd=1, colour = colour_fill) +
    geom_line(aes(x=t, y=lower), linetype = "solid", alpha = 0.2, colour = colour_fill) +
    geom_line(aes(x=t, y=upper), linetype = "solid", alpha = 0.2, colour = colour_fill) +
    geom_vline(xintercept = max(control$years), colour = "gray30", linetype = "dashed")#%.>%
  
  #haz_plot_single_arm  
  haz_plot_single_arm <- haz_plot_single_arm %.>%
    population_data_plot(plot = ., dataset = dataset[i]) +
    xlab("Time (Years)") + 
    scale_y_continuous("Hazard", limits = c(0, 0.75))
  
  
  assign(paste0("surv_plot_single_arm", i), surv_plot_single_arm)
  assign(paste0("haz_plot_single_arm", i), haz_plot_single_arm)
  
  plot_surv_haz_single_arm <- surv_plot_single_arm + 
    theme(plot.title = element_text(hjust = hjust_vec[i],
                                    size=8, face="bold",
                                    lineheight = 1.2),
          plot.tag = element_text(size = 8))+
    labs(title = title_vec[i])+
    haz_plot_single_arm +
    plot_layout(nrow = 1, ncol = 2)
  
  #plot_surv_haz_single_arm 
  assign(paste0("plot_surv_haz_single_arm", i), plot_surv_haz_single_arm)
  
}

plot_all_single_arm <- plot_grid(
  plot_surv_haz_single_arm1,
  NULL,
  plot_surv_haz_single_arm2,
  NULL,
  plot_surv_haz_single_arm3,
  align = "hv",
  rel_heights = c(1,-0.05,1,-0.05,1),
  ncol = 1)

plot_all_single_arm

tiff(file = "figure_survival_hazard_single_arm.tiff",   
     width = 5.8, 
     height = 5.3,
     units = 'in',  
     res = 300, 
     compression = "lzw")
print(plot_all_single_arm)
dev.off()




################################################################################
## Two-arm modelling
## Just considering models with no extra knots, and models with extra knots at 10, 15 and 25 years
if(run_two_arm_models){
  two_arm_models <- list()
  models_cetux <- list()
  models_cont <- list()
  
  # We define the default model and investigate 1-way sensitivity analyses to this set-up
  df_base <- 10 # degrees of freedom
  ek_base <- NULL # extra knots after end of trial follow-up
  l_base <- 1 # rate for sigma prior 
  prior_hsd_base <- p_gamma(2, l_base) # changing prior for sigma parameter
  mspline_base <- mspline_spec(Surv(years, d) ~ 1, data=cetux, df=df_base, add_knots=ek_base)
  external_data_base <- "Trial only" # no background rates and no SEER data for default
  external_base <- NULL
  backhaz_base <- NULL
  
  
  # Prior for hazard ratio for treatment effects
  # Weak prior. Prior median of 1 and upper 95% confidence interval of 50
  prior_loghr_base <- p_hr(median=1, upper=50)
  prior_hr(prior_loghr_base)
  
  # Prior for hazard ratio variability tau in NPH model
  # Gamma(2,3) represents judgment that the HR could vary about 10 fold over 20 years
  prior_hrsd_base <- p_gamma(2, 3)
  
  # PH model, Non-PH model and Separate arm modelling
  # extra knots and different data sources
  extra_knots <- list(NULL, 25, c(10, 25), c(10, 15, 25))
  external_data <- c("Trial only", "Trial + Population rates", "Trial + Population rates + Registry")
  for(ek in extra_knots){
    print(paste("extra knots=",paste(ek, collapse=",")))
    for(ed in external_data){
      print(paste("external data=",ed))
      # Spline specification
      mspline <- mspline_spec(Surv(years, d) ~ treat, data=cetux, df=df_base, add_knots=ek)
      
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
      
      # PH modelling
      print(paste0("PH model"))
      two_arm_models[[paste0("model=PH; ", "df=",df_base,"; extra_knots=", 
                             paste(ek, collapse=","),"; prior_rate=", l_base,
                             "; external_data=", ed)]] <-  
        survextrap_mem(Surv(years, d) ~ treat, data=cetux, mspline=mspline,
                       chains=chains, iter=iter,
                       prior_hscale=prior_hscale, prior_hsd = prior_hsd_base,
                       prior_loghr=prior_loghr_base,
                       external=external, 
                       backhaz=backhaz,
                       smooth_model = smooth_model
        )
      
      # NPH modelling
      print(paste0("NPH model"))
      two_arm_models[[paste0("model=NPH; ", "df=",df_base,"; extra_knots=", 
                             paste(ek, collapse=","),"; prior_rate=", l_base,
                             "; external_data=", ed)]] <-  
        survextrap_mem(Surv(years, d) ~ treat, data=cetux, mspline=mspline,
                       chains=chains, iter=iter,
                       nonprop = TRUE,
                       prior_hscale=prior_hscale, prior_hsd = prior_hsd_base,
                       prior_loghr=prior_loghr_base,
                       prior_hrsd=prior_hrsd_base,
                       external=external, 
                       backhaz=backhaz,
                       smooth_model = smooth_model
        )
      
      # Separate arms (Control arm only)
      print(paste0("Separate arms model - control arm"))
      models_cont[[paste0("df=",df_base,"; extra_knots=", 
                          paste(ek, collapse=","),"; prior_rate=", l_base,
                          "; external_data=", ed)]] <- 
        survextrap_mem(Surv(years, d) ~ 1, data=control, mspline=mspline,
                       chains=chains, iter=iter,
                       prior_hscale=prior_hscale, prior_hsd = prior_hsd_base,
                       external=external, 
                       backhaz=backhaz,
                       smooth_model = smooth_model)
      
      
      # Separate arms (Cetuximab arm only)
      # There is no registry data for the cetux arm
      print(paste0("Separate arms model - cetux arm"))
      models_cetux[[paste0("df=",df_base,"; extra_knots=", 
                           paste(ek, collapse=","),"; prior_rate=", l_base,
                           "; external_data=", ed)]] <-  
        survextrap_mem(Surv(years, d) ~ 1, data=cetuximab, mspline=mspline,
                       chains=chains, iter=iter,
                       prior_hscale=prior_hscale, prior_hsd = prior_hsd_base,
                       external=external_base,   # this must be NULL to reflect no external data for cetux arm
                       backhaz=backhaz,
                       smooth_model = smooth_model
        )
    }
  }
  
  #########################################
  ## RMST and difference in RMST at 40 years
  # rmst metrics at 40 years
  rmst_two_arms <- irmst_two_arms <- data.frame()
  rmst_cont <- data.frame()
  rmst_cetux <- data.frame()
  
  # Separate arms modelling (control arm)
  for(k in seq_along(models_cont)){
    print(k)
    name <- names(models_cont)[k]
    # Use a regular expression to extract 
    df <- as.numeric(str_extract(name, "(?<=df=)[^;]+"))
    ek <- str_extract(name, "(?<=extra_knots=)[^;]+")
    prior_rate <- as.numeric(str_extract(name, "(?<=prior_rate=)[^;]+"))
    ed <- str_extract(name, "(?<=external_data=)[^;]+")
    
    rmst_cont <- rbind(rmst_cont,
                       data.frame(treat = "Control", df = df, ek = ek, prior_rate = prior_rate, ed = ed,
                                  rmst_mem(models_cont[[k]], t=40, niter = 1000) %>% mutate(wane=0))
    )
  }
  save(rmst_cont, file =  "rmst_cont.Rdata") # save rmst stats to disk
  # Separate arms modelling (cetux arm)
  for(k in seq_along(models_cetux)){
    print(k)
    name <- names(models_cetux)[k]
    # Use a regular expression to extract 
    df <- as.numeric(str_extract(name, "(?<=df=)[^;]+"))
    ek <- str_extract(name, "(?<=extra_knots=)[^;]+")
    prior_rate <- as.numeric(str_extract(name, "(?<=prior_rate=)[^;]+"))
    ed <- str_extract(name, "(?<=external_data=)[^;]+")
    
    rmst_cetux <- rbind(rmst_cetux,
                        data.frame(treat = "Cextuximab", df = df, ek = ek, prior_rate = prior_rate, ed = ed,
                                   rmst_mem(models_cetux[[k]], t=40, niter = 1000) %>% mutate(wane=0))
    )
    
  }
  save(rmst_cetux, file =  "rmst_cetux.Rdata") # save rmst stats to disk
  # Separate arms modelling (irmst)
  for(k in seq_along(models_cetux)){
    print(k)
    name <- names(models_cetux)[k]
    # Use a regular expression to extract 
    df <- as.numeric(str_extract(name, "(?<=df=)[^;]+"))
    ek <- str_extract(name, "(?<=extra_knots=)[^;]+")
    prior_rate <- as.numeric(str_extract(name, "(?<=prior_rate=)[^;]+"))
    ed <- str_extract(name, "(?<=external_data=)[^;]+")
    
    irmst_sam <- rmst_mem(models_cetux[[k]], t=40, niter = 1000, sample = TRUE) - 
      rmst_mem(models_cont[[k]], t=40, niter = 1000, sample = TRUE)
    res <- survextrap:::summarise_output(irmst_sam, t=40, summ_fns = NULL, newdata = NULL,
                                         summ_name = "irmst", sample = FALSE)
    irmst_two_arms <- rbind(irmst_two_arms,
                            data.frame(model = "Separate fits", df = df, ek = ek, prior_rate = prior_rate, ed = ed,
                                       res %>% mutate(wane=0))
    )
    
  }
  
  # PH and NPH models
  for(k in seq_along(two_arm_models)){
    print(k)
    name <- names(two_arm_models)[k]
    # Use a regular expression to extract 
    model <- str_extract(name, "(?<=model=)[^;]+")
    df <- as.numeric(str_extract(name, "(?<=df=)[^;]+"))
    ek <- str_extract(name, "(?<=extra_knots=)[^;]+")
    prior_rate <- as.numeric(str_extract(name, "(?<=prior_rate=)[^;]+"))
    ed <- str_extract(name, "(?<=external_data=)[^;]+")
    
    rmst_two_arms <- rbind(rmst_two_arms,
                           data.frame(model = model, df = df, ek = ek, prior_rate = prior_rate, ed = ed,
                                      rmst_mem(two_arm_models[[k]], t=40, niter = 1000) %>% mutate(wane=0))
    )
    irmst_two_arms <- rbind(irmst_two_arms,
                            data.frame(model = model, df = df, ek = ek, prior_rate = prior_rate, ed = ed,
                                       irmst_mem(two_arm_models[[k]], t=40, niter = 1000) %>% mutate(wane=0))
    )
    # Waning at 6, 10 and 20 years for PH model only 
    if(model == "PH"){
      nd <- data.frame(treat=c("Cetuximab"))
      nd0 <- data.frame(treat=c("Control"))
      ind <- data.frame(treat=c("Control", "Cetuximab"))
      ind0 <- data.frame(treat=c("Control", "Control"))
      for(w in c(6, 10, 20)){
        print(paste0("waning =", w))
        rmst_two_arms <- rbind(rmst_two_arms,
                               data.frame(model = model, df = df, ek = ek, prior_rate = prior_rate, ed = ed,
                                          rmst_mem(two_arm_models[[k]], t=40, niter = 1000, 
                                                   wane_period = c(tmax, w), newdata = nd, newdata0 = nd0) %>% 
                                            mutate(wane=w))
        )
        
        irmst_two_arms <- rbind(irmst_two_arms,
                                data.frame(model = model, df = df, ek = ek, prior_rate = prior_rate, ed = ed,
                                           irmst_mem(two_arm_models[[k]], t=40, niter = 1000,
                                                     wane_period = c(tmax, w), newdata = ind, newdata0 = ind0) %>%
                                             mutate(wane=w))
        )
      }
    }
  }    
  save(rmst_two_arms, file =  "rmst_two_arms.Rdata") # save rmst stats to disk
  save(irmst_two_arms, file =  "irmst_two_arms.Rdata") # save rmst stats to disk
  
  
  ## Survival plots
  surv_plots_two_arms <- data.frame()
  for(k in seq_along(two_arm_models)){
    print(k)
    name <- names(two_arm_models)[k]
    # Use regular expressions
    model <- str_extract(name, "(?<=model=)[^;]+")
    df <- as.numeric(str_extract(name, "(?<=df=)[^;]+"))
    ek <- str_extract(name, "(?<=extra_knots=)[^;]+")
    prior_rate <- as.numeric(str_extract(name, "(?<=prior_rate=)[^;]+"))
    ed <- str_extract(name, "(?<=external_data=)[^;]+")
    
    surv_plots_two_arms <- rbind(surv_plots_two_arms,
                                 data.frame(model = model, df = df, ek = ek, prior_rate = prior_rate, ed = ed,
                                            survival_mem(two_arm_models[[k]], tmax=40))
    )
    
  }
  surv_plots_two_arms$df <- factor(surv_plots_two_arms$df)
  surv_plots_two_arms$`Extra knots` <- factor(ifelse(is.na(surv_plots_two_arms$ek), "None", surv_plots_two_arms$ek))
  surv_plots_two_arms$`Prior for sigma` <- factor(surv_plots_two_arms$prior_rate)
  surv_plots_two_arms$`External data` <- factor(surv_plots_two_arms$ed, levels = c("Trial only", "Trial + Population rates", "Trial + Population rates + Registry"))
  save(surv_plots_two_arms, file =  "surv_plots_two_arms.Rdata") # save survival stats to disk
  
  # Cetux arm
  surv_plots_cetux <- data.frame()
  for(k in seq_along(models_cetux)){
    print(k)
    name <- names(models_cetux)[k]
    # Use regular expressions
    model <- str_extract(name, "(?<=model=)[^;]+")
    df <- as.numeric(str_extract(name, "(?<=df=)[^;]+"))
    ek <- str_extract(name, "(?<=extra_knots=)[^;]+")
    prior_rate <- as.numeric(str_extract(name, "(?<=prior_rate=)[^;]+"))
    ed <- str_extract(name, "(?<=external_data=)[^;]+")
    
    surv_plots_cetux <- rbind(surv_plots_cetux,
                                 data.frame(df = df, ek = ek, prior_rate = prior_rate, ed = ed,
                                            survival_mem(models_cetux[[k]], tmax=40))
    )
    
  }
  surv_plots_cetux$df <- factor(surv_plots_cetux$df)
  surv_plots_cetux$`Extra knots` <- factor(ifelse(is.na(surv_plots_cetux$ek), "None", surv_plots_cetux$ek))
  surv_plots_cetux$`Prior for sigma` <- factor(surv_plots_cetux$prior_rate)
  surv_plots_cetux$`External data` <- factor(surv_plots_cetux$ed, levels = c("Trial only", "Trial + Population rates", "Trial + Population rates + Registry"))
  save(surv_plots_cetux, file =  "surv_plots_cetux.Rdata") # save survival stats to disk

  #########################################
  ## Hazard plots
  haz_plots_two_arms <- data.frame()
  for(k in seq_along(two_arm_models)){
    print(k)
    name <- names(two_arm_models)[k]
    # Use regular expressions
    model <- str_extract(name, "(?<=model=)[^;]+")
    df <- as.numeric(str_extract(name, "(?<=df=)[^;]+"))
    ek <- str_extract(name, "(?<=extra_knots=)[^;]+")
    prior_rate <- as.numeric(str_extract(name, "(?<=prior_rate=)[^;]+"))
    ed <- str_extract(name, "(?<=external_data=)[^;]+")
    
    haz_plots_two_arms <- rbind(haz_plots_two_arms,
                                 data.frame(model = model, df = df, ek = ek, prior_rate = prior_rate, ed = ed,
                                            hazard_mem(two_arm_models[[k]], tmax=40))
    )
    
  }
  haz_plots_two_arms$df <- factor(haz_plots_two_arms$df)
  haz_plots_two_arms$`Extra knots` <- factor(ifelse(is.na(haz_plots_two_arms$ek), "None", haz_plots_two_arms$ek))
  haz_plots_two_arms$`Prior for sigma` <- factor(haz_plots_two_arms$prior_rate)
  haz_plots_two_arms$`External data` <- factor(haz_plots_two_arms$ed, levels = c("Trial only", "Trial + Population rates", "Trial + Population rates + Registry"))
  save(haz_plots_two_arms, file =  "haz_plots_two_arms.Rdata") # save hazard stats to disk
  
  # Cetux arm
  haz_plots_cetux <- data.frame()
  for(k in seq_along(models_cetux)){
    print(k)
    name <- names(models_cetux)[k]
    # Use regular expressions
    model <- str_extract(name, "(?<=model=)[^;]+")
    df <- as.numeric(str_extract(name, "(?<=df=)[^;]+"))
    ek <- str_extract(name, "(?<=extra_knots=)[^;]+")
    prior_rate <- as.numeric(str_extract(name, "(?<=prior_rate=)[^;]+"))
    ed <- str_extract(name, "(?<=external_data=)[^;]+")
    
    haz_plots_cetux <- rbind(haz_plots_cetux,
                              data.frame(df = df, ek = ek, prior_rate = prior_rate, ed = ed,
                                         hazard_mem(models_cetux[[k]], tmax=40))
    )
    
  }
  haz_plots_cetux$df <- factor(haz_plots_cetux$df)
  haz_plots_cetux$`Extra knots` <- factor(ifelse(is.na(haz_plots_cetux$ek), "None", haz_plots_cetux$ek))
  haz_plots_cetux$`Prior for sigma` <- factor(haz_plots_cetux$prior_rate)
  haz_plots_cetux$`External data` <- factor(haz_plots_cetux$ed, levels = c("Trial only", "Trial + Population rates", "Trial + Population rates + Registry"))
  save(haz_plots_cetux, file =  "haz_plots_cetux.Rdata") # save hazard stats to disk
  
  #########################################
  ## Hazard Ratio plots
  haz_ratio_plots <- data.frame()
  # PH and NPH models
  for(k in seq_along(two_arm_models)){
    print(k)
    name <- names(two_arm_models)[k]
    # Use regular expressions
    model <- str_extract(name, "(?<=model=)[^;]+")
    df <- as.numeric(str_extract(name, "(?<=df=)[^;]+"))
    ek <- str_extract(name, "(?<=extra_knots=)[^;]+")
    prior_rate <- as.numeric(str_extract(name, "(?<=prior_rate=)[^;]+"))
    ed <- str_extract(name, "(?<=external_data=)[^;]+")
    
    haz_ratio_plots <- rbind(haz_ratio_plots,
                                data.frame(model = model, df = df, ek = ek, prior_rate = prior_rate, ed = ed,
                                           hazard_ratio_mem(two_arm_models[[k]], tmax=40))
    )
    
  }
  # Hazard Ratio separate arm modelling
  for(k in seq_along(models_cetux)){
    print(k)
    name <- names(models_cetux)[k]
    # Use a regular expression to extract 
    df <- as.numeric(str_extract(name, "(?<=df=)[^;]+"))
    ek <- str_extract(name, "(?<=extra_knots=)[^;]+")
    prior_rate <- as.numeric(str_extract(name, "(?<=prior_rate=)[^;]+"))
    ed <- str_extract(name, "(?<=external_data=)[^;]+")
    
    hr_sam <- hazard_mem(models_cetux[[k]], tmax=40, niter = 1000, sample = TRUE)/hazard_mem(models_cont[[k]], tmax=40, niter = 1000, sample = TRUE)
    res <- survextrap:::summarise_output(hr_sam, t = seq(0,40, length=101), summ_fns = NULL, newdata = NULL,
                                         sample = FALSE)
    haz_ratio_plots <- rbind(haz_ratio_plots,
                             data.frame(model = "Separate fits", df = df, ek = ek, prior_rate = prior_rate, ed = ed,
                                        res)
    )
  }
  haz_ratio_plots$df <- factor(haz_ratio_plots$df)
  haz_ratio_plots$`Extra knots` <- factor(ifelse(is.na(haz_ratio_plots$ek), "None", haz_ratio_plots$ek))
  haz_ratio_plots$`Prior for sigma` <- factor(haz_ratio_plots$prior_rate)
  haz_ratio_plots$`External data` <- factor(haz_ratio_plots$ed, levels = c("Trial only", "Trial + Population rates", "Trial + Population rates + Registry"))
  save(haz_ratio_plots, file =  "haz_ratio_plots.Rdata") # save hazard ratio stats to disk
  
    
} else{
  # Load from single arms too.
  load("surv_plots.Rdata") # load survival stats
  load("haz_plots.Rdata") # load survival stats
  
  
  load("surv_plots_two_arms.Rdata") # load survival stats
  load("surv_plots_cetux.Rdata") # load survival stats
  load("haz_plots_two_arms.Rdata") # load hazard stats
  load("haz_plots_cetux.Rdata") # load hazard stats
  load("haz_ratio_plots.Rdata") # load hazard stats
  load("rmst_cont.Rdata") 
  load("rmst_cetux.Rdata")
  load("rmst_two_arms.Rdata") # load rmst stats
  load("irmst_two_arms.Rdata") # load rmst stats
}

# Survival varying extra knots across trial only, trial + BH, trial + BH + seer
surv_plots_combined <-
  rbind(surv_plots_two_arms,
        surv_plots %>% mutate(model = "Separate fits", treat="Control"),
        surv_plots_cetux %>% mutate(model = "Separate fits", treat="Cetuximab")                           
) %>%
  mutate(model = factor(model, levels = c("PH", "NPH", "Separate fits"),
                        labels = c("Proportional \nhazards", 
                                   "Non-proportional \nhazards",
                                   "Separate arms"))) %>%
  mutate("External data" = factor(ed, levels = c("Trial only", "Trial + Population rates", "Trial + Population rates + Registry"),
                     labels = c("Trial only", "Trial + \nPopulation rates", "Trial + \nPopulation rates + \nRegistry")))
                          
s5 <- ggplot(surv_plots_combined %>% filter(df == 10 & prior_rate==1 & `Extra knots`!="10,25")) +
  geom_line(aes(x=t, y=median, 
                col=`Extra knots`, linetype = treat), lwd=1, alpha = 0.7) +
  facet_grid(model ~`External data`) +
  geom_step(data=km_control, aes(x=time, y=surv), lwd=1, linetype = "dashed",
            inherit.aes = FALSE, alpha = 0.7) +
  geom_step(data=km_cetux, aes(x=time, y=surv), lwd=1,
            inherit.aes = FALSE, alpha = 0.6) +
  xlab("Time (Years)") + 
  scale_y_continuous("Overall Survival", limits = c(0,1), labels = scales::percent) +
  scale_linetype_discrete("Treatment") +
  geom_vline(xintercept = max(control$years), colour = "gray30", linetype = "dashed") +
  theme_classic()+
  theme_paper2()+
  theme(plot.title = element_text(hjust = 0))+
  theme(axis.title.x = element_text(size =8),
        strip.text.x = element_text(size = 8),
        strip.text.y = element_text(size = 8))+  
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "grey20"),
        legend.position = "bottom",
        legend.title=element_text(size=10,
                                  margin = unit(c(0.0,0.2,0.0,0), "cm")),
         legend.text = element_text(size=8),
         legend.spacing.x = unit(0.1, 'cm'),
         legend.margin = margin(unit(4*c(1,1,1,1), "cm")),
          legend.box.margin = unit(c(0.0,0.2,0.0,0), "cm"), 
         legend.background = element_rect(linetype = 1, 
                                          linewidth = 0.5, 
                                          colour = "gray30"))+
  guides(linetype  = guide_legend(order = 1))

s5

tiff(file = "figure_survival_cetux_two_arm.tiff",   
     width = 5.4, 
     height = 4.6,
     units = 'in',  
     res = 300, 
     compression = "lzw")
print(s5)
dev.off()


# Hazard varying extra knots across trial only, trial + BH, trial + BH + seer
haz_plots_combined <- rbind(haz_plots_two_arms,
                            haz_plots %>% mutate(model = "Separate fits", treat="Control"),
                             haz_plots_cetux %>% mutate(model = "Separate fits", treat="Cetuximab")
) %>%
  mutate(model = factor(model, levels = c("PH", "NPH", "Separate fits"),
                        labels = c("Proportional \nhazards", 
                                   "Non-proportional \nhazards",
                                   "Separate arms"))) %>%
  mutate("External data" = factor(ed, levels = c("Trial only", "Trial + Population rates", "Trial + Population rates + Registry"),
                                  labels = c("Trial only", "Trial + \nPopulation rates", "Trial + \nPopulation rates + \nRegistry")))


cetux_seer_plot_two_arm <- cbind(cetux_seer, `External data` =  factor("Trial + Population rates + Registry", 
                                                                       levels = c("Trial only", "Trial + Population rates", "Trial + Population rates + Registry"),
                                                                       labels = c("Trial only", "Trial + \nPopulation rates", "Trial + \nPopulation rates + \nRegistry")))

cetux_bh_plot_two_arm <- rbind(
  cbind(cetux_bh, `External data` =  factor("Trial + Population rates", 
                                            levels = c("Trial only", "Trial + Population rates", "Trial + Population rates + Registry"),
                                            labels = c("Trial only", "Trial + \nPopulation rates", "Trial + \nPopulation rates + \nRegistry"))),
  cbind(cetux_bh, `External data` =  factor("Trial + Population rates + Registry", 
                                            levels = c("Trial only", "Trial + Population rates", "Trial + Population rates + Registry"),
                                            labels = c("Trial only", "Trial + \nPopulation rates", "Trial + \nPopulation rates + \nRegistry")))
)

#cetux_bh_plot_two_arm
bh_text_two_arm <- tibble(x=35, y=0.06, lab = "Population \nmortality",
                  `External data` =  factor(c("Trial + Population rates", "Trial + Population rates + Registry"), 
                                            levels = c("Trial only", "Trial + Population rates", "Trial + Population rates + Registry"),
                                            labels = c("Trial only", "Trial + \nPopulation rates", "Trial + \nPopulation rates + \nRegistry")))
#bh_text_two_arm
seer_text_two_arm <- tibble(x=20, y=0.3, lab = "Registry data",
                    `External data` =  factor("Trial + Population rates + Registry", 
                                              levels = c("Trial only", "Trial + Population rates", "Trial + Population rates + Registry"),
                                              labels = c("Trial only", "Trial + \nPopulation rates", "Trial + \nPopulation rates + \nRegistry")))


h5 <- ggplot(haz_plots_combined %>% filter(df == 10 & prior_rate==1 & `Extra knots`!="5,10,25")) +
  geom_line(aes(x=t, y=median, 
                col=`Extra knots`, linetype = treat), lwd=1) +
  geom_step(data=cetux_seer_plot_two_arm, 
            aes(x=start, y=haz), inherit.aes = FALSE, colour= "brown4") +
  geom_step(data=cetux_bh_plot_two_arm %>% filter(time < 40), 
            aes(x=time, y=hazard), inherit.aes = FALSE, colour="gray30") +
  geom_text(data = bh_text_two_arm, aes(x=x, y=y, label = lab), 
            size = 6/.pt, colour="gray30", lineheight = 0.85) +
  geom_text(data = seer_text_two_arm, aes(x=x, y=y, label = lab),
            size = 6/.pt, colour= "brown4", lineheight = 0.85 ) +
  facet_grid(model~`External data`) +
  xlab("Time (Years)") + 
  scale_y_continuous("Hazard") +
  scale_linetype_discrete("Treatment") +
  geom_vline(xintercept = max(control$years), colour = "gray30", linetype = "dashed") +
  # ggtitle("Varying extra knot locations after trial follow-up")+
  theme_classic()+
  theme_paper2()+
  theme(plot.title = element_text(hjust = 0))+
  theme(axis.title.x = element_text(size =8),
        strip.text.x = element_text(size = 8),
        strip.text.y = element_text(size = 8))+  
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "grey20"),
        legend.position = "bottom",
        legend.title=element_text(size=10,
                                  margin = unit(c(0.0,0.2,0.0,0), "cm")),
        legend.text = element_text(size=8),
        legend.spacing.x = unit(0.1, 'cm'),
        legend.margin = margin(unit(4*c(1,1,1,1), "cm")),
        legend.box.margin = unit(c(0.0,0.2,0.0,0), "cm"), 
        legend.background = element_rect(linetype = 1, 
                                         linewidth = 0.5, 
                                         colour = "gray30"))+
  guides(linetype  = guide_legend(order = 1))

h5

tiff(file = "figure_hazard_cetux_two_arm.tiff",   
     width = 5.4, 
     height = 4.6,
     units = 'in',  
     res = 300, 
     compression = "lzw")
print(h5)
dev.off()


# Hazard Ratio plots
hr1 <- haz_ratio_plots %>%
  mutate(model = factor(model, levels = c("PH", "NPH", "Separate fits"),
                        labels = c("Proportional \nhazards", 
                                   "Non-proportional \nhazards",
                                   "Separate arms"))) %>%
  mutate("External data" = factor(
    ed, 
    levels = c("Trial only", "Trial + Population rates", "Trial + Population rates + Registry"),
    labels = c("Trial only", "Trial + \nPopulation rates", "Trial + \nPopulation rates + \nRegistry"))) %>%
  filter(t > 0)  %>%
  ggplot() +
  geom_line(aes(x=t, y=median, 
                col=`Extra knots`), lwd=1) +
  geom_line(aes(x=t, y=upper, col=`Extra knots`), linetype = "dashed") +
  geom_line(aes(x=t, y=lower, col=`Extra knots`), linetype = "dashed") +
  facet_grid(model~`External data`) +
  xlab("Years after diagnosis") + 
  ylab("Hazard ratio") +
  geom_vline(xintercept = max(control$years)) +
  geom_vline(xintercept = c(10,25), linetype = "dotted" , colour = "gray30") +
  geom_hline(yintercept = 1) +
  xlab("Time (Years)") + 
  scale_y_continuous(breaks = c(0.25, 0.5, 1, 2, 4)) +
  coord_trans(y = "log10", ylim=c(0.2, 4)) +
  theme_classic()+
  theme_paper2()+
  theme(plot.title = element_text(hjust = 0))+
  theme(axis.title.x = element_text(size =8),
        strip.text.x = element_text(size = 8),
        strip.text.y = element_text(size = 8))+  
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "grey20"),
        legend.position = "bottom",
        legend.title=element_text(size=10,
                                  margin = unit(c(0.0,0.2,0.0,0), "cm")),
        legend.text = element_text(size=8),
        legend.spacing.x = unit(0.1, 'cm'),
        legend.margin = margin(unit(4*c(1,1,1,1), "cm")),
        legend.box.margin = unit(c(0.0,0.2,0.0,0), "cm"), 
        legend.background = element_rect(linetype = 1, 
                                         linewidth = 0.5, 
                                         colour = "gray30"))

hr1
tiff(file = "figure_hr_cetux_two_arm.tiff",   
     width = 5.4, 
     height = 4.6,
     units = 'in',  
     res = 300, 
     compression = "lzw")
print(hr1)
dev.off()



## Difference in RMST table 

irmst_two_arms %>%
  filter(wane == 0, model == "PH")

irmst_two_arms %>% 
  mutate(irmst = sprintf("%.2f (%.2f, %.2f)", round(median,2), 
                 round(lower,2), round(upper,2))) %>%
  select(-c(df, prior_rate, variable, t, median, lower, upper)) %>%
  pivot_wider(names_from = ed, values_from = irmst) %>%
  mutate(model = factor(model, levels = c("PH", "NPH", "Separate fits"))) %>%
  arrange(wane, model) %>%
  mutate(wane = ifelse(wane==0, NA, wane)) %>%
  rename(Model = model) %>%
  rename(`Extra knots` = ek) %>%
  rename(`Waning time (years)` = wane) %>%
  flextable() %>% 
  add_header_row(colwidths = c(3, 3),
                 values = c("", "Difference in 40-year RMST (years)"),
                 
  ) %>% 
  # perhaps avoid saving binary files to repository.
  # save_as_docx(path = "aim2-tables/cetux_case_study_irmst.docx")
  save_as_docx(path = "cetux_case_study_irmst.docx")




#########################################################
# Two arms, manuscript plot with both survival and HRs.
#########################################################


two_arm_margins1 <- unit(c(0.0,0.0,0.0,0), "cm")
two_arm_margins2 <- unit(c(0.0,0.4,0.0,0), "cm")
two_arm_hjust_vec <- c(-1.9/.pt, 1.5/.pt,-1.1/.pt)

two_arm_models <- c("Proportional \nhazards", 
                    "Non-proportional \nhazards",
                    "Separate arms")

two_arm_dataset <- c("Trial only", 
                     "Trial + \nPopulation rates",
                     "Trial + \nPopulation rates + \nRegistry")

two_arm_title_vec <- c("(a) Proportional excess hazards",
               "(b) Flexible non-proportional hazards           ",
               "(c) Separate arms modelling")

two_arm_colour_fill <- "darkgreen"

for(i in 1:3){

  #i <- 1
  legend_text <- element_text(size=8)
  legend_key <- unit(0.3, "lines")
  legend_space_y <- unit(1.5, "pt")

  surv_plot_two_arm <- surv_plots_combined %>% 
    filter(df == 10, prior_rate==1, model == two_arm_models[i])  %>%
    filter(`Extra knots` == "10,15,25") %>%
    filter(`External data` == "Trial + \nPopulation rates + \nRegistry") %>%
    mutate(treat = factor(treat, levels = c("Cetuximab", "Control"),
                          labels = c("Cetuximab+\nRadiotherapy", "Radiotherapy"))) %>%
    mutate(treat_dummy = treat) %>%
    ggplot(aes(colour = treat)) +
    theme_classic()+
    theme(axis.text.y=element_text(size = 6),
          axis.title.y=element_text(size = 8),
          axis.text.x = element_text(size = 6),
          axis.title.x = element_text(size = 8),
          strip.text.x = element_text(size = 4),
          plot.margin = two_arm_margins1,
          legend.position = "inside" ,
          legend.position.inside = c(0.72, 0.75),
          legend.text=legend_text,
          legend.key.size=legend_key,
          legend.key.spacing.y = legend_space_y)+
    geom_ribbon(aes(x = t, ymin = lower, ymax = upper, group = treat, fill = treat), alpha = 0.1, colour = NA)+
    geom_line(aes(x=t, y=median), lwd=0.85) +
    geom_line(aes(x=t, y=lower, group = treat, alpha = treat), linetype = "solid") +
    geom_line(aes(x=t, y=upper, group = treat, alpha = treat), linetype = "solid") +
    geom_vline(xintercept = max(control$years), colour = "gray30", linetype = "dashed") +
    scale_colour_manual(NULL, values = c("chocolate3","royalblue1")) +
    scale_fill_manual(NULL, values =  c("chocolate3","royalblue1")) +
    scale_alpha_manual(NULL, values = c(0.25, 0.5))+
    xlab("Time (Years)") + 
    scale_y_continuous("Overall Survival", limits = c(0,1), labels = scales::percent) 
  
  surv_plot_two_arm
  
  haz_ratio_plot_two_arm <- haz_ratio_plots %>%
    mutate(model = factor(model, levels = c("PH", "NPH", "Separate fits"),
                          labels = c("Proportional \nhazards", 
                                     "Non-proportional \nhazards",
                                     "Separate arms"))) %>%
    mutate("External data" = factor(
      ed, 
      levels = c("Trial only", "Trial + Population rates", "Trial + Population rates + Registry"),
      labels = c("Trial only", "Trial + \nPopulation rates", "Trial + \nPopulation rates + \nRegistry"))) %>%
    filter(df == 10, prior_rate==1, model == two_arm_models[i])  %>%
    filter(`Extra knots` == "10,15,25") %>%
    filter(`External data` == "Trial + \nPopulation rates + \nRegistry") %>%
    filter(t>0) %>%
    ggplot() +
    theme_classic()+
    theme(axis.text.y=element_text(size = 6),
          axis.title.y=element_text(size = 8),
          axis.text.x = element_text(size = 6),
          axis.title.x = element_text(size = 8),
          strip.text.x = element_text(size = 4),
          plot.margin = two_arm_margins2)+
    geom_ribbon(aes(x = t, ymin = lower, ymax = upper), fill= two_arm_colour_fill, alpha = 0.12)+
    geom_line(aes(x=t, y=median), lwd=1, colour = two_arm_colour_fill) +
    geom_line(aes(x=t, y=lower), linetype = "solid", alpha = 0.2, colour = two_arm_colour_fill) +
    geom_line(aes(x=t, y=upper), linetype = "solid", alpha = 0.2, colour = two_arm_colour_fill) +
    geom_vline(xintercept = max(control$years), colour = "gray30", linetype = "dashed") +
    geom_hline(yintercept = 1, linewidth = 1, 
               alpha = 0.55, colour = "gray40",linetype = "solid")+
    xlab("Time (Years)") + 
    scale_y_continuous("Hazard ratio", breaks = c(0.25, 0.5, 1, 2, 4)) +
    coord_trans(y = "log10", ylim=c(0.2, 4))
  
 
  assign(paste0("surv_plot_two_arm", i), surv_plot_two_arm)
  assign(paste0(" haz_ratio_plot_two_arm", i),  haz_ratio_plot_two_arm )
  
  plot_surv_hr_two_arm <- surv_plot_two_arm + 
    theme(plot.title = element_text(hjust = two_arm_hjust_vec[i],
                                    size=8, face="bold",
                                    lineheight = 1.2),
          plot.tag = element_text(size = 8))+
    labs(title = two_arm_title_vec[i])+
    haz_ratio_plot_two_arm +
    plot_layout(nrow = 1, ncol = 2)
  
  plot_surv_hr_two_arm
  assign(paste0("plot_surv_hr_two_arm", i), plot_surv_hr_two_arm)
  
}

plot_all_two_arm <- plot_grid(
  plot_surv_hr_two_arm1,
  NULL,
  plot_surv_hr_two_arm2,
  NULL,
  plot_surv_hr_two_arm3,
  align = "hv",
  rel_heights = c(1,-0.05,1,-0.05,1),
  ncol = 1)

print(plot_all_two_arm)

tiff(file = "figure_survival_hr_two_arm.tiff",   
     width = 5.8, 
     height = 5.3,
     units = 'in',  
     res = 300, 
     compression = "lzw")
print(plot_all_two_arm)
dev.off()


#########################################################
# Two arms, manuscript table with RMST and iRMST.
#########################################################

summary(as.factor(irmst_two_arms$ek))

control_waning <- rmst_two_arms %>% 
  filter(ek == "10,15,25") %>%
  filter(model == "PH") %>%
  filter(treat == "Control") %>%
  select(-wane) %>%
  expand_grid(wane = c(6, 10, 20))

summary(as.factor(irmst_two_arms$ed))

two_arm_table <- irmst_two_arms %>% 
  mutate(treat = "Contrast") %>%
  bind_rows(rmst_two_arms) %>%
  bind_rows(rmst_cont %>% mutate(model = "Separate fits", treat = "Control")) %>%
  bind_rows(rmst_cetux %>% mutate(model = "Separate fits", treat = "Cetuximab")) %>%
  bind_rows(control_waning) %>%
  mutate(treat = factor(treat, levels = c("Cetuximab", "Control", "Contrast"))) %>%
  filter(ek == "10,15,25") %>%
  mutate(estimand = sprintf("%.2f (%.2f, %.2f)", round(median,2), 
                         round(lower,2), round(upper,2))) %>%
  mutate(ed = factor(ed, 
                     levels = c("Trial only", "Trial + Population rates",  "Trial + Population rates + Registry"),
                     labels = c("trial_only", "trial_gpms", "trial_gpm_registry")))  %>%
  select(-c(df, ek, prior_rate, variable, t, median, lower, upper))%>%
  filter(treat %in% c("Control", "Contrast")) %>%
  arrange(ed) %>%
  arrange(treat) %>%
  #pivot_wider(names_from = treat, values_from = c("trial_only", "trial_gpms", "trial_gpm_registry"))# %>%
  pivot_wider(names_from = c(ed, treat), values_from = estimand) %>%
  mutate(model = factor(model, levels = c("PH", "NPH", "Separate fits"))) %>%
  arrange(wane, model) %>%
  mutate(wane = ifelse(wane==0, NA, wane)) %>%
  rename(Model = model) %>%
  rename(`Waning time (years)` = wane)

two_arm_table
write_csv(two_arm_table, "two_arm_irmst_table.csv")



