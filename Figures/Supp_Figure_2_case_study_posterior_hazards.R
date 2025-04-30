library(survextrap)
library(reshape2)
library(dplyr)
library(tibble)
library(tidyr)


wd <- "/projects/aa/statistical_innovation/survextrap/cetux_case_study"
setwd(wd)

df_base <- 10 # degrees of freedom
ek_base <- NULL # extra knots after end of trial follow-up
l_base <- 1 # rate for sigma prior 
prior_hsd_base <- p_gamma(2, l_base) # changing prior for sigma parameter
mspline_base <- mspline_spec(Surv(years, d) ~ 1, data=cetux, df=df_base, add_knots=ek_base)
external_data_base <- "Trial only" # no background rates and no SEER data for default
external_base <- NULL
backhaz_base <- NULL
#

mspline_base <- mspline_spec(Surv(years, d) ~ 1, data=cetux, df=10)
prior_hscale <- p_meansurv(median=25, upper=100, mspline=mspline_base) # prior for the baseline log hazard scale parameter. We'll leave this untouched once calculated (from 10df mspline)
prior_haz_const(mspline_base, prior_hscale = prior_hscale)

chains <- 4; iter <- 2000
smooth_model <- "exchangeable"

extra_knots <- list(NULL, 25, c(10, 25), c(10, 15, 25))

models_control <- NULL

for(ek in extra_knots){
  print(paste("extra knots=",paste(ek, collapse=",")))
  # Spline specification
  mspline <- mspline_spec(Surv(years, d) ~ 1, data=cetux, df=df_base, add_knots=ek)
  
  external <- NULL
  backhaz <- NULL
  
  models_control[[paste0("df=",df_base,"; extra_knots=", paste(ek, collapse=","),"; prior_rate=", l_base,
                         "; external_data=", external_data_base)]] <-  
    survextrap(Surv(years, d) ~ 1, data=control, mspline=mspline,
               chains=chains, iter=iter,
               prior_hscale=prior_hscale, prior_hsd = prior_hsd_base,
               external=external, 
               backhaz=backhaz,
               smooth_model = smooth_model
    )
}

set.seed(17473)

linewidth <- 0.8
alpha <- 0.5
ylim <- c(0,0.78)


# No extra knots
test <- hazard(models_control[[1]], t=c(seq(0,6,length=100),7:40), sample = TRUE)
# sample 20 of the posteriors
samp <- sample(1:dim(test)[2], 100, replace = FALSE)
samp
test2 <- as.data.frame(test[,samp,1]) %>% mutate(t = c(seq(0,6,length=100),7:40))
# test4 <- test2 %>%
#   pivot_longer(cols = starts_with("V"), names_to = "sample", values_to = "hazard")
# head(test4)
# View(test4)  
# head(test3)
test3 <- reshape2::melt(test2, id.vars = "t", variable.name = "sample", value.name = "hazard")
plot1 <-ggplot(test3) +
  geom_line(aes(x=t, y=hazard, 
                col=sample), 
            lwd=linewidth, alpha =alpha) +
  xlab("Time (Years)") + 
  ylab("Hazard") +
  ylim(ylim) +
  geom_vline(xintercept = max(control$years)) +
  theme_classic()+
  theme_paper2()+
  theme(plot.title = element_text(hjust = 0)) +
  theme(legend.position="none")


# extra knots at 25 years
test <- hazard(models_control[[2]], t=c(seq(0,6,length=100),7:40), sample = TRUE)
# sample 20 of the posteriors
samp <- sample(1:dim(test)[2], 20, replace = FALSE)
test2 <- as.data.frame(test[,samp,1]) %>% mutate(t = c(seq(0,6,length=100),7:40))
test3 <- reshape2::melt(test2, id.vars = "t", variable.name = "sample", value.name = "hazard")
plot2 <- ggplot(test3) +
  geom_line(aes(x=t, y=hazard, 
                col=sample),  
            lwd=linewidth, alpha =alpha) +
  xlab("Time (Years)") + 
  ylab("Hazard") +
  ylim(ylim) +
  geom_vline(xintercept = max(control$years)) +
  theme_classic()+
  theme_paper2()+
  theme(plot.title = element_text(hjust = 0)) +
  theme(legend.position="none")


# extra knots at 10 and 25 years
test <- hazard(models_control[[3]], t=c(seq(0,6,length=100),7:40), sample = TRUE)
# sample 20 of the posteriors
samp <- sample(1:dim(test)[2], 20, replace = FALSE)
test2 <- as.data.frame(test[,samp,1]) %>% mutate(t = c(seq(0,6,length=100),7:40))
test3 <- reshape2::melt(test2, id.vars = "t", variable.name = "sample", value.name = "hazard")
plot3 <- ggplot(test3) +
  geom_line(aes(x=t, y=hazard, 
                col=sample),  
            lwd=linewidth, alpha =alpha) +
  xlab("Time (Years)") + 
  ylab("Hazard") +
  ylim(ylim) +
  geom_vline(xintercept = max(control$years)) +
  theme_classic()+
  theme_paper2()+
  theme(plot.title = element_text(hjust = 0)) +
  theme(legend.position="none")

# extra knots at 10, 15 and 25 years
test <- hazard(models_control[[4]], t=c(seq(0,6,length=100),7:40), sample = TRUE)
# sample 20 of the posteriors
samp <- sample(1:dim(test)[2], 20, replace = FALSE)
test2 <- as.data.frame(test[,samp,1]) %>% mutate(t = c(seq(0,6,length=100),7:40))
test3 <- reshape2::melt(test2, id.vars = "t", variable.name = "sample", value.name = "hazard")
plot4 <- ggplot(test3) +
  geom_line(aes(x=t, y=hazard, 
                col=sample),  
            lwd=linewidth, alpha =alpha) +
  xlab("Time (Years)") + 
  ylab("Hazard") +
  ylim(ylim) +
  geom_vline(xintercept = max(control$years)) +
  theme_classic()+
  theme_paper2()+
  theme(plot.title = element_text(hjust = 0)) +
  theme(legend.position="none")


add_theme <-  function(){
  theme(plot.title = element_text(hjust = 0,
                                  size = 10),
        plot.margin = unit(c(.2,0,0,0), "cm"))
  
  }
plot1_label <- plot1+
  theme(legend.position="none",
        plot.title = element_text(size=10, face="bold",
                                  lineheight = 1.2))+
  ggtitle("(a) No extra knots")+
  add_theme()
plot2_label <- plot2+
  theme(legend.position="none",
        plot.title = element_text(size=10, face="bold",
                                  lineheight = 1.2))+
  ggtitle("(b) Extra knot at 25 years")+
  add_theme()
plot3_label <- plot3+
  theme(legend.position="none",
        plot.title = element_text(size=10, face="bold",
                                  lineheight = 1.2))+
  ggtitle("(c) Extra knots at 10 and 25 years")+
  add_theme()
plot4_label <- plot4+
  theme(legend.position="none",
        plot.title = element_text(size=10, face="bold",
                                  lineheight = 1.2))+
  ggtitle("(b) Extra knots at 10, 15 and 25 years")+
  add_theme()

plot_all <- plot_grid(plot1_label, plot2_label, plot3_label, plot4_label,
          ncol = 1)


tiff(file = "figure_posterior_hazard_cetux.tiff",   
     width = 5.0, 
     height = 7.2,
     units = 'in',  
     res = 300, 
     compression = "lzw")
print(plot_all)
dev.off()

