## This R script is written by Xiuhan Zhang (xiuhan.zhang17@imperial.ac.uk) 
## for journal article: 
## "The gibbon umbrella: Acoustic monitoring demonstrates conservation spill-over for tropical rainforest birds"
## Please do not re-distribute without permission.

## When using the script, install R packages in advance.

##Set console language
Sys.setenv(LANG = "en")

##Set working directory
setwd("E:/Coursework/RA/海南鸟类数据/Publication")

##Load packages required
library(readr)
library(dplyr)
library(tidyr)
library(raster)
library(sf)
library("PerformanceAnalytics")
library("FactoMineR")
library("factoextra")
library(vegan)
library(ggplot2)
library(ape)
library(lme4)
library(lmerTest)
library(rsq)
library(caret)
library(emmeans)
library(ggeffects)
library(ggbeeswarm)
library(DHARMa)
library(glmmTMB)
library(MuMIn)
library(car)
library(performance)

# Load CSV dataset created in last section
hainan_all_hour_diet <- read.csv("Data/hainan_all_hour_diet.csv")

# Calculate taxa/foraging guild activity at per recording hour and month at each site
# (for plotting points on figures)
hainan_all_sub_all <- hainan_all_hour_diet%>%#filter(Total>0)%>%
  group_by(Site,Hour,Month)%>%
  summarise(Total=mean(Total),n=mean(n),Rapt=mean(Rapt),Inse=mean(Inse),
            Omni=mean(Omni),Alt=unique(Alt),Slope=unique(Slope),
            Habitat=unique(Habitat),Cover=unique(Cover))

#####
# 1. Spatial+temporal patterns of hourly gibbon activity
zig_gibbon_env <- glmmTMB(n~Habitat+log(Alt)+as.factor(Hour)+as.factor(Month)+
                            (Cover)+as.factor(Slope)+
                            (1|Date)+(1|Group), ziformula=~1, family = poisson, 
                          data = hainan_all_hour_diet)
# Calculate model AIC
AIC(zig_gibbon_env)
# Plot VIF graph of model variables
plot(check_collinearity(zig_gibbon_env,component = c('all')))
# Test for model zero-inflation
testZeroInflation(zig_gibbon_env)
# Model summary
summary(zig_gibbon_env)
# Calculate model R-squared value
r.squaredGLMM(zig_gibbon_env, adj = TRUE)
# Compute significance of model variables
car::Anova(zig_gibbon_env,test.statistic="Chisq")

#####
# 1.1. Plot model predictions
## 1.1.1. Hour
# Compute model predictions for recording hour
pred_hour <- ggpredict(zig_gibbon_env, terms = "Hour")
# Plot predictions
ggplot(pred_hour, aes(x = x, y = predicted)) +
  # Slightly scatter points for presentation
  geom_jitter(data = hainan_all_sub_all,aes(x = as.factor(Hour), y = n),
              alpha = 0.2,width = 0.25)+
  # Zoom in on main density of data
  coord_cartesian(ylim = c(0, 2)) +
  # Plot mean model prediction
  geom_point(data = pred_hour,aes(x = as.factor(x), y = predicted),size = 3) +
  # Plot model prediction confidence intervals
  geom_errorbar(data = pred_hour,
                aes(x = as.factor(x), ymin = conf.low, ymax = conf.high),
                width = 0.15,size = 0.9) +
  # Classic plotting theme and font size of axis labels
  theme_classic() + theme(axis.text=element_text(size=15))+
  # Axis titles
  labs(x = "Hour", y = "Hourly gibbon activity")

## 1.1.2. Month
# Compute model predictions for recording month
# (check 1.1.1 annotations for functions of each row)
pred_month <- ggpredict(zig_gibbon_env, terms = "Month")
ggplot(pred_month, aes(x = x, y = predicted)) +
  geom_jitter(data = hainan_all_sub_all,aes(x = as.factor(Month), y = n),
              alpha = 0.2,width = 0.25)+
  coord_cartesian(ylim = c(0, 2)) + 
  geom_point(data = pred_month,aes(x = as.factor(x), y = predicted),size = 3) +
  geom_errorbar(data = pred_month,
                aes(x = as.factor(x), ymin = conf.low, ymax = conf.high),
                width = 0.15,size = 0.9) +
  theme_classic() +theme(axis.text=element_text(size=15))+
  labs(x = "Month", y = "Hourly gibbon activity")

## 1.1.3. Slope angle
# Compute model predictions for slope angle
# (check 1.1.1 annotations for functions of each row)
pred_slope <- ggpredict(zig_gibbon_env, terms = "Slope")
ggplot() +
  geom_jitter(data = hainan_all_sub_all,aes(x = as.factor(Slope), y = n),
              alpha = 0.2,width = 0.25)+
  coord_cartesian(ylim = c(0, 2)) + 
  geom_point(data = pred_slope,aes(x = x, y = predicted),size = 3) +
  geom_errorbar(data = pred_slope,
                aes(x = x, ymin = conf.low, ymax = conf.high),
                width = 0.15,size = 0.9) +
  theme_classic() +theme(axis.text=element_text(size=15))+  
  # Convert numerical values (1,2,3,4) into angle categories
  scale_x_discrete(labels = c("0-5°", "5-15°", "15-25°","25-35°"), 
                   breaks = c(1, 2, 3,4))+
  labs(x = "Slope angle category", y = "Hourly gibbon activity")

## 1.1.4. Forest habitat type
# Compute model predictions for forest type
# (check 1.1.1 annotations for functions of each row)
pred_habitat <- ggpredict(zig_gibbon_env, terms = "Habitat")
ggplot(pred_habitat, aes(x = x, y = predicted)) +
  geom_jitter(data = hainan_all_sub_all,aes(x = as.factor(Habitat), y = n),
              alpha = 0.2,width = 0.25)+
  coord_cartesian(ylim = c(0, 2)) +
  geom_point(data = pred_habitat,aes(x = as.factor(x), y = predicted),size = 3) +
  geom_errorbar(data = pred_habitat,
                aes(x = as.factor(x), ymin = conf.low, ymax = conf.high),
                width = 0.15,size = 0.9) +
  theme_classic() + theme(axis.text=element_text(size=15))+
  labs(x = "Type of forest habitat", y = "Hourly gibbon activity")

## 1.1.5. Canopy cover
# Compute model predictions for canopy cover
pred_cover <- ggpredict(zig_gibbon_env, terms = "Cover")
ggplot() +
  # Slightly scatter points for presentation
  geom_jitter(data = hainan_all_sub_all,aes(x = Cover, y = n),alpha = 0.2,
              width=0.01)+
  # Plot confidence band for model prediction
  geom_ribbon(data = pred_cover,
              aes(x = x, ymin = conf.low, ymax = conf.high),alpha = 0.20) +
  # Plot line for mean model prediction 
  geom_line(data = pred_cover,aes(x = x, y = predicted),linewidth = 1) +
  # Zoom in on main density of data
  coord_cartesian(ylim = c(0, 2)) +
  # Classic plotting theme and font size of axis labels
  theme_classic() +theme(axis.text = element_text(size = 15))+
  # Axis titles
  labs(x = "Canopy cover rate", y = "Hourly gibbon activity")

## 1.1.6. Month
# Compute model predictions for recording hour
# (check 1.1.5 annotations for functions of each row)
pred_alt <- ggpredict(zig_gibbon_env, terms = "Alt")
ggplot() +
  geom_jitter(data = hainan_all_sub_all,aes(x = Alt, y = n),alpha = 0.2,
              width=3)+
  geom_ribbon(data = pred_alt,
              aes(x = x, ymin = conf.low, ymax = conf.high),alpha = 0.20) +
  geom_line(data = pred_alt,aes(x = x, y = predicted),linewidth = 1) +
  coord_cartesian(ylim = c(0, 2)) +
  theme_classic() +theme(axis.text = element_text(size = 15))+
  labs(x = "Altitude (m)", y = "Hourly gibbon activity")

#####
# 2. Spatial+temporal patterns of hourly bird-community activity
zig_bird_env_all = glmmTMB(Total~Habitat+log(Alt)+log(Cover)+as.factor(Hour)+
                             as.factor(Month)+as.factor(Slope)+
                             (1|Date)+(1|Group),ziformula=~1, family = poisson, 
                           data = hainan_all_hour_diet)
# Calculate model AIC
AIC(zig_bird_env_all)
# Plot VIF graph of model variables
plot(check_collinearity(zig_bird_env_all,component = c('all')))
# Test for model zero-inflation
testZeroInflation(zig_bird_env_all)
# Model summary
summary(zig_bird_env_all)
# Calculate model R-squared value
r.squaredGLMM(zig_bird_env_all, adj = TRUE)
# Compute significance of model variables
car::Anova(zig_bird_env_all,test.statistic="Chisq")

#####
# 2.1. Plot model predictions
## 2.1.1. Hour
# Compute model predictions for recording hour
# (check 1.1.1 annotations for functions of each row)
pred_hour <- ggpredict(zig_bird_env_all, terms = "Hour")
ggplot(pred_hour, aes(x = x, y = predicted)) +
  geom_jitter(data = hainan_all_sub_all,aes(x = as.factor(Hour), y = Total),
              alpha = 0.2,width = 0.25)+
  coord_cartesian(ylim = c(0, 21)) +  
  geom_point(data = pred_hour,aes(x = as.factor(x), y = predicted),size = 3) +
  geom_errorbar(data = pred_hour,
                aes(x = as.factor(x), ymin = conf.low, ymax = conf.high),
                width = 0.15,size = 0.9) +
  theme_classic() + theme(axis.text=element_text(size=15))+
  labs(x = "Hour", y = "Hourly bird activity")

## 2.1.2. Month
# Compute model predictions for recording month
# (check 1.1.1 annotations for functions of each row)
pred_month <- ggpredict(zig_bird_env_all, terms = "Month")
ggplot(pred_month, aes(x = x, y = predicted)) +
  geom_jitter(data = hainan_all_sub_all,aes(x = as.factor(Month), y = Total),
              alpha = 0.2,width = 0.25)+
  coord_cartesian(ylim = c(0, 21)) + 
  geom_point(data = pred_month,aes(x = as.factor(x), y = predicted),size = 3) +
  geom_errorbar(data = pred_month,
                aes(x = as.factor(x), ymin = conf.low, ymax = conf.high),
                width = 0.15,size = 0.9) +
  theme_classic() +theme(axis.text=element_text(size=15))+
  labs(x = "Month", y = "Hourly bird activity")

## 2.1.3. Slope angle
# Compute model predictions for slope angle
# (check 1.1.1 annotations for functions of each row)
pred_slope <- ggpredict(zig_bird_env_all, terms = "Slope")
ggplot() +
  geom_jitter(data = hainan_all_sub_all,aes(x = as.factor(Slope), y = Total),
              alpha = 0.2,width = 0.25)+
  coord_cartesian(ylim = c(0, 21)) + 
  geom_point(data = pred_slope,aes(x = x, y = predicted),size = 3) +
  geom_errorbar(data = pred_slope,
                aes(x = x, ymin = conf.low, ymax = conf.high),
                width = 0.15,size = 0.9) +
  theme_classic() +
  theme(axis.text=element_text(size=15))+
  # Convert numerical values (1,2,3,4) into angle categories
  scale_x_discrete(labels = c("0-5°", "5-15°", "15-25°","25-35°"), 
                   breaks = c(1, 2, 3,4))+  
  labs(x = "Slope angle category", y = "Hourly bird activity")

## 2.1.4. Forest habitat type
# Compute model predictions for forest type
# (check 1.1.1 annotations for functions of each row)
pred_habitat <- ggpredict(zig_bird_env_all, terms = "Habitat")
ggplot(pred_habitat, aes(x = x, y = predicted)) +
  geom_jitter(data = hainan_all_sub_all,aes(x = as.factor(Habitat), y = Total),
              alpha = 0.2,width = 0.25)+
  coord_cartesian(ylim = c(0, 21)) +
  geom_point(data = pred_habitat,aes(x = as.factor(x), y = predicted),size = 3) +
  geom_errorbar(data = pred_habitat,
                aes(x = as.factor(x), ymin = conf.low, ymax = conf.high),
                width = 0.15,size = 0.9) +
  theme_classic() +theme(axis.text=element_text(size=15))+
  labs(x = "Type of forest habitat", y = "Hourly bird activity")

## 2.1.5. Canopy cover
# Compute model predictions for canopy cover
# (check 1.1.5 annotations for functions of each row)
pred_cover <- ggpredict(zig_bird_env_all, terms = "Cover")
ggplot() +
  geom_jitter(data = hainan_all_sub_all,aes(x = Cover, y = Total),alpha = 0.2,
              width=0.01) +
  geom_ribbon(data = pred_cover,aes(x = x, ymin = conf.low, ymax = conf.high),
              alpha = 0.20) +
  geom_line(data = pred_cover,aes(x = x, y = predicted),linewidth = 1) +
  coord_cartesian(ylim = c(0, 21)) +
  theme_classic() +theme(axis.text = element_text(size = 15))+
  labs(x = "Canopy cover rate", y = "Hourly bird activity")

## 2.1.6. Altitude
# Compute model predictions for altitude
# (check 1.1.5 annotations for functions of each row)
pred_alt <- ggpredict(zig_bird_env_all, terms = "Alt")
ggplot() +
  geom_jitter(data = hainan_all_sub_all,aes(x = Alt, y = Total),alpha = 0.2,
              width=3)+
  geom_ribbon(data = pred_alt,aes(x = x, ymin = conf.low, ymax = conf.high),
              alpha = 0.20) +
  geom_line(data = pred_alt,aes(x = x, y = predicted),linewidth = 1) +
  coord_cartesian(ylim = c(0, 21)) +
  theme_classic() +theme(axis.text = element_text(size = 15))+
  labs(x = "Altitude (m)", y = "Hourly bird activity") 


#####
# 3. Spatial+temporal patterns of hourly raptor activity
zig_bird_env_rapt = glmmTMB(Rapt~Habitat+log(Alt)+as.factor(Hour)+
                              as.factor(Month)+log(Cover)+as.factor(Slope)
                            +(1|Date)+(1|Group),ziformula=~1, family = poisson, 
                            data = hainan_all_hour_diet)
# Calculate model AIC
AIC(zig_bird_env_rapt)
# Plot VIF graph of model variables
plot(check_collinearity(zig_bird_env_rapt,component = c('all')))
# Test for model zero-inflation
testZeroInflation(zig_bird_env_rapt)
# Model summary
summary(zig_bird_env_rapt)
# Calculate model R-squared value
r.squaredGLMM(zig_bird_env_rapt, adj = TRUE)
# Compute significance of model variables
car::Anova(zig_bird_env_rapt,test.statistic="Chisq")

#####
# 3.1. Plot model predictions
## 3.1.1. Hour
# Compute model predictions for recording hour
# (check 1.1.1 annotations for functions of each row)
pred_hour <- ggpredict(zig_bird_env_rapt, terms = "Hour")
ggplot(pred_hour, aes(x = x, y = predicted)) +
  geom_jitter(data = hainan_all_sub_all,aes(x = as.factor(Hour), y = Rapt),
              alpha = 0.2,width = 0.25)+
  coord_cartesian(ylim = c(0, 20)) + 
  geom_point(data = pred_hour,aes(x = as.factor(x), y = predicted),size = 3) +
  geom_errorbar(data = pred_hour,
                aes(x = as.factor(x), ymin = conf.low, ymax = conf.high),
                width = 0.15,size = 0.9) +
  theme_classic() + theme(axis.text=element_text(size=15))+
  labs(x = "Hour", y = "Hourly raptor activity")

## 3.1.2. Month
# Compute model predictions for recording month
# (check 1.1.1 annotations for functions of each row)
pred_month <- ggpredict(zig_bird_env_rapt, terms = "Month")
ggplot(pred_month, aes(x = x, y = predicted)) +
  geom_jitter(data = hainan_all_sub_all,aes(x = as.factor(Month), y = Rapt),
              alpha = 0.2,width = 0.25)+
  coord_cartesian(ylim = c(0, 30)) + 
  geom_point(data = pred_month,aes(x = as.factor(x), y = predicted),size = 3) +
  geom_errorbar(data = pred_month,
                aes(x = as.factor(x), ymin = conf.low, ymax = conf.high),
                width = 0.15,size = 0.9) +
  theme_classic() + theme(axis.text=element_text(size=15))+
  labs(x = "Month", y = "Hourly raptor activity")

## 3.1.3. Slope
# Compute model predictions for slope angle
# (check 1.1.1 annotations for functions of each row)
pred_slope <- ggpredict(zig_bird_env_rapt, terms = "Slope")
ggplot() +
  geom_jitter(data = hainan_all_sub_all,aes(x = as.factor(Slope), y = Rapt),
              alpha = 0.2,width = 0.25)+
  coord_cartesian(ylim = c(0, 20)) +
  geom_point(data = pred_slope,aes(x = x, y = predicted),size = 3) +
  geom_errorbar(data = pred_slope,
                aes(x = x, ymin = conf.low, ymax = conf.high),
                width = 0.15,size = 0.9) +
  theme_classic() + theme(axis.text=element_text(size=15))+  
  scale_x_discrete(labels = c("0-5°", "5-15°", "15-25°","25-35°"), breaks = c(1, 2, 3,4))+
  labs(x = "Slope angle category", y = "Hourly raptor activity")

## 3.1.4. Forest habitat type
# Compute model predictions for forest type
# (check 1.1.1 annotations for functions of each row)
pred_habitat <- ggpredict(zig_bird_env_rapt, terms = "Habitat")
ggplot(pred_habitat, aes(x = x, y = predicted)) +
  geom_jitter(data = hainan_all_sub_all,aes(x = as.factor(Habitat), y = Rapt),
              alpha = 0.2,width = 0.25)+
  coord_cartesian(ylim = c(0, 19)) +  
  geom_point(data = pred_habitat,aes(x = as.factor(x), y = predicted),size = 3) +
  geom_errorbar(data = pred_habitat,
                aes(x = as.factor(x), ymin = conf.low, ymax = conf.high),
                width = 0.15,size = 0.9) +
  theme_classic() + theme(axis.text=element_text(size=15))+
  labs(x = "Type of forest habitat", y = "Hourly raptor activity")

## 3.1.5. Canopy cover
# Compute model predictions for canopy cover
# (check 1.1.5 annotations for functions of each row)
pred_cover <- ggpredict(zig_bird_env_rapt, terms = "Cover")
ggplot() +
  geom_jitter(data = hainan_all_sub_all,aes(x = Cover, y = Rapt),alpha = 0.2,
              width=0.01)+
  geom_ribbon(data = pred_cover,
              aes(x = x, ymin = conf.low, ymax = conf.high),alpha = 0.20) +
  geom_line(data = pred_cover,aes(x = x, y = predicted),linewidth = 1) +
  coord_cartesian(ylim = c(0, 70)) +
  theme_classic() +theme(axis.text = element_text(size = 15))+
  labs(x = "Canopy cover rate", y = "Hourly raptor activity") 

## 3.1.6. Altitude
# Compute model predictions for altitude
# (check 1.1.5 annotations for functions of each row)
pred_alt <- ggpredict(zig_bird_env_rapt, terms = "Alt")
ggplot() +
  geom_jitter(data = hainan_all_sub_all,aes(x = Alt, y = Rapt),alpha = 0.2,
              width=3)+
  geom_ribbon(data = pred_alt,
              aes(x = x, ymin = conf.low, ymax = conf.high),alpha = 0.20) +
  geom_line(data = pred_alt,aes(x = x, y = predicted),linewidth = 1) +
  coord_cartesian(ylim = c(0, 20)) +
  theme_classic() +theme(axis.text = element_text(size = 15))+
  labs(x = "Altitude (m)", y = "Hourly raptor activity") 

#####
# 4. Spatial+temporal patterns of hourly insectivore activity
zig_bird_env_inse = glmmTMB(Inse~Habitat+log(Alt)+as.factor(Hour)+
                              as.factor(Month)+(Cover)+as.factor(Slope)+
                              (1|Date)+(1|Group),ziformula=~1, family = poisson, 
                            data = hainan_all_hour_diet)
# Calculate model AIC
AIC(zig_bird_env_inse)
# Plot VIF graph of model variables
plot(check_collinearity(zig_bird_env_inse,component = c('all')))
# Test for model zero-inflation
testZeroInflation(zig_bird_env_inse)
# Model summary
summary(zig_bird_env_inse)
# Calculate model R-squared value
r.squaredGLMM(zig_bird_env_inse, adj = TRUE)
# Compute significance of model variables
car::Anova(zig_bird_env_inse,test.statistic="Chisq")

#####
# 4.1. Plot model predictions
## 4.1.1. Hour
# Compute model predictions for recording hour
# (check 1.1.1 annotations for functions of each row)
pred_hour <- ggpredict(zig_bird_env_inse, terms = "Hour")
ggplot(pred_hour, aes(x = x, y = predicted)) +
  geom_jitter(data = hainan_all_sub_all,aes(x = as.factor(Hour), y = Inse),
              alpha = 0.2,width = 0.25)+
  coord_cartesian(ylim = c(0, 11)) +  
  geom_point(data = pred_hour,aes(x = as.factor(x), y = predicted),size = 3) +
  geom_errorbar(data = pred_hour,
                aes(x = as.factor(x), ymin = conf.low, ymax = conf.high),
                width = 0.15,size = 0.9) +
  theme_classic() + theme(axis.text=element_text(size=15))+
  labs(x = "Hour", y = "Hourly insectivore activity")

## 4.1.2. Month
# Compute model predictions for recording month
# (check 1.1.1 annotations for functions of each row)
pred_month <- ggpredict(zig_bird_env_inse, terms = "Month")
ggplot(pred_month, aes(x = x, y = predicted)) +
  geom_jitter(data = hainan_all_sub_all,aes(x = as.factor(Month), y = Inse),
              alpha = 0.2,width = 0.25)+
  coord_cartesian(ylim = c(0, 11)) +  
  geom_point(data = pred_month,aes(x = as.factor(x), y = predicted),size = 3) +
  geom_errorbar(data = pred_month,
                aes(x = as.factor(x), ymin = conf.low, ymax = conf.high),
                width = 0.15,size = 0.9) +
  theme_classic() + theme(axis.text=element_text(size=15))+
  labs(x = "Month", y = "Hourly insectivore activity")

## 4.1.3. Slope
# Compute model predictions for slope angle
# (check 1.1.1 annotations for functions of each row)
pred_slope <- ggpredict(zig_bird_env_inse, terms = "Slope")
ggplot() +
  geom_jitter(data = hainan_all_sub_all,aes(x = as.factor(Slope), y = Inse),
              alpha = 0.2,width = 0.25)+
  coord_cartesian(ylim = c(0, 11)) + 
  geom_point(data = pred_slope,aes(x = x, y = predicted),size = 3) +
  geom_errorbar(data = pred_slope,
                aes(x = x, ymin = conf.low, ymax = conf.high),
                width = 0.15,size = 0.9) +
  theme_classic() + theme(axis.text=element_text(size=15))+  
  scale_x_discrete(labels = c("0-5°", "5-15°", "15-25°","25-35°"), 
                   breaks = c(1, 2, 3,4))+
  labs(x = "Slope angle category", y = "Hourly insectivore activity")

## 4.1.4. Forest habitat type
# Compute model predictions for forest type
# (check 1.1.1 annotations for functions of each row)
pred_habitat <- ggpredict(zig_bird_env_inse, terms = "Habitat")
ggplot(pred_habitat, aes(x = x, y = predicted)) +
  geom_jitter(data = hainan_all_sub_all,aes(x = as.factor(Habitat), y = Inse),
              alpha = 0.2,width = 0.25)+
  coord_cartesian(ylim = c(0, 9)) +   
  geom_point(data = pred_habitat,aes(x = as.factor(x), y = predicted),size = 3) +
  geom_errorbar(data = pred_habitat,
                aes(x = as.factor(x), ymin = conf.low, ymax = conf.high),
                width = 0.15,size = 0.9) +
  theme_classic() + theme(axis.text=element_text(size=15))+
  labs(x = "Type of forest habitat", y = "Hourly insectivore activity")

## 4.1.5. Canopy cover
# Compute model predictions for canopy cover
# (check 1.1.5 annotations for functions of each row)
pred_cover <- ggpredict(zig_bird_env_inse, terms = "Cover")
ggplot() +
  geom_jitter(data = hainan_all_sub_all,aes(x = Cover, y = Inse),alpha = 0.2,
              width=0.01)+
  geom_ribbon(data = pred_cover,aes(x = x, ymin = conf.low, ymax = conf.high),alpha = 0.20) +
  geom_line(data = pred_cover,aes(x = x, y = predicted),linewidth = 1) +
  coord_cartesian(ylim = c(0, 15)) +
  theme_classic() + theme(axis.text = element_text(size = 15))+
  labs(x = "Canopy cover rate", y = "Hourly insectivore activity") 

## 4.1.6. Altitude
# Compute model predictions for altitude
# (check 1.1.5 annotations for functions of each row)
pred_alt <- ggpredict(zig_bird_env_inse, terms = "Alt")
ggplot() +
  geom_jitter(data = hainan_all_sub_all,aes(x = Alt, y = Inse),alpha = 0.2,
              width=3)+
  geom_ribbon(data = pred_alt,
              aes(x = x, ymin = conf.low, ymax = conf.high),alpha = 0.20) +
  geom_line(data = pred_alt,aes(x = x, y = predicted),linewidth = 1) +
  coord_cartesian(ylim = c(0, 13)) +
  theme_classic() +theme(axis.text = element_text(size = 15))+
  labs(x = "Altitude (m)", y = "Hourly insectivore activity")

#####
# 5. Spatial+temporal patterns of hourly omnivore activity
zig_bird_env_omni = glmmTMB(Omni~Habitat+log(Alt)+as.factor(Hour)+
                              as.factor(Month)+(Cover)+as.factor(Slope)+
                              (1|Date)+(1|Group),ziformula=~1, family = poisson, 
                            data = hainan_all_hour_diet)
# Calculate model AIC
AIC(zig_bird_env_omni)
# Plot VIF graph of model variables
plot(check_collinearity(zig_bird_env_omni,component = c('all')))
# Test for model zero-inflation
testZeroInflation(zig_bird_env_omni)
# Model summary
summary(zig_bird_env_omni)
# Calculate model R-squared value
r.squaredGLMM(zig_bird_env_omni, adj = TRUE)
# Compute significance of model variables
car::Anova(zig_bird_env_omni,test.statistic="Chisq")

#####
# 5.1. Plot model predictions
## 5.1.1. Hour
# Compute model predictions for recording hour
# (check 1.1.1 annotations for functions of each row)
pred_hour <- ggpredict(zig_bird_env_omni, terms = "Hour")
ggplot(pred_hour, aes(x = x, y = predicted)) +
  geom_jitter(data = hainan_all_sub_all,aes(x = as.factor(Hour), y = Omni),
              alpha = 0.2,width = 0.25)+
  coord_cartesian(ylim = c(0, 21)) +   
  geom_point(data = pred_hour,aes(x = as.factor(x), y = predicted),size = 3) +
  geom_errorbar(data = pred_hour,
                aes(x = as.factor(x), ymin = conf.low, ymax = conf.high),
                width = 0.15,size = 0.9) +
  theme_classic() +
  theme(axis.text=element_text(size=15))+
  labs(x = "Hour", y = "Hourly omnivore activity")

## 5.1.2. Month
# Compute model predictions for recording month
# (check 1.1.1 annotations for functions of each row)
pred_month <- ggpredict(zig_bird_env_omni, terms = "Month")
ggplot(pred_month, aes(x = x, y = predicted)) +
  geom_jitter(data = hainan_all_sub_all,aes(x = as.factor(Month), y = Omni),
              alpha = 0.2,width = 0.25)+
  coord_cartesian(ylim = c(0, 21)) +   
  geom_point(data = pred_month,aes(x = as.factor(x), y = predicted),size = 3) +
  geom_errorbar(data = pred_month,
                aes(x = as.factor(x), ymin = conf.low, ymax = conf.high),
                width = 0.15,size = 0.9) +
  theme_classic() +theme(axis.text=element_text(size=15))+
  labs(x = "Month", y = "Hourly omnivore activity")

## 5.1.3. Slope angle
# Compute model predictions for slope angle
# (check 1.1.1 annotations for functions of each row)
pred_slope <- ggpredict(zig_bird_env_omni, terms = "Slope")
ggplot() +
  geom_jitter(data = hainan_all_sub_all,aes(x = as.factor(Slope), y = Omni),
              alpha = 0.2,width = 0.25)+
  coord_cartesian(ylim = c(0, 21)) +   
  geom_point(data = pred_slope,aes(x = x, y = predicted),size = 3) +
  geom_errorbar(data = pred_slope,
                aes(x = x, ymin = conf.low, ymax = conf.high),
                width = 0.15,size = 0.9) +
  theme_classic() + theme(axis.text=element_text(size=15))+  
  scale_x_discrete(labels = c("0-5°", "5-15°", "15-25°","25-35°"), 
                   breaks = c(1, 2, 3,4))+
  labs(x = "Slope angle category", y = "Hourly omnivore activity")

## 5.1.4. Forest habitat type
# Compute model predictions for forest type
# (check 1.1.1 annotations for functions of each row)
pred_habitat <- ggpredict(zig_bird_env_omni, terms = "Habitat")
ggplot(pred_habitat, aes(x = x, y = predicted)) +
  geom_jitter(data = hainan_all_sub_all,aes(x = as.factor(Habitat), y = Omni),
              alpha = 0.2,width = 0.25)+
  coord_cartesian(ylim = c(0, 21)) +
  geom_point(data = pred_habitat,aes(x = as.factor(x), y = predicted),size = 3) +
  geom_errorbar(data = pred_habitat,
                aes(x = as.factor(x), ymin = conf.low, ymax = conf.high),
                width = 0.15,size = 0.9) +
  theme_classic() +theme(axis.text=element_text(size=15))+
  labs(x = "Type of forest habitat", y = "Hourly omnivore activity")

## 5.1.5. Canopy cover
# Compute model predictions for canopy cover
# (check 1.1.5 annotations for functions of each row)
pred_cover <- ggpredict(zig_bird_env_omni, terms = "Cover")
ggplot() +
  geom_jitter(data = hainan_all_sub_all,aes(x = Cover, y = Omni),alpha = 0.2,
              width=0.01)+
  geom_ribbon(data = pred_cover,
              aes(x = x, ymin = conf.low, ymax = conf.high),alpha = 0.20) +
  geom_line(data = pred_cover,aes(x = x, y = predicted),linewidth = 1) +
  coord_cartesian(ylim = c(0, 21)) +
  theme_classic() +theme(axis.text = element_text(size = 15))+
  labs(x = "Canopy cover rate", y = "Hourly omnivore activity") 

## 5.1.6. Altitude
# Compute model predictions for altitude
# (check 1.1.5 annotations for functions of each row)
pred_alt <- ggpredict(zig_bird_env_omni, terms = "Alt")
ggplot() +
  geom_jitter(data = hainan_all_sub_all,aes(x = Alt, y = Omni),alpha = 0.2,
              width=3)+
  geom_ribbon(data = pred_alt,
              aes(x = x, ymin = conf.low, ymax = conf.high),alpha = 0.20) +
  geom_line(data = pred_alt,aes(x = x, y = predicted),linewidth = 1) +
  coord_cartesian(ylim = c(0, 11)) +
  theme_classic() +theme(axis.text = element_text(size = 15))+
  labs(x = "Altitude (m)", y = "Hourly omnivore activity")
