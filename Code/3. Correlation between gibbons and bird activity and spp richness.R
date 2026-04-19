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

#####
# 1. Load CSV dataset created in last section
hainan_all_hour_diet <- read.csv("Data/hainan_all_hour_diet.csv")

#####
# 2. Prepare site-level average and peak activity of taxa/guilds
# Calculate daily total, mean and maximum activity of each site of each taxa/guild
hainan_bird_day <- hainan_all_hour_diet%>%group_by(Site,Date)%>%
  # Only keep 6am-6pm and June-September data
  filter(Month!="202110" & (Hour >5 & Hour <18))%>%
  # Total activity per day per site
  summarise(RaptT = sum(Rapt),InseT=sum(Inse),OmniT=sum(Omni),FlorT=sum(Flor),
            TotalT=sum(Total),GibbonT=sum(n),count=n(),
            # Mean hourly activity per site
            RaptD = mean(Rapt),InseD=mean(Inse),OmniD=mean(Omni),FlorD=mean(Flor),
            TotalD=mean(Total),GibbonD=mean(n),
            # Daily peak activity per site
            RaptM = max(Rapt),InseM=max(Inse),OmniM=max(Omni),FlorM=max(Flor),
            TotalM=max(Total),GibbonM=max(n),
            # Include environmental information per site
            Alt=mean(Alt),Slope=mean(Slope),Habitat=unique(Habitat),
            Cover=mean(Cover),Month=unique(Month))

# Calculate site-level total and average activity of each taxa/guild
hainan_bird_site_bysite <- hainan_all_hour_diet%>%
  # Only keep 6am-6pm and June-September data, group by site
  filter(Month!="202110" & (Hour >5 & Hour <18))%>%group_by(Site)%>%
  # Total activity per site
  summarise(RaptT = sum(Rapt),InseT=sum(Inse),OmniT=sum(Omni),FlorT=sum(Flor),
            TotalT=sum(Total),GibbonT=sum(n),count=n(),
            # Average hourly activity per site
            RaptD = mean(Rapt),InseD=mean(Inse),OmniD=mean(Omni),FlorD=mean(Flor),
            TotalD=mean(Total),GibbonD=mean(n),
            # Include environmental information per site
            Alt=mean(Alt),Slope=mean(Slope),Habitat=unique(Habitat),Cover=mean(Cover))

# Calculate site-level peak activity of each taxa/guild
hainan_bird_site_byday <- hainan_bird_day%>%
  filter(Month!="202110")%>%group_by(Site)%>%
  # Maximum daily peak activity
  summarise(RaptMaxMax = max(RaptM),InseMaxMax=max(InseM),OmniMaxMax=max(OmniM),
            FlorMaxMax=max(FlorM),TotalMaxMax=max(TotalM),GibbonMaxMax=max(GibbonM),
            # Average daily peak activity
            RaptAMax = mean(RaptM),InseAMax=mean(InseM),OmniAMax=mean(OmniM),
            FlorAMax=mean(FlorM),TotalAMax=mean(TotalM),GibbonAMax=mean(GibbonM),
            # Maximum daily average activity
            RaptMaxA = max(RaptD),InseMaxA=max(InseD),OmniMaxA=max(OmniD),
            FlorMaxA=max(FlorD),TotalMaxA=max(TotalD),GibbonMaxA=max(GibbonD))

# Merge two datasets
hainan_bird_site <- merge(hainan_bird_site_bysite,hainan_bird_site_byday,by="Site")

#####
# 3. Calculate total and foraging-guild-specific species richness of each site
# Calculate total and guild-specific species richness
hainan_bird_diversity <- hainan_bird_hour_select%>%
  group_by(Site,DietClass,Spp1) %>% summarise(count=n())
# Extract raptor species richness
hainan_diversity_rapt <- hainan_bird_diversity %>% filter(DietClass=="猛禽") %>% 
  group_by(Site) %>% summarise(RaptCount=n())
# Extract insectivorous bird species richness
hainan_diversity_inse <- hainan_bird_diversity %>% filter(DietClass=="食虫") %>% 
  group_by(Site) %>% summarise(InseCount=n())
# Extract omnivorous bird species richness
hainan_diversity_omni <- hainan_bird_diversity %>% filter(DietClass=="杂食性") %>% 
  group_by(Site) %>% summarise(OmniCount=n())
# Extract total bird species richness
hainan_diversity_allbird <- hainan_bird_hour_select %>% group_by(Site,Spp1)%>%
  summarise(count=n()) %>% group_by(Site)%>%summarise(SppCount=n())

# Merge all species richness measures into one dataset
hainan_diversity_total <- merge(hainan_diversity_allbird,hainan_diversity_rapt,by="Site")
hainan_diversity_total <- merge(hainan_diversity_total,hainan_diversity_inse,by="Site")
hainan_diversity_total <- merge(hainan_diversity_total,hainan_diversity_omni,by="Site")
# Compute group of sample point clusters
hainan_diversity_total$Group <- substr(hainan_diversity_total$Site,start=1,stop=1)
# Merge bird activity and species richness dataset
hainan_bird_site <- merge(hainan_bird_site, hainan_diversity_total,by="Site")

#####
# 4. Correlation between gibbon activity and bird total activity and species richness
## 4.1. Model specification and stats
# Correlation between average daily peak activity of gibbons and birds
zig_bg_all_Amax =  lm(log(TotalAMax)~log(GibbonAMax), data = hainan_bird_site)
# Model summary
summary(zig_bg_all_Amax)

# Correlation between average hourly activity of gibbons and birds
zig_bg_all_mean =  lm(log(TotalD)~log(GibbonD), data = hainan_bird_site)
# Model summary
summary(zig_bg_all_mean)

# Correlation between average gibbon daily peak activity and bird spp. richness
zig_bg_all_SppA =  glm(SppCount~log(GibbonAMax), family="quasipoisson",data = hainan_bird_site)
# Test for model dispersion (for poisson distribution)
# (underdispersion found, so switched to quasipoisson distribution)
#simRes <- simulateResiduals(zig_bg_all_SppA)
#testDispersion(simRes)
# Model summary
summary(zig_bg_all_SppA)
# Check for significance of fixed effect
car::Anova(zig_bg_all_SppA,test.statistic="F")
# Calculate model R-squared value
r2(zig_bg_all_SppA)

# Correlation between average gibbon hourly activity and bird spp. richness
# (check annotations above for function of each row)
zig_bg_all_SppM =  glm(SppCount~log(GibbonD), family="quasipoisson", data = hainan_bird_site)
#simRes <- simulateResiduals(zig_bg_all_SppM)
#testDispersion(simRes)
summary(zig_bg_all_SppM)
car::Anova(zig_bg_all_SppM,test.statistic="F")
r2(zig_bg_all_SppM)

#####
## 4.2. Plot model predictions
# 4.2.1. Average of daily max hourly activity in gibbons vs birds
# Create a sequence of predictor values
new_data <- data.frame(GibbonAMax = seq(min(hainan_bird_site$GibbonAMax),
                                        max(hainan_bird_site$GibbonAMax),
                                        length.out = 100))
# Get log-scale predictions + confidence intervals
pred_log_ci <- predict(zig_bg_all_Amax,newdata = new_data,interval = "confidence")
# Convert to a data frame and bind to new_data
pred_log_ci <- as.data.frame(pred_log_ci)
new_data <- cbind(new_data, pred_log_ci)
# Bias-corrected back-transformation
sigma2 <- summary(zig_bg_all_Amax)$sigma^2
# Final log-transformed predictions for plotting
new_data <- transform(new_data,
                      fit_bc = exp(fit + sigma2 / 2),
                      lwr_bc = exp(lwr + sigma2 / 2),
                      upr_bc = exp(upr + sigma2 / 2))
# Plot model predictions
ggplot() +
  # Plot raw data points
  geom_point(data = hainan_bird_site,aes(x = GibbonAMax, y = TotalAMax),alpha = 0.6) +
  # Plot confidence band for model predictions
  geom_ribbon(data = new_data,
              aes(x = GibbonAMax, ymin = lwr_bc, ymax = upr_bc),
              fill = "black", alpha = 0.15) +
  # Plot line for mean model predictions
  geom_line(data = new_data,
            aes(x = GibbonAMax, y = fit_bc),color = "black", size = 1.2) +
  # Classic plot theme and font size of axis labels
  theme_classic()+theme(axis.text=element_text(size=15))+
  # Axis titles
  labs(x = "Mean daily peak gibbon activity",y = "Mean daily peak bird total activity") 

# 4.2.2. Average of hourly activity in gibbons vs birds
# (check 4.2.1 annotations for functions of each row)
new_data <- data.frame(GibbonD = seq(min(hainan_bird_site$GibbonD),
                                     max(hainan_bird_site$GibbonD),
                                     length.out = 100))
pred_log_ci <- predict(zig_bg_all_mean,newdata = new_data,interval = "confidence")
pred_log_ci <- as.data.frame(pred_log_ci)
new_data <- cbind(new_data, pred_log_ci)
sigma2 <- summary(zig_bg_all_mean)$sigma^2
new_data <- transform(new_data,
                      fit_bc = exp(fit + sigma2 / 2),
                      lwr_bc = exp(lwr + sigma2 / 2),
                      upr_bc = exp(upr + sigma2 / 2))
ggplot() +
  geom_point(data = hainan_bird_site,aes(x = GibbonD, y = TotalD),alpha = 0.6) +
  geom_ribbon(data = new_data,
              aes(x = GibbonD, ymin = lwr_bc, ymax = upr_bc),
              fill = "black", alpha = 0.15) +
  geom_line(data = new_data,
            aes(x = GibbonD, y = fit_bc),color = "black", size = 1.2) +
  coord_cartesian(xlim = c(0.05, 0.3))+
  theme_classic()+theme(axis.text=element_text(size=15))+
  labs(x = "Mean hourly gibbon activity",y = "Mean hourly bird total activity")

# 4.2.3. Average of daily max hourly activity in gibbons vs bird spp richness
# (check 4.2.1 annotations for functions of each row)
new_data <- data.frame(GibbonAMax = seq(min(hainan_bird_site$GibbonAMax),
                                        max(hainan_bird_site$GibbonAMax),
                                        length.out = 100))
pred <- predict(zig_bg_all_SppA,new_data,type = "link",se.fit = T)
new_data$fit <- exp(pred$fit)
new_data$lwr <- exp(pred$fit - 1.96 * pred$se.fit)
new_data$upr <- exp(pred$fit + 1.96 * pred$se.fit)
ggplot(hainan_bird_site, aes(x = GibbonAMax, y = SppCount)) +
  geom_point(alpha = 0.6) +
  geom_line(data = new_data,
            aes(x = GibbonAMax, y = fit),colour = "black",size = 1.2) +
  geom_ribbon(data = new_data,
              aes(x = GibbonAMax, ymin = lwr, ymax = upr),
              alpha = 0.15,fill = "black",inherit.aes = FALSE) +
  theme_classic() + theme(axis.text=element_text(size=15))+
  labs(x = "Mean daily peak gibbon activity",y = "Bird total species richness")  

# 4.2.4. Average of hourly activity in gibbons vs bird spp richness
# (check 4.2.1 annotations for functions of each row)
new_data <- data.frame(GibbonD = seq(min(hainan_bird_site$GibbonD),
                                     max(hainan_bird_site$GibbonD),
                                     length.out = 100))
pred <- predict(zig_bg_all_SppM,new_data,type = "link",se.fit = T)
new_data$fit <- exp(pred$fit)
new_data$lwr <- exp(pred$fit - 1.96 * pred$se.fit)
new_data$upr <- exp(pred$fit + 1.96 * pred$se.fit)
ggplot(hainan_bird_site, aes(x = GibbonD, y = SppCount)) +
  geom_point(alpha = 0.6) +
  geom_line(data = new_data,
            aes(x = GibbonD, y = fit),colour = "black",size = 1.2) +
  geom_ribbon(data = new_data,
              aes(x = GibbonD, ymin = lwr, ymax = upr),
              alpha = 0.15,fill = "black",inherit.aes = FALSE) +
  theme_classic() +theme(axis.text=element_text(size=15))+
  coord_cartesian(xlim = c(0.05, 0.3))+
  labs(x = "Mean hourly gibbon activity",y = "Bird total species richness") 

#####
# 5. Correlation between gibbon activity and raptor activity and species richness
## 5.1. Model specification and stats
# (check 4.1. annotations for functions of each row)
# Correlation between average daily peak activity of gibbons and raptors
zig_bg_rapt_AMax =  lm((RaptAMax)~log(GibbonAMax), data = hainan_bird_site)
summary(zig_bg_rapt_AMax)

# Correlation between average hourly activity of gibbons and raptors
zig_bg_rapt_mean =  lm((RaptD)~log(GibbonD), data = hainan_bird_site)
summary(zig_bg_rapt_mean)

# Correlation between average gibbon daily peak activity and raptor species richness
zig_bg_rapt_SppA =  glm(RaptCount~log(GibbonAMax), family="quasipoisson",data = hainan_bird_site)
#simRes <- simulateResiduals(zig_bg_rapt_SppA)
#testDispersion(simRes)
summary(zig_bg_rapt_SppA)
car::Anova(zig_bg_rapt_SppA,test.statistic="F")
r2(zig_bg_rapt_SppA)

# Correlation between average gibbon hourly activity and raptor species richness
zig_bg_rapt_SppM =  glm(RaptCount~log(GibbonD), family="quasipoisson",data = hainan_bird_site)
#simRes <- simulateResiduals(zig_bg_rapt_SppM)
#testDispersion(simRes)
summary(zig_bg_rapt_SppM)
car::Anova(zig_bg_rapt_SppM,test.statistic="F")
r2(zig_bg_rapt_SppM)

#####
## 5.2. Plot model predictions
# 5.2.1. Average of daily max hourly activity in gibbons vs raptors
# Create a sequence of predictor values
new_data <- data.frame(GibbonAMax = seq(min(hainan_bird_site$GibbonAMax),
                                        max(hainan_bird_site$GibbonAMax),
                                        length.out = 100))
# Predict model confidence intervals
pred <- predict(zig_bg_rapt_AMax,new_data,interval = "confidence")
# Get average and lower/upper confidence intervals of model prediction
new_data$fit <- pred[, "fit"]
new_data$lwr <- pred[, "lwr"]
new_data$upr <- pred[, "upr"]
# Plot model predictions
ggplot(hainan_bird_site, aes(x = GibbonAMax, y = RaptAMax)) +
  # Plot raw data points
  geom_point(alpha = 0.6) +
  # Plot line for mean model predictions
  geom_line(data = new_data,
            aes(x = GibbonAMax, y = fit),colour = "black",size = 1.2) +
  # Plot confidence band for model predictions
  geom_ribbon(data = new_data,
              aes(x = GibbonAMax, ymin = lwr, ymax = upr),
              alpha = 0.15,fill = "black",inherit.aes = FALSE) +
  # Classic plot theme and axis label font size
  theme_classic() + theme(axis.text=element_text(size=15))+
  # Axis titles
  labs(x = "Mean daily peak gibbon activity",y = "Mean daily peak raptor activity")

# 5.2.2. Average of hourly activity in gibbons vs raptors
# (check 5.2.1 annotations for functions of each row)
new_data <- data.frame(GibbonD = seq(min(hainan_bird_site$GibbonD),
                                     max(hainan_bird_site$GibbonD),
                                     length.out = 100))
pred <- predict(zig_bg_rapt_mean,new_data,interval = "confidence")
new_data$fit <- pred[, "fit"]
new_data$lwr <- pred[, "lwr"]
new_data$upr <- pred[, "upr"]
ggplot(hainan_bird_site, aes(x = GibbonD, y = RaptD)) +
  geom_point(alpha = 0.6) +
  geom_line(data = new_data,
            aes(x = GibbonD, y = fit),colour = "black",size = 1.2) +
  geom_ribbon(data = new_data,
              aes(x = GibbonD, ymin = lwr, ymax = upr),
              alpha = 0.15,fill = "black",inherit.aes = FALSE) +
  coord_cartesian(xlim = c(0.05, 0.3))+
  theme_classic() + theme(axis.text=element_text(size=15))+
  labs(x = "Mean hourly gibbon activity",y = "Mean hourly raptor activity")

# 5.2.3. Average of daily max hourly activity in gibbons vs raptor spp richness
# (check 5.2.1 annotations for functions of each row)
new_data <- data.frame(GibbonAMax = seq(min(hainan_bird_site$GibbonAMax),
                                        max(hainan_bird_site$GibbonAMax),
                                        length.out = 100))
pred <- predict(zig_bg_rapt_SppA,new_data,type = "link",se.fit = T)
new_data$fit <- exp(pred$fit)
new_data$lwr <- exp(pred$fit - 1.96 * pred$se.fit)
new_data$upr <- exp(pred$fit + 1.96 * pred$se.fit)
ggplot(hainan_bird_site, aes(x = GibbonAMax, y = RaptCount)) +
  geom_point(alpha = 0.6) +
  geom_line(data = new_data,
            aes(x = GibbonAMax, y = fit),colour = "black",size = 1.2) +
  geom_ribbon(data = new_data,
              aes(x = GibbonAMax, ymin = lwr, ymax = upr),
              alpha = 0.15,fill = "black",inherit.aes = FALSE) +
  theme_classic() + theme(axis.text=element_text(size=15)) +
  labs(x = "",y = "") 

# 5.2.4. Average of hourly activity in gibbons vs raptor spp richness
# (check 5.2.1 annotations for functions of each row)
new_data <- data.frame(GibbonD = seq(min(hainan_bird_site$GibbonD),
                                     max(hainan_bird_site$GibbonD),
                                     length.out = 100))
pred <- predict(zig_bg_rapt_SppM,new_data,type = "link",se.fit = T)
new_data$fit <- exp(pred$fit)
new_data$lwr <- exp(pred$fit - 1.96 * pred$se.fit)
new_data$upr <- exp(pred$fit + 1.96 * pred$se.fit)
ggplot(hainan_bird_site, aes(x = GibbonD, y = RaptCount)) +
  geom_point(alpha = 0.6) +
  geom_line(data = new_data,
            aes(x = GibbonD, y = fit),colour = "black",size = 1.2) +
  geom_ribbon(data = new_data,
              aes(x = GibbonD, ymin = lwr, ymax = upr),
              alpha = 0.15,fill = "black",inherit.aes = FALSE) +
  theme_classic() + theme(axis.text=element_text(size=15)) +
  coord_cartesian(xlim = c(0.05, 0.3))+
  labs(x = "Mean hourly gibbon activity",y = "Raptor bird species richness")

#####
# 6. Correlation between gibbon activity and omnivorous bird activity and species richness
## 6.1. Model specification and stats
# (check 4.1. annotations for functions of each row)
# Correlation between average daily peak activity of gibbons and omnivorous birds
zig_bg_omni_AMax =  lm(log(OmniAMax)~log(GibbonAMax), data = hainan_bird_site)
summary(zig_bg_omni_AMax)

# Correlation between average hourly activity of gibbons and omnivorous birds
zig_bg_omni_mean =  lm(log(OmniD)~log(GibbonD), data = hainan_bird_site)
summary(zig_bg_omni_mean)

# Correlation between average gibbon daily peak activity and omnivorous bird species richness
zig_bg_omni_SppA =  glm(OmniCount~log(GibbonAMax), family="quasipoisson",data = hainan_bird_site)
#simRes <- simulateResiduals(zig_bg_omni_SppA)
#testDispersion(simRes)
summary(zig_bg_omni_SppA)
car::Anova(zig_bg_omni_SppA,test.statistic="F")
r2(zig_bg_omni_SppA)

# Correlation between average gibbon hourly activity and omnivorous bird species richness
zig_bg_omni_SppM =  glm(OmniCount~log(GibbonD), family="quasipoisson", data = hainan_bird_site)
#simRes <- simulateResiduals(zig_bg_omni_SppM)
#testDispersion(simRes)
summary(zig_bg_omni_SppM)
car::Anova(zig_bg_omni_SppM,test.statistic="F")
r2(zig_bg_omni_SppM)

#####
## 6.2. Plot model predictions
# 6.2.1. Average of daily peak activity in gibbons vs omnivores
# Create a sequence of predictor values
new_data <- data.frame(GibbonAMax = seq(min(hainan_bird_site$GibbonAMax),
                                        max(hainan_bird_site$GibbonAMax),
                                        length.out = 100))
# Predict log-transformed confidence intervals of model predictions
pred_log_ci <- predict(zig_bg_omni_AMax,newdata = new_data,interval = "confidence")
# Transform into data frame
pred_log_ci <- as.data.frame(pred_log_ci)
# Back-transform into un-logged forms (mean predictions and lower/upper intervals)
new_data <- cbind(new_data, pred_log_ci)
sigma2 <- summary(zig_bg_omni_AMax)$sigma^2
new_data <- transform(new_data,fit_bc = exp(fit + sigma2 / 2),
                      lwr_bc = exp(lwr + sigma2 / 2),
                      upr_bc = exp(upr + sigma2 / 2))
# Plot model predictions
ggplot() +
  # Plot raw data points
  geom_point(data = hainan_bird_site,
             aes(x = GibbonAMax, y = OmniAMax),alpha = 0.6) +
  # Plot confidence band of model predictions
  geom_ribbon(data = new_data,
              aes(x = GibbonAMax, ymin = lwr_bc, ymax = upr_bc),
              fill = "black", alpha = 0.15) +
  # Plot line of mean model predictions
  geom_line(data = new_data,
            aes(x = GibbonAMax, y = fit_bc),color = "black", size = 1.2) + 
  # Classic plot theme and axis label font size
  theme_classic() + theme(axis.text=element_text(size=15))+
  labs(x = "Mean daily peak gibbon activity",y = "Mean daily peak omnivore activity") 

# 6.2.2. Average of hourly activity in gibbons vs omnivores
# (check 6.2.1. annotations for functions of each row)
new_data <- data.frame(GibbonD = seq(min(hainan_bird_site$GibbonD),
                                     max(hainan_bird_site$GibbonD),
                                     length.out = 100))
pred_log_ci <- predict(zig_bg_omni_mean,newdata = new_data,interval = "confidence")
pred_log_ci <- as.data.frame(pred_log_ci)
new_data <- cbind(new_data, pred_log_ci)
sigma2 <- summary(zig_bg_omni_mean)$sigma^2
new_data <- transform(new_data,fit_bc = exp(fit + sigma2 / 2),
                      lwr_bc = exp(lwr + sigma2 / 2),
                      upr_bc = exp(upr + sigma2 / 2))
ggplot() +
  geom_point(data = hainan_bird_site,
             aes(x = GibbonD, y = OmniD),alpha = 0.6) +
  geom_ribbon(data = new_data,
              aes(x = GibbonD, ymin = lwr_bc, ymax = upr_bc),
              fill = "black", alpha = 0.15) +
  geom_line(data = new_data,
            aes(x = GibbonD, y = fit_bc),color = "black", size = 1.2)+ 
  coord_cartesian(xlim = c(0.05, 0.3))+
  theme_classic()+theme(axis.text=element_text(size=15)) +
  labs(x = "Mean hourly gibbon activity",y = "Mean hourly omnivore activity") 

# 6.2.3. Average of daily peak activity in gibbons vs omnivore spp richness
# (check 6.2.1. annotations for functions of each row)
new_data <- data.frame(GibbonAMax = seq(min(hainan_bird_site$GibbonAMax),
                                        max(hainan_bird_site$GibbonAMax),
                                        length.out = 100))
pred <- predict(zig_bg_omni_SppA,new_data,type = "link",se.fit = T)
new_data$fit <- exp(pred$fit)
new_data$lwr <- exp(pred$fit - 1.96 * pred$se.fit)
new_data$upr <- exp(pred$fit + 1.96 * pred$se.fit)
ggplot(hainan_bird_site, aes(x = GibbonAMax, y = OmniCount)) +
  geom_point(alpha = 0.6) +
  geom_line(data = new_data,
            aes(x = GibbonAMax, y = fit),colour = "black",size = 1.2) +
  geom_ribbon(data = new_data,
              aes(x = GibbonAMax, ymin = lwr, ymax = upr),
              alpha = 0.15,fill = "black",inherit.aes = FALSE) +
  scale_y_continuous(labels = scales::label_number(accuracy = 1))+
  theme_classic() + theme(axis.text=element_text(size=15))+
  labs(x = "Mean daily peak gibbon activity",y = "Omnivorous bird species richness") 

# 6.2.4. Average of hourly activity in gibbons vs omnivore spp richness
# (check 6.2.1. annotations for functions of each row)
new_data <- data.frame(GibbonD = seq(min(hainan_bird_site$GibbonD),
                                     max(hainan_bird_site$GibbonD),
                                     length.out = 100))
pred <- predict(zig_bg_omni_SppM,new_data,type = "link",se.fit = T)
new_data$fit <- exp(pred$fit)
new_data$lwr <- exp(pred$fit - 1.96 * pred$se.fit)
new_data$upr <- exp(pred$fit + 1.96 * pred$se.fit)
ggplot(hainan_bird_site, aes(x = GibbonD, y = OmniCount)) +
  geom_point(alpha = 0.6) +
  geom_line(data = new_data,
            aes(x = GibbonD, y = fit),colour = "black",size = 1.2) +
  geom_ribbon(data = new_data,
              aes(x = GibbonD, ymin = lwr, ymax = upr),
              alpha = 0.15,fill = "black",inherit.aes = FALSE) +
  coord_cartesian(xlim = c(0.05, 0.3))+
  theme_classic() + theme(axis.text=element_text(size=15)) +
  labs(x = "Mean hourly gibbon activity",y = "Omnivorous bird species richness")

#####
# 7. Correlation between gibbon activity and insectivorous bird activity and species richness
## 7.1. Model specification and stats
# (check 4.1. annotations for functions of each row)
# Correlation between average daily peak activity of gibbons and insectivorous birds
zig_bg_inse_AMax =  lm((InseAMax)~(GibbonAMax), data = hainan_bird_site)
summary(zig_bg_inse_AMax)

# Correlation between average daily peak activity of gibbons and insectivorous birds
zig_bg_inse_mean =  lm((InseD)~(GibbonD), data = hainan_bird_site)
summary(zig_bg_inse_mean)

# Correlation between average gibbon daily peak activity and insectivorous bird species richness
zig_bg_inse_SppA =  glm(InseCount~(GibbonAMax), family="quasipoisson",data = hainan_bird_site)
#simRes <- simulateResiduals(zig_bg_inse_SppA)
#testDispersion(simRes)
summary(zig_bg_inse_SppA)
car::Anova(zig_bg_inse_SppA,test.statistic="F")
r2(zig_bg_inse_SppA)

# Correlation between average gibbon hourly activity and insectivorous bird species richness
zig_bg_inse_SppM =  glm(InseCount~(GibbonD), family="quasipoisson", data = hainan_bird_site)
#simRes <- simulateResiduals(zig_bg_inse_SppM)
#testDispersion(simRes)
summary(zig_bg_inse_SppM)
car::Anova(zig_bg_inse_SppM,test.statistic="F")
r2(zig_bg_inse_SppM)

#####
## 7.2. Plot model predictions
# 7.2.1. Average of daily peak activity in gibbons vs insectivores
# Create a sequence of predictor values
new_data <- data.frame(GibbonAMax = seq(min(hainan_bird_site$GibbonAMax),
                                        max(hainan_bird_site$GibbonAMax),
                                        length.out = 100))
# Compute model predictions
pred <- predict(zig_bg_inse_AMax,new_data,interval = "confidence")
# Extract mean and confidence intervals of model predictions
new_data$fit <- pred[, "fit"]
new_data$lwr <- pred[, "lwr"]
new_data$upr <- pred[, "upr"]
# Plot model predictions
ggplot() +
  # Plot raw data points
  geom_point(data = hainan_bird_site,
             aes(x = GibbonAMax, y = InseAMax),alpha = 0.6) +
  # Plot line of mean model predictions
  geom_line(data = new_data,
            aes(x = GibbonAMax, y = fit),color = "black", size = 1.2) +
  # Plot confidence band for model predictions
  geom_ribbon(data = new_data,
              aes(x = GibbonAMax, ymin = lwr, ymax = upr),
              fill = "black", alpha = 0.15) + 
  # Classic plot theme and font of axis labels
  theme_classic() + theme(axis.text=element_text(size=15))+
  # Axis titles
  labs(x = "Mean daily peak gibbon activity",y = "Mean daily peak insectivorous bird activity") 

# 7.2.2. Average of hourly activity in gibbons vs insectivores
# (check 7.2.1. annotations for functions of each row)
new_data <- data.frame(GibbonD = seq(min(hainan_bird_site$GibbonD),
                                     max(hainan_bird_site$GibbonD),
                                     length.out = 100))
pred <- predict(zig_bg_inse_mean,new_data,interval = "confidence")
new_data$fit <- pred[, "fit"]
new_data$lwr <- pred[, "lwr"]
new_data$upr <- pred[, "upr"]
ggplot() +
  geom_point(data = hainan_bird_site,
             aes(x = GibbonD, y = InseD),alpha = 0.6) +
  geom_line(data = new_data,
            aes(x = GibbonD, y = fit),color = "black", size = 1.2) +
  geom_ribbon(data = new_data,
              aes(x = GibbonD, ymin = lwr, ymax = upr),
              fill = "black", alpha = 0.15) +
  coord_cartesian(xlim = c(0.05, 0.3))+ 
  theme_classic()+theme(axis.text=element_text(size=15)) +
  labs(x = "Mean hourly gibbon activity",y = "Mean hourly insectivorous activity") 

# 7.2.3. Average of daily max hourly activity in gibbons vs insectivore spp richness
# (check 7.2.1. annotations for functions of each row)
new_data <- data.frame(GibbonAMax = seq(min(hainan_bird_site$GibbonAMax),
                                        max(hainan_bird_site$GibbonAMax),
                                        length.out = 100))
pred <- predict(zig_bg_inse_SppA,new_data,type = "link",se.fit = TRUE)
new_data$fit <- exp(pred$fit)
new_data$lwr <- exp(pred$fit - 1.96 * pred$se.fit)
new_data$upr <- exp(pred$fit + 1.96 * pred$se.fit)
ggplot() +
  geom_point(data = hainan_bird_site,
             aes(x = GibbonAMax, y = InseCount),alpha = 0.6) +
  geom_ribbon(data = new_data,
              aes(x = GibbonAMax, ymin = lwr, ymax = upr),
              fill = "black", alpha = 0.15) +
  geom_line(data = new_data,
            aes(x = GibbonAMax, y = fit),color = "black", size = 1.2) + 
  theme_classic() + theme(axis.text=element_text(size=15))+
  labs(x = "Mean daily peak gibbon activity",y = "Insectivorous bird species richness") 

# 7.2.4. Average of hourly activity in gibbons vs insectivore spp richness
# (check 7.2.1. annotations for functions of each row)
new_data <- data.frame(GibbonD = seq(min(hainan_bird_site$GibbonD),
                                     max(hainan_bird_site$GibbonD),
                                     length.out = 100))
pred <- predict(zig_bg_inse_SppM,new_data,type = "link",se.fit = TRUE)
new_data$fit <- exp(pred$fit)
new_data$lwr <- exp(pred$fit - 1.96 * pred$se.fit)
new_data$upr <- exp(pred$fit + 1.96 * pred$se.fit)
ggplot() +
  geom_point(data = hainan_bird_site,
             aes(x = GibbonD, y = InseCount),alpha = 0.6) +
  geom_ribbon(data = new_data,
              aes(x = GibbonD, ymin = lwr, ymax = upr),
              fill = "black", alpha = 0.15) +
  geom_line(data = new_data,
            aes(x = GibbonD, y = fit),color = "black", size = 1.2) +
  coord_cartesian(xlim = c(0.05, 0.3)) + 
  theme_classic() + theme(axis.text=element_text(size=15))+
  labs(x = "Mean hourly gibbon activity",y = "Insectivorous bird species richness")

#####
# 8. Write dataset into CSV file for next section
write.csv(hainan_bird_site,"Data/hainan_bird_site.csv",row.names=F)
