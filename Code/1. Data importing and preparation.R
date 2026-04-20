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
# 1. Load raw data
# Import soundscape indices data
indices_hainan <- read.csv("indices_hainan.csv")
# Import bird detection file
hainan_bird_raw <- read.csv("Data/bird_detection.csv")
# Import gibbon detection file
hainan_gibbons <- read.csv("Data/gibbon_detection.csv")
# Import site environmental variables
hainan_env <- read.csv("Data/site_env_vars.csv")
# Import list of bird species with high error rates
hainan_bird_error <- read.csv("Data/bird_species_higherror.csv")
# Import list of bird species with >5 detections with >90% confidence
hainan_bird_filter <- read.csv("Data/bird_species_abundant.csv")
# Import diet of bird species in China (Wang et al, 2021)
hainan_bird_diet <- read.csv("Data/bird_species_diet.csv")

# Standardise column name of bird diet dataset
colnames(hainan_bird_diet) <- c("Spp1","Latin.Name","Diet")
# Remove NAs in dataset
hainan_bird_diet <- na.omit(hainan_bird_diet)


#####
# 2. Gibbon call processing
# Compute site from gibbon detection directory paths
hainan_gibbons$Site <- sub("\\).*", "", sub(".*\\(", "", hainan_gibbons$path))
# Remove ".wav" at end of file name
hainan_gibbons$Filename <- substr(hainan_gibbons$file, start = 1, stop = 24)
# Separate file name into machine ID, date and time
hainan_gibbons <- hainan_gibbons %>% 
  separate(Filename, c("Machine","Date","Time"), "_", remove=F)
# Compute hour of recording from time
hainan_gibbons$Hour <- as.numeric(substr(hainan_gibbons$Time, start = 1, stop = 2))
# Select gibbon detections above the rhythm score threshold
hainan_gibbon_list <- hainan_gibbons%>%filter(ribbit>=0.275)
# Count number of reliable gibbon detections per day per hour at each site
hainan_gibbon_hour = hainan_gibbon_list %>% group_by(Site, Date, Hour)%>%count()


#####
# 3. Bird call processing
# Filter for bird detections above 90% confidence threshold
hainan_bird_hour <- hainan_bird_raw %>% filter(Prob1>=0.9) %>% 
  group_by(Site, Date, Hour, Spp1)%>%count()
# Remove bird species with high error rates
hainan_bird_hour <- hainan_bird_hour[!hainan_bird_hour$Spp1 %in% 
                                       unique(hainan_bird_error$Chinese.Name),]
# Only keep birds with >5 calls
hainan_bird_hour <- hainan_bird_hour[hainan_bird_hour$Spp1 %in% 
                                       unique(hainan_bird_filter$Chinese.Name),]
# Merge bird diet dataset with bird detections
hainan_bird_hour_select <- merge(hainan_bird_hour,hainan_bird_diet,by="Spp1")
# Merge carnivorous birds with raptors
hainan_bird_hour_select$DietClass <- ifelse(hainan_bird_hour_select$Diet=="食虫，食肉"|
                                              hainan_bird_hour_select$Diet=="食肉"|
                                              hainan_bird_hour_select$Diet=="食肉，食虫","猛禽",
                                            hainan_bird_hour_select$Diet)

# Count number of bird species detections per day per hour at each site
hainan_bird_hour_species <- hainan_bird_hour_select%>%dplyr::select(Spp1:n)
# Count number of bird foraging guild detections per day per hour at each site
hainan_bird_hour_diet <- hainan_bird_hour_select%>%
  group_by(Site,Date, Hour,DietClass)%>% summarise(num=sum(n))

#hainan_timecount<- indices_hainan%>%group_by(Site,Date)%>%count()

###
# Make species and data frames wider (each species/foraging guild a column)
hainan_bird_hour2_sp <- pivot_wider(hainan_bird_hour_species, names_from="Spp1",
                                    values_from="n")
hainan_bird_hour2_diet <- pivot_wider(hainan_bird_hour_diet, names_from="DietClass",
                                      values_from="num")

# Some machines have 5am data recorded but not properly identified, 
# this chunk is to separate them
hainan_bird_hour2_sp$Hour <- as.numeric(hainan_bird_hour2_sp$Hour)
hainan_bird_hour2_sp$Hour[is.na(hainan_bird_hour2_sp$Hour)] <- 5
hainan_bird_hour2_diet$Hour <- as.numeric(hainan_bird_hour2_diet$Hour)
hainan_bird_hour2_diet$Hour[is.na(hainan_bird_hour2_diet$Hour)] <- 5

# Convert NAs (no detection) into 0-values
hainan_bird_hour2_sp[is.na(hainan_bird_hour2_sp)] <- 0
hainan_bird_hour2_diet[is.na(hainan_bird_hour2_diet)] <- 0
# Rename columns into English abbreviations
# (raptor, insectivore, omnivore and florivore)
hainan_bird_hour3_diet <- hainan_bird_hour2_diet
colnames(hainan_bird_hour3_diet)<-c(colnames(hainan_bird_hour3_diet)[1:3],
                                    "Rapt","Inse","Omni","Flor")
# Calculate total bird activities per hour at each site
hainan_bird_hour3_diet$Total <- hainan_bird_hour3_diet$Rapt+
  hainan_bird_hour3_diet$Inse+hainan_bird_hour3_diet$Omni+
  hainan_bird_hour3_diet$Flor

#####
# 4. Prepare dataset for full analysis
# Merge bird feeding guild detections with gibbon detections
hainan_all_hour_diet <- merge(hainan_bird_hour3_diet,hainan_gibbon_hour,
                              by=c("Site", "Date", "Hour"),all=TRUE)
hainan_all_hour_diet[is.na(hainan_all_hour_diet)]=0
# Merge with soundscape indices (to include files with no bird/gibbons detected)
hainan_all_hour_diet <- merge(hainan_all_hour_diet, indices_hainan, 
                              by=c("Site", "Date", "Hour"),all=TRUE)

# Merge dataset with environmental variables of each site
hainan_all_hour_diet <- merge(hainan_all_hour_diet,hainan_env,by=c("Site"))
# Convert hours with no detections (NAs) into 0-values
hainan_all_hour_diet$Rapt[is.na(hainan_all_hour_diet$Rapt)]=0
hainan_all_hour_diet$Inse[is.na(hainan_all_hour_diet$Inse)]=0
hainan_all_hour_diet$Omni[is.na(hainan_all_hour_diet$Omni)]=0
hainan_all_hour_diet$Flor[is.na(hainan_all_hour_diet$Flor)]=0
hainan_all_hour_diet$Total[is.na(hainan_all_hour_diet$Total)]=0
hainan_all_hour_diet$n[is.na(hainan_all_hour_diet$n)]=0
# Only keep recordings between 6am and 6pm
hainan_all_hour_diet <- hainan_all_hour_diet%>%filter(Month!="202110" & (Hour >5 & Hour <18))
# Compute groups (A,B,C,D,E) of sample points
hainan_all_hour_diet$Group <- substr(hainan_all_hour_diet$Site,start=1,stop=1)

# Write full dataset into a CSV file
write.csv(hainan_all_hour_diet,"Data/hainan_all_hour_diet.csv",row.names=F)
# This data frame has Chinese letters, so it needs to be written using a different function
write_excel_csv(hainan_bird_hour2_sp,file="Data/hainan_bird_hour_sp.csv")
