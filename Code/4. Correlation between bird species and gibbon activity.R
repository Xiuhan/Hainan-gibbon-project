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
library(ggpubr)

#####
# 1. Load CSV dataset created in last section
hainan_bird_hour_sp <- read.csv("Data/hainan_bird_hour_sp.csv")
hainan_bird_site <- read.csv("Data/hainan_bird_site.csv")

#####
# 2. Prepare dataset for correlation test between gibbons and bird species
# Compute month of recording (for excluding abnormal data)
hainan_bird_hour_sp$Month=substr(hainan_bird_hour_sp$Date,1,6)
# Calculate total and mean activity of bird species at each site
hainan_bird_hour_spsum <- hainan_bird_hour_sp%>%
  # Remove recordings in abnormal months and time of the day
  filter(Month!="202110" & (Hour >5 & Hour <18))%>%
  # Group by site
  group_by(Site)%>%
  # Count the number of activities detected for each species at each site
  summarise(白腹凤鹛T = sum(白腹凤鹛),斑头鸺鹠T=sum(斑头鸺鹠),叉尾太阳鸟T=sum(叉尾太阳鸟),
            橙腹叶鹎T=sum(橙腹叶鹎),大鹰鹃T=sum(大鹰鹃),发冠卷尾T=sum(发冠卷尾),
            海南蓝仙鹟T=sum(海南蓝仙鹟),黑短脚鹎T=sum(黑短脚鹎),黑喉噪鹛T=sum(黑喉噪鹛),
            红头穗鹛T=sum(红头穗鹛),红头咬鹃T=sum(红头咬鹃),红胸啄花鸟T=sum(红胸啄花鸟),
            红嘴蓝鹊T=sum(红嘴蓝鹊),黄颊山雀T=sum(黄颊山雀),黄嘴角鸮T=sum(黄嘴角鸮),
            黄嘴栗啄木鸟T=sum(黄嘴栗啄木鸟),灰树鹊T=sum(灰树鹊),栗背短脚鹎T=sum(栗背短脚鹎),
            领角鸮T=sum(领角鸮),领鸺鹠T=sum(领鸺鹠),绿翅短脚鹎T=sum(绿翅短脚鹎),
            蛇雕T=sum(蛇雕),四声杜鹃T=sum(四声杜鹃),乌鹃T=sum(乌鹃),朱背啄花鸟T=sum(朱背啄花鸟),
            棕颈钩嘴鹛T=sum(棕颈钩嘴鹛),棕脸鹟莺T=sum(棕脸鹟莺),
            # Calculate average hourly activity for each species at each site
            白腹凤鹛D = mean(白腹凤鹛),斑头鸺鹠D=mean(斑头鸺鹠),叉尾太阳鸟D=mean(叉尾太阳鸟),
            橙腹叶鹎D=mean(橙腹叶鹎),大鹰鹃D=mean(大鹰鹃),发冠卷尾D=mean(发冠卷尾),
            海南蓝仙鹟D=mean(海南蓝仙鹟),黑短脚鹎D=mean(黑短脚鹎),黑喉噪鹛D=mean(黑喉噪鹛),
            红头穗鹛D=mean(红头穗鹛),红头咬鹃D=mean(红头咬鹃),红胸啄花鸟D=mean(红胸啄花鸟),
            红嘴蓝鹊D=mean(红嘴蓝鹊),黄颊山雀D=mean(黄颊山雀),黄嘴角鸮D=mean(黄嘴角鸮),
            黄嘴栗啄木鸟D=mean(黄嘴栗啄木鸟),灰树鹊D=mean(灰树鹊),栗背短脚鹎D=mean(栗背短脚鹎),
            领角鸮D=mean(领角鸮),领鸺鹠D=mean(领鸺鹠),绿翅短脚鹎D=mean(绿翅短脚鹎),
            蛇雕D=mean(蛇雕),四声杜鹃D=mean(四声杜鹃),乌鹃D=mean(乌鹃),
            朱背啄花鸟D=mean(朱背啄花鸟),棕颈钩嘴鹛D=mean(棕颈钩嘴鹛),棕脸鹟莺D=mean(棕脸鹟莺)
  )

# Calculate daily peak activity of each species at each site
hainan_bird_day_spsum <- hainan_bird_hour_sp%>%
  # Remove recordings in abnormal months and time of the day
  filter(Month!="202110" & (Hour >5 & Hour <18))%>%
  # Group by site and day
  group_by(Site,Date)%>%
  # Count the maximum hourly activity per species at each site
  summarise(白腹凤鹛Max = max(白腹凤鹛),斑头鸺鹠Max=max(斑头鸺鹠),
            叉尾太阳鸟Max=max(叉尾太阳鸟),橙腹叶鹎Max=max(橙腹叶鹎),
            大鹰鹃Max=max(大鹰鹃),发冠卷尾Max=max(发冠卷尾),
            海南蓝仙鹟Max=max(海南蓝仙鹟),黑短脚鹎Max=max(黑短脚鹎),
            黑喉噪鹛Max=max(黑喉噪鹛),红头穗鹛Max=max(红头穗鹛),
            红头咬鹃Max=max(红头咬鹃),红胸啄花鸟Max=max(红胸啄花鸟),
            红嘴蓝鹊Max=max(红嘴蓝鹊),黄颊山雀Max=max(黄颊山雀),
            黄嘴角鸮Max=max(黄嘴角鸮),黄嘴栗啄木鸟Max=max(黄嘴栗啄木鸟),
            灰树鹊Max=max(灰树鹊),栗背短脚鹎Max=max(栗背短脚鹎),
            领角鸮Max=max(领角鸮),领鸺鹠Max=max(领鸺鹠),绿翅短脚鹎Max=max(绿翅短脚鹎),
            蛇雕Max=max(蛇雕),四声杜鹃Max=max(四声杜鹃),乌鹃Max=max(乌鹃),
            朱背啄花鸟Max=max(朱背啄花鸟),棕颈钩嘴鹛Max=max(棕颈钩嘴鹛),
            棕脸鹟莺Max=max(棕脸鹟莺)
            
  )

# Calculate mean daily peak activity per species at each site
# Group by site
hainan_bird_day_meanmax <- hainan_bird_day_spsum%>%group_by(Site)%>%
  # Calculate mean daily peak activity per species
  summarise(白腹凤鹛AMax = mean(白腹凤鹛Max),斑头鸺鹠AMax=mean(斑头鸺鹠Max),
            叉尾太阳鸟AMax=mean(叉尾太阳鸟Max),橙腹叶鹎AMax=mean(橙腹叶鹎Max),
            大鹰鹃AMax=mean(大鹰鹃Max),发冠卷尾AMax=mean(发冠卷尾Max),
            海南蓝仙鹟AMax=mean(海南蓝仙鹟Max),黑短脚鹎AMax=mean(黑短脚鹎Max),
            黑喉噪鹛AMax=mean(黑喉噪鹛Max),红头穗鹛AMax=mean(红头穗鹛Max),
            红头咬鹃AMax=mean(红头咬鹃Max),红胸啄花鸟AMax=mean(红胸啄花鸟Max),
            红嘴蓝鹊AMax=mean(红嘴蓝鹊Max),黄颊山雀AMax=mean(黄颊山雀Max),
            黄嘴角鸮AMax=mean(黄嘴角鸮Max),黄嘴栗啄木鸟AMax=mean(黄嘴栗啄木鸟Max),
            灰树鹊AMax=mean(灰树鹊Max),栗背短脚鹎AMax=mean(栗背短脚鹎Max),
            领角鸮AMax=mean(领角鸮Max),领鸺鹠AMax=mean(领鸺鹠Max),
            绿翅短脚鹎AMax=mean(绿翅短脚鹎Max),蛇雕AMax=mean(蛇雕Max),
            四声杜鹃AMax=mean(四声杜鹃Max),乌鹃AMax=mean(乌鹃Max),
            朱背啄花鸟AMax=mean(朱背啄花鸟Max),棕颈钩嘴鹛AMax=mean(棕颈钩嘴鹛Max),
            棕脸鹟莺AMax=mean(棕脸鹟莺Max)
  )

# Extract average daily peak and hourly activity of gibbons
hainan_gibbon_activity <- hainan_bird_site %>% dplyr::select (Site, GibbonAMax, GibbonD)
# Merge three datasets together
hainan_all_abundance <- merge (hainan_bird_hour_spsum,hainan_gibbon_activity,by="Site")
hainan_all_abundance <-merge(hainan_all_abundance,hainan_bird_day_meanmax,by="Site")

#####
# 3. Pearson's correlation test between bird species and gibbon activity
## 3.1. Average daily peak activity
cor.test(hainan_all_abundance$白腹凤鹛AMax,hainan_all_abundance$GibbonAMax )
cor.test(hainan_all_abundance$斑头鸺鹠AMax,hainan_all_abundance$GibbonAMax ) 
cor.test(hainan_all_abundance$叉尾太阳鸟AMax,hainan_all_abundance$GibbonAMax ) 
cor.test(hainan_all_abundance$橙腹叶鹎AMax,hainan_all_abundance$GibbonAMax )
cor.test(hainan_all_abundance$大鹰鹃AMax,hainan_all_abundance$GibbonAMax )
cor.test(hainan_all_abundance$发冠卷尾AMax,hainan_all_abundance$GibbonAMax )
cor.test(hainan_all_abundance$海南蓝仙鹟AMax,hainan_all_abundance$GibbonAMax )
cor.test(hainan_all_abundance$黑短脚鹎AMax,hainan_all_abundance$GibbonAMax )
cor.test(hainan_all_abundance$黑喉噪鹛AMax,hainan_all_abundance$GibbonAMax )
cor.test(hainan_all_abundance$红头穗鹛AMax,hainan_all_abundance$GibbonAMax ) 
cor.test(hainan_all_abundance$红头咬鹃AMax,hainan_all_abundance$GibbonAMax ) # p=0.07
cor.test(hainan_all_abundance$红胸啄花鸟AMax,hainan_all_abundance$GibbonAMax )
cor.test(hainan_all_abundance$红嘴蓝鹊AMax,hainan_all_abundance$GibbonAMax )
cor.test(hainan_all_abundance$黄颊山雀AMax,hainan_all_abundance$GibbonAMax )
cor.test(hainan_all_abundance$黄嘴角鸮AMax,hainan_all_abundance$GibbonAMax )
cor.test(hainan_all_abundance$黄嘴栗啄木鸟AMax,hainan_all_abundance$GibbonAMax )
cor.test(hainan_all_abundance$灰树鹊AMax,hainan_all_abundance$GibbonAMax )
cor.test(hainan_all_abundance$栗背短脚鹎AMax,hainan_all_abundance$GibbonAMax )
cor.test(hainan_all_abundance$领角鸮AMax,hainan_all_abundance$GibbonAMax ) # p<0.01, but too few
cor.test(hainan_all_abundance$领鸺鹠AMax,hainan_all_abundance$GibbonAMax )
cor.test(hainan_all_abundance$绿翅短脚鹎AMax,hainan_all_abundance$GibbonAMax )
cor.test(hainan_all_abundance$蛇雕AMax,hainan_all_abundance$GibbonAMax )
cor.test(hainan_all_abundance$四声杜鹃AMax,hainan_all_abundance$GibbonAMax )
cor.test(hainan_all_abundance$乌鹃AMax,hainan_all_abundance$GibbonAMax ) # p=0.03
cor.test(hainan_all_abundance$朱背啄花鸟AMax,hainan_all_abundance$GibbonAMax )
cor.test(hainan_all_abundance$棕颈钩嘴鹛AMax,hainan_all_abundance$GibbonAMax )
cor.test(hainan_all_abundance$棕脸鹟莺AMax,hainan_all_abundance$GibbonAMax )

## 3.1. Average hourly activity
cor.test(hainan_all_abundance$白腹凤鹛D,hainan_all_abundance$GibbonD )
cor.test(hainan_all_abundance$斑头鸺鹠D,hainan_all_abundance$GibbonD ) # p=0.08
cor.test(hainan_all_abundance$叉尾太阳鸟D,hainan_all_abundance$GibbonD ) 
cor.test(hainan_all_abundance$橙腹叶鹎D,hainan_all_abundance$GibbonD )
cor.test(hainan_all_abundance$大鹰鹃D,hainan_all_abundance$GibbonD )
cor.test(hainan_all_abundance$发冠卷尾D,hainan_all_abundance$GibbonD )
cor.test(hainan_all_abundance$海南蓝仙鹟D,hainan_all_abundance$GibbonD )
cor.test(hainan_all_abundance$黑短脚鹎D,hainan_all_abundance$GibbonD )
cor.test(hainan_all_abundance$黑喉噪鹛D,hainan_all_abundance$GibbonD )
cor.test(hainan_all_abundance$红头穗鹛D,hainan_all_abundance$GibbonD ) 
cor.test(hainan_all_abundance$红头咬鹃D,hainan_all_abundance$GibbonD ) # p=0.10
cor.test(hainan_all_abundance$红胸啄花鸟D,hainan_all_abundance$GibbonD )
cor.test(hainan_all_abundance$红嘴蓝鹊D,hainan_all_abundance$GibbonD )
cor.test(hainan_all_abundance$黄颊山雀D,hainan_all_abundance$GibbonD )
cor.test(hainan_all_abundance$黄嘴角鸮D,hainan_all_abundance$GibbonD )
cor.test(hainan_all_abundance$黄嘴栗啄木鸟D,hainan_all_abundance$GibbonD )
cor.test(hainan_all_abundance$灰树鹊D,hainan_all_abundance$GibbonD )
cor.test(hainan_all_abundance$栗背短脚鹎D,hainan_all_abundance$GibbonD )
cor.test(hainan_all_abundance$领角鸮D,hainan_all_abundance$GibbonD ) ## p=0.04, but too few
cor.test(hainan_all_abundance$领鸺鹠D,hainan_all_abundance$GibbonD ) 
cor.test(hainan_all_abundance$绿翅短脚鹎D,hainan_all_abundance$GibbonD )
cor.test(hainan_all_abundance$蛇雕D,hainan_all_abundance$GibbonD ) # p=0.09
cor.test(hainan_all_abundance$四声杜鹃D,hainan_all_abundance$GibbonD )
cor.test(hainan_all_abundance$乌鹃D,hainan_all_abundance$GibbonD )
cor.test(hainan_all_abundance$朱背啄花鸟D,hainan_all_abundance$GibbonD )
cor.test(hainan_all_abundance$棕颈钩嘴鹛D,hainan_all_abundance$GibbonD )
cor.test(hainan_all_abundance$棕脸鹟莺D,hainan_all_abundance$GibbonD )


#####
# 4. Plot correlation for species with p<0.1
# 4.1. Otus lettia (significant but too few, so not reported in the paper)
# Average daily peak activity
ggscatter(hainan_all_abundance, x = "GibbonAMax", y = "领角鸮AMax", 
          # Add a regression line and confidence band
          add = "reg.line", conf.int = TRUE,
          # Add colours to regression line and confidence band
          add.params = list(color = "black", fill = "gray", alpha = 0.5), 
          cor.method = "pearson")+ 
  # Classic plot theme and axis label font size
  theme_classic() + theme(axis.text=element_text(size=15))+
  labs(x = "Mean daily peak gibbon activity",
       y = "Mean daily peak Otus lettia activity") 
# Average daily peak activity
ggscatter(hainan_all_abundance, x = "GibbonD", y = "领角鸮D", 
          add = "reg.line", conf.int = TRUE,
          add.params = list(color = "black", fill = "gray", alpha = 0.5), 
          cor.method = "pearson")+
  coord_cartesian(xlim = c(0.05, 0.3)) + 
  theme_classic() + theme(axis.text=element_text(size=15))+
  labs(x = "Mean hourly gibbon activity",
       y = "Mean hourly Otus lettia activity") 

# 4.2. Surniculus lugubris
# (check 4.1. annotations for functions of each row)
ggscatter(hainan_all_abundance, x = "GibbonAMax", y = "乌鹃AMax", 
          add = "reg.line", conf.int = TRUE,
          add.params = list(color = "black", fill = "gray", alpha = 0.5), 
          cor.method = "pearson") + 
  theme_classic()+theme(axis.text=element_text(size=15))+
  labs(x = "Mean daily peak gibbon activity",
       y = "Mean daily peak Surniculus lugubris activity") 
ggscatter(hainan_all_abundance, x = "GibbonD", y = "乌鹃D", 
          add = "reg.line", conf.int = TRUE,
          add.params = list(color = "black", fill = "gray", alpha = 0.5), 
          cor.method = "pearson")+
  coord_cartesian(xlim = c(0.05, 0.3))+ 
  theme_classic()+theme(axis.text=element_text(size=15))+
  labs(x = "Mean hourly gibbon activity",
       y = "Mean hourly Surniculus lugubris activity") 

# 4.3. Harpactes erythrocephalus
# (check 4.1. annotations for functions of each row)
ggscatter(hainan_all_abundance, x = "GibbonAMax", y = "红头咬鹃AMax", 
          add = "reg.line", conf.int = TRUE,
          add.params = list(color = "black", fill = "gray", alpha = 0.5), 
          cor.method = "pearson") + 
  theme_classic()+theme(axis.text=element_text(size=15))+
  labs(x = "Mean daily peak gibbon activity",
       y = "Mean daily peak Harpactes erythrocephalus activity") 
ggscatter(hainan_all_abundance, x = "GibbonD", y = "红头咬鹃D", 
          add = "reg.line", conf.int = TRUE,
          add.params = list(color = "black", fill = "gray", alpha = 0.5), 
          cor.method = "pearson")+
  coord_cartesian(xlim = c(0.05, 0.3)) + 
  theme_classic()+theme(axis.text=element_text(size=15))+
  labs(x = "Mean hourly gibbon activity",
       y = "Mean hourly Harpactes erythrocephalus activity") 

# 4.4. Spilornis cheela
# (check 4.1. annotations for functions of each row)
ggscatter(hainan_all_abundance, x = "GibbonAMax", y = "蛇雕AMax", 
          add = "reg.line", conf.int = TRUE,
          add.params = list(color = "black", fill = "gray", alpha = 0.5), 
          cor.method = "pearson") + 
  theme_classic()+theme(axis.text=element_text(size=15))+
  labs(x = "Mean daily peak gibbon activity",
       y = "Mean daily peak Spilornis cheela activity") 
ggscatter(hainan_all_abundance, x = "GibbonD", y = "蛇雕D", 
          add = "reg.line", conf.int = TRUE,
          add.params = list(color = "black", fill = "gray", alpha = 0.5), 
          cor.method = "pearson")+
  coord_cartesian(xlim = c(0.05, 0.3)) + 
  theme_classic()+theme(axis.text=element_text(size=15))+
  labs(x = "Mean hourly gibbon activity",
       y = "Mean hourly Spilornis cheela activity") 

# 4.5. Glaucidium cuculoides
# (check 4.1. annotations for functions of each row)
ggscatter(hainan_all_abundance, x = "GibbonAMax", y = "斑头鸺鹠AMax", 
          add = "reg.line", conf.int = TRUE,
          add.params = list(color = "black", fill = "gray", alpha = 0.5), 
          cor.method = "pearson") + 
  theme_classic()+theme(axis.text=element_text(size=15))+
  labs(x = "Mean daily peak gibbon activity",
       y = "Mean daily peak Glaucidium cuculoides activity") 
ggscatter(hainan_all_abundance, x = "GibbonD", y = "斑头鸺鹠D", 
          add = "reg.line", conf.int = TRUE,
          add.params = list(color = "black", fill = "gray", alpha = 0.5), 
          cor.method = "pearson")+
  coord_cartesian(xlim = c(0.05, 0.3))+
  theme_classic()+theme(axis.text=element_text(size=15))+
labs(x = "Mean hourly gibbon activity",
     y = "Mean hourly Glaucidium cuculoides activity") 
