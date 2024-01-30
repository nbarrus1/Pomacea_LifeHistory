#the purpose of this code is it to perform summary statistics on the variables used for modelling egg clutchs
#and standardize all the variables in preparation for the modelling effort.

#note the code from 001 needs to be ran before this code is run

#remove everything from the environment but the data we will be using for modelling

rm(list = setdiff(ls(),"summ.LILA.egg"))

#-----------------------------------------------
#####libraries#####
#-----------------------------------------------

library(plot3D)          #packaged to create 3D regression plane
library(MASS)            #model fitting package contains
library(lme4)            #model fitting package for mixed effect models GLMM models 
library(lmerTest)        #package that gets model summary information more readable
library(qpcR)            #            
library(lemon)           #package for creating more readable plot pannels
library(cowplot)         #package that puts multiple plots together
library(ggpubr)          #adds important ggplot features neccesary for publishing
library(mgcv)

#-----------------------------------------------
####check correlations####
#-----------------------------------------------

cordata <- summ.LILA.egg %>% 
  ungroup() %>% 
  dplyr::select(ave.depth.cm,depth.change,sd.depth.cm,Season,pom.dens,ave.temp.c, photoperiod) %>% 
  mutate(Season = if_else(Season == "dry",
                          true = 0,
                          false = 1))


cordata <- as.matrix(cordata)
cormatrix <- cor(cordata)
cormatrix

pairs(cordata)

#------------------------------------------------------------------------------------------
####Summary Stats####
#------------------------------------------------------------------------------------------

#explanatory variables#

#environmental variables

summary_stat <- summ.LILA.egg 

summary_stat_env <- summary_stat %>% 
  ungroup() %>% 
  group_by(Pomacea_Species) %>% 
  summarise(n = n(),
            ave.depth = mean(ave.depth.cm, na.rm = T),
            sd.depth = sd(ave.depth.cm, na.rm = T),
            max.depth = max(ave.depth.cm, na.rm = T),
            min.depth = min(ave.depth.cm, na.rm = T),
            ave.dens = mean(pom.dens, na.rm = T))
mutate(se.depth = sd.depth/sqrt(n)),
            sd.dens = sd(pom.dens, na.rm = T),
            max.dens = max(pom.dens, na.rm = T),
            min.dens = min(pom.dens, na.rm = T),
            ave.temp = mean(ave.temp.c, na.rm = T),
            sd.temp = sd(ave.temp.c, na.rm = T),
            max.temp = max(ave.temp.c, na.rm = T),
            min.temp = min(ave.temp.c, na.rm = T),
            ave.sddepth = mean(sd.depth.cm, na.rm = T),
            sd.sddepth = sd(sd.depth.cm, na.rm = T),
            max.sddepth = max(sd.depth.cm, na.rm = T),
            min.sddepth = min(sd.depth.cm, na.rm = T),
            ave.change = mean(depth.change, na.rm = T),
            sd.change = sd(depth.change, na.rm = T),
            max.change = max(depth.change, na.rm = T),
            min.change = min(depth.change, na.rm = T),
            ave.photo = mean(photoperiod, na.rm = T),
            min.photo = min(photoperiod, na.rm = T),
            max.photo = max(photoperiod, na.rm = T),
            sd.photo = sd(photoperiod, na.rm = T)) 

summary_stat_env

#apple snail densities


summary_stat <- summ.LILA.egg 

summary_stat_dens <- summary_stat%>% 
  dplyr::select(pom.dens, Pomacea_Species, Cell, year) %>% 
  ungroup() %>% 
  group_by(Pomacea_Species, Cell, year) %>% 
  summarise(ave.dens = mean(pom.dens, na.rm = T))

summary_stat_dens

####egg clutch summary####

summary_stat_counts <- summ.LILA.egg %>% 
  mutate(month = month(Date)) %>% 
  ungroup() %>% 
  group_by(Pomacea_Species, year, Cell) %>% 
  summarise(ave.count = mean(count, na.rm = T),
            sum.count = sum(count, na.rm = T),
            max.count = max(count, na.rm = T),
            n = n()) %>% 
  mutate(spp.count = sum(sum.count))

summary_stat_counts

#---------------------------------------------------------------------------------------
#####standardize all explanotory variables variables (for use in mixed effect modeling)#####
#---------------------------------------------------------------------------------------
summ.LILA.egg <- summ.LILA.egg %>% 
  ungroup() %>% 
  mutate(zavedepth = (ave.depth.cm - mean(ave.depth.cm, na.rm = TRUE))/ sd(ave.depth.cm, na.rm = T),
         zsddepth = (sd.depth.cm - mean(sd.depth.cm, na.rm = TRUE))/ sd(sd.depth.cm, na.rm = T),
         zavetemp = (ave.temp.c - mean(ave.temp.c, na.rm = TRUE))/ sd(ave.temp.c, na.rm = T),
         zdeltadepth = (depth.change - mean(depth.change, na.rm = TRUE))/ sd(depth.change, na.rm = T),
         zphotoperiod = (photoperiod - mean(photoperiod, na.rm = TRUE))/ sd(photoperiod, na.rm = T))

