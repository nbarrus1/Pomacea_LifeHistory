#the purpose of this code is to model the egg clutch counts using a GLMM negative binomial and compare models using AICc

#remove everything from the environment but the data we will be using for modelling

rm(list = setdiff(ls(),"summ.LILA.egg"))

#-----------------------------
#####load in libraries######
#-----------------------------


install.packages("tidyverse")

library(tidyverse)
library(plot3D)          #packaged to create 3D regression plane
library(MASS)            #model fitting package contains
library(lme4)            #model fitting package for mixed effect models GLMM models 
library(lmerTest)        #package that gets model summary information more readable
library(qpcR)            #            
library(lemon)           #package for creating more readable plot pannels
library(cowplot)         #package that puts multiple plots together
library(ggpubr)          #adds important ggplot features neccesary for publishing


#--------------------------------------------
####Zscores####
#--------------------------------------------
summ.LILA.egg <- summ.LILA.egg %>% 
  ungroup() %>% 
  mutate(zavedepth = (ave.depth.cm - mean(ave.depth.cm, na.rm = TRUE))/ sd(ave.depth.cm, na.rm = T),
         zsddepth = (sd.depth.cm - mean(sd.depth.cm, na.rm = TRUE))/ sd(sd.depth.cm, na.rm = T),
         zavetemp = (ave.temp.c - mean(ave.temp.c, na.rm = TRUE))/ sd(ave.temp.c, na.rm = T),
         zdeltadepth = (depth.change - mean(depth.change, na.rm = TRUE))/ sd(depth.change, na.rm = T),
         zphotoperiod = (photoperiod - mean(photoperiod, na.rm = TRUE))/ sd(photoperiod, na.rm = T))

mean(summ.LILA.egg$ave.depth.cm, na.rm = TRUE)
mean(summ.LILA.egg$ave.temp.c,na.rm = TRUE)
sd(summ.LILA.egg$ave.depth.cm, na.rm = TRUE)
sd(summ.LILA.egg$ave.temp.c,na.rm = TRUE)
max(summ.LILA.egg$ave.depth.cm)
max(summ.LILA.egg$ave.temp.c)
min(summ.LILA.egg$ave.depth.cm)
min(summ.LILA.egg$ave.temp.c)
#-------------------------------------------
#####model construction for P. paludosa#####
#-------------------------------------------

#####filter to get just the paludosa data
ppal.egg <- summ.LILA.egg %>% 
  ungroup() %>% 
  filter(Pomacea_Species == "Paludosa") %>% 
  mutate(zpaldens = (pom.dens - mean(pom.dens, na.rm = TRUE))/ sd(pom.dens, na.rm = T))

ppal.egg[,14:ncol(ppal.egg)]

#create a list of models for comparison

#avedepth only
pal.mefit1 <- glmer.nb(formula = count ~ zavedepth + (1|Cell), data = ppal.egg)

#deltadepth only
pal.mefit2 <- glmer.nb(formula = count ~ zdeltadepth+ (1|Cell), data = ppal.egg)

#sddepth only
pal.mefit3 <- glmer.nb(formula = count ~ zsddepth+ (1|Cell), data = ppal.egg)

#avedepth and deltadepth
pal.mefit4 <- glmer.nb(formula = count ~ zavedepth + zdeltadepth+ (1|Cell), data = ppal.egg)

#sd depth and deltadepth
pal.mefit5 <- glmer.nb(formula = count ~ zsddepth + zdeltadepth+ (1|Cell), data = ppal.egg)

#avedepth and sddepth
pal.mefit6 <- glmer.nb(formula = count ~ zavedepth + zsddepth+ (1|Cell), data = ppal.egg)

#ave depth sd depth, deltadepth
pal.mefit7 <- glmer.nb(formula = count ~ zavedepth + zsddepth + zdeltadepth+ (1|Cell), data = ppal.egg)

#average depth binomial 
pal.mefit8 <- glmer.nb(formula = count ~ zavedepth + I(zavedepth^2)+ (1|Cell), data = ppal.egg)   ####singularity issues

tt <- getME(pal.mefit8,"theta")   ###check singularity
ll <- getME(pal.mefit8,"lower")  
min(tt[ll==0])                    ###shoun't have a problem with singularity
derivs1 <- pal.mefit8@optinfo$derivs                       ####check gradient
sc_grad1 <- with(derivs1,solve(Hessian,gradient))
max(abs(sc_grad1))                                         ####shouldnt have issues try more interations
pal.mefit8.2 <- update(pal.mefit8,control=glmerControl(optCtrl=list(maxfun=2e4)))  ###fixed the problem
ss <- getME(pal.mefit8,c("theta","fixef"))                                                          ###start with previous theta
pal.mefit8.3 <- update(pal.mefit8,start=ss,control=glmerControl(optCtrl=list(maxfun=2e4)))         ###still has issues
pal.mefit8.4<- update(pal.mefit8,start=ss,control=glmerControl(optimizer="bobyqa",                 ###different optimzer works
                                                               optCtrl=list(maxfun=2e5))) 

#season only
pal.mefit9 <- glmer.nb(formula = count ~ Season + (1|Cell), data = ppal.egg)

#depth season interaction
pal.mefit10 <- glmer.nb(formula = count ~ zavedepth*Season + (1|Cell), data = ppal.egg)

#deltadepth, season interaction, this should not be included
pal.mefit11 <- glmer.nb(formula = count ~ zdeltadepth*Season + (1|Cell), data = ppal.egg)

#sdseason interaction
pal.mefit12 <- glmer.nb(formula = count ~ zsddepth*Season + (1|Cell), data = ppal.egg)

#season depth binomial interaction
pal.mefit13 <- glmer.nb(formula = count ~ zavedepth*Season + I(zavedepth^2)*Season+ (1|Cell), data = ppal.egg)

#density only
pal.mefit14 <- glmer.nb(formula = count ~ zpaldens + (1|Cell), data = ppal.egg)

tt <- getME(pal.mefit14,"theta")   ###check singularity
ll <- getME(pal.mefit14,"lower")  
min(tt[ll==0])                    ###shoun't have a problem with singularity
derivs1 <- pal.mefit14@optinfo$derivs                       ####check gradient
sc_grad1 <- with(derivs1,solve(Hessian,gradient))
max(abs(sc_grad1))                                         ####shouldnt have issues try more interations

pal.mefit14.2 <- update(pal.mefit14,control=glmerControl(optCtrl=list(maxfun=2e4)))  ###still issues try starting at previous theta
ss <- getME(pal.mefit14,c("theta","fixef"))                                                          ###start with previous theta
pal.mefit14.3 <- update(pal.mefit14,start=ss,control=glmerControl(optCtrl=list(maxfun=2e4)))         ###still has issues
pal.mefit14.4<- update(pal.mefit14,start=ss,control=glmerControl(optimizer="bobyqa",                 ###different optimzer works
                                                                 optCtrl=list(maxfun=2e5)))  

#depth density
pal.mefit15 <- glmer.nb(formula = count ~ zavedepth + zpaldens + (1|Cell), data = ppal.egg)

#depth density season
pal.mefit16 <- glmer.nb(formula = count ~ zavedepth + zpaldens + Season + (1|Cell), data = ppal.egg)

tt <- getME(pal.mefit16,"theta")   ###check singularity
ll <- getME(pal.mefit16,"lower")  
min(tt[ll==0])                    ###shoun't have a problem with singularity
derivs1 <- pal.mefit16@optinfo$derivs                       ####check gradient
sc_grad1 <- with(derivs1,solve(Hessian,gradient))
max(abs(sc_grad1))                                         ####shouldnt have issues try more interations

pal.mefit16.2 <- update(pal.mefit16,control=glmerControl(optCtrl=list(maxfun=2e4)))  ###still issues try starting at previous theta
ss <- getME(pal.mefit16,c("theta","fixef"))                                                          ###start with previous theta
pal.mefit16.3 <- update(pal.mefit16,start=ss,control=glmerControl(optCtrl=list(maxfun=2e4)))         ###worked
pal.mefit16.4 <- update(pal.mefit16,start=ss,control=glmerControl(optimizer="bobyqa",                 ###different optimzer works
                                                                  optCtrl=list(maxfun=2e5)))  
pal.mefit16.5 <- update(pal.mefit16,start=ss,control=glmerControl(optimizer="Nelder_Mead",                 ###different optimzer works
                                                                  optCtrl=list(maxfun=2e5)))  
pal.mefit16.6 <- update(pal.mefit16,start=ss,control=glmerControl(optimizer="nlminbwrap",                 ###different optimzer works
                                                                  optCtrl=list(maxfun=2e5)))  
pal.mefit16.7 <- update(pal.mefit16,start=ss,control=glmerControl(optimizer="nloptwrap",                 ###different optimzer works
                                                                  optCtrl=list(maxfun=2e5,
                                                                               ftol_abs=1e-4,
                                                                               xtol_abs=1e-4)))         ###still has issues remove from list of models

#density season
pal.mefit17 <- glmer.nb(formula = count ~ zpaldens + Season + (1|Cell), data = ppal.egg)

tt <- getME(pal.mefit17,"theta")   ###check singularity
ll <- getME(pal.mefit17,"lower")  
min(tt[ll==0])                    ###shoun't have a problem with singularity
derivs1 <- pal.mefit17@optinfo$derivs                       ####check gradient
sc_grad1 <- with(derivs1,solve(Hessian,gradient))
max(abs(sc_grad1))                                         ####shouldnt have issues try more interations

pal.mefit17.2 <- update(pal.mefit17,control=glmerControl(optCtrl=list(maxfun=2e4))) #good

#season density interaction
pal.mefit18 <- glmer.nb(formula = count ~ Season*zpaldens+ (1|Cell), data = ppal.egg)

#density, binomial depth
pal.mefit19 <- glmer.nb(formula = count ~ zpaldens + zavedepth + I(zavedepth^2)+ (1|Cell), data = ppal.egg)

#season binomial depth
pal.mefit20 <- glmer.nb(formula = count ~ Season + zavedepth + I(zavedepth^2)+ (1|Cell), data = ppal.egg)

#density, season, binomial depth
pal.mefit21 <- glmer.nb(formula = count ~ zpaldens + Season + zavedepth + I(zavedepth^2)+ (1|Cell), data = ppal.egg)

#density, season x binomial depth 
pal.mefit22 <- glmer.nb(formula = count ~ zpaldens + Season + zavedepth + I(zavedepth^2) +
                          Season*zavedepth + Season* I(zavedepth^2)+ (1|Cell), data = ppal.egg)

#temp alone
pal.mefit23 <- glmer.nb(formula = count ~ zavetemp + (1|Cell), data = ppal.egg)

#temp depth
pal.mefit24 <- glmer.nb(formula = count ~ zavetemp + zavedepth + (1|Cell), data = ppal.egg)

#temp x depth interaction
pal.mefit25 <- glmer.nb(formula = count ~ zavetemp + zavedepth +zavetemp*zavedepth+ (1|Cell), data = ppal.egg)

summary(pal.mefit25)
#temp, depth ,dens
pal.mefit26 <- glmer.nb(formula = count ~ zavetemp + zavedepth + zpaldens + (1|Cell), data = ppal.egg)

#temp, dens
pal.mefit27 <- glmer.nb(formula = count ~ zavetemp + zpaldens+ (1|Cell), data = ppal.egg)

#temp, binomial depth
pal.mefit28 <- glmer.nb(formula = count ~ zavetemp + zavedepth + I(zavedepth^2)+ (1|Cell), data = ppal.egg)

#temp, binomial depth, density
pal.mefit29 <- glmer.nb(formula = count ~ zavetemp + zpaldens + zavedepth + I(zavedepth^2)+ (1|Cell), data = ppal.egg)

#density, temp x binomial depth
pal.mefit30 <- glmer.nb(formula = count ~ zpaldens + zavetemp + zavedepth + I(zavedepth^2) +                          
                          zavetemp*zavedepth + zavetemp* I(zavedepth^2)+ (1|Cell), data = ppal.egg)  ###convergence issues needs fixed
ss <- getME(pal.mefit30,c("theta","fixef"))                                                          ###start with previous theta
pal.mefit30.2 <- update(pal.mefit30,start=ss,control=glmerControl(optCtrl=list(maxfun=2e4)))         ###still has issues
pal.mefit30.3<- update(pal.mefit30,start=ss,control=glmerControl(optimizer="bobyqa",                 ###different optimzer works
                                                                 optCtrl=list(maxfun=2e5)))  

#temp x binomial depth
pal.mefit31 <- glmer.nb(formula = count ~ zavedepth*zavetemp + I(zavedepth^2)*zavetemp+ (1|Cell), data = ppal.egg)

summary(pal.mefit31)
#photo
pal.mefit32 <- glmer.nb(formula = count ~ zphotoperiod+ (1|Cell), data = ppal.egg)

#dens, photoperiod x binomial depth
pal.mefit33 <- glmer.nb(formula = count ~ zpaldens + zphotoperiod + zavedepth + I(zavedepth^2) +                                 #best model
                          zphotoperiod*zavedepth + zphotoperiod* I(zavedepth^2)+ (1|Cell), data = ppal.egg)

tt <- getME(pal.mefit33,"theta")   ###check singularity
ll <- getME(pal.mefit33,"lower")  
min(tt[ll==0])                    ###shoun't have a problem with singularity
derivs1 <- pal.mefit33@optinfo$derivs                       ####check gradient
sc_grad1 <- with(derivs1,solve(Hessian,gradient))
max(abs(sc_grad1))                                         ####shouldnt have issues try more interations

pal.mefit33.2 <- update(pal.mefit33,control=glmerControl(optCtrl=list(maxfun=2e4)))  ###still issues try starting at previous theta
ss <- getME(pal.mefit33,c("theta","fixef"))                                                          ###start with previous theta
pal.mefit33.3 <- update(pal.mefit33,start=ss,control=glmerControl(optCtrl=list(maxfun=2e4)))         ###worked

#photoperiod x binomial depth
pal.mefit34 <- glmer.nb(formula = count ~ zphotoperiod + zavedepth + I(zavedepth^2) +                                 
                          zphotoperiod*zavedepth + zphotoperiod* I(zavedepth^2)+ (1|Cell), data = ppal.egg )

#density, photoperiod, binomial depth
pal.mefit35 <-  glmer.nb(formula = count ~zpaldens + zphotoperiod + zavedepth + I(zavedepth^2)+ (1|Cell), data = ppal.egg)

#pca season
pal.mefit36<- glmer.nb(formula = count ~ season.pca + (1|Cell), data = ppal.egg)

#density, season.pca x binomial depth 
pal.mefit37 <- glmer.nb(formula = count ~ zpaldens + season.pca + zavedepth + I(zavedepth^2) +                                 #best model
                          season.pca*zavedepth + season.pca* I(zavedepth^2)+ (1|Cell), data = ppal.egg)

tt <- getME(pal.mefit37,"theta")   ###check singularity
ll <- getME(pal.mefit37,"lower")  
min(tt[ll==0])                    ###shoun't have a problem with singularity
derivs1 <- pal.mefit37@optinfo$derivs                       ####check gradient
sc_grad1 <- with(derivs1,solve(Hessian,gradient))
max(abs(sc_grad1))                                         ####shouldnt have issues try more interations

pal.mefit37.2 <- update(pal.mefit37,control=glmerControl(optCtrl=list(maxfun=2e4)))  ###still issues try starting at previous theta
ss <- getME(pal.mefit37,c("theta","fixef"))                                                          ###start with previous theta
pal.mefit37.3 <- update(pal.mefit37,start=ss,control=glmerControl(optCtrl=list(maxfun=2e4)))         ###worked
pal.mefit37.4<- update(pal.mefit37,start=ss,control=glmerControl(optimizer="bobyqa",                 ###different optimzer works
                                                                 optCtrl=list(maxfun=2e5)))  #it worked!

#season.pca x binomial depth
pal.mefit38 <- glmer.nb(formula = count ~ season.pca + zavedepth + I(zavedepth^2) +                                 
                          season.pca*zavedepth + season.pca* I(zavedepth^2)+ (1|Cell), data = ppal.egg )

#compare all the models using AICc values

pal.mes <- c(pal.mefit1,pal.mefit2,pal.mefit3,pal.mefit4,pal.mefit5,pal.mefit6,pal.mefit7,pal.mefit8.2,pal.mefit9,pal.mefit10,
             pal.mefit12,pal.mefit13,pal.mefit14.4,pal.mefit15,pal.mefit17.2,pal.mefit18,pal.mefit19,pal.mefit20,
             pal.mefit21,pal.mefit22,pal.mefit23,pal.mefit24,pal.mefit25,pal.mefit26,pal.mefit27,pal.mefit28,pal.mefit29,pal.mefit30.3,
             pal.mefit31,pal.mefit32, pal.mefit33.3, pal.mefit34, pal.mefit35,pal.mefit36,pal.mefit37.4,pal.mefit38)

pal.AIC <- tibble(model = c(1:10,12:15,17:38),
                  aic = sapply(pal.mes, AIC)) %>% 
  mutate(deltaAIC = aic-min(aic),
         w = akaike.weights(aic)$weights)

pal.AIC %>% 
  filter(deltaAIC < 4)
#model 30 was the best < 2 AIC

#-------------------------------------------
#####model construction for P. maculata#####
#-------------------------------------------

#filter by only the p. maculata data
pmac.egg <- summ.LILA.egg %>% 
  ungroup() %>% 
  filter(Pomacea_Species == "Maculata") %>% 
  mutate(zmacdens = (pom.dens - mean(pom.dens, na.rm = TRUE))/ sd(pom.dens, na.rm = T))

#create a list of models for comparison

#avedepth only
mac.mefit1 <- glmer.nb(formula = count ~ zavedepth + (1|Cell), data = pmac.egg)

#deltadepth only
mac.mefit2 <- glmer.nb(formula = count ~ zdeltadepth+ (1|Cell), data = pmac.egg,
                       control = glmerControl(optCtrl=list(maxfun=10000000)))  ###not converging must remove from list

#sddepth only
mac.mefit3 <- glmer.nb(formula = count ~ zsddepth+ (1|Cell), data = pmac.egg)

#avedepth and deltadepth
mac.mefit4 <- glmer.nb(formula = count ~ zavedepth + zdeltadepth+ (1|Cell), data = pmac.egg)

#sd depth and deltadepth
mac.mefit5 <- glmer.nb(formula = count ~ zsddepth + zdeltadepth+ (1|Cell), data = pmac.egg,
                       control = glmerControl(optCtrl=list(maxfun=10000000))) ####not converging must remove from list

#avedepth and sddepth
mac.mefit6 <- glmer.nb(formula = count ~ zavedepth + zsddepth+ (1|Cell), data = pmac.egg)

#ave depth sd depth, deltadepth
mac.mefit7 <- glmer.nb(formula = count ~ zavedepth + zsddepth + zdeltadepth+ (1|Cell), data = pmac.egg)

#average depth binomial 
mac.mefit8<- glmer.nb(formula = count ~ zavedepth + I(zavedepth^2)+ (1|Cell), data = pmac.egg)

#season only
mac.mefit9 <- glmer.nb(formula = count ~ Season + (1|Cell), data = pmac.egg)

#depth season interaction
mac.mefit10 <- glmer.nb(formula = count ~ zavedepth*Season + (1|Cell), data = pmac.egg)

#deltadepth, season interaction, this should not be included
mac.mefit11 <- glmer.nb(formula = count ~ zdeltadepth*Season + (1|Cell), data = pmac.egg)

#sd season interaction
mac.mefit12 <- glmer.nb(formula = count ~ zsddepth*Season + (1|Cell), data = pmac.egg)

#season depth binomial interaction
mac.mefit13 <- glmer.nb(formula = count ~ zavedepth*Season + I(zavedepth^2)*Season+ (1|Cell), data = pmac.egg)

#depth density
mac.mefit14 <- glmer.nb(formula = count ~ zmacdens + (1|Cell), data = pmac.egg)

#depth density season
mac.mefit15 <- glmer.nb(formula = count ~ zavedepth + zmacdens + (1|Cell), data = pmac.egg)

#depth density season
mac.mefit16 <- glmer.nb(formula = count ~ zavedepth + zmacdens + Season + (1|Cell), data = pmac.egg)

#density season
mac.mefit17 <- glmer.nb(formula = count ~ zmacdens + Season + (1|Cell), data = pmac.egg)

#season density interaction
mac.mefit18 <- glmer.nb(formula = count ~ Season*zmacdens+ (1|Cell), data = pmac.egg)

#density, binomial depth
mac.mefit19 <- glmer.nb(formula = count ~ zmacdens + zavedepth + I(zavedepth^2)+ (1|Cell), data = pmac.egg)

#season binomial depth
mac.mefit20 <- glmer.nb(formula = count ~ Season + zavedepth + I(zavedepth^2)+ (1|Cell), data = pmac.egg)

#density, season, binomial depth
mac.mefit21 <- glmer.nb(formula = count ~ zmacdens + Season + zavedepth + I(zavedepth^2)+ (1|Cell), data = pmac.egg)

#density, season x binomial depth 
mac.mefit22 <- glmer.nb(formula = count ~ zmacdens + Season + zavedepth + I(zavedepth^2) +
                          Season*zavedepth + Season* I(zavedepth^2)+ (1|Cell), data = pmac.egg)

tt <- getME(mac.mefit22,"theta")   ###check singularity
ll <- getME(mac.mefit22,"lower")  
min(tt[ll==0])                    ###shoun't have a problem with singularity
derivs1 <- mac.mefit22@optinfo$derivs                       ####check gradient
sc_grad1 <- with(derivs1,solve(Hessian,gradient))
max(abs(sc_grad1))                                         ####shouldnt have issues try more interations

mac.mefit22.2 <- update(mac.mefit22,control=glmerControl(optCtrl=list(maxfun=2e4)))  ###fixed it

#temp alone
mac.mefit23 <- glmer.nb(formula = count ~ zavetemp + (1|Cell), data = pmac.egg)

#temp depth
mac.mefit24 <- glmer.nb(formula = count ~ zavetemp + zavedepth + (1|Cell), data = pmac.egg)

#temp x depth interaction
mac.mefit25 <- glmer.nb(formula = count ~ zavetemp + zavedepth +zavetemp*zavedepth+ (1|Cell), data = pmac.egg)

#temp, depth ,dens
mac.mefit26 <- glmer.nb(formula = count ~ zavetemp + zavedepth + zmacdens + (1|Cell), data = pmac.egg)

#temp, dens
mac.mefit27 <- glmer.nb(formula = count ~ zavetemp + zmacdens+ (1|Cell), data = pmac.egg)

#temp, binomial depth
mac.mefit28 <- glmer.nb(formula = count ~ zavetemp + zavedepth + I(zavedepth^2)+ (1|Cell), data = pmac.egg)

#temp, binomial depth, density
mac.mefit29 <- glmer.nb(formula = count ~ zavetemp + zmacdens + zavedepth + I(zavedepth^2)+ (1|Cell), data = pmac.egg)

#density, temp x binomial depth
mac.mefit30 <- glmer.nb(formula = count ~ zmacdens + zavetemp + zavedepth + I(zavedepth^2) +
                          zavetemp*zavedepth + zavetemp* I(zavedepth^2)+ (1|Cell), data = pmac.egg)
tt <- getME(mac.mefit30,"theta")   ###check singularity
ll <- getME(mac.mefit30,"lower")  
min(tt[ll==0])                    ###shoun't have a problem with singularity
derivs1 <- mac.mefit30@optinfo$derivs                       ####check gradient
sc_grad1 <- with(derivs1,solve(Hessian,gradient))
max(abs(sc_grad1))                                         ####shouldnt have issues try more interations

mac.mefit30.2 <- update(mac.mefit30,control=glmerControl(optCtrl=list(maxfun=2e4)))  ###still issues try starting at previous theta
ss <- getME(mac.mefit30,c("theta","fixef"))                                                          ###start with previous theta
mac.mefit30.3 <- update(mac.mefit30,start=ss,control=glmerControl(optCtrl=list(maxfun=2e4)))         ###still has issues
mac.mefit30.4 <- update(mac.mefit30,start=ss,control=glmerControl(optimizer="bobyqa",                 ###different optimzer works
                                                                  optCtrl=list(maxfun=2e5)))  

#temp x binomial depth
mac.mefit31 <- glmer.nb(formula = count ~ zavedepth*zavetemp + I(zavedepth^2)*zavetemp+ (1|Cell), data = pmac.egg)

#photo
mac.mefit32 <- glmer.nb(formula = count ~ zphotoperiod+ (1|Cell), data = pmac.egg)

#dens, photoperiod x binomial depth
mac.mefit33 <- glmer.nb(formula = count ~ zmacdens + zphotoperiod + zavedepth + I(zavedepth^2) +                                 #best model
                          zphotoperiod*zavedepth + zphotoperiod* I(zavedepth^2)+ (1|Cell), data = pmac.egg)

#photoperiod x binomial depth
mac.mefit34 <- glmer.nb(formula = count ~ zphotoperiod + zavedepth + I(zavedepth^2) +                                 
                          zphotoperiod*zavedepth + zphotoperiod* I(zavedepth^2)+ (1|Cell), data = pmac.egg )

#density, photoperiod, binomial depth
mac.mefit35 <-  glmer.nb(formula = count ~zmacdens + zphotoperiod + zavedepth + I(zavedepth^2)+ (1|Cell), data = pmac.egg)

#pca season
mac.mefit36<- glmer.nb(formula = count ~ season.pca + (1|Cell), data = pmac.egg)

#density, season.pca x binomial depth 
mac.mefit37 <- glmer.nb(formula = count ~ zmacdens + season.pca + zavedepth + I(zavedepth^2) +                                 #best model
                          season.pca*zavedepth + season.pca* I(zavedepth^2)+ (1|Cell), data = pmac.egg)

#season.pca x binomial depth
mac.mefit38 <- glmer.nb(formula = count ~ season.pca + zavedepth + I(zavedepth^2) +                                 
                          season.pca*zavedepth + season.pca* I(zavedepth^2)+ (1|Cell), data = pmac.egg )

#compare all the models using AICc values

mac.mes <- c(mac.mefit1,mac.mefit3,mac.mefit4,mac.mefit6,mac.mefit7,mac.mefit8,mac.mefit9,mac.mefit10,
             mac.mefit11,mac.mefit12,mac.mefit13,mac.mefit14,mac.mefit15,mac.mefit16,mac.mefit17,mac.mefit18,mac.mefit19,mac.mefit20,
             mac.mefit21,mac.mefit22.2,mac.mefit23,mac.mefit24,mac.mefit25,mac.mefit26,mac.mefit27,mac.mefit28,mac.mefit29,mac.mefit30.4,
             mac.mefit31, mac.mefit32, mac.mefit33, mac.mefit34, mac.mefit35,mac.mefit36,mac.mefit37,mac.mefit38)

mac.AIC <- tibble(model = c(1,3,4,6,7:38),
                  aic = sapply(mac.mes, AIC)) 
mac.AIC <- mac.AIC %>% 
  mutate(deltaAIC = aic-min(aic),
         w = akaike.weights(mac.AIC$aic)$weights)

mac.AIC %>% 
  filter(deltaAIC < 4)

#model 33 and 35 were the best