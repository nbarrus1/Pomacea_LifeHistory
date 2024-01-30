#the purpose of this code is it to a Principal Components Analysis on the three "seasonal" type variables
#used for modelling egg clutch and add this to the this variables that we will use in the modelling effort.

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

#------------------------------------------------
#####PCA prep#####
#------------------------------------------------

#only take the zscore data for the three "seasonal" variables (i.e., deltaDepth, Temp, Photoperiod)
eggmass.seasonal.matrix <- as.matrix(summ.LILA.egg[,21:23])

#get the number of obs
egg.obs <- nrow(eggmass.seasonal.matrix)
#get the number of variables
egg.var <- ncol(eggmass.seasonal.matrix)

#get the correlation matrix
egg.cormatrix <- cor(eggmass.seasonal.matrix)

#print the correlation matrix
egg.cormatrix

#view the plot for the correlation matrix
pairs(eggmass.seasonal.matrix)

#-------------------------------------------------
####PCA analysis####
#-------------------------------------------------

#run the PCA
fit.pca <- princomp(eggmass.seasonal.matrix, cor = T)

#look at the different metrics for determining how many components to retain
pca_lam <- fit.pca$sdev^2
pca_lam
summary(fit.pca)               #the first component explains 78.8% of the total variance
plot(fit.pca,type = "lines")   #scree plot suggest only one component should be retained 

######3from the metrics only retain the first component!!!!
#now lets look at the loadings on PCA1

fit.pca$loadings[,1]

#the three variables are loading positively on this component 
#(i.e., as PCA1 increases so do the three variables)


#save the results of component 1 into our tibble in preperation for the modelling effort
summ.LILA.egg$season.pca <- fit.pca$scores[,1]
