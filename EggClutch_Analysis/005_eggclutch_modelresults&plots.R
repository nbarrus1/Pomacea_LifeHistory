#the prupose of this code is to explor the results of the best model

#-----------------------------
#####load in libraries######
#-----------------------------

library(plot3D)          #packaged to create 3D regression plane
library(MASS)            #model fitting package contains
library(lme4)            #model fitting package for mixed effect models GLMM models 
library(lmerTest)        #package that gets model summary information more readable
library(qpcR)            #            
library(lemon)           #package for creating more readable plot pannels
library(cowplot)         #package that puts multiple plots together
library(ggpubr)          #adds important ggplot features neccesary for publishing


#--------------------------------------------------------------------------------------------
#######Best Model Exploration#########
#--------------------------------------------------------------------------------------------


###p. paludosa####
pal.mefit30 <- glmer.nb(formula = count ~ zpaldens + zavetemp + zavedepth + I(zavedepth^2) +                                 #best model
                          zavetemp*zavedepth + zavetemp* I(zavedepth^2)+ (1|Cell), data = ppal.egg)

summary(pal.mefit30)
summary(pal.mefit30.3)

pal.meNULL <- glmer.nb(formula = count~ 1 + (1|Cell), data = ppal.egg)

pal.pseudoR2 <- MuMIn::r.squaredGLMM(pal.mefit30.3, pal.meNULL)  ###trigamma is recommended
pal.pseudoR2 


ppal.egg <- ppal.egg %>% 
  mutate(predict.count = exp(predict(pal.mefit30.3)))

#make a 3d version of the interaction to see if that helps interpretation



grid.lines = 26
x.pred <- seq(min(ppal.egg$zavedepth), max(ppal.egg$zavedepth), length.out = grid.lines)
y.pred <- seq(min(ppal.egg$zavetemp), max(ppal.egg$zavetemp), length.out = grid.lines)
xy <- expand.grid( zavedepth = x.pred,  zavetemp = y.pred)
xy$Cell <- rep("M2", times = length(xy$zavedepth))
xy$zpaldens <- rep(0, times = length(xy$zavedepth))
xy$count <- exp(predict(pal.mefit30.3, newdata = xy, type = "link", allow.new.levels = TRUE))
table(is.na(xy))
z.pred <- matrix(exp(predict(pal.mefit30.3, newdata = xy, type = "link")), 
                 nrow = grid.lines, ncol = grid.lines)

(-0.6018598*sd(summ.LILA.egg$ave.depth.cm))+mean(summ.LILA.egg$ave.depth.cm)

(-2.2182418*sd(summ.LILA.egg$ave.temp.c))+mean(summ.LILA.egg$ave.temp.c)

(x.pred[19]*sd(summ.LILA.egg$ave.depth.cm))+mean(summ.LILA.egg$ave.depth.cm)
(y.pred[26]*sd(summ.LILA.egg$ave.temp.c))+mean(summ.LILA.egg$ave.temp.c)

sum(z.pred[4:14,1:10])/sum(z.pred)
(x.pred[4:14]*sd(summ.LILA.egg$ave.depth.cm))+mean(summ.LILA.egg$ave.depth.cm)
(y.pred[1:10]*sd(summ.LILA.egg$ave.temp.c))+mean(summ.LILA.egg$ave.temp.c)
# scatter plot with regression plane
par(mfrow = c(1, 1))
par(mai = c(0.4,0.4,0.4,0.4),
    cex.axis = 0.75)


depthtempinter3d <- scatter3D(xy$zavedepth, xy$zavetemp, xy$count+100, pch = NULL, cex = .0001, 
                              theta = 160, phi = 20, ticktype = "detailed",col ="black",
                              NAcol = "blue",colvar = xy$count,
                              zlim = c(0,80), bty = "b2", 
                              colkey = list(plot = FALSE), clim = c(0,80),
                              xlab = "Depth (cm)", ylab = "Temp", zlab = "Count",
                              surf = list(x = x.pred, y = y.pred, z = z.pred, 
                                          facets = NA))
points3D(x = ppal.egg$zavedepth,y = ppal.egg$zavetemp, z = ppal.egg$count,
         type = "h", add = TRUE, pch = 16, cex = 1.4,
         col = ramp.col(c("yellowgreen","darkgreen", "midnight blue")))

ggsave("depthtempinter3d.png", plot = depthtempinter3d, device = "png", width = 6.67, height = 7,
       units = "in")

par(mfrow = c(3, 3))
for (alph in c(0.25, 0.75))
  image2D(volcano, alpha = alph,
          main = paste("jet.col, alpha = ", alph))
image2D(volcano, main = "jet.col")
image2D(volcano, col = jet2.col(100), main = "jet2.col")
image2D(volcano, col = gg.col(100), main = "gg.col")
image2D(volcano, col = gg2.col(100), main = "gg2.col")
image2D(volcano, col = rainbow(100), main = "rainbow")
image2D(volcano, col = terrain.colors(100), main = "terrain.colors")
image2D(volcano, col = ramp.col(c("blue", "yellow", "green", "red")),
        main = "ramp.col")

mean(ppal.egg$ave.depth.cm)
sd(ppal.egg$ave.depth.cm)
mean(ppal.egg$ave.temp.c)
sd(ppal.egg$ave.temp.c)


####predicted versus actuals

ppal.egg %>%
  ggplot(aes(x = dayofyear, y = count))+
  facet_grid(Cell~as.character(year))+
  theme_classic()+
  geom_point()+
  geom_line(aes(x = dayofyear, y = count), color = "black")+
  geom_point(aes(x =dayofyear, y = predict.count), shape = 17, color = "#666666")+
  geom_line(aes(x=dayofyear, y = predict.count), color = "#666666", linetype = "dashed")+
  scale_x_continuous(labels = function(x) format(as.Date(as.character(x), "%j"), "%d-%b"))

#####best model for P. maculata#####

mac.mefit35 <-  glmer.nb(formula = count ~zmacdens + zphotoperiod + zavedepth + I(zavedepth^2)+ (1|Cell), data = pmac.egg)
summary(mac.mefit35)

mac.meNULL <-glmer.nb(formula = count ~ 1 + (1|Cell), data = pmac.egg)

mac.pseudoR2<-MuMIn::r.squaredGLMM(mac.mefit35, mac.meNULL) ###trigamma is recommended
mac.pseudoR2

####depth effect on P. maculata

depthinterdata <- data.frame(
  zavedepth = seq(from = min(pmac.egg$zavedepth), to = max(pmac.egg$zavedepth),
                  length.out = 1000),
  Cell = rep("M4", times = 1000),
  zmacdens = rep(0, times = 1000),
  zphotoperiod = rep (0, times = 1000))


depthinterdata <- cbind(depthinterdata, predict(mac.mefit35, depthinterdata, type = "link",
                                                allow.new.levels = TRUE)) 
depthinterdata$lncount.predict <- depthinterdata$`predict(mac.mefit35, depthinterdata, type = "link", allow.new.levels = TRUE)`
depthinterdata$count.predict <- exp(depthinterdata$lncount.predict)

depthinterdata.low<- data.frame(
  zavedepth.low = seq(from = min(pmac.egg$zavedepth), to = max(pmac.egg$zavedepth),
                      length.out = 1000),
  Cell.low = rep("M4", times = 1000),
  zmacdens.low = rep(-1.6338, times = 1000),
  zphotoperiod.low = rep (0, times = 1000))


depthinterdata.low <- cbind(depthinterdata, predict(mac.mefit35, depthinterdata, type = "link",
                                                    allow.new.levels = TRUE)) 
depthinterdata.low$lncount.predict.low <- depthinterdata$`predict(mac.mefit35, depthinterdata, type = "link", allow.new.levels = TRUE)`
depthinterdata.low$count.predict.low <- exp(depthinterdata$lncount.predict)

depthinterdata %>% 
  ggplot(aes(x = zavedepth, y = count.predict)) +
  geom_line(size = 1.5)+
  geom_line(aes(x = zavedepth.low, y = count.predict.low),data = depthinterdata.low, color = "red")+
  geom_point(data = pmac.egg, aes(x = zavedepth, y = count))+
  theme_classic()+
  labs()+
  scale_x_continuous(breaks = c(-2, -1.5,-1,-0.5,0,0.5,1,1.5,2))+
  scale_y_continuous(limits = c(0,150),
                     breaks = c(0,30,60,90,120,50))

pmac.egg <- pmac.egg %>% 
  mutate(predict.count = exp(predict(mac.mefit35)))

sum(depthinterdata$count.predict)

sum(depthinterdata$count.predict[depthinterdata$zavedepth > -0.1 & depthinterdata$zavedepth < 1.1])/
  sum(depthinterdata$count.predict)

sort((depthinterdata$zavedepth[depthinterdata$zavedepth > -0.1 & depthinterdata$zavedepth < 1.1]*
        sd(summ.LILA.egg$ave.depth.cm, na.rm = T))+mean(summ.LILA.egg$ave.depth.cm, na.rm = T))

(depthinterdata$zavedepth[depthinterdata$count.predict == max(depthinterdata$count.predict)]*
    sd(summ.LILA.egg$ave.depth.cm, na.rm = T))+mean(summ.LILA.egg$ave.depth.cm, na.rm = T)
####predicted versus actuals



pmac.egg %>%
  ggplot(aes(x = dayofyear, y = count))+
  facet_grid(Cell~as.character(year))+
  theme_classic()+
  geom_point()+
  geom_line(aes(x = dayofyear, y = count), color = "black")+
  geom_point(aes(x =dayofyear, y = predict.count), shape = 17, color = "#666666")+
  geom_line(aes(x=dayofyear, y = predict.count), color = "#666666", linetype = "dashed")+
  scale_x_continuous(labels = function(x) format(as.Date(as.character(x), "%j"), "%d-%b"))


pmac.egg %>%
  ggplot(aes(x = dayofyear, y = zavetemp))+
  facet_grid(Cell~as.character(year))+
  theme_classic()+
  geom_point()+
  geom_line(aes(x = dayofyear, y = zavetemp), color = "black")+
  geom_point(aes(x =dayofyear, y = zavedepth), shape = 17, color = "#666666")+
  geom_line(aes(x=dayofyear, y = zavedepth), color = "#666666")+
  scale_x_continuous(labels = function(x) format(as.Date(as.character(x), "%j"), "%d-%b"))+
  labs(y = "Z-score")

###leave one out cross validation####

#work flow

###set up data
ppal.egg <- ppal.egg %>% 
  mutate(predict.count = 0,
         obs = 1:length(ppal.egg$Pomacea_Species),
         upp = 0,
         low = 0)


##try with first row of data

train <- ppal.egg %>% 
  filter(obs == 1)
model <- ppal.egg %>% 
  filter(obs != 1)

pal.mefit30 <- glmer.nb(formula = count ~ zpaldens + zavetemp + zavedepth + I(zavedepth^2) +                          
                          zavetemp*zavedepth + zavetemp* I(zavedepth^2)+ (1|Cell), data = model)  ###convergence issues needs fixed
ss <- getME(pal.mefit30,c("theta","fixef"))                                                          ###start with previous theta
pal.mefit30.2 <- update(pal.mefit30,start=ss,control=glmerControl(optCtrl=list(maxfun=2e4)))         ###still has issues
pal.mefit30.3<- update(pal.mefit30,start=ss,control=glmerControl(optimizer="bobyqa",                 ###different optimzer works
                                                                 optCtrl=list(maxfun=2e5))) 
ppal.egg$predict.count[51] <- exp(predict(pal.mefit30.3, train, type = "link",
                                          allow.new.levels = TRUE))

b <- bootMer(pal.mefit30.3, nsim = 100, function(x) predict(x, newdata = train))
ppal.egg$upp[1] = exp(quantile(b$t, probs = 0.975))
ppal.egg$low[1] = exp(quantile(b$t, probs = 0.025))

ppal.egg[,23:27]  ###it worked

####now create loop
for (i in 1:length(ppal.egg$predict.count)) {
  train <- ppal.egg %>% 
    filter(obs == i)
  model <- ppal.egg %>% 
    filter(obs != i)
  
  pal.mefit30 <- glmer.nb(formula = count ~ zpaldens + zavetemp + zavedepth + I(zavedepth^2) +                          
                            zavetemp*zavedepth + zavetemp* I(zavedepth^2)+ (1|Cell), data = model)  ###convergence issues needs fixed
  ss <- getME(pal.mefit30,c("theta","fixef"))                                                          ###start with previous theta
  pal.mefit30.2 <- update(pal.mefit30,start=ss,control=glmerControl(optCtrl=list(maxfun=2e4)))         ###still has issues
  pal.mefit30.3<- update(pal.mefit30,start=ss,control=glmerControl(optimizer="bobyqa",                 ###different optimzer works
                                                                   optCtrl=list(maxfun=2e5))) 
  ppal.egg$predict.count[i] <- exp(predict(pal.mefit30.3, train, type = "link",
                                           allow.new.levels = TRUE))
  
  b <- bootMer(pal.mefit30.3, nsim = 100, function(x) predict(x, newdata = train))
  ppal.egg$upp[i] = exp(quantile(b$t, probs = 0.975))
  ppal.egg$low[i] = exp(quantile(b$t, probs = 0.025))
}

#create loop for p. maculata

pmac.egg <- pmac.egg %>% 
  mutate(predict.count = 0,
         upp = 0,
         low = 0,
         obs = 1:length(pmac.egg$Pomacea_Species))

for (i in 1:length(pmac.egg$predict.count)) {
  train <- pmac.egg %>% 
    filter(obs == i)
  model <- pmac.egg %>% 
    filter(obs != i)
  mac.mefit35 <-  glmer.nb(formula = count ~zmacdens + zphotoperiod + zavedepth + I(zavedepth^2)+ (1|Cell), data = pmac.egg) 
  pmac.egg$predict.count[i] <- exp(predict(mac.mefit35, train, type = "link",
                                           allow.new.levels = TRUE))
  
  b <- bootMer(mac.mefit35, nsim = 100, function(x) predict(x, newdata = train))
  pmac.egg$upp[i] = exp(quantile(b$t, probs = 0.975))
  pmac.egg$low[i] = exp(quantile(b$t, probs = 0.025))
}


###ppal.egg plot with 95% prediction intervals####
ppal.egg %>%
  ggplot(aes(x = dayofyear, y = count))+
  facet_rep_grid(Cell~as.character(year))+
  theme_classic()+
  geom_pointrange(aes(x =dayofyear, y = predict.count, ymin = low, ymax = upp,
                      color = "predict.count", shape = "predict.count"),
                  show.legend = F)+
  geom_point(aes(color = "count", shape = "count"),
             show.legend = T)+
  geom_line(aes(x = dayofyear, y = count,
                color = "count"), show.legend = F)+
  scale_x_continuous(labels = function(x) format(as.Date(as.character(x), "%j"), "%b"))+
  coord_cartesian(ylim = c(0,105))+
  theme(axis.text.x = element_text(vjust = 0.5),
        axis.title = element_text(size = 12, hjust = 0.5),
        strip.background = element_rect(fill = "white", linetype = 0),
        strip.text = element_text(size = 12, face = "bold", color = "black"))+
  labs(x = NULL, y = "Total Egg Clutches in Transects")+
  scale_y_continuous(breaks = c(0,20,40,60,80,100))+
  scale_color_manual(name = NULL,
                     values = c( "predict.count" = "#666666", "count" = "black"),
                     labels = c("Observed", "Predicted"))+
  scale_shape_manual(name = NULL,
                     values = c("predict.count" = 17,"count" = 16),
                     labels = c("Observed", "Predicted"))

ggsave("C:/Users/Nathan Barrus/Documents/FAU/Masters Thesis/AppleSnail_LifeHistory_manuscript/Figures&Tables/Figure3_Paludosa_avp.pdf",
       plot = last_plot(), device = "pdf", units = "in", width = 7.5, height = 6)

####predicted versus actuals p. maculata

pmac.egg %>%
  ggplot(aes(x = dayofyear, y = count))+
  facet_rep_grid(Cell~as.character(year))+
  theme_classic()+
  geom_pointrange(aes(x =dayofyear, y = predict.count, ymin = low, ymax = upp,
                      color = "predict.count", shape = "predict.count"), show.legend = F)+
  geom_point(aes(y = count, color = "count", shape = "count"),size = 2, show.legend = T)+
  geom_line(aes(x = dayofyear, y = count, color = "count"),show.legend = F)+
  scale_x_continuous(labels = function(x) format(as.Date(as.character(x), "%j"), "%b"))+
  labs(x = NULL, y = "Total Egg Masses in Transects" ,
       color = "legend")+
  theme(axis.text.x = element_text(vjust = 0.5),
        axis.title = element_text(size = 12, hjust = 0.5),
        strip.background = element_rect(fill = "white", linetype = 0),
        strip.text = element_text(size = 12, face = "bold", color = "black"))+
  scale_y_continuous(breaks = c(0,30,60,90,120,150))+
  scale_color_manual(name = NULL,
                     values = c( "predict.count" = "#666666", "count" = "black"),
                     labels = c("Observed", "Predicted"))+
  scale_shape_manual(name = NULL,
                     values = c("predict.count" = 17,"count" = 16),
                     labels = c("Observed", "Predicted"))+
  coord_cartesian(ylim = c(0,150))



ggsave("C:/Users/Nathan Barrus/Documents/FAU/Masters Thesis/AppleSnail_LifeHistory_manuscript/Figures&Tables/Figure4_Maculata_avp.pdf",
       plot = last_plot(), device = "pdf", units = "in", width = 10, height = 8)


####model diagnostics####
library(DHARMa)


ppal.egg <- ppal.egg %>% 
  mutate(predict.count = exp(predict(pal.mefit30)),
         RQR.residual = qres.binom(pal.mefit30))
qresiduals(pal.mefit30)


pal.mefit30.3 <- glmer.nb(formula = count ~ zpaldens + zavetemp + zavedepth + I(zavedepth^2) +                                 #best model
                            zavetemp*zavedepth + zavetemp* I(zavedepth^2)+ (1|Cell), data = ppal.egg)

####Get Randomized Quantile Residuals####
testDispersion(pal.mefit30.3)
testOutliers(pal.mefit30.3)
simulationOutput <- simulateResiduals(fittedModel = pal.mefit30.3, plot = F)
simulationOutput
residuals(simulationOutput)

ppal.egg <- ppal.egg %>% 
  mutate(RQR = residuals(simulationOutput))

plot(simulationOutput)

testDispersion(mac.mefit35)
testOutliers(mac.mefit35)
simulationOutput <- simulateResiduals(fittedModel = mac.mefit35, plot = F)
simulationOutput
residuals(simulationOutput)

pmac.egg <- pmac.egg %>% 
  mutate(RQR = residuals(simulationOutput))
