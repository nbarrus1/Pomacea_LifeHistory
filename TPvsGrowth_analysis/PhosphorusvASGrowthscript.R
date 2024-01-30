rm(list = ls())

#---------------------------------------------------------------
#####libraries####
#---------------------------------------------------------------

library(tidyverse)
library(readxl)
library(cowplot)
library(gridExtra)
library(lmerTest)

#--------------------------------------------------------------
####read in and datawrangling data####
#-------------------------------------------------------------

####GrowthData####

growthdata <- read_excel("FieldGrowthExperiments-2019-2020.xlsx")

growth.summ <- growthdata %>% 
  group_by(year,cell,treatment,cage,species) %>% 
  summarise(fl_mean = mean(length_sl, na.rm = T),
            fm_mean = mean(mass, na.rm = T)) %>% 
  mutate(cage = as.character(cage))

####obtain initial masses 

#create regression to get dry weights
#P.paludosa length weight regression
ppal.lengthweight <- function(SL) {
  ln.tot.drymass <- 2.64*log(SL) - 1.84
}

#P.maculata length weight regression

pmac.lengthweight <- function(SL) {
  ln.tot.drymass <- 2.70*log(SL) - 2.21
}

####get initial lengths
initialdata <- tibble(year = c(2019,2019,2020,2020),
                      species = c("P.paludosa","P.maculata","P.paludosa","P.maculata"),
                      il_mean = c(3.56,3.61,4.40,4.60)) %>% 
  mutate(im_mean = pmac.lengthweight(SL = il_mean))
  
growth.summ <- growth.summ %>% 
  left_join(initialdata, by = c("year", "species"))

####Nutrient Data####

nutrientdata <- read_excel("RawTP_2019&2020.xlsx", sheet = 5)

nutrientdata <- nutrientdata %>% 
  select(Year, Cell, Cage, 'TP_ug/g') %>% 
  rename(cell = Cell,
         cage = Cage,
         TP = 'TP_ug/g',
         year = Year)

head(nutrientdata)

growth.summ <- growth.summ %>% 
  left_join(nutrientdata, by = c("year", "cell", "cage")) 

growth.summ <-growth.summ%>% 
  mutate(SGR_mass = (log(fm_mean)-log(im_mean))/35) %>% 
  mutate(Growth_mass = (fm_mean - im_mean)/35,
         SGR_length = (log(fl_mean)-log(il_mean))/35,
         Growth_length = (fl_mean - il_mean)/35) 

TP_mean <- mean(growth.summ$TP, na.rm = T)
TP_sd <- sd(growth.summ$TP, na.rm = T)

growth.summ <- growth.summ %>% 
  mutate(z_TP = (TP-TP_mean)/TP_sd)

##survival data

survival_data <- read_excel("survival.xlsx") %>% 
  mutate(cage = as.character(cage))

growth.summ <- growth.summ %>% 
  mutate(cage = as.character(cage)) %>% 
  left_join(survival_data, by = c("year", "cell", "treatment", "cage", "species"))

head(growth.summ[10:16])
tail(growth.summ[10:16])

#------------------------------------------------------
####Plot Growth Rates vs TP####
#------------------------------------------------------
#growth with outlier

growth.summ %>% 
  filter(species == "P.paludosa" | species == "P.maculata") %>% 
  ggplot(aes(x = TP, y = Growth_mass, color = species))+
  geom_point()+
  #geom_smooth(aes(fill = species), method = "lm")+
  theme_classic()+
  scale_color_manual(values = c("black", "#666666"))+
  scale_fill_manual(values = c("#666666", "#999999"))+
  labs(x = "Total Phosphorus (ug/g)",
       y = "Growth (mg/day)")

#growth without outlier

absolute_mass <- growth.summ %>% 
  filter(species == "P.paludosa" | species == "P.maculata") %>% 
  filter(TP<1000) %>% 
  ggplot(aes(x = TP, y = Growth_mass, color = species, shape = species))+
  geom_smooth(aes(fill = species), method = "lm", show.legend = F)+
  geom_point(size = 3, show.legend = F)+
  theme_classic()+
  scale_color_manual(values = c("black", "#666666"))+
  scale_fill_manual(values = c("#666666", "#999999"))+
  labs(x = "Total Phosphorus (\U00B5g/g)",
       y = "Absolute Mass",
       title =  "A.")+
  theme(title = element_text(size = 24),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 16))

absolute_mass
#SGR with outlier

growth.summ %>% 
  filter(species == "P.paludosa" | species == "P.maculata") %>% 
  #filter(TP<1000) %>% 
  ggplot(aes(x = TP, y = SGR_mass, color = species))+
  geom_point(aes(shape = species))+
  #geom_smooth(aes(fill = species), method = "lm")+
  theme_classic()+
  scale_color_manual(values = c("black", "#666666"))+
  scale_fill_manual(values = c("#666666", "#999999"))+
  labs(x = "Total Phosphorus (ug/g)",
       y = "SGR (daily porportional growth)")

#SGR without outlier

SGR_mass <- growth.summ %>% 
  filter(species == "P.paludosa" | species == "P.maculata") %>% 
  filter(TP<1000) %>% 
  ggplot(aes(x = TP, y = SGR_mass, color = species))+
  geom_point(aes(shape = species), size = 3, show.legend = F)+
  geom_smooth(aes(fill = species), method = "lm",show.legend = F)+
  theme_classic()+
  scale_color_manual(values = c("black", "#666666"))+
  scale_fill_manual(values = c("#666666", "#999999"))+
  labs(x = "Total Phosphorus (\U00B5g/g)",
       y = "SGR Mass ",
       title =  "A.")+
  scale_y_continuous(breaks = c(0.08,0.11,0.14,0.17,0.20))+
  theme(title = element_text(size = 24),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 16))

SGR_mass

#length with outlier

growth.summ %>% 
  filter(species == "P.paludosa" | species == "P.maculata") %>% 
  #filter(TP<1000) %>% 
  ggplot(aes(x = TP, y = SGR_length, color = species))+
  geom_point(aes(shape = species))+
  #geom_smooth(aes(fill = species), method = "lm")+
  theme_classic()+
  scale_color_manual(values = c("black", "#666666"))+
  scale_fill_manual(values = c("#666666", "#999999"))+
  labs(x = "Total Phosphorus (ug/g)",
       y = "SGR (daily porportional growth)")


#lengths without outlier

SGR_length <- growth.summ %>% 
  filter(species == "P.paludosa" | species == "P.maculata") %>% 
  filter(TP<1000) %>% 
  ggplot(aes(x = TP, y = SGR_length, color = species))+
  geom_point(aes(shape = species), size = 3, show.legend = F)+
  geom_smooth(aes(fill = species), method = "lm", show.legend = F)+
  theme_classic()+
  scale_color_manual(values = c("black", "#666666"))+
  scale_fill_manual(values = c("#666666", "#999999"))+
  labs(x = "Total Phosphorus (\U00B5g/g)",
       y = "SGR Length",
       title = "B.")+
  theme(title = element_text(size = 24),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 16))

SGR_length

abs_length <- growth.summ %>% 
  filter(species == "P.paludosa" | species == "P.maculata") %>% 
  filter(TP<1000) %>% 
  ggplot(aes(x = TP, y = Growth_length, color = species))+
  geom_point(aes(shape = species), size = 3, show.legend = F)+
  geom_smooth(aes(fill = species), method = "lm", show.legend = F)+
  theme_classic()+
  scale_color_manual(values = c("black", "#666666"))+
  scale_fill_manual(values = c("#666666", "#999999"))+
  labs(x = "Total Phosphorus (\U00B5g/g)",
       y = "Absolute Length",
       title = "C.")+
  theme(title = element_text(size = 24),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 16))

abs_length

##---------------------------------------------------------
#####grwoth analysis####
##---------------------------------------------------------
c <- growth.summ %>% 
  drop_na(TP) 

cor(c$year,c$TP)

fit <- lm(SGR_mass~z_TP + species + z_TP:species, data = growth.summ[growth.summ$TP < 1000,]) 
fit1 <- lm(SGR_mass~species + as.character(year) + z_TP + as.character(year):z_TP +
             as.character(year):species + species:z_TP + species:z_TP:as.character(year),
           data = growth.summ[growth.summ$TP < 1000,])
summary(fit1)

fit1 <- lm(SGR_mass~species + as.character(year) + z_TP + as.character(year):z_TP +
             as.character(year):species + species:z_TP,
           data = growth.summ[growth.summ$TP < 1000,])
summary(fit1)
fit1 <- lm(SGR_mass~species + as.character(year) + z_TP + as.character(year):z_TP +
             species:z_TP,
           data = growth.summ[growth.summ$TP < 1000,])
summary(fit1)
fit1 <- lm(SGR_mass~species + as.character(year) + z_TP + as.character(year):z_TP +
             species:z_TP,
           data = growth.summ[growth.summ$TP < 1000,])
summary(fit1)
fitme <- lmer(formula = SGR_mass ~ species + z_TP +species*z_TP,
              data = growth.summ[growth.summ$TP < 1000,]) 


co <- c %>% 
  filter(TP < 1000) %>% 
  select(SGR_mass, species, z_TP, year)
fitme <- lmer(formula = SGR_mass ~ species + z_TP + species:z_TP + (1| year),
              data = co) 

summary(fitme)

growth.summ[growth.summ$TP < 1000,] %>% 
  drop_na(z_TP) %>% 
  ggplot(aes(y = SGR_mass, x = z_TP, col = as.character(year)))+
  geom_point()+
  facet_wrap(~species)+
  geom_smooth(method = "lm")


summary(fit)
summary(fit1)

plot(fit)

###----------------------------------------------------------
######Survival vs TP plot####
###----------------------------------------------------------
survivalvTP <- growth.summ %>% 
  filter(species == "P.paludosa" | species == "P.maculata") %>% 
  filter(TP<1000) %>% 
  ggplot(aes(x = TP, y = Survival, color = species))+
  geom_point(aes(shape = species),show.legend = F, size = 3)+
  scale_color_manual(values = c("black", "#666666"))+
  theme_classic()+
  scale_y_continuous(limits = c(0,1))+
  labs(x = "Total Phosphorus (\U00B5g/g)",
       y = "Survival ",
       title = "B.")+
  theme(title = element_text(size = 24),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 16))

survivalvAFDM <- growth.summ %>% 
  filter(species == "P.paludosa" | species == "P.maculata") %>% 
  ggplot(aes(x = AFDM, y = Survival, color = species))+
  geom_point(aes(shape = species),show.legend = F)+
  scale_color_manual(values = c("black", "#666666"))+
  theme_classic()+
  scale_y_continuous(limits = c(0,1))+
  labs(x = "AFDM (mg/cm^2)",
       y = NULL,
       title = "E.")

survivalvperAFDM <- growth.summ %>% 
  filter(species == "P.paludosa" | species == "P.maculata") %>% 
  filter(percent_AFDM>20) %>% 
  ggplot(aes(x = percent_AFDM, y = Survival, color = species))+
  geom_point(aes(shape = species),show.legend = F)+
  scale_color_manual(values = c("black", "#666666"))+
  theme_classic()+
  scale_y_continuous(limits = c(0,1))+
  labs(x = "% Carbon",
       y = NULL,
       title = "F.")

legend <- growth.summ %>% 
  filter(TP<1000) %>% 
  ggplot(aes(x = TP, y = SGR_mass, color = species,shape = species))+
  geom_point(size = 3)+
  geom_smooth(aes(fill = species), method = "lm", show.legend = F)+
  theme_classic()+
  scale_color_manual(values = c("black", "#666666"))+
  scale_fill_manual(values = c("#333333", "#999999"))+
  labs(x = "Total Phosphorus (ug/g)",
       y = NULL)+
  theme(legend.text = element_text(size = 12),
        legend.title = element_blank())

legend_a <- cowplot::get_legend(
  legend+
    theme(legend.text = element_text(face = "italic"))
)


library(gridExtra)

Fig5 <- grid.arrange(SGR_mass,survivalvTP,legend_a,
          ncol = 7, nrow = 1,
          layout_matrix = rbind(c(1,1,1,2,2,2,3)))

ggsave("C:/Users/Nathan Barrus/Documents/FAU/Masters Thesis/AppleSnail_LifeHistory_manuscript/Figures&Tables/Figure6_nutvgrowthsurv.pdf",
       plot = Fig5, device = "pdf", width = 10, height = 4, units = "in",
       limitsize = FALSE)


summ <- growth.summ %>% 
  filter(TP<1000) %>% 
  filter(percent_AFDM > 20) %>% 
  group_by(species) %>% 
  summarise(ncage = n(),
         ave.TP = mean(TP, na.rm = T),
         min.TP = min(TP, na.rm = T),
         max.TP = max(TP, na.rm = T),
         ave.AFDM = mean(AFDM, na.rm = T),
         min.AFDM = min(AFDM, na.rm = T),
         max.AFDM = max(AFDM, na.rm = T),
         ave.perAFDM = mean(percent_AFDM, na.rm = T),
         min.perAFDM = min(percent_AFDM, na.rm = T),
         max.perAFDM = max(percent_AFDM, na.rm = T),
         ave.SGR.m = mean(SGR_mass, na.rm = T),
         min.SGR.m = min(SGR_mass, na.rm = T),
         max.SGR.m = max(SGR_mass, na.rm = T),
         sd.SGR.m = sd(SGR_mass, na.rm = T),
         ave.survival = mean(Survival, na.rm = T),
         min.survival = min(Survival, na.rm = T),
         max.survival = max(Survival, na.rm = T),
         sd.survival = sd(Survival, na.rm = T),
         ave.SGR.l = mean(SGR_length, na.rm = T),
         min.SGR.l = min(SGR_length, na.rm = T),
         max.SGR.l = max(SGR_length, na.rm = T),
         sd.SGR.l = sd(SGR_length, na.rm = T))
summ

summ <- growth.summ %>% 
  filter(TP<1000) %>% 
  filter(percent_AFDM > 20) %>% 
  group_by(species) %>%
  summarise(ncages = n(),
            ave.SGR.l = mean(SGR_length, na.rm = T),
            min.SGR.l = min(SGR_length, na.rm = T),
            max.SGR.l = max(SGR_length, na.rm = T),
            sd.SGR.l = sd(SGR_length, na.rm = T),
            ave.SGR.m = mean(SGR_mass, na.rm = T),
            min.SGR.m = min(SGR_mass, na.rm = T),
            max.SGR.m = max(SGR_mass, na.rm = T),
            sd.SGR.m = sd(SGR_mass, na.rm = T),
            ave.m = mean(Growth_mass, na.rm = T),
            min.m = min(Growth_mass, na.rm = T),
            max.m = max(Growth_mass, na.rm = T),
            sd.m = sd(Growth_mass, na.rm = T),
            ave.l = mean(Growth_length, na.rm = T),
            min.l = min(Growth_length, na.rm = T),
            max.l = max(Growth_length, na.rm = T),
            sd.l = sd(Growth_length, na.rm = T))

####length plots supplementary figure####

#SGR with treatment
TPvSGR <- growth.summ %>% 
  filter(species == "P.paludosa" | species == "P.maculata") %>% 
  filter(TP<1000) %>% 
  ggplot(aes(x = TP, y = SGR_length, color = species))+
  geom_point(aes(shape = species), show.legend = F)+
  geom_smooth(aes(fill = species), method = "lm",show.legend = F, linetype = "dashed")+
  theme_classic()+
  scale_color_manual(values = c("black", "#666666"))+
  scale_fill_manual(values = c("#333333", "#999999"))+
  labs(x = "Total Phosphorus (ug/g)",
       y = "SGR (daily propotional growth)",
       title = "A.")+
  scale_y_continuous(limits = c(0.025,0.07),
                     breaks = c(0.03,0.04,0.05,0.06,0.07))

fit <- lm(SGR_length~TP + species + TP:species, data = growth.summ[growth.summ$TP < 1000,]) #nonsignificant interaction so pool
summary(fit)

plot(fit)


####AFDM data ####

AFDMvSGR <- growth.summ %>% 
  filter(species == "P.paludosa" | species == "P.maculata") %>% 
  drop_na(AFDM) %>% 
  ggplot(aes(x = AFDM, y = SGR_length, color = species))+
  geom_point(aes(shape = species),show.legend = F)+
  geom_line(aes(y = AFDM.predict$fit), linetype = "dashed", size = 1,
            show.legend = F)+
  geom_ribbon(aes(ymin = AFDM.predict$lwr, ymax = AFDM.predict$upr,
                  color = NULL, fill = species), alpha = 0.4,
              show.legend = F)+
  theme_classic()+
  scale_color_manual(values = c("black", "#666666"))+
  scale_fill_manual(values = c("#333333", "#999999"))+
  labs(x = "AFDM (mg/cm^2)",
       y = NULL,
       title = "B.")+
  scale_y_continuous(limits = c(0.025,0.07),
                     breaks = c(0.03,0.04,0.05,0.06,0.07))

fit <- lm(SGR_length ~ AFDM + species + AFDM:species, data = growth.summ)
summary(fit)

AFDM.predict <- growth.summ %>% 
  drop_na(AFDM)

AFDM_fit <- lm(SGR_length ~ AFDM + species, data = AFDM.predict)
summary(AFDM_fit)

AFDM.predict <- data.frame(predict(AFDM_fit, interval = "confidence"))

####percent AFDM #####
perAFDMvSGR <- growth.summ %>% 
  filter(species == "P.paludosa" | species == "P.maculata") %>% 
  filter(percent_AFDM > 20) %>% 
  drop_na(percent_AFDM) %>% 
  ggplot(aes(x = percent_AFDM, y = SGR_length, color = species))+
  geom_point(aes(shape = species),show.legend = F)+
  geom_line(aes(y = percAFDM_predict$fit), linetype = "dashed", size = 1, show.legend = F)+
  geom_ribbon(aes(ymin = percAFDM_predict$lwr, ymax = percAFDM_predict$upr,
                  color = NULL, fill = species), alpha = 0.4,
              show.legend = F)+
  theme_classic()+
  scale_color_manual(values = c("black", "#666666"))+
  scale_fill_manual(values = c("#333333", "#999999"))+
  labs(x = "% Carbon",
       y = NULL,
       title = "C.")+
  scale_y_continuous(limits = c(0.025,0.07),
                     breaks = c(0.03,0.04,0.05,0.06,0.07))



fit <- lm(SGR_length ~ percent_AFDM + species + percent_AFDM:species, data = growth.summ)
summary(fit)

percAFDM_fit <- lm(SGR_length ~ percent_AFDM + species, data = growth.summ[growth.summ$percent_AFDM > 20,])
summary(percAFDM_fit)
plot(fit)

percAFDM_predict <- data.frame(predict(percAFDM_fit, interval = "confidence"))


legend <- growth.summ %>% 
  filter(TP<1000) %>% 
  ggplot(aes(x = TP, y = SGR_mass, color = species,shape = species))+
  geom_point()+
  geom_smooth(aes(fill = species), method = "lm", linetype = "dashed")+
  theme_classic()+
  scale_color_manual(values = c("black", "#666666"))+
  scale_fill_manual(values = c("#333333", "#999999"))+
  labs(x = "Total Phosphorus (ug/g)",
       y = NULL)

legend_a <- cowplot::get_legend(
  legend+
    theme(legend.text = element_text(face = "italic"))
)

SuppFig1 <- grid.arrange(absolute_mass,SGR_length,abs_length,legend_a, 
                     ncol = 10, nrow = 1,
                     layout_matrix = rbind(c(1,1,1,2,2,2,3,3,3,4)))

ggsave("C:/Users/Nathan Barrus/Documents/FAU/Masters Thesis/AppleSnail_LifeHistory_manuscript/Figures&Tables/SuppFig1_nutvgrowthlength.pdf",
       plot = SuppFig1, device = "pdf", width = 14, height = 5, units = "in",
       limitsize = FALSE)
