#the purpose of this code it to read in data and to calculate the
#environmental variables used for our reproduction (egg clutch laying) analysis
#for P. paludosa and P. maculata.  


rm(list = ls())  #clears global environment

####Important!!! packages that I have below will need to be installed if not installed 
####previously

#------------------------------------------------
#####libraries#####
#------------------------------------------------

library(tidyverse)       #package for languange and syntax that makes code a little more friendly                        #includes a host of other packages like dplyr, ggplot2, etc.
library(readxl)          #package used to read excel formatted data
library(lubridate)       #library used to make working with date easier
library(geosphere)       #package to find day lengths (photoperiod)
library(here)

#-------------------------------------------------
####load in Data#####
#-------------------------------------------------

####egg clutch data

LILA.egg.counts <- read_excel(here("Pomacea/LifeHistory_manuscript/data","EggMassData_v1.12_nb.xlsx"), sheet = 2) %>% 
  mutate(year = year(Date),
         month = months(Date),
         dayofyear = as.numeric(format(Date, "%j"))) 

####transect data
#note) this data is included in the metadeta sheet in "EggMassData_v1.12_nb.xlsx"

transectarea <- tibble(Cell = c("M4","M4","M4","M4","M3","M3","M3","M3","M2","M2","M2","M2","M1","M1","M1","M1"),
                       Transect = c("T1","T2","T3","T4","T1","T2","T3","T4","T1","T2","T3","T4","T1","T2","T3","T4"),
                       Area = c(1408,1268,1240,2156,1368,1264,1188,2136,1340,1276,1280,2100,1320,1292,1216,2124))

####depth data

depth.raw <- read_excel(here("Pomacea/LifeHistory_manuscript/data","LILA_EnvironmentalData_060118-081021.xlsx"), sheet = 2, skip = 6) %>% 
  rename(Cell = Station,
         Date = 'Daily Date',
         Depth.ft = 'Data Value') 

#rename and calculate depths in cm

depth.raw <- depth.raw %>%  
  mutate(Cell = if_else(Cell == "LILA1O",
                        true = "M1",
                        false = if_else(Cell == "LILA2O",
                                        true = "M2",
                                        false = if_else(Cell == "LILA3O",
                                                        true = "M3",
                                                        false = "M4"))),
         depth.cm = (Depth.ft - 13.5)*12*2.54,
         treatment = if_else(Cell == "M1" | Cell == "M3",
                             true = "unconstrained",
                             false = "constrained"))
#####Temperature data

temp.raw <- read_excel(here("Pomacea/LifeHistory_manuscript/data","LILA_EnvironmentalData_060118-081021.xlsx"), sheet = 4,
                       skip = 3) %>% 
  rename(Date = 'Daily Date',
         Temp.c = 'Data Value') 

##density data from Throw Trapping

TTdata <- read_excel(here("Pomacea/LifeHistory_manuscript/data","LILA_Stage_Contrast_TTdata_WS2018-DS2021.xlsx"), sheet = 4)


#--------------------------------------
#####check the data#####
#--------------------------------------

####check data

#egg clutch data

head(LILA.egg.counts)                    #look at first 6 rows of data
tail(LILA.egg.counts)                    #look at last 6 rows of data
table(LILA.egg.counts$Cell)              #check the Cell row for misspelled data lines
table(LILA.egg.counts$Transect)          #check the transect rows for misspelling
table(LILA.egg.counts$Habitat)           #check habitat for misspelling
table(LILA.egg.counts$Pomacea_Species)   #check species for misspelling
table(LILA.egg.counts$year)              #check years
table(LILA.egg.counts$month)             #check months
hist(LILA.egg.counts$Count)              #look at distribution of counts (definitely non-normal)

#transect data

head(transectarea)                    #look at first 6 rows of data
tail(transectarea)                    #look at last 6 rows of data
table(transectarea$Cell)              #check the Cell row for misspelled data lines
table(transectarea$Transect)          #check the transect rows for misspelling
hist(transectarea$Area)              #look at distribution of counts (definitely non-normal)

#depth data

head(depth.raw)                    #look at first 6 rows of data
tail(depth.raw)                    #look at last 6 rows of data
table(depth.raw$Cell)              #check the Cell row for misspelled data lines
table(depth.raw$Depth.ft)          #check the transect rows for misspelling
table(depth.raw$Season)           #check habitat for misspelling
hist(depth.raw$Depth.ft)  
hist(depth.raw$depth.cm)

#temp data

head(temp.raw)                    #look at first 6 rows of data
tail(temp.raw)                    #look at last 6 rows of data
hist(temp.raw$Temp.c) 

#throwtrap data

head(TTdata)                    #look at first 6 rows of data
tail(TTdata)                    #look at last 6 rows of data
hist(TTdata$pommac) 
hist(TTdata$pompal) 
table(TTdata$Macrocosm)              #check the Cell row for misspelled data lines
table(TTdata$Throw_Trap)          #check the transect rows for misspelling
table(TTdata$Location) 
table(TTdata$Session)


#--------------------------------------------------------
####Merge LILA.egg.counts with Transect area####
#--------------------------------------------------------

LILA.egg.counts <- LILA.egg.counts %>%                  
  left_join(transectarea, by = c("Cell", "Transect"))   ## add transect area indexed by Cell and Transect

head(LILA.egg.counts[,5:11])           #check that areas got looking at beginning of data
tail(LILA.egg.counts[,5:11])           #check that areas got looking at end of data     

#------------------------------------------------------------------------------------------
#####Combine Ecotone Transect Counts and Areas, as well as calculate densities####
#------------------------------------------------------------------------------------------

#####Get Counts for Modelling#####

#lines 149 -> create a new data frame for our modelling data 
#lines 150 -> only use transect one and two as these are the deep slough-ridge transects
#lines 151 -> group everything together to combine transects
#lines 152-153 -> summ the two tansects counts and areas together
#lines 154 -> reformat our date from Poisix to Date structure, 
#lines 155 -> convert counts to density (meters squared)
#lines 157-159 -> create week.date for use in merging weekly covariates

summ.LILA.egg <- LILA.egg.counts %>% 
  filter(Transect == "T1" | Transect == "T2") %>% 
  group_by(Date, Pomacea_Species, year, month, Cell, dayofyear) %>% 
  summarise(count = sum(Count),
            area = sum(Area)) %>% 
  mutate(Date = as_date(Date),
         density = count/area,
         week.date = floor_date(Date,
                                unit = "week",
                                week_start  = getOption("lubridate.week.start", 1)))


##check our data

table(summ.LILA.egg$Pomacea_Species)     #looks good we omit all zero counts from M4
table(summ.LILA.egg$year)                #We had differing amounts of surveys each year
table(summ.LILA.egg$month)               #all the monts we have transects
table(summ.LILA.egg$Cell)                #Balanced number of surveys by M1-M3, except M4
#we omited p. paludosa from m4 (extirpated)
table(summ.LILA.egg$dayofyear)           #some of our counts occured on the same day of year
hist(summ.LILA.egg$count)                #again the distribution is non-normal (poission or negbinom)
hist(summ.LILA.egg$density)              #same for desnity

##lets check out our count distribution assumptions poisson has mean and variance that are equal
mean(summ.LILA.egg$density)
var(summ.LILA.egg$density) #mean does not equal variance so not poisson distribuiton 
#needs a dispersion parameter aka negative binomial

#--------------------------------------------------------------------
####calculate covariates for modelling####
#--------------------------------------------------------------------

#####average depths and sd for weekly water levels####

#save all the depth data into a new data
depth.summ <- depth.raw 

#lines 228 -> save our changes
#lines 229 -> Group our depth data by Cell
#lines 230-232 -> create a new variable to have weekly dates starting on monday
#lines 233 -> regroup by Cell and week date 
#lines 234 -> get weekly averages 
#lines 235 -> get weekly standard deviations

depth.summ <- depth.summ %>% 
  group_by(Cell) %>% 
  mutate(week.date = ceiling_date(Date,
                                  unit = "week",
                                  week_start  = getOption("lubridate.week.start", 1))) %>% 
  group_by(Cell,week.date) %>% 
  summarise(ave.depth.cm = mean(depth.cm, na.rm = T),
            sd.depth.cm = sd(depth.cm, na.rm = T))

#check the data
hist(depth.summ$ave.depth.cm)     #looks good
hist(depth.summ$sd.depth.cm)  #interesting

####average weekly temperature####

#lines 260 -> save our changes to temperature summ
#lines 261-263 -> create a week date by monday
#lines 264 -> get weekly averages 
#lines 265 -> get average temperature

temp.summ <- temp.raw%>% 
  mutate(week.date = ceiling_date(Date,
                                  unit = "week",
                                  week_start  = getOption("lubridate.week.start", 1))) %>% 
  group_by(week.date) %>% 
  summarise(ave.temp.c = mean(Temp.c, na.rm = T))

#look at the yearly summaries

temp.summ %>% 
  mutate(year = year(week.date)) %>% 
  group_by(year) %>% 
  summarise(ave = mean(ave.temp.c),
            n = n())

#check data
hist(temp.summ$ave.temp.c)      #interesting our temperate is schewed

#####Get Change in Depths####

#Separate depths by wetland (lines 231-253) in order to use the diff function correctly

M1.depth <- depth.raw %>% 
  filter(Cell == "M1") %>% 
  mutate(year = year(Date),
         month = months(Date),
         dayofyear = as.numeric(format(Date, "%j")))

M2.depth <- depth.raw %>% 
  filter(Cell == "M2") %>% 
  mutate(year = year(Date),
         month = months(Date),
         dayofyear = as.numeric(format(Date, "%j")))

M3.depth <- depth.raw %>% 
  filter(Cell == "M3") %>% 
  mutate(year = year(Date),
         month = months(Date),
         dayofyear = as.numeric(format(Date, "%j")))

M4.depth <- depth.raw %>% 
  filter(Cell == "M4") %>% 
  mutate(year = year(Date),
         month = months(Date),
         dayofyear = as.numeric(format(Date, "%j")))

#create a tibble of our change in depths using the diff function, including dates, and cell info

change <- tibble(depth.change = c(diff(M1.depth$depth.cm, lag = 20),
                                  diff(M2.depth$depth.cm, lag = 20),
                                  diff(M3.depth$depth.cm, lag = 20),
                                  diff(M4.depth$depth.cm, lag = 20)),
                 week.date = c(M1.depth$Date[21:length(M1.depth$depth.cm)],
                               M2.depth$Date[21:length(M1.depth$depth.cm)],
                               M3.depth$Date[21:length(M1.depth$depth.cm)],
                               M4.depth$Date[21:length(M1.depth$depth.cm)]),
                 Cell = c(rep("M1", times = length(M1.depth$Date[21:length(M1.depth$depth.cm)])),
                          rep("M2", times = length(M1.depth$Date[21:length(M1.depth$depth.cm)])),
                          rep("M3", times = length(M1.depth$Date[21:length(M1.depth$depth.cm)])),
                          rep("M4", times = length(M1.depth$Date[21:length(M1.depth$depth.cm)]))))

###check it
hist(change$depth.change)      #looks pretty good


####Get Apple Snail Densities####

#reformat the TTdata to change the data structure from wide to long,
#and rename some of the variables to match the egg clutch format, 
#add calendar year for merging, calculate apple snail densities, remove Paludosa in M4

summ.TTdata <- TTdata %>% 
  filter(Location != "CR") %>% 
  group_by(wateryr, Macrocosm) %>% 
  summarise(n.TT = n(),
            pompal = sum(pompal),
            pommac = sum(pommac)) %>% 
  gather(key = "Pomacea_Species", value = pom.count, 4:5) %>% 
  mutate(pom.dens = pom.count/n.TT,
         year = wateryr,
         Pomacea_Species = if_else(Pomacea_Species == "pompal", true = "Paludosa",
                                   false = "Maculata"),
         year = if_else(year == 2018,
                        true = 2019,
                        false = if_else(year == 2019,
                                        true = 2020,
                                        false = 2021))) %>% 
  rename(Cell = Macrocosm) %>% 
  ungroup() %>% 
  mutate(Cell = if_else(condition = Cell == 1,
                        true = "M1",
                        false = if_else(condition = Cell == 2,
                                        true = "M2",
                                        false = if_else(condition = Cell == 3,
                                                        true = "M3",
                                                        false = "M4"))),
         pom.count = if_else(condition = Cell == "M4" & Pomacea_Species == "Paludosa",
                             true = 100000,
                             false = pom.count)) %>% 
  filter(pom.count < 1000) %>% 
  mutate(pom.dens = if_else(pom.dens > 0,
                            true = pom.dens,
                            false = 0.01)) %>% 
  dplyr::select( Cell, Pomacea_Species, year, pom.dens,n.TT)

#save as a table for our supplementary data

write_csv(summ.TTdata, file = "SuppleTable4.csv")


###calculate photo period###

#with the daylength function obtain photoperiods with latitudes and dates

summ.LILA.egg <- summ.LILA.egg %>% 
  mutate(photoperiod = daylength(lat = 26.4993,
                                 doy = Date))

####Get seasons####

seasons <- depth.raw %>% 
  filter(Date %in% depth.summ$week.date) %>% 
  dplyr::select(Cell, Date, Season) %>% 
  rename(week.date = Date)


#-----------------------------------
####merge the covariates to the LILA counts #####
#-----------------------------------

#merge the depth data

summ.LILA.egg <- summ.LILA.egg %>% 
  left_join(depth.summ, by = c("week.date","Cell"))

#merge the temperature data

summ.LILA.egg <- summ.LILA.egg %>% 
  left_join(temp.summ, by = "week.date")

#merge the changing depth  data

summ.LILA.egg <- summ.LILA.egg %>% 
  left_join(change, by = c("week.date", "Cell"))

#merge the seasons data

summ.LILA.egg <- summ.LILA.egg %>% 
  left_join(seasons, by = c("week.date", "Cell"))

#merge the snail densities

summ.LILA.egg <- summ.LILA.egg %>% 
  left_join(summ.TTdata, by = c("Cell", "Pomacea_Species", "year"))

#-----------------------------------------------------------------------------------------
###### create hydrograph #####
#-----------------------------------------------------------------------------------------

#save our survey period dates for shading in the plot

rect_data  <- tibble(start.date = as_date(c("2019-03-18","2020-02-25","2021-01-25")),
                     end.date = as_date(c("2019-07-22","2020-08-12", "2021-08-02")))

#make the hydrogragh with average weekly depths

depth.raw %>% 
  mutate(week.date = ceiling_date(Date,
                                  unit = "week",
                                  week_start  = getOption("lubridate.week.start", 1))) %>% 
  group_by(week.date, Cell) %>% 
  summarise(depth.cm = mean(depth.cm, na.rm = T)) %>% 
  ggplot(aes(x = week.date, y = depth.cm, color = Cell))+
  theme_classic()+
  geom_rect(data= rect_data, inherit.aes = FALSE,
            aes(xmin=as.POSIXct(start.date), xmax=as.POSIXct(end.date), ymin=-Inf, ymax=+Inf), 
            fill='Grey', alpha=0.3)+
  geom_line(size = 1.25)+
  labs(y = "Water Depth (cm)", x = NULL,
       title = " ")+
  theme(axis.text.x = element_text(angle = 60, vjust = 0.5, size = 18),
        axis.text.y = element_text(size = 12),
        title = element_text(size = 36),
        axis.title = element_text(size = 20),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12))+
  scale_x_datetime(date_breaks = "3 months",
                   date_labels = "%b-%y")+
  scale_color_manual(values = c("Black","#666666","#333333", "#999999"),
                     name = NULL)
