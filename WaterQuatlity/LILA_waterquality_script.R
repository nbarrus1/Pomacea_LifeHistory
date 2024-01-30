rm(list = ls())

###libraries####

library(tidyverse)
library(lubridate)
library(readxl)

###load in data###

waterdata.raw <- read_excel("LILA_waterchemistry_all.xlsx", sheet = 2)

####data clean####

unique(waterdata.raw$Site)

waterdata.clean <- waterdata.raw %>% 
  mutate(cell = if_else(Site == "M1E4" | Site == "M1E5"| Site == "M1ESW"|
                          Site == "M1W8DEEP"| Site == "M1W5" | Site == "M1W9"|
                          Site == "M1Eb"| Site == "M1Wb" |Site == "M1EB" | Site == "M1WB",
                        true = "M1", 
                        false = if_else(Site ==  "M2E4" | Site == "M2E5"| Site == "M2ESW"|
                                          Site == "M2W5" | Site == "M2W9" |
                                          Site == "M2Eb"| Site == "M2Wb" |Site == "M2EB" | Site == "M2WB",
                                        true = "M2",
                                        false = if_else(Site == "M3E4" | Site == "M3E5"| Site == "M3ESW"|
                                                          Site == "M3E8DEEP"| Site == "M3W5" | Site == "M3W9"|
                                                          Site == "M3Eb"| Site == "M3Wb" |Site == "M3EB" | Site == "M3WB",
                                                        true = "M3",
                                                        false = "M4"))),
         year = year(Date)) %>% 
  select(cell, pH, 'Ca2+', TP) %>% 
  gather(key = "measure", value = value, ,pH, 'Ca2+',TP) %>% 
  drop_na()

waterdata.summ.comb <- waterdata.clean %>% 
  group_by(measure) %>% 
  summarise(n = n(),
            mean = round(mean(value, na.rm = T), digits = 2),
            sd = round(sd(value, na.rm = T), digits = 2),
            min = round(min(value, na.rm = T), digits = 2),
            max = round(max(value, na.rm = T), digits = 2)) %>% 
  mutate(cell = "combined") 

waterdata.summ.comb

waterdata.summ.cell <- waterdata.clean %>% 
  group_by(cell, measure) %>% 
  summarise(n = n(),
            mean = round(mean(value, na.rm = T), digits = 2),
            sd = round(sd(value, na.rm = T), digits = 2),
            min = round(min(value, na.rm = T), digits = 2),
            max = round(max(value, na.rm = T), digits = 2))  

waterdata.summ <- waterdata.summ.cell %>% 
  bind_rows(waterdata.summ.comb)

write_csv(waterdata.summ, file = "summarytable.csv")
