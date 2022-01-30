# Summarizing co2 data from sensors. Sensor data reads in as a text file with
# measurements every 5 minutes.

# load libraries
library(tidyverse); library(lubridate); library(here)

## Freshwater site ####

# Read in individual chamber data for FW site
amb <- read.table(here("supp_data", "FW_CO2_amb3.TXT"), sep = ",", stringsAsFactors = F)
elv2 <- read.table(here("supp_data", "FW_CO2_elv2.TXT"), sep = ",", stringsAsFactors = F)
elv5 <- read.table(here("supp_data", "FW_CO2_elv5.TXT"), sep = ",", stringsAsFactors = F)
elv6 <- read.table(here("supp_data", "FW_CO2_elv6.TXT"), sep = ",", stringsAsFactors = F)

column_names <- c("unixtime", "timestamp", "id", "co2", "duty1", "co2ave1", "co2ave5",
                  "dutyave5", "NA", "temp", "pressure", "humidity", "gas", "alt")

colnames(amb) <- column_names
colnames(elv2) <- column_names
colnames(elv5) <- column_names
colnames(elv6) <- column_names

# make a new column for chamber
amb$chamber <- rep("amb", nrow(amb))
elv2$chamber <- rep("elv2", nrow(elv2))
elv5$chamber <- rep("elv5", nrow(elv5))
elv6$chamber <- rep("elv6", nrow(elv6))

# bind together data for one data.frame
full <- rbind(amb, elv2, elv5, elv6)
full$timestamp <- ymd_hms(full$timestamp)

# Summary statistics. 
hours_on <- seq(7,19,1)/24

full %>% 
  mutate(hour = hour(ymd_hms(timestamp)),
         co2ave1 = as.numeric(co2ave1)) %>% 
  # Filter out date to be after 22 June Filter data to be
  # between during the day (0700 - 1900 in program)
  filter(timestamp > "2019/6/22 0:00:00" & hour %in% 7:19) %>% 
  group_by(chamber) %>% 
  summarize(quant_25 = quantile(co2ave1, 0.25, na.rm = T),
            mean = mean(co2ave1, na.rm = T),
            quant_75 = quantile(co2ave1, 0.75, na.rm = T)) -> FW_summary


## Brackish site ####
gcrew <- read.csv(here("supp_data", "phragCR23X_final_storage_1 (1).dat"), header = F)

gcrew %>% filter(V1 == 103) %>% 
  select(V3, V4, V12, V13, V14, V15) %>% 
  # Filter out days
  filter(V3 %in% seq(173,284)) %>% 
  # Filter times that CO2 was on
  filter(V4 %in% seq(630, 2000))-> co2_g

colnames(co2_g) <- c("day","time", "g1", "g2", "g3", "g4")

# Summarize by chamber
co2_g %>% 
  gather(key = chamber, value = co2, g1:g4) %>% 
  group_by(chamber) %>% 
  summarize(lower_25 = quantile(co2, 0.25),
            mean = mean(co2),
            upper_75 = quantile(co2, 0.75)) -> GC_summary

FW_summary
GC_summary
