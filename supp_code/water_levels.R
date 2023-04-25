fw_1 <- read_csv("~/Downloads/processed/BG19_fw_WL_processed_02222019-08092019.csv", skip = 1)
fw_2 <- read_csv("~/Downloads/processed/BG19_fw_WL_processed_08092019-10312019.csv", skip = 1)
gc1 <- read_csv("~/Downloads/processed/BG19_gc_WL_processed_07122019-08162019.csv", skip = 1)
gc2 <- read_csv("~/Downloads/processed/BG19_gc_WL_processed_08162019-10312019.csv", skip = 1)

library(lubridate)

fw_1 %>% 
  mutate(date = mdy_hms(`Date Time, GMT-05:00`),
         month = month(date),
         day = day(date),
         yday = yday(date)) %>% 
  filter(yday > 162) %>% 
  select(month, day, yday, water = `Water Level, meters (LGR S/N: 1265315)`) %>% 
  drop_na() %>% 
  mutate(order = 1:2803) -> fw1_data

fw_2 %>% 
  mutate(date = mdy_hms(`Date Time, GMT-04:00`),
         month = month(date),
         day = day(date),
         yday = yday(date)) %>% 
  filter(yday > 221 & yday < 285) %>% 
  select(month, day, yday, water = `Water Level, meters (LGR S/N: 1265315)`) %>% 
  drop_na() %>% 
  mutate(order = 1:3024+2803) -> fw2_data
  

rbind(fw1_data, fw2_data) %>% 
  ggplot(aes(x = order, y = water)) +
  geom_line()
  
