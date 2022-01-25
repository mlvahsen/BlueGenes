# Plot water grab nutrient samples through time from Blue Genes 2019 experiment

library(tidyverse)
nut <- read_csv("chp2/data/salinity_nutrients/BlueGeneNutrientResults.csv")
names(nut)
head(nut)

theme_set(theme_classic())

nut %>% 
  mutate(site_day = substr(Sample, nchar(Sample) - 4, nchar(Sample))) %>% 
  mutate(site = substr(str_squish(site_day), 1, 2)) %>% 
  mutate(sample_day = parse_number(site_day)) %>% 
  select(site,
         sample_day,
         p_dissolved = `Total Dissolved Phosphorus (mg-P/L)`,
         p = `Total Phosphorus (mg-P/L)`,
         nh3 = `Ammonia (mg-N/L)`,
         no2_no3 = `Nitrate+Nitrite (mg-N/L)`) -> nut_tab

nut_tab %>% 
  ggplot(aes(x = sample_day, y = p_dissolved, color = site)) + 
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = 1:10) +
  xlab("sample day") +
  ylab("total dissolved phosphorus (mg-P/L)")-> p_dissolved_graph

nut_tab %>% 
  ggplot(aes(x = sample_day, y = p, color = site)) + 
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = 1:10) +
  xlab("sample day") +
  ylab("total phosphorus (mg-P/L)")-> p_graph

nut_tab %>% 
  ggplot(aes(x = sample_day, y = nh3, color = site)) + 
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = 1:10) +
  xlab("sample day") +
  ylab("ammonia (mg-N/L)")-> nh3_graph

nut_tab %>% 
  ggplot(aes(x = sample_day, y = no2_no3, color = site)) + 
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = 1:10) +
  xlab("sample day") +
  ylab("nitrate+nitrite (mg-N/L)")-> no2_no3_graph

png("chp2/figs/nutrients_through_time.png", height = 7, width = 10, res = 300, units = "in")
cowplot::plot_grid(p_dissolved_graph, p_graph,
                   nh3_graph, no2_no3_graph)
dev.off()