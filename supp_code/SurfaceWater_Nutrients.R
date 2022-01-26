# Plot water grab nutrient samples through time from Blue Genes 2019 experiment

library(tidyverse); library(geomtextpath)
nut <- read_csv(here("supp_data", "BlueGeneNutrientResults.csv"))

nut %>% 
  mutate(site_day = substr(Sample, nchar(Sample) - 4, nchar(Sample))) %>% 
  mutate(site = substr(str_squish(site_day), 1, 2)) %>% 
  mutate(sample_day = rep(as.Date(c("2019-07-22", "2019-07-30", "2019-08-07",
                                "2019-08-13", "2019-08-20", "2019-08-26",
                                "2019-09-10", "2019-09-19", "2019-09-25",
                                "2019-10-02")), each = 2)) %>% 
  dplyr::select(site,
         sample_day,
         p_dissolved = `Total Dissolved Phosphorus (mg-P/L)`,
         p = `Total Phosphorus (mg-P/L)`,
         nh3 = `Ammonia (mg-N/L)`,
         no2_no3 = `Nitrate+Nitrite (mg-N/L)`) -> nut_tab

# nut_tab %>% 
#   ggplot(aes(x = sample_day, y = p_dissolved, color = site)) + 
#   geom_line() +
#   geom_point() +
#   xlab("sample day") +
#   ylab("total dissolved phosphorus (mg-P/L)")-> p_dissolved_graph

salinity_colors <- c("#ff7f00", "#fdbf6f")

nut_tab %>% 
  ggplot(aes(x = sample_day, y = p, color = site)) + 
  geom_line() +
  geom_point(size = 3) +
  xlab("") +
  ylab("total phosphorus (mg-P/L)") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14)) +
  scale_color_manual(values = salinity_colors)+
  scale_x_date(date_breaks = "1 month", date_labels = c("Aug", "Sep", "Oct")) -> phosphorus_plot

nut_tab %>% 
  ggplot(aes(x = sample_day, y = nh3, color = site)) + 
  geom_line() +
  geom_point(size = 3) +
  xlab("") +
  ylab("ammonia (mg-N/L)") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14)) +
  scale_color_manual(values = salinity_colors) +
  scale_x_date(date_breaks = "1 month", date_labels = c("Aug", "Sep", "Oct")) -> ammonia_plot

# nut_tab %>% 
#   ggplot(aes(x = sample_day, y = no2_no3, color = site)) + 
#   geom_line() +
#   geom_point() +
#   scale_x_continuous(breaks = 1:10) +
#   xlab("sample day") +
#   ylab("nitrate+nitrite (mg-N/L)")-> no2_no3_graph

# Biggest differences between sites are are total phosphorus and ammonia so save
# those plots for exp design figure

ggsave(here("figs", "phosphorus.png"), phosphorus_plot, height = 3.4, width = 4.2, units = "in")
ggsave(here("figs", "ammonia.png"), ammonia_plot, height = 3.4, width = 4.2, units = "in")
