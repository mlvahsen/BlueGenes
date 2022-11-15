# Plot water grab nutrient samples through time from Blue Genes 2019 experiment
# as shown in Fig 2

library(tidyverse); library(geomtextpath)
# Read in most nutrient data
nut <- read_csv(here("supp_data", "BlueGeneNutrientResults.csv"))
# Read in particulate N
part_n <- read_csv(here("supp_data", "BG_EA_CN_2019.csv"))
# Format it to get put together with the rest of the nutrients
part_n %>% 
  mutate(site = case_when(site == "fresh" ~ "FW",
                          T ~ "GC")) %>% 
  mutate(sample_day = date) %>% 
  select(site, sample_day, n_mgL = N_mg.L, c_mgL = C_mg.L) -> part_n_sub

# Format rest of nutrient data
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

# Merge two datasets together
right_join(nut_tab, part_n_sub) -> all_nut

salinity_colors <- c("#ff7f00", "#fdbf6f")

all_nut %>% 
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

# Particulate N
all_nut %>% 
  ggplot(aes(x = sample_day, y = n_mgL, color = site)) + 
  geom_line() +
  geom_point(size = 3) +
  xlab("") +
  ylab("particulate N (mg-N/L)") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14)) +
  scale_color_manual(values = salinity_colors) +
  scale_x_date(date_breaks = "1 month", date_labels = c("Aug", "Sep", "Oct"))

# Particulate C
all_nut %>% 
  ggplot(aes(x = sample_day, y = c_mgL, color = site)) + 
  geom_line() +
  geom_point(size = 3) +
  xlab("") +
  ylab("particulate C (mg-N/L)") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14)) +
  scale_color_manual(values = salinity_colors) +
  scale_x_date(date_breaks = "1 month", date_labels = c("Aug", "Sep", "Oct"))

# Particulate N/C
all_nut %>% 
  ggplot(aes(x = sample_day, y = n_mgL/c_mgL, color = site)) + 
  geom_line() +
  geom_point(size = 3) +
  xlab("") +
  ylab("particulate N/C response ratio") +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 14),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14)) +
  scale_color_manual(values = salinity_colors) +
  scale_x_date(date_breaks = "1 month", date_labels = c("Aug", "Sep", "Oct"))

# Biggest differences between sites are are total phosphorus and ammonia so save
# those plots for exp design figure

ggsave(here("figs", "Fig2_phosphorus.png"), phosphorus_plot, height = 3.4, width = 4.2, units = "in")
ggsave(here("figs", "Fig2_ammonia.png"), ammonia_plot, height = 3.4, width = 4.2, units = "in")
