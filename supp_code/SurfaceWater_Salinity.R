# This code generates the salinity data plot in Fig 2 and calculates overall
# salinity differences between the two common garden sites

library(tidyverse); library(wql);
library(lubridate); library(oce);
library(gt); library(here);library(geomtextpath)

sal_fw1 <- read_csv(here("supp_data", "FW_cond_6.24.2019-08.09.2019_corr.csv"), skip = 1)
sal_fw2 <- read_csv(here("supp_data", "FW_cond_08.09.2019-10.31.2019_corr.csv"), skip = 1)
sal_fw <- rbind(sal_fw1, sal_fw2)

sal_gc1 <- read_csv(here("supp_data", "GCREW_cond_6.28.2019-8.16.2019_corr.csv"), skip = 1)
sal_gc2 <- read_csv(here("supp_data", "GCREW_cond_08.16.2019-10.31.2019_corr.csv"), skip = 1)
names(sal_gc1) <- names(sal_gc2)
sal_gc <- rbind(sal_gc1, sal_gc2)

# Clean up an convert conductivities to mS/cm -- use swSCTp() to calculate
# salinity
sal_fw %>%
  dplyr::select(date_time = `Date Time, GMT-04:00`,
         conductivity = `Full Range, μS/cm (LGR S/N: 20625411, SEN S/N: 20625411)`,
         temperature = `Temp, °C (LGR S/N: 20625411, SEN S/N: 20625411)`) %>%
  mutate(date_time = mdy_hms(date_time)) %>%
  mutate(conductivity = conductivity / 1000) %>%
  mutate(
    salinity = swSCTp(
      conductivity = conductivity,
      temperature = temperature,
      pressure = NULL,
      conductivityUnit = "mS/cm"
    )
  ) %>%
  mutate(site = "fw") -> sal_tab_swSCTp

# Now calculate salinity for GCREW
sal_gc %>%
  dplyr::select(date_time = `Date Time, GMT-05:00`,
         conductivity = `High Range, μS/cm (LGR S/N: 20625481, SEN S/N: 20625481)`,
         temperature = `Temp, °C (LGR S/N: 20625481, SEN S/N: 20625481)`) %>%
  mutate(date_time = mdy_hms(date_time)) %>%
  mutate(conductivity = conductivity / 1000) %>%
  mutate(
    salinity = swSCTp(
      conductivity = conductivity,
      temperature = temperature,
      pressure = NULL,
      conductivityUnit = "mS/cm"
    )
  ) %>%
  mutate(site = "gc") -> sal_gc_tab_swSCTp

# Bind together data freshwater and gcrew datasets
salinity_all <- rbind(sal_tab_swSCTp, sal_gc_tab_swSCTp)

# Plot both sites together
salinity_colors <- c("#ff7f00", "#fdbf6f")

salinity_all %>%
  filter(date_time < "2019-10-12") %>% 
  mutate(site = case_when(site == "fw" ~ "freshwater site",
                          T ~ "brackish site")) %>% 
  ggplot(aes(x = date_time, y = salinity, color = site)) +
  geom_line(alpha = 0.2) +
  geom_labelsmooth(aes(label = site), text_smoothing = 30,
                   method = "loess", formula = y ~ x,
                   size = 4, linewidth = 1, boxlinewidth = 0.3,
                   hjust = 0.8) +
  theme_classic() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.text.y = element_text(size = 12)) +
  xlab("") +
  ylab("salinity (ppt)") +
  scale_color_manual(values = salinity_colors) -> salinity_plot

# Figure out which points are from which months
salinity_all$date_time_num <- as.numeric(salinity_all$date_time)
salinity_all$month <- month(salinity_all$date_time)
salinity_all$day <- day(salinity_all$date_time)

# Calculate averages across sites
salinity_all %>% 
  filter(month %in% 6:9 | month + day < 15) %>% 
  group_by(site) %>% 
  summarize(mean = round(mean(salinity)))
  
# Save plot
ggsave(here("figs", "Fig2_salinity.png"), salinity_plot, height = 2.8, width = 3.4, units = "in")
