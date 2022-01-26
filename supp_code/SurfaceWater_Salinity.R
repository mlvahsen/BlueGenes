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
  
  
# Calculate total average and averages by month for each site
predicted_values <- ggplot_build(salinity_plot)$data[[2]]

# Figure out which points are from which months
salinity_all$date_time_num <- as.numeric(salinity_all$date_time)
salinity_all$month <- month(salinity_all$date_time)

july_range <- range(salinity_all[salinity_all$month == 7,"date_time_num"])
august_range <- range(salinity_all[salinity_all$month == 8,"date_time_num"])
sept_range <- range(salinity_all[salinity_all$month == 9,"date_time_num"])
oct_range <- range(salinity_all[salinity_all$month == 10,"date_time_num"])

# Calculate salinity averages by site and month
predicted_values %>% 
  mutate(month = ifelse(x > july_range[1] & x < july_range[2], "07",
                        ifelse(x > august_range[1] & x < august_range[2], "08",
                               ifelse(x > sept_range[1] & x < sept_range[2], "09",
                                      ifelse(x > oct_range[1] & x < oct_range[2], "10", NA))))) %>% 
  filter(complete.cases(month)) %>% 
  mutate(group = ifelse(group == 1, "freshwater", "gcrew"))-> predicted_values

predicted_values %>% 
  group_by(group, month) %>% 
  summarize(`salinity (ppt)` = round(mean(y),2))

# Calculate overall averages by site for july - october
predicted_values %>% 
  dplyr::select(site = group,
         y = y) %>% 
  group_by(site) %>% 
  summarize(salinity = round(mean(y),2))

# site       salinity
# freshwater     6.79
# gcrew          4.49

# Save plot
ggsave(here("figs", "salinity.png"), salinity_plot, height = 2.8, width = 3.4, units = "in")
  
