# Plot water grab nutrient samples through time from Blue Genes 2019 experiment
# as shown in Fig 2

library(tidyverse); library(geomtextpath); library(here); library(patchwork)
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
         # This is actually ammonium
         nh4 = `Ammonia (mg-N/L)`,
         no2_no3 = `Nitrate+Nitrite (mg-N/L)`) -> nut_tab

# Merge two datasets together
right_join(nut_tab, part_n_sub) -> all_nut

salinity_colors <- c("#fdbf6f", "#ff7f00")

# Change labeling of sites
nut_tab %>% 
  mutate(salinity = case_when(site == "FW" ~ "freshwater site (4 ppt)",
                              T ~ "brackish site (6 ppt)")) -> nut_tab

## TOTAL PHOSPHORUS #### 
nut_tab %>% 
  # Covert to micrograms
  mutate(p_mu = p*1000) %>% 
  ggplot(aes(x = sample_day, y = p_mu, color = salinity)) + 
  geom_line() +
  geom_point(size = 3) +
  xlab("") +
  ylab(expression(paste("total phosphorus (", mu, "g P ", L^-1, ")"))) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14, angle=45, hjust = 1),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14)) +
  scale_color_manual(values = rev(salinity_colors))+
  scale_x_date(breaks = unique(nut_tab$sample_day), date_labels = "%b %d") -> phosphorus_plot

## DISSOLVED PHOSPHORUS #### 
nut_tab %>% 
  # Convert to micrograms
  mutate(p_dissolved_mu = p_dissolved*1000) %>% 
  ggplot(aes(x = sample_day, y = p_dissolved_mu, color = salinity)) + 
  geom_line() +
  geom_point(size = 3) +
  xlab("") +
  ylab(expression(paste("dissolved phosphorus (", mu, "g P ", L^-1, ")"))) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14, angle=45, hjust = 1),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14)) +
  scale_color_manual(values = rev(salinity_colors))+
  scale_y_continuous(breaks = seq(0,350,50), limits = c(0,300))+
  scale_x_date(breaks = unique(nut_tab$sample_day), date_labels = "%b %d") -> dissolved_phosphorus_plot

## DISSOLVED NITROGEN #### 
nut_tab %>% 
  # Convert to micrograms
  mutate(din = nh4*1000 + no2_no3*1000) %>% 
  ggplot(aes(x = sample_day, y = din, color = salinity)) + 
  geom_line() +
  geom_point(size = 3) +
  xlab("") +
  ylab(expression(paste("dissolved nitrogen (", mu, "g N ", L^-1, ")"))) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14, angle=45, hjust = 1),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14)) +
  scale_color_manual(values = rev(salinity_colors)) +
  scale_y_continuous(breaks = seq(0,350,50), limits = c(0,300)) +
  scale_x_date(breaks = unique(nut_tab$sample_day), date_labels = "%b %d") -> din_plot

## TOTAL NITROGEN #### 
all_nut %>% 
  mutate(salinity = case_when(site == "FW" ~ "freshwater site (4 ppt)",
                              T ~ "brackish site (6 ppt)")) -> all_nut

# Subset nut_tab data to get dissolved nitrogen on sampling dates in all_nut
nut_tab %>% 
  filter(sample_day %in% all_nut$sample_day) %>% 
  select(sample_day, nh4, no2_no3, salinity) -> dN_5

# Create column for total N
merge(all_nut, dN_5) %>% 
  mutate(totalN = nh4 + no2_no3 + n_mgL) %>% 
  # Convert to micrograms
  mutate(n_mu = totalN*1000) %>% 
  ggplot(aes(x = sample_day, y = n_mu, color = salinity)) + 
  geom_line() +
  geom_point(size = 3) +
  xlab("") +
  ylab(expression(paste("total nitrogen (", mu, "g N ", L^-1, ")"))) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14, angle=45, hjust = 1),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14)) +
  scale_color_manual(values = rev(salinity_colors)) +
  scale_x_date(breaks = unique(all_nut$sample_day), date_labels = "%b %d")-> totalN_plot

## PARTICULATE CARBON #### 
all_nut %>% 
  # Convert to micrograms 
  mutate(c_mu = c_mgL*1000) %>% 
  ggplot(aes(x = sample_day, y = c_mu, color = salinity)) + 
  geom_line() +
  geom_point(size = 3) +
  xlab("") +
  ylab(expression(paste("particulate carbon (", mu, "g C ", L^-1, ")"))) +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14, angle=45, hjust = 1),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14)) +
  scale_color_manual(values = rev(salinity_colors)) +
  scale_x_date(breaks = unique(all_nut$sample_day), date_labels = "%b %d") -> particulateC_plot

## PARTICULATE NITROGEN:CARBON #### 
all_nut %>% 
  ggplot(aes(x = sample_day, y = n_mgL/c_mgL, color = salinity)) + 
  geom_line() +
  geom_point(size = 3) +
  xlab("") +
  ylab("particulate nitrogen:carbon ratio") +
  theme_classic() +
  theme(axis.text.x = element_text(size = 14, angle=45, hjust = 1),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 14)) +
  scale_color_manual(values = rev(salinity_colors)) +
  scale_x_date(breaks = unique(all_nut$sample_day), date_labels = "%b %d") -> particulateNC_plot

## Bring plots together ####

(phosphorus_plot + dissolved_phosphorus_plot + totalN_plot + din_plot +
  particulateC_plot + particulateNC_plot)  +
  plot_layout(nrow = 3, guides = "collect") & theme(legend.position = "top", legend.text = element_text(size = 14),
                                                    legend.title = element_text(size = 14))-> nutrient_plot

ggsave(here("figs", "FigS1_nutrients.png"), nutrient_plot, height = 12.5, width = 8.5, units = "in")
ggsave(here("figs", "Fig2_phosphorus.png"), dissolved_phosphorus_plot, height = 4, width = 4.2, units = "in")
ggsave(here("figs", "Fig2_nitrogen.png"), din_plot, height = 4, width = 4.2, units = "in")

## Create linear models to test which factors are the most different ####

# Total phosphorus
mod_phosphorus <- lm(p ~ sample_day + site, data = nut_tab)
car::Anova(mod_phosphorus)# highly significant effect of site (***)

# Dissolved phosphorus
mod_dissolved_phosphorus <- lm(p_dissolved ~ sample_day + site, data = nut_tab)
car::Anova(mod_dissolved_phosphorus)# no significant effect of site

# Ammonium
mod_din <- lm(din ~ sample_day + site, data = nut_tab %>% mutate(din = no2_no3 + nh4))
car::Anova(mod_din) # significant effect of site (*)

# Particulate N
mod_partN <- lm(n_mgL ~ sample_day + site, data = all_nut)
car::Anova(mod_partN) # significant effect of site (*)

# Particulate C
mod_partC <- lm(c_mgL ~ sample_day + site, data = all_nut)
car::Anova(mod_partC) # marginally significant effect of site (.)

# N/C response ratio
mod_RR <- lm(n_mgL/c_mgL ~ sample_day + site, data = all_nut)
car::Anova(mod_RR) # no significant effect of site
