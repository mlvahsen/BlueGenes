# Figure 4. Eco-evo effects across age cohorts for AGB and R:S
library(tidyverse); library(sjPlot)

## Bring in trait data and model ####
source("supp_code/CompileTraitData.R")

bg_full %>% 
  filter(comp == 0 & location != "blackwater" & level < 5 &
           agb_scam > 0) -> traits_nocomp

# Center and scale elevation values
traits_nocomp %>% 
  mutate(elevation_sc = scale(elevation)[,1],
         elevation_sc2 = scale(elevation^2)[,1]) -> traits_nocomp

traits_nocomp %>% 
  mutate(rs = total_bg / agb_scam) -> traits_nocomp

traits_nocomp %>% 
  filter(rs < 6 & pot_no !=165 & pot_no !=176) -> traits_nocomp_rs

agb_evo_model <- readRDS("derived_data/agb_model.rda")
rs_evo_model <- readRDS("derived_data/rs_model.rda")

## Get tidal data to plot flooding on the x-axis ####

# Read in tidal data
tidal_1 <- read_csv(here("supp_data", "tidal_1.csv"))
tidal_2 <- read_csv(here("supp_data", "tidal_2.csv"))
tidal_3 <- read_csv(here("supp_data", "tidal_3.csv"))
tidal_4 <- read_csv(here("supp_data", "tidal_4.csv"))

# Bind all together
tidal_all <- rbind(tidal_1, tidal_2, tidal_3, tidal_4)

# Change column names
mean_tide <- mean(tidal_all$`Verified (m)`)

## AGB ~ age x co2 x salinity x elevation ####
agb_plot <- plot_model(agb_evo_model, terms = c("elevation_sc[all]", "co2", "age", "salinity"), type = "emm")

# Collect data from plot_model object
plot_agb_data <- tibble(elevation_sc = c(agb_plot[[1]]$data$x, agb_plot[[2]]$data$x),
                        co2 = c(agb_plot[[1]]$data$group, agb_plot[[2]]$data$group),
                        age = c(agb_plot[[1]]$data$facet, agb_plot[[2]]$data$facet),
                        salinity = c(agb_plot[[1]]$data$panel, agb_plot[[2]]$data$panel),
                        agb_scam = c(agb_plot[[1]]$data$predicted, agb_plot[[2]]$data$predicted),
                        lower_ci = c(agb_plot[[1]]$data$conf.low, agb_plot[[2]]$data$conf.low),
                        upper_ci = c(agb_plot[[1]]$data$conf.high, agb_plot[[2]]$data$conf.high)) %>% 
  # unscale elevation_sc data for plotting
  mutate(elevation = elevation_sc * sd(traits_nocomp$elevation) + mean(traits_nocomp$elevation),
         flooding = mean_tide - elevation)

# Change the labels in raw data for plots
traits_nocomp %>% 
  mutate(age = case_when(age == "modern" ~ "descendant cohort (2000-2020)",
                         T ~ "ancestral cohort (1900-1950)")) %>% 
  mutate(salinity = case_when(salinity == "fresh" ~ "freshwater site (4 ppt)",
                              T ~ "brackish site (6 ppt)")) %>% 
  mutate(flooding = mean_tide - elevation) -> traits_nocomp_plot

# Create plot
plot_agb_data %>% 
  mutate(age = case_when(age == "modern" ~ "descendant cohort (2000-2020)",
                         T ~ "ancestral cohort (1900-1950)")) %>% 
  mutate(salinity = case_when(salinity == "fresh" ~ "freshwater site (4 ppt)",
                              T ~ "brackish site (6 ppt)")) %>% 
  ggplot(aes(x = flooding, y = agb_scam, color = co2)) +
  geom_point(data = traits_nocomp_plot, shape = 1, alpha = 0.6, stroke = 0.7) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = co2), alpha = 0.2, color = NA) +
  geom_line(size = 1.2) +
  facet_grid(salinity ~ age) +
  scale_color_manual(values = c("#b2df8a","#33a02c"), labels = c(expression(paste("ambient ", CO[2])),
                                                                  expression(paste("elevated ", CO[2])))) +
  scale_fill_manual(values = c("#b2df8a","#33a02c"), labels = c(expression(paste("ambient ", CO[2])),
                                                                 expression(paste("elevated ", CO[2])))) +
  ylab("aboveground biomass (g)") +
  xlab("average inundation (m)") +
  theme_bw(base_size = 14) + theme(legend.position = "top") + labs(color = NULL, fill = NULL) -> agb_plot

## RS ~ age x co2 x elevation ####

rs_plot <- summary(emmeans(rs_evo_model, ~elevation_sc:co2:age,
                           at = list(elevation_sc = seq(min(traits_nocomp$elevation_sc),
                                                        max(traits_nocomp$elevation_sc), length.out = 50),
                                     weights = "proportional"), type = "response"))

# Collect data from plot_model object
plot_rs_data <- tibble(elevation_sc = rs_plot$elevation_sc,
                       co2 = rs_plot$co2,
                       age = rs_plot$age,
                       rs = rs_plot$response,
                       lower_ci = rs_plot$lower.CL,
                       upper_ci = rs_plot$upper.CL) %>% 
  # unscale elevation data
  mutate(elevation = elevation_sc * sd(traits_nocomp_rs$elevation) + mean(traits_nocomp_rs$elevation),
         flooding = mean_tide - elevation)

# Change the labels in raw data for plots
traits_nocomp_rs %>% 
  mutate(age = case_when(age == "modern" ~ "descendant cohort (2000-2020)",
                         T ~ "ancestral cohort (1900-1950)")) %>% 
  mutate(salinity = case_when(salinity == "fresh" ~ "freshwater site (4 ppt)",
                              T ~ "brackish site (6 ppt)")) %>% 
  mutate(flooding = mean_tide - elevation) -> traits_nocomp_rs_plot

plot_rs_data %>% 
  mutate(age = case_when(age == "modern" ~ "descendant cohort (2000-2020)",
                         T ~ "ancestral cohort (1900-1950)")) %>% 
  ggplot(aes(x = flooding, y = rs, color = co2)) +
  geom_point(data = traits_nocomp_rs_plot, shape = 1, alpha = 0.6, stroke = 0.7) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = co2), alpha = 0.2, color = NA) +
  geom_line(size = 1.2) +
  facet_wrap(~age) +
  scale_color_manual(values = c("#b2df8a","#33a02c"), labels = c(expression(paste("ambient ", CO[2])),
                                                                 expression(paste("elevated ", CO[2])))) +
  scale_fill_manual(values = c("#b2df8a","#33a02c"), labels = c(expression(paste("ambient ", CO[2])),
                                                                expression(paste("elevated ", CO[2])))) +
  ylab("root-to-shoot ratio") +
  xlab("average inundation (m)") +
  theme_bw(base_size = 14) +
  theme(legend.position = "top") + 
  labs(color = NULL, fill = NULL) -> rs_plot

## Bring plots together ####

Fig4 <- agb_plot / rs_plot + plot_layout(guides = 'collect', heights = c(3,2)) + plot_annotation(tag_levels = 'a') & theme(legend.position = 'top')
ggsave(here("figs", "Fig4.png"), Fig4, height = 8.5, width = 6, units = "in")
