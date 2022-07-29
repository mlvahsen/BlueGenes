# Figure 3. Age effects

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
  mutate(elevation = elevation_sc * sd(traits_nocomp$elevation) + mean(traits_nocomp$elevation))

# Change the labels in raw data for plots
traits_nocomp %>% 
  mutate(age = case_when(age == "modern" ~ "descendant cohort (2000-2020)",
                         T ~ "ancestral cohort (1900-1950)")) %>% 
  mutate(salinity = case_when(salinity == "fresh" ~ "freshwater site (4ppt)",
                              T ~ "brackish site (6ppt)")) -> traits_nocomp_plot

# Create plot
plot_agb_data %>% 
  mutate(age = case_when(age == "modern" ~ "descendant cohort (2000-2020)",
                         T ~ "ancestral cohort (1900-1950)")) %>% 
  mutate(salinity = case_when(salinity == "fresh" ~ "freshwater site (4ppt)",
                              T ~ "brackish site (6ppt)")) %>% 
  ggplot(aes(x = elevation, y = agb_scam, color = co2)) +
  geom_point(data = traits_nocomp_plot, shape = 1, alpha = 0.6, stroke = 0.7) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = co2), alpha = 0.2, color = NA) +
  geom_line(aes(label = co2), size = 1.2) +
  facet_grid(salinity ~ age) +
  scale_color_manual(values = c("#b2df8a","#33a02c"), labels = c(expression(paste("ambient ", CO[2])),
                                                                  expression(paste("elevated ", CO[2])))) +
  scale_fill_manual(values = c("#b2df8a","#33a02c"), labels = c(expression(paste("ambient ", CO[2])),
                                                                 expression(paste("elevated ", CO[2])))) +
  ylab("aboveground biomass (g)") +
  xlab("elevation (m NAVD88)") +
  theme_bw(base_size = 14) + theme(legend.position = "top") + labs(color = NULL, fill = NULL) -> agb_plot

## RS ~ age x co2 x elevation ####

rs_plot <- summary(emmeans(rs_evo_model, ~elevation_sc:co2:age,
                           at = list(elevation_sc = seq(min(traits_nocomp$elevation_sc),
                                                        max(traits_nocomp$elevation_sc), length.out = 50)), type = "response"))

# Collect data from plot_model object
plot_rs_data <- tibble(elevation_sc = rs_plot$elevation_sc,
                       co2 = rs_plot$co2,
                       age = rs_plot$age,
                       rs = rs_plot$response,
                       lower_ci = rs_plot$lower.CL,
                       upper_ci = rs_plot$upper.CL) %>% 
  # unscale elevation data
  mutate(elevation = elevation_sc * sd(traits_nocomp_rs$elevation) + mean(traits_nocomp_rs$elevation))

# Change the labels in raw data for plots
traits_nocomp_rs %>% 
  mutate(age = case_when(age == "modern" ~ "descendant cohort (2000-2020)",
                         T ~ "ancestral cohort (1900-1950)")) %>% 
  mutate(salinity = case_when(salinity == "fresh" ~ "freshwater site (4ppt)",
                              T ~ "brackish site (6ppt)")) -> traits_nocomp_rs_plot

plot_rs_data %>% 
  mutate(age = case_when(age == "modern" ~ "descendant cohort (2000-2020)",
                         T ~ "ancestral cohort (1900-1950)")) %>% 
  ggplot(aes(x = elevation, y = rs, color = co2)) +
  geom_point(data = traits_nocomp_rs_plot, shape = 1, alpha = 0.6, stroke = 0.7) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = co2), alpha = 0.2, color = NA) +
  geom_line(aes(label = co2), size = 1.2) +
  facet_wrap(~age) +
  scale_color_manual(values = c("#b2df8a","#33a02c"), labels = c(expression(paste("ambient ", CO[2])),
                                                                 expression(paste("elevated ", CO[2])))) +
  scale_fill_manual(values = c("#b2df8a","#33a02c"), labels = c(expression(paste("ambient ", CO[2])),
                                                                expression(paste("elevated ", CO[2])))) +
  ylab("root-to-shoot ratio") +
  xlab("elevation (m NAVD88)") +
  theme_bw(base_size = 14) +
  theme(legend.position = "top") + 
  labs(color = NULL, fill = NULL) -> rs_plot

## Bring plots together ####

Fig3 <- agb_plot / rs_plot + plot_layout(guides = 'collect', heights = c(3,2)) + plot_annotation(tag_levels = 'a') & theme(legend.position = 'top')
 
ggsave(here("figs", "Fig3.png"), Fig3, height = 8.5, width = 6, units = "in")
