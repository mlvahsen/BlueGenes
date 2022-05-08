# Figure 4. Age effects

## AGB ~ age x co2 x salinity x elevation ####
agb_plot <- plot_model(agb_evo_model, terms = c("elevation[all]", "co2", "age", "salinity"), type = "emm")

# Collect data from plot_model object
plot_agb_data <- tibble(elevation = c(agb_plot[[1]]$data$x, agb_plot[[2]]$data$x),
                        co2 = c(agb_plot[[1]]$data$group, agb_plot[[2]]$data$group),
                        age = c(agb_plot[[1]]$data$facet, agb_plot[[2]]$data$facet),
                        salinity = c(agb_plot[[1]]$data$panel, agb_plot[[2]]$data$panel),
                        agb_scam = c(agb_plot[[1]]$data$predicted, agb_plot[[2]]$data$predicted),
                        lower_ci = c(agb_plot[[1]]$data$conf.low, agb_plot[[2]]$data$conf.low),
                        upper_ci = c(agb_plot[[1]]$data$conf.high, agb_plot[[2]]$data$conf.high))

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
  theme_bw(base_size = 13) + theme(legend.position = "top") + labs(color = NULL, fill = NULL) -> agb_plot

## RS ~ age x co2 x elevation ####

rs_plot <- summary(emmeans(rs_evo_model, ~elevation:co2:age,
                           at = list(elevation = seq(0.156, 0.544, length.out = 50)), type = "response"))

# Collect data from plot_model object
plot_rs_data <- tibble(elevation = rs_plot$elevation,
                       co2 = rs_plot$co2,
                       age = rs_plot$age,
                       rs = rs_plot$response,
                       lower_ci = rs_plot$lower.CL,
                       upper_ci = rs_plot$upper.CL)

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
  theme_bw(base_size = 13) +
  theme(legend.position = "top") + 
  labs(color = NULL, fill = NULL) -> rs_plot

## Bring plots together ####

Fig4 <- agb_plot / rs_plot + plot_layout(guides = 'collect', heights = c(3,2)) + plot_annotation(tag_levels = 'a') & theme(legend.position = 'top')
 
ggsave(here("figs", "Fig4.png"), Fig4, height = 8, width = 6, units = "in")
