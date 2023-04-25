# Figure 5 - Provenance by eco interactions

## Get tidal data to plot flooding on the x-axis for stem width ####

# Read in tidal data
tidal_1 <- read_csv(here("supp_data", "tidal_1.csv"))
tidal_2 <- read_csv(here("supp_data", "tidal_2.csv"))
tidal_3 <- read_csv(here("supp_data", "tidal_3.csv"))
tidal_4 <- read_csv(here("supp_data", "tidal_4.csv"))

# Bind all together
tidal_all <- rbind(tidal_1, tidal_2, tidal_3, tidal_4)

# Change column names
mean_tide <- mean(tidal_all$`Verified (m)`)

## Stem width (mm) ####

# Get emmeans values for plots
width_EL_plot <- summary(emmeans::emmeans(width_evo_model, ~elevation_sc:location,
                                          weights = "proportional",
                                          at = list(elevation_sc = seq(min(traits_nocomp$elevation_sc), max(traits_nocomp$elevation_sc), length.out = 50))))
width_SL_plot <- summary(emmeans::emmeans(width_evo_model, ~salinity:location,
                                          weights = "proportional"))

# Create column in raw data for average inundation
traits_nocomp %>% 
  mutate(flooding = mean_tide - elevation,
         provenance = location) -> traits_nocomp

tibble(width_scam_mid = width_EL_plot$emmean,
       elevation_sc = width_EL_plot$elevation_sc,
       lower = width_EL_plot$lower.CL,
       upper = width_EL_plot$upper.CL,
       provenance = rep(c("corn", "kirkpatrick"), each = 50)) %>% 
  mutate(elevation = elevation_sc * sd(traits_nocomp$elevation) + mean(traits_nocomp$elevation),
         flooding = mean_tide - elevation) %>% 
  ggplot(aes(x = flooding, y = width_scam_mid, color = provenance)) +
  geom_point(data = traits_nocomp, aes(x = flooding, y = width_scam_mid, color = provenance), shape = 1,
             stroke = 0.8, alpha = 0.6, size = 1.2) +
  geom_line(size = 1.2) +
  theme_bw(base_size = 14) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = provenance), alpha = 0.2, linetype = "dashed", color = NA) +
  ylab("mean stem width (mm)") +
  xlab("average inundation (m)") + 
  scale_color_manual(values = c("#cab2d6", "#6a3d9a")) +
  scale_fill_manual(values = c("#cab2d6", "#6a3d9a")) +
  theme(legend.position = "none")-> width_provenance

width_SL_plot <- summary(emmeans::emmeans(width_evo_model, ~salinity:location, weights = "proportional"))

tibble(width_scam_mid = width_SL_plot$emmean,
       lower = width_SL_plot$lower.CL,
       upper = width_SL_plot$upper.CL,
       location = rep(c("corn", "kirkpatrick"), each = 2),
       salinity = c("fresh", "salt", "fresh", "salt")) %>% 
  ggplot(aes(x = salinity, y = width_scam_mid, color = location, group = location)) +
  theme_bw(base_size = 14) +
  geom_jitter(data = traits_nocomp, aes(x = salinity, y = width_scam_mid, color = location), shape = 1,
              stroke = 0.8, alpha = 0.6, size = 1.2, height = 0, width = 0.1,show.legend = F) +
  geom_line(position=position_dodge(width=0.3), linetype = "dashed", color = "black") +
  geom_errorbar(aes(ymin = lower, ymax = upper), position=position_dodge(width=0.3), width = 0.2, color = "black") +
  geom_point(size = 3, position=position_dodge(width=0.3), show.legend = F) +
  geom_point(size = 3, position=position_dodge(width=0.3), shape = 1, stroke = 0.8, color = "black") +
  ylab("mean stem width (mm)") +
  scale_color_manual(values = c("#cab2d6", "#6a3d9a")) +
  labs(color = "provenance") +
  theme(legend.position = "none",
        plot.margin = margin(0, 0, 0, 0)) +
  scale_x_discrete(breaks=c("fresh", "salt"),
                   labels=c("freshwater site", "brackish site"))-> width_provenance2



## Root distribution parameter ####
# Get emmeans values for plots
beta_SL_plot <- summary(emmeans::emmeans(beta_evo_model, ~salinity:location,
                                         at = list(elevation = seq(0.156, 0.544, length.out = 50)),
                                         weights = "proportional"))

tibble(beta = beta_SL_plot$emmean,
       lower = beta_SL_plot$lower.CL,
       upper = beta_SL_plot$upper.CL,
       location = rep(c("corn", "kirkpatrick"), each = 2),
       salinity = c("fresh", "salt", "fresh", "salt")) %>% 
  ggplot(aes(x = salinity, y = beta, color = location, group = location)) +
  theme_bw(base_size = 14) +
  geom_jitter(data = traits_nocomp_rs, aes(x = salinity, y = beta, color = location), shape = 1,
              stroke = 0.8, alpha = 0.6, size = 1.2, height = 0, width = 0.1, show.legend = F) +
  geom_line(position=position_dodge(width=0.3), linetype = "dashed", color = "black") +
  geom_errorbar(aes(ymin = lower, ymax = upper), position=position_dodge(width=0.3), width = 0.2, color = "black") +
  geom_point(size = 3, position=position_dodge(width=0.3), show.legend = F) +
  geom_point(size = 3, position=position_dodge(width=0.3), shape = 1, stroke = 0.8, color = "black") +
  ylab("root parameter") +
  scale_color_manual(values = c("#cab2d6", "#6a3d9a")) +
  labs(color = "provenance") +
  theme(legend.position = "none",
        plot.margin = margin(0, 0, 0, 0)) +
  scale_x_discrete(breaks=c("fresh", "salt"),
                   labels=c("freshwater site", "brackish site"))-> beta_provenance

## Stem height ####

height_SL_plot <- summary(emmeans::emmeans(height_evo_model, ~salinity:location,
                                           weights = "proportional"))

tibble(height_scam_tot = height_SL_plot$emmean,
       lower = height_SL_plot$lower.CL,
       upper = height_SL_plot$upper.CL,
       location = rep(c("corn", "kirkpatrick"), each = 2),
       salinity = c("fresh", "salt", "fresh", "salt")) %>% 
  ggplot(aes(x = salinity, y = height_scam_tot, color = location, group = location)) +
  theme_bw(base_size = 14) +
  geom_jitter(data = traits_nocomp_height, aes(x = salinity, y = height_scam_tot, color = location), shape = 1,
              stroke = 0.8, alpha = 0.6, size = 1.2, height = 0, width = 0.1, show.legend = F) +
  geom_line(position=position_dodge(width=0.3), linetype = "dashed", color = "black") +
  geom_errorbar(aes(ymin = lower, ymax = upper), position=position_dodge(width=0.3), width = 0.2, color = "black") +
  geom_point(size = 3, position=position_dodge(width=0.3), show.legend = F) +
  geom_point(size = 3, position=position_dodge(width=0.3), shape = 1, stroke = 0.8, color = "black") +
  ylab("mean stem height (cm)") +
  scale_color_manual(values = c("#cab2d6", "#6a3d9a")) +
  labs(color = "provenance") +
  theme(legend.position = "none",
        plot.margin = margin(0, 0, 0, 0)) +
  scale_x_discrete(breaks=c("fresh", "salt"),
                   labels=c("freshwater site", "brackish site"))-> height_provenance

## Bring plots together ####

Fig5 <- guide_area() + (width_provenance + width_provenance2 +
                          height_provenance + beta_provenance) + 
  plot_layout(guides = "collect",
              nrow = 2, heights = c(0.5,10)) +
  plot_annotation(tag_levels = 'a') & theme(legend.position = "top") 
  

png("figs/Fig5.png", height = 6, width = 7, res = 300, units = "in")
Fig5
dev.off()
