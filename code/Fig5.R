## Stem width (mm) ####

# Get emmeans values for plots
width_EL_plot <- summary(emmeans::emmeans(width_evo_model, ~elevation:location, at = list(elevation = seq(0.156, 0.544, length.out = 50))))
width_SL_plot <- summary(emmeans::emmeans(width_evo_model, ~salinity:location))

tibble(width_scam_mid = width_EL_plot$emmean,
       elevation = width_EL_plot$elevation,
       lower = width_EL_plot$lower.CL,
       upper = width_EL_plot$upper.CL,
       location = rep(c("corn", "kirkpatrick"), each = 50)) %>% 
  ggplot(aes(x = elevation, y = width_scam_mid, color = location)) +
  geom_point(data = traits_nocomp, aes(x = elevation, y = width_scam_mid, color = location), shape = 1,
             stroke = 0.8, alpha = 0.6, size = 0.8) +
  geom_line(size = 1.2) +
  theme_bw(base_size = 13) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = location), alpha = 0.2, linetype = "dashed", color = NA) +
  ylab("mean stem width (mm)") +
  xlab("elevation (m NAVD88)") + 
  scale_color_manual(values = c("#cab2d6", "#6a3d9a")) +
  scale_fill_manual(values = c("#cab2d6", "#6a3d9a")) +
  theme(legend.position = "none")  -> width_provenance

width_SL_plot <- summary(emmeans::emmeans(width_evo_model, ~salinity:location))

tibble(width_scam_mid = width_SL_plot$emmean,
       lower = width_SL_plot$lower.CL,
       upper = width_SL_plot$upper.CL,
       location = rep(c("corn", "kirkpatrick"), each = 2),
       salinity = c("fresh", "salt", "fresh", "salt")) %>% 
  ggplot(aes(x = salinity, y = width_scam_mid, color = location, group = location)) +
  theme_bw(base_size = 13) +
  geom_jitter(data = traits_nocomp, aes(x = salinity, y = width_scam_mid, color = location), shape = 1,
              stroke = 0.8, alpha = 0.6, size = 0.8, height = 0, width = 0.1) +
  geom_line(position=position_dodge(width=0.3), linetype = "dashed", color = "black") +
  geom_errorbar(aes(ymin = lower, ymax = upper), position=position_dodge(width=0.3), width = 0.2, color = "black") +
  geom_point(size = 3, position=position_dodge(width=0.3)) +
  geom_point(size = 3, position=position_dodge(width=0.3), shape = 1, stroke = 0.8, color = "black") +
  ylab("stem width (mm)") +
  scale_color_manual(values = c("#cab2d6", "#6a3d9a")) +
  labs(color = "provenance") +
  theme(legend.position = "none",
        plot.margin = margin(0, 0, 0, 0)) +
  scale_x_discrete(breaks=c("fresh", "salt"),
                   labels=c("freshwater site", "brackish site"))-> width_provenance2



## Root distribution parameter ####
# Get emmeans values for plots
beta_SL_plot <- summary(emmeans::emmeans(beta_evo_model, ~salinity:location, at = list(elevation = seq(0.156, 0.544, length.out = 50))))

tibble(beta = beta_SL_plot$emmean,
       lower = beta_SL_plot$lower.CL,
       upper = beta_SL_plot$upper.CL,
       location = rep(c("corn", "kirkpatrick"), each = 2),
       salinity = c("fresh", "salt", "fresh", "salt")) %>% 
  ggplot(aes(x = salinity, y = beta, color = location, group = location)) +
  theme_bw(base_size = 13) +
  geom_jitter(data = traits_nocomp_rs, aes(x = salinity, y = beta, color = location), shape = 1,
              stroke = 0.8, alpha = 0.6, size = 0.8, height = 0, width = 0.1) +
  geom_line(position=position_dodge(width=0.3), linetype = "dashed", color = "black") +
  geom_errorbar(aes(ymin = lower, ymax = upper), position=position_dodge(width=0.3), width = 0.2, color = "black") +
  geom_point(size = 3, position=position_dodge(width=0.3)) +
  geom_point(size = 3, position=position_dodge(width=0.3), shape = 1, stroke = 0.8, color = "black") +
  ylab("root parameter") +
  xlab("") +
  scale_color_manual(values = c("#cab2d6", "#6a3d9a")) +
  labs(color = "provenance") +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        plot.margin = margin(0, 0, 0, 0)) -> beta_provenance

## Stem height ####

height_SL_plot <- summary(emmeans::emmeans(height_evo_model, ~salinity:location))

tibble(height_scam_tot = height_SL_plot$emmean,
       lower = height_SL_plot$lower.CL,
       upper = height_SL_plot$upper.CL,
       location = rep(c("corn", "kirkpatrick"), each = 2),
       salinity = c("fresh", "salt", "fresh", "salt")) %>% 
  ggplot(aes(x = salinity, y = height_scam_tot, color = location, group = location)) +
  theme_bw(base_size = 13) +
  geom_jitter(data = traits_nocomp_height, aes(x = salinity, y = height_scam_tot, color = location), shape = 1,
              stroke = 0.8, alpha = 0.6, size = 0.8, height = 0, width = 0.1) +
  geom_line(position=position_dodge(width=0.3), linetype = "dashed", color = "black") +
  geom_errorbar(aes(ymin = lower, ymax = upper), position=position_dodge(width=0.3), width = 0.2, color = "black") +
  geom_point(size = 3, position=position_dodge(width=0.3)) +
  geom_point(size = 3, position=position_dodge(width=0.3), shape = 1, stroke = 0.8, color = "black") +
  ylab("stem height (cm)") +
  xlab("") +
  scale_color_manual(values = c("#cab2d6", "#6a3d9a")) +
  labs(color = "provenance") +
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        plot.margin = margin(0, 0, 0, 0)) -> height_provenance

## Bring plots together ####


Fig5 <- width_provenance + (beta_provenance / height_provenance / width_provenance2) +
  plot_layout(widths = c(3,2)) + plot_annotation(tag_levels = 'a')

png("figs/Fig5.png", height = 6.1, width = 6.7, res = 300, units = "in")
Fig5
dev.off()
