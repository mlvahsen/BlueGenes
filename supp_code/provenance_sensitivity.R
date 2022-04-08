## Stem width (mm) ####

# Get emmeans values for plots
width_EL_plot <- summary(emmeans::emmeans(width_step_mod, ~elevation:location, at = list(elevation = seq(0.156, 0.544, length.out = 50))))
width_SL_plot <- summary(emmeans::emmeans(width_step_mod, ~salinity:location))

tibble(width_scam_mid = width_EL_plot$emmean,
       elevation = width_EL_plot$elevation,
       lower = width_EL_plot$lower.CL,
       upper = width_EL_plot$upper.CL,
       location = rep(c("corn", "kirkpatrick"), each = 50)) %>% 
  ggplot(aes(x = elevation, y = width_scam_mid, color = location)) +
  geom_point(data = traits_nocomp, aes(x = elevation, y = width_scam_mid, color = location), shape = 1,
             stroke = 0.8, alpha = 0.6, size = 0.8) +
  geom_line(size = 1.2) +
  theme_bw() +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = location), alpha = 0.2, linetype = "dashed", color = NA) +
  ylab("mean stem width (mm)") +
  xlab("elevation (m NAVD88)") + 
  scale_color_manual(values = c("#cab2d6", "#6a3d9a")) +
  scale_fill_manual(values = c("#cab2d6", "#6a3d9a")) +
  theme(legend.position = "none")  -> width_provenance

width_SL_plot <- summary(emmeans::emmeans(width_step_mod, ~salinity:location))

tibble(width_scam_mid = width_SL_plot$emmean,
       lower = width_SL_plot$lower.CL,
       upper = width_SL_plot$upper.CL,
       location = rep(c("corn", "kirkpatrick"), each = 2),
       salinity = c("fresh", "salt", "fresh", "salt")) %>% 
  ggplot(aes(x = salinity, y = width_scam_mid, color = location, group = location)) +
  theme_bw() +
  geom_jitter(data = traits_nocomp, aes(x = salinity, y = width_scam_mid, color = location), shape = 1,
              stroke = 0.8, alpha = 0.6, size = 0.8, height = 0, width = 0.1) +
  geom_line(position=position_dodge(width=0.3), linetype = "dashed", color = "black") +
  geom_errorbar(aes(ymin = lower, ymax = upper), position=position_dodge(width=0.3), width = 0.2, color = "black") +
  geom_point(size = 3, position=position_dodge(width=0.3)) +
  geom_point(size = 3, position=position_dodge(width=0.3), shape = 1, stroke = 0.8, color = "black") +
  ylab("mean stem width (mm)") +
  scale_color_manual(values = c("#cab2d6", "#6a3d9a")) +
  labs(color = "provenance") +
  theme(legend.position = "none") -> width_provenance2


## Root-to-shoot ####

# Collect means
rs_EL_plot <- summary(emmeans::emmeans(rs_step_mod, ~elevation:location, at = list(elevation = seq(0.156, 0.544, length.out = 50))))

tibble(rs = exp(rs_EL_plot$emmean),
       elevation = rs_EL_plot$elevation,
       lower = exp(rs_EL_plot$lower.CL),
       upper = exp(rs_EL_plot$upper.CL),
       location = rep(c("corn", "kirkpatrick"), each = 50)) %>% 
  ggplot(aes(x = elevation, y = rs, color = location)) +
  geom_point(data = traits_nocomp_rs, aes(x = elevation, y = rs, color = location), shape = 1,
             stroke = 0.8, alpha = 0.6, size = 0.8) +
  geom_line(size = 1.2) +
  theme_bw() +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = location), alpha = 0.2, linetype = "dashed", color = NA) +
  ylab("root-to-shoot ratio") +
  xlab("") + 
  scale_color_manual(values = c("#cab2d6", "#6a3d9a")) +
  scale_fill_manual(values = c("#cab2d6", "#6a3d9a")) +
  theme(legend.position = "none",
        axis.text.x = element_blank()) +
  annotate("text", x = 0.25, y = 4.5, label = "Corn Island", fontface = 2, color = "#cab2d6") +
  annotate("text", x = 0.29, y = 4.15, label = "Kirkpatrick Marsh", fontface = 2, color = "#6a3d9a") +
  geom_point(aes(x = 0.165, y = 4.5), color = "#cab2d6", size = 2) +
  geom_point(aes(x = 0.165, y = 4.15), color = "#6a3d9a", size = 2)-> rs_provenance

## Root distribution parameter ####
# Get emmeans values for plots
beta_EL_plot <- summary(emmeans::emmeans(beta_mod1, ~elevation:location, at = list(elevation = seq(0.156, 0.544, length.out = 50))))

tibble(beta = beta_EL_plot$emmean,
       elevation = beta_EL_plot$elevation,
       lower = beta_EL_plot$lower.CL,
       upper = beta_EL_plot$upper.CL,
       location = rep(c("corn", "kirkpatrick"), each = 50)) %>% 
  ggplot(aes(x = elevation, y = beta, color = location)) +
  geom_point(data = traits_nocomp, aes(x = elevation, y = beta, color = location), shape = 1,
             stroke = 0.8, alpha = 0.6, size = 0.8) +
  geom_line(size = 1.2) +
  theme_bw() +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = location), alpha = 0.2, linetype = "dashed", color = NA) +
  ylab("root distribution parameter") +
  xlab("") + 
  scale_color_manual(values = c("#cab2d6", "#6a3d9a")) +
  scale_fill_manual(values = c("#cab2d6", "#6a3d9a")) +
  theme(legend.position = "none",
        axis.text.x = element_blank()) -> beta_provenance

## Stem height ####

height_SL_plot <- summary(emmeans::emmeans(height_step_mod, ~salinity:location))

tibble(height_scam_tot = height_SL_plot$emmean,
       lower = height_SL_plot$lower.CL,
       upper = height_SL_plot$upper.CL,
       location = rep(c("corn", "kirkpatrick"), each = 2),
       salinity = c("fresh", "salt", "fresh", "salt")) %>% 
  ggplot(aes(x = salinity, y = height_scam_tot, color = location, group = location)) +
  theme_bw() +
  geom_jitter(data = traits_nocomp_noOut, aes(x = salinity, y = height_scam_tot, color = location), shape = 1,
              stroke = 0.8, alpha = 0.6, size = 0.8, height = 0, width = 0.1) +
  geom_line(position=position_dodge(width=0.3), linetype = "dashed", color = "black") +
  geom_errorbar(aes(ymin = lower, ymax = upper), position=position_dodge(width=0.3), width = 0.2, color = "black") +
  geom_point(size = 3, position=position_dodge(width=0.3)) +
  geom_point(size = 3, position=position_dodge(width=0.3), shape = 1, stroke = 0.8, color = "black") +
  ylab("mean stem height (cm)") +
  xlab("") +
  scale_color_manual(values = c("#cab2d6", "#6a3d9a")) +
  labs(color = "provenance") +
  theme(legend.position = "none",
        axis.text.x = element_blank()) -> height_provenance

## Trait sensitivities ####
tibble(trait = c("rs", "beta", "width", " width", "height"),
       environment = c("sea-level rise", "sea-level rise", "sea-level rise", "salinity", "salinity"),
       change_slope_mean = c(abs(fixef(rs_step_mod)[8]),
                             abs(fixef(beta_mod1)[6]),
                             fixef(width_step_mod)[8],
                             fixef(width_step_mod)[7],
                             fixef(height_step_mod)[5]),
       change_slope_lower = c(abs(confint(rs_step_mod)[10,2]),
                              abs(confint(beta_mod1)[10,2]),
                              confint(width_step_mod)[10,1],
                              confint(width_step_mod)[9,1],
                              confint(height_step_mod)[7,1]),
       change_slope_upper = c(abs(confint(rs_step_mod)[10,1]),
                              abs(confint(beta_mod1)[10,1]),
                              confint(width_step_mod)[10,2],
                              confint(width_step_mod)[9,2],
                              confint(height_step_mod)[7,2])) -> change_in_slopes

change_in_slopes %>% 
  mutate(trait = factor(trait, levels = c(" width", "height", "width", "beta", "rs"))) %>% 
  ggplot(aes(x = trait, y = change_slope_mean, color = environment)) +
  geom_point() +
  geom_errorbar(aes(ymin = change_slope_lower, ymax = change_slope_upper), width=0) +
  coord_flip() +
  theme(legend.position = "none") +
  ylab(expression(paste("| ",Delta,slope[(Corn - Kirkpatrick)], " |"))) +
  scale_color_manual(values = c("orange", "dodgerblue")) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", color = "gray47") +
  annotate("text", x = 1.5, y = 2, label = "salinity", color = "orange", fontface = 3) +
  annotate("text", x = 4, y = 2, label = "sea-level rise", color = "dodgerblue", fontface = 3) -> sensitivity_plot

## Bring plots together ####

(rs_provenance / beta_provenance / width_provenance + plot_annotation(tag_levels = c("a"))) & theme(plot.margin = margin(t = 0,r = 5,b = 5,l = 5)) -> first_plot
(height_provenance / width_provenance2/ sensitivity_plot + plot_layout(tag_level = "new") + plot_annotation(tag_levels = list(c("d", "e", "f")))) & theme(plot.margin = margin(t = 0,r = 5,b = 5,l = 5)) -> second_plot 

png("figs/FigX_provenance_sensitivity.png", height = 7.78, width = 6.64, units = "in", res = 300)
cowplot::plot_grid(first_plot, second_plot)
dev.off()


