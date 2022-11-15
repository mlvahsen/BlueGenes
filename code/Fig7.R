# Figure 7 - Competition by elevation influence on aboveground biomass

# Create plot based on model object 
plot_agb_data_corn <- summary(emmeans::emmeans(agb_corn_evo_model, ~comp:elevation_sc,
                                               weights = "proportional",
                                               at = list(elevation_sc = seq(min(traits_corn$elevation_sc),
                                                                            max(traits_corn$elevation_sc), length.out = 50))))

plot_agb_data_corn %>%
  mutate(agb_scam = emmean^2,
         elevation = elevation_sc*sd(traits_corn$elevation) + mean(traits_corn$elevation)) %>% 
  ggplot(aes(x = elevation, y = agb_scam, group = comp, linetype = as.factor(comp))) +
  geom_line(size = 0.7) +
  geom_point(data = traits_corn, aes(x = elevation, y = agb_scam, shape = as.factor(comp)), alpha = 0.5, size = 2) +
  ylab(expression(paste(italic("S. americanus")," aboveground biomass (g)"))) +
  xlab("elevation (m NAVD88)") +
  labs(linetype = "competition:", shape = "competition:") +
  theme_bw(base_size = 14) +
  theme(legend.position = "top", legend.key.width = unit(1.5,"cm"),
        legend.text.align = 0) +
  scale_linetype_manual(labels=expression(paste("without ", italic("S. patens")), paste("with ", italic("S. patens"))), values = c("dashed", "solid")) +
  scale_shape_manual(labels=expression(paste("without ", italic("S. patens")), paste("with ", italic("S. patens"))), values = c(1,16)) +
  guides(shape = guide_legend(override.aes = list(alpha = 1))) -> agb_comp_plot

png("figs/Fig7.png", height = 5.5, width = 5.5, res = 300, units = "in")
agb_comp_plot
dev.off()
