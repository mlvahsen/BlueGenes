# Figure 7 - Competition by elevation influence on aboveground biomass

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

## Create AGB competition plot ####

# Create flooding column in raw data
traits_corn %>% 
  mutate(flooding = mean_tide - elevation) -> traits_corn

# Create plot based on model object 
plot_agb_data_corn <- summary(emmeans::emmeans(agb_corn_evo_model, ~comp:elevation_sc,
                                               weights = "proportional",
                                               at = list(elevation_sc = seq(min(traits_corn$elevation_sc),
                                                                            max(traits_corn$elevation_sc), length.out = 50)))) 

plot_agb_data_corn %>%
  mutate(agb_scam = emmean^2,
         elevation = elevation_sc*sd(traits_corn$elevation) + mean(traits_corn$elevation),
         flooding = mean_tide - elevation) %>% 
  ggplot(aes(x = flooding, y = agb_scam, group = comp, linetype = as.factor(comp))) +
  geom_line(size = 0.7) +
  geom_point(data = traits_corn, aes(x = flooding, y = agb_scam, shape = as.factor(comp)), alpha = 0.5, size = 2) +
  ylab(expression(paste(italic("S. americanus")," aboveground biomass (g)"))) +
  xlab("average inundation (m)") +
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
