##
# Corn only (competition analysis)
##


## Data manipulation ####
bg_full %>% 
  filter(location == "corn") -> corn_only

corn_only %>% 
  filter(level < 5) %>% 
  filter(agb_scam > 0) -> traits_corn

# Count to see how many we have for each genotype in this data set
traits_corn %>% 
  group_by(genotype) %>% 
  summarize(n = n()) %>% 
  pull(n) %>% range()
# There are between 5 and 26 reps per genotype in this dataset

## Aboveground biomass ####

# Most complex model
agb_mod_corn <- lmer(sqrt(agb_scam) ~ weight_init + date_cloned_grp + origin_lab +
                       (age + comp + co2 + salinity + elevation)^5 + I(elevation^2) +
                       (1+co2*elevation + co2*comp + co2*salinity + elevation*comp + elevation*salinity + comp*salinity|genotype) + (1|site_frame), data = traits_corn)

# Try lmerTest::step
agb_corn_evo_model <- get_model(lmerTest::step(agb_mod_corn)) 

# Check for convergence issues
summary(agb_corn_evo_model)
# All good

# Check assumptions
plot_model(agb_corn_evo_model, type = "diag")
# Looks good

# Create plots 
plot_agb_data_corn <- summary(emmeans::emmeans(agb_corn_evo_model, ~comp:elevation, at = list(elevation = seq(0.156, 0.544, length.out = 50))))

plot_agb_data_corn %>%
  mutate(agb_scam = emmean^2) %>% 
  ggplot(aes(x = elevation, y = agb_scam, group = comp, linetype = as.factor(comp))) +
  geom_line(size = 0.7) +
  geom_point(data = traits_corn, aes(x = elevation, y = agb_scam, shape = as.factor(comp)), alpha = 0.5, size = 2) +
  ylab("aboveground biomass (g)") +
  xlab("elevation (m NAVD88)") +
  labs(linetype = "competition:", shape = "competition:") +
  theme_bw() +
  theme(legend.position = "top", legend.key.width = unit(1.5,"cm"),
        legend.text.align = 0) +
  scale_linetype_manual(labels=expression(paste("without ", italic("S. patens")), paste("with ", italic("S. patens"))), values = c("dashed", "solid")) +
  scale_shape_manual(labels=expression(paste("without ", italic("S. patens")), paste("with ", italic("S. patens"))), values = c(1,16)) +
  guides(shape = guide_legend(override.aes = list(alpha = 1))) -> agb_comp_plot

png("figs/FigX_comp.png", height = 4.2, width = 4.7, res = 300, units = "in")
agb_comp_plot
dev.off()

##
# Corn only
##

## Stem height ####
# Need to remove pot 1830 because most of the stems were eaten by caterpillars
traits_corn %>% filter(pot_no !=1830 ) ->traits_corn_height

height_corn_mod <- lmer(height_scam_tot ~ weight_init + date_cloned_grp + origin_lab +
                                (age + comp + co2 + salinity + elevation)^5 + I(elevation^2) +
                                (1+co2*elevation + co2*comp + co2*salinity + elevation*comp + elevation*salinity + comp*salinity|genotype) + (1|site_frame), data = traits_corn_sub)

height_corn_evo_model <- get_model(step(height_corn_mod))

# Check for convergence errors
summary(height_corn_evo_model)
# All good

# Check assumptions
plot_model(height_corn_evo_model, type = "diag")

## Stem width ####

width_mod_corn <- lmer(width_scam_mid ~ weight_init + date_cloned_grp + origin_lab +
                         (age + comp + co2 + salinity + elevation)^5 + I(elevation^2) +
                         (1|site_frame) + (1+co2*elevation + co2*comp + co2*salinity + elevation*comp + elevation*salinity + comp*salinity|genotype), data = traits_corn)

# Get stepwise model
width_corn_evo_model <- get_model(lmerTest::step(width_mod_corn))

# Check assumptions.
plot_model(width_corn_evo_model, type = "diag") # All looks good!

## Stem density ####

# Most complex model
density_mod_corn <- lmer(dens_scam_live ~ weight_init + date_cloned_grp + origin_lab +
                           (age + comp + co2 + salinity + elevation)^5 + I(elevation^2) +
                           (1|site_frame) + (1+co2*elevation + co2*comp + co2*salinity + elevation*comp + elevation*salinity + comp*salinity|genotype), data = traits_corn)

# Get stepwise model
density_corn_evo_model <- get_model(lmerTest::step(density_mod_corn))

# Check assumptions.
plot_model(density_corn_evo_model, type = "diag") # All looks good!

## Fig S8: effect of SPPA competition on stem height and width ####
# Get emmeans values for plots
width_comp_plot <- summary(emmeans::emmeans(width_corn_evo_model, ~comp))

tibble(width_scam_mid = width_comp_plot$emmean,
       lower = width_comp_plot$lower.CL,
       upper = width_comp_plot$upper.CL,
       comp = c(0, 1)) %>% 
  ggplot(aes(x = as.factor(comp), y = width_scam_mid, color = as.factor(comp))) +
  geom_jitter(data = traits_corn, aes(x = as.factor(comp), y = width_scam_mid, color = as.factor(comp)), shape = 1,
              stroke = 0.8, alpha = 0.5, size = 0.8)+
  theme_bw() +
  geom_errorbar(aes(ymin = lower, ymax = upper), position=position_dodge(width=0.3), width = 0.2, color = "black") +
  geom_point(size = 3, position=position_dodge(width=0.3)) +
  geom_point(size = 3, position=position_dodge(width=0.3), shape = 1, stroke = 0.8, color = "black") +
  ylab("mean stem width (mm)") +
  scale_color_manual(values = c("gray47", "black")) +
  scale_fill_manual(values = c("gray47", "black")) +
  theme(legend.position = "none") + 
  scale_x_discrete(labels=expression(paste("without ", italic("S. patens")), paste("with ", italic("S. patens")))) +
  xlab("") -> width_comp

# Repeat for stem height
height_comp_plot <- summary(emmeans::emmeans(height_corn_evo_model, ~comp))

tibble(height_scam_tot = height_comp_plot$emmean,
       lower = height_comp_plot$lower.CL,
       upper = height_comp_plot$upper.CL,
       comp = c(0, 1)) %>% 
  ggplot(aes(x = as.factor(comp), y = height_scam_tot, color = as.factor(comp))) +
  geom_jitter(data = traits_corn, aes(x = as.factor(comp), y = height_scam_tot, color = as.factor(comp)), shape = 1,
              stroke = 0.8, alpha = 0.5, size = 0.8)+
  theme_bw() +
  geom_errorbar(aes(ymin = lower, ymax = upper), position=position_dodge(width=0.3), width = 0.2, color = "black") +
  geom_point(size = 3, position=position_dodge(width=0.3)) +
  geom_point(size = 3, position=position_dodge(width=0.3), shape = 1, stroke = 0.8, color = "black") +
  ylab("mean stem height (cm)") +
  scale_color_manual(values = c("gray47", "black")) +
  scale_fill_manual(values = c("gray47", "black")) +
  theme(legend.position = "none") + 
  scale_x_discrete(labels=expression(paste("without ", italic("S. patens")), paste("with ", italic("S. patens")))) +
  xlab("") -> height_comp

FigS8 <- height_comp + width_comp + plot_annotation(tag_levels = "a")

ggsave(here("figs", "FigS8_hw_comp.png"), FigS8, height = 3.5, width = 5.5, units = "in")

