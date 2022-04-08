##
# Corn only (competition analysis)
##


## Data manipulation ####

corn_only %>% 
  filter(level < 5) %>% 
  filter(agb_scam > 0) %>% 
  # Also scale elevation to help with convergence
  mutate(elevation_sc = scale(elevation)[,1]) -> traits_corn

# Count to see how many we have for each genotype in this data set
traits_corn %>% 
  group_by(genotype) %>% 
  summarize(n = n()) %>% 
  pull(n) %>% range()
# There are between 5 and 26 reps per genotype in this dataset

## Fit model and model selection (fixed effects) ####

# Most complex model
agb_mod_corn <- lmer(sqrt(agb_scam) ~ weight_init + date_cloned_grp + origin_lab +
                       (age + comp + co2 + salinity + elevation)^5 + I(elevation^2) +
                       (1|genotype) + (1|site_frame), data = traits_corn)

# Try lmerTest::step
get_model(lmerTest::step(agb_mod_corn), keep = keep_model_terms_corn) 

agb_corn_step_mod <- lmer(sqrt(agb_scam) ~ weight_init + date_cloned_grp + age + comp +  
                            co2 + salinity + elevation + I(elevation^2) + (1 | genotype) +  
                            age:comp + age:co2 + age:salinity + age:elevation + comp:co2 +  
                            comp:salinity + comp:elevation + co2:salinity + co2:elevation +  
                            salinity:elevation + age:comp:co2 + age:comp:salinity + age:comp:elevation +  
                            age:co2:salinity + age:co2:elevation + age:salinity:elevation +  
                            comp:co2:salinity + comp:co2:elevation + comp:salinity:elevation +  
                            co2:salinity:elevation + age:comp:co2:salinity + age:comp:co2:elevation +  
                            age:comp:salinity:elevation + age:co2:salinity:elevation +  
                            comp:co2:salinity:elevation + age:comp:co2:salinity:elevation, data = traits_corn)

# This is the final FIXED effects model. Check assumptions.
plot_model(agb_corn_step_mod, type = "diag") # All looks good!

## Model selection (random slopes) ####

# Now we are going to try and add up to 2-way interactions as random slopes.

# Here are the possibilities:
# co2:elevation
# co2:salinity
# salinity:elevation
# co2 + elevation
# co2 + salinity
# salinity + elevation
# co2
# elevation
# salinity

agb_mod_corn1 <- update(agb_corn_step_mod, .~.-(1|genotype) + (1+co2*elevation|genotype))
agb_mod_corn2 <- update(agb_corn_step_mod, .~.-(1|genotype) + (1+co2*salinity|genotype))
agb_mod_corn3 <- update(agb_corn_step_mod, .~.-(1|genotype) + (1+elevation*salinity|genotype))
agb_mod_corn4 <- update(agb_corn_step_mod, .~.-(1|genotype) + (1+comp*salinity|genotype))
agb_mod_corn5 <- update(agb_corn_step_mod, .~.-(1|genotype) + (1+comp*elevation|genotype))
agb_mod_corn6 <- update(agb_corn_step_mod, .~.-(1|genotype) + (1+comp*co2|genotype))
agb_mod_corn7 <- update(agb_corn_step_mod, .~.-(1|genotype) + (1+co2+elevation|genotype))
agb_mod_corn8 <- update(agb_corn_step_mod, .~.-(1|genotype) + (1+elevation + salinity|genotype))
agb_mod_corn9 <- update(agb_corn_step_mod, .~.-(1|genotype) + (1+salinity + co2|genotype))
agb_mod_corn10 <- update(agb_corn_step_mod, .~.-(1|genotype) + (1+comp+salinity|genotype))
agb_mod_corn11 <- update(agb_corn_step_mod, .~.-(1|genotype) + (1+comp+elevation|genotype))
agb_mod_corn12 <- update(agb_corn_step_mod, .~.-(1|genotype) + (1+comp+co2|genotype))
agb_mod_corn13 <- update(agb_corn_step_mod, .~.-(1|genotype) + (1+co2|genotype))
agb_mod_corn14<- update(agb_corn_step_mod, .~.-(1|genotype) + (1+elevation|genotype))#*
agb_mod_corn15 <- update(agb_corn_step_mod, .~.-(1|genotype) + (1+salinity|genotype))
agb_mod_corn16 <- update(agb_corn_step_mod, .~.-(1|genotype) + (1+comp|genotype))

# Compare to base model (elevation is significant)
anova(agb_corn_step_mod, agb_mod_corn14, refit = F)

# Create plots 
plot_agb_data_corn <- summary(emmeans::emmeans(agb_corn_step_mod, ~elevation:co2:comp:salinity:age, at = list(elevation = seq(0.156, 0.544, length.out = 50))))

traits_corn %>% 
  mutate(age = case_when(age == "modern" ~ "descendant cohort (2000-2020)",
                         T ~ "ancestral cohort (1900-1950)")) %>% 
  mutate(salinity = case_when(salinity == "fresh" ~ "freshwater site (4ppt)",
                              T ~ "brackish site (6ppt)")) -> traits_corn_plot

plot_agb_data_corn %>% 
  filter(comp == 0) %>% 
  mutate(agb_scam = emmean^2) %>% 
  mutate(age = case_when(age == "modern" ~ "descendant cohort (2000-2020)",
                         T ~ "ancestral cohort (1900-1950)")) %>% 
  mutate(salinity = case_when(salinity == "fresh" ~ "freshwater site (4ppt)",
                              T ~ "brackish site (6ppt)")) %>% 
  ggplot(aes(x = elevation, y = agb_scam, color = co2)) +
  geom_point(data = traits_corn_plot %>% filter(comp == 0), shape = 1, alpha = 0.6, stroke = 0.7) +
  geom_ribbon(aes(ymin = lower.CL^2, ymax = upper.CL^2, fill = co2), alpha = 0.2, color = NA) +
  geom_textline(aes(label = co2), linewidth = 1.2, fontface = 2, hjust = 0.2) +
  facet_grid(salinity ~ age) +
  scale_color_manual(values = c("#b2df8a","#33a02c")) +
  scale_fill_manual(values = c("#b2df8a","#33a02c")) +
  ylab("aboveground biomass (g)") +
  xlab("elevation (m NAVD88)") +
  theme_bw() + ylim(0,15) +
  theme(legend.position = "none") -> agb_plot_nocomp

plot_agb_data_corn %>% 
  filter(comp == 1) %>% 
  mutate(agb_scam = emmean^2) %>% 
  mutate(age = case_when(age == "modern" ~ "descendant cohort (2000-2020)",
                         T ~ "ancestral cohort (1900-1950)")) %>% 
  mutate(salinity = case_when(salinity == "fresh" ~ "freshwater site (4ppt)",
                              T ~ "brackish site (6ppt)")) %>% 
  ggplot(aes(x = elevation, y = agb_scam, color = co2)) +
  geom_point(data = traits_corn_plot %>% filter(comp == 1), shape = 1, alpha = 0.6, stroke = 0.7) +
  geom_ribbon(aes(ymin = lower.CL^2, ymax = upper.CL^2, fill = co2), alpha = 0.2, color = NA) +
  geom_textline(aes(label = co2), linewidth = 1.2, fontface = 2, hjust = 0.2) +
  facet_grid(salinity ~ age) +
  scale_color_manual(values = c("#b2df8a","#33a02c")) +
  scale_fill_manual(values = c("#b2df8a","#33a02c")) +
  ylab("") +
  xlab("elevation (m NAVD88)") +
  theme_bw() + ylim(0,15) +
  theme(legend.position = "none") -> agb_plot_comp

agb_plot_nocomp + agb_plot_comp + plot_annotation(tag_levels = "a") -> agb_comp_comparison

ggsave(here("figs", "SuppFig_AgbCompComparison.png"), agb_comp_comparison, height = 4, width = 9, units = "in")

# Repeat this process for the competition pots

##
# Corn only
##

## Fit model and model selection (fixed effects) ####

# Most complex model
height_mod_corn <- lmer(height_scam_tot ~ weight_init + date_cloned_grp + origin_lab +
                          (age + comp + co2 + salinity + elevation)^5 + I(elevation^2) +
                          (1|genotype) + (1|site_frame), data = traits_corn)

which(resid(height_mod_corn) < -15)

# Need to remove pot 1830 because most of the stems were cut

traits_corn %>% filter(pot_no !=1830 ) ->traits_corn_sub
height_corn_noOut_mod <- lmer(height_scam_tot ~ weight_init + date_cloned_grp + origin_lab +
                                (age + comp + co2 + salinity + elevation)^5 + I(elevation^2) +
                                (1|genotype) + (1|site_frame), data = traits_corn_sub)

get_model(step(height_corn_noOut_mod, keep = keep_model_terms_corn))

height_corn_noOut_mod_step <- lmer(height_scam_tot ~ weight_init + age + comp + co2 + salinity +  
                                     elevation + (1 | genotype) + co2:elevation, data = traits_corn_sub)


# Check assumptions
plot_model(height_corn_noOut_mod_step, type = "diag")

## Model selection (random slopes) ####

# Now we are going to try and add up to 2-way interactions as random slopes.

height_mod_corn1 <- update(height_corn_noOut_mod_step, .~.-(1|genotype) + (1+co2*elevation|genotype))
height_mod_corn2 <- update(height_corn_noOut_mod_step, .~.-(1|genotype) + (1+co2+elevation|genotype))#*
height_mod_corn3 <- update(height_corn_noOut_mod_step, .~.-(1|genotype) + (1+co2+comp|genotype))
height_mod_corn4 <- update(height_corn_noOut_mod_step, .~.-(1|genotype) + (1+comp+elevation|genotype))
height_mod_corn5 <- update(height_corn_noOut_mod_step, .~.-(1|genotype) + (1+comp+salinity|genotype))
height_mod_corn6 <- update(height_corn_noOut_mod_step, .~.-(1|genotype) + (1+elevation + salinity|genotype))
height_mod_corn7 <- update(height_corn_noOut_mod_step, .~.-(1|genotype) + (1+salinity + co2|genotype))
height_mod_corn8 <- update(height_corn_noOut_mod_step, .~.-(1|genotype) + (1+co2|genotype))#*
height_mod_corn9<- update(height_corn_noOut_mod_step, .~.-(1|genotype) + (1+elevation|genotype))#*
height_mod_corn10 <- update(height_corn_noOut_mod_step, .~.-(1|genotype) + (1+salinity|genotype))
height_mod_corn11 <- update(height_corn_noOut_mod_step, .~.-(1|genotype) + (1+comp|genotype))

# Compare to base model (all ns)

anova(height_corn_noOut_mod_step, height_mod_corn2)
anova(height_corn_noOut_mod_step, height_mod_corn8)
anova(height_corn_noOut_mod_step, height_mod_corn9)

anova(height_corn_noOut_mod_step)

height_corn_comp_plot <- summary(emmeans::emmeans(height_corn_noOut_mod_step, ~comp))

tibble(height_scam_tot = height_corn_EC_plot$emmean,
       co2 = height_corn_EC_plot$co2,
       elevation = height_corn_EC_plot$elevation,
       lower = height_corn_EC_plot$lower.CL,
       upper = height_corn_EC_plot$upper.CL) %>% 
  ggplot(aes(x = elevation, y = height_scam_tot, color = co2)) +
  geom_textline(label = height_corn_EC_plot$co2, hjust = 0.9) +
  theme_bw() +
  geom_point(data = traits_corn_sub, aes(x = elevation, y = height_scam_tot, color = co2), shape = 1, stroke = 0.8, alpha = 0.7, size = 0.8) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = co2), alpha = 0.2, color = NA) +
  ylab("mean stem height (cm)") +
  xlab("elevation (m NAVD88)") +
  scale_color_manual(values = c("#b2df8a", "#33a02c")) +
  scale_fill_manual(values = c("#b2df8a", "#33a02c"))


# Repeat this process for the competition pots

##
# Corn only
##

## Fit model and model selection (fixed effects) ####

# Most complex model
width_mod_corn <- lmer(width_scam_mid ~ weight_init + date_cloned_grp + origin_lab +
                         (age + comp + co2 + salinity + elevation)^5 + I(elevation^2) +
                         (1|site_frame) + (1|genotype), data = traits_corn)

# Try step mod
get_model(lmerTest::step(width_mod_corn, keep = keep_model_terms_corn))

# Step mod
width_corn_step_mod <- lmer(width_scam_mid ~ weight_init + age + comp + co2 + salinity +  
                              elevation + (1 | genotype) + age:co2 + age:salinity + age:elevation +  
                              co2:salinity + co2:elevation + salinity:elevation + age:co2:salinity + age:salinity:elevation, data = traits_corn)

# This is the final FIXED effects model. Check assumptions.
plot_model(width_corn_step_mod, type = "diag") # All looks good!

## Model selection (random slopes) ####

# Now we are going to try and add up to 2-way interactions as random slopes.

# Here are the possibilities:
# co2:elevation
# co2:salinity
# salinity:elevation
# co2 + elevation
# co2 + salinity
# salinity + elevation
# co2
# elevation
# salinity

width_mod_corn1 <- update(width_corn_step_mod, .~.-(1|genotype) + (1+co2*elevation|genotype))
width_mod_corn2 <- update(width_corn_step_mod, .~.-(1|genotype) + (1+co2*salinity|genotype))
width_mod_corn3 <- update(width_corn_step_mod, .~.-(1|genotype) + (1+elevation*salinity|genotype))
width_mod_corn4 <- update(width_corn_step_mod, .~.-(1|genotype) + (1+comp*salinity|genotype))#*
width_mod_corn5 <- update(width_corn_step_mod, .~.-(1|genotype) + (1+elevation*comp|genotype))
width_mod_corn6 <- update(width_corn_step_mod, .~.-(1|genotype) + (1+comp*co2|genotype))
width_mod_corn7 <- update(width_corn_step_mod, .~.-(1|genotype) + (1+co2+elevation|genotype))
width_mod_corn8 <- update(width_corn_step_mod, .~.-(1|genotype) + (1+elevation + salinity|genotype))
width_mod_corn9 <- update(width_corn_step_mod, .~.-(1|genotype) + (1+salinity + co2|genotype))
width_mod_corn10 <- update(width_corn_step_mod, .~.-(1|genotype) + (1+comp+salinity|genotype))
width_mod_corn11 <- update(width_corn_step_mod, .~.-(1|genotype) + (1+elevation+comp|genotype))
width_mod_corn12<- update(width_corn_step_mod, .~.-(1|genotype) + (1+comp+co2|genotype))
width_mod_corn13<- update(width_corn_step_mod, .~.-(1|genotype) + (1+co2|genotype))
width_mod_corn14<- update(width_corn_step_mod, .~.-(1|genotype) + (1+comp|genotype))
width_mod_corn15<- update(width_corn_step_mod, .~.-(1|genotype) + (1+elevation|genotype))
width_mod_corn16 <- update(width_corn_step_mod, .~.-(1|genotype) + (1+salinity|genotype))#*

# Compare to base model (all ns)
anova(width_corn_step_mod, width_mod_corn4)
anova(width_corn_step_mod, width_mod_corn16)

# Competition model for density

## Fit model and model selection (fixed effects) ####

# Most complex model
density_mod_corn <- lmer(dens_scam_live ~ weight_init + date_cloned_grp + origin_lab +
                           (age + comp + co2 + salinity + elevation)^5 + I(elevation^2) +
                           (1|site_frame) + (1|genotype), data = traits_corn)

# Try step mod
get_model(lmerTest::step(density_mod_corn, keep = keep_model_terms_corn))

# Step mod
density_corn_step_mod <- lmer(dens_scam_live ~ weight_init + date_cloned_grp +
                                age + comp + co2 + salinity + elevation + I(elevation^2) + (1 | genotype) +  
                                comp:elevation, data = traits_corn)

# This is the final FIXED effects model. Check assumptions.
plot_model(density_corn_step_mod, type = "diag") # All looks good!

## Model selection (random slopes) ####

# Now we are going to try and add up to 2-way interactions as random slopes.

density_mod_corn1 <- update(density_corn_step_mod, .~.-(1|genotype) + (1+co2*elevation|genotype))
density_mod_corn2 <- update(density_corn_step_mod, .~.-(1|genotype) + (1+co2*salinity|genotype))
density_mod_corn3 <- update(density_corn_step_mod, .~.-(1|genotype) + (1+elevation*salinity|genotype))
density_mod_corn4 <- update(density_corn_step_mod, .~.-(1|genotype) + (1+comp*salinity|genotype))#*
density_mod_corn5 <- update(density_corn_step_mod, .~.-(1|genotype) + (1+elevation*comp|genotype))
density_mod_corn6 <- update(density_corn_step_mod, .~.-(1|genotype) + (1+comp*co2|genotype))
density_mod_corn7 <- update(density_corn_step_mod, .~.-(1|genotype) + (1+co2+elevation|genotype))
density_mod_corn8 <- update(density_corn_step_mod, .~.-(1|genotype) + (1+elevation + salinity|genotype))
density_mod_corn9 <- update(density_corn_step_mod, .~.-(1|genotype) + (1+salinity + co2|genotype))
density_mod_corn10 <- update(density_corn_step_mod, .~.-(1|genotype) + (1+comp+salinity|genotype))
density_mod_corn11<- update(density_corn_step_mod, .~.-(1|genotype) + (1+elevation+comp|genotype))
density_mod_corn12<- update(density_corn_step_mod, .~.-(1|genotype) + (1+comp+co2|genotype))
density_mod_corn13<- update(density_corn_step_mod, .~.-(1|genotype) + (1+co2|genotype))
density_mod_corn14<- update(density_corn_step_mod, .~.-(1|genotype) + (1+comp|genotype))
density_mod_corn15<- update(density_corn_step_mod, .~.-(1|genotype) + (1+elevation|genotype))
density_mod_corn16<- update(density_corn_step_mod, .~.-(1|genotype) + (1+salinity|genotype))#*
# nothing fits

plot_model(density_corn_step_mod, type = "diag") # All looks good!
