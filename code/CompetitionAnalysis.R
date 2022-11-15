# Fit linear models to subset of data to assess if S. patens presence influenced
# S. americanus traits and their response to environmental factors.

# Read in libraries 
library(here); library(tidyverse); library(emmeans);
library(lme4); library(lmerTest); library(blme);
library(brglm); library(lmtest); library(sjPlot)

# Read in compiled trait data
source(here("supp_code", "CompileTraitData.R"))

## Data manipulation ####

# Filter for only Corn Island genotypes, within levels 1-4, for extant plants 
bg_full %>% 
  filter(location == "corn" & level < 5 & agb_scam > 0) -> traits_corn

# Count to see how many we have for each genotype in this data set
traits_corn %>% 
  group_by(genotype) %>% 
  summarize(n = n()) %>% 
  pull(n) %>% range()
# There are between 5 and 26 reps per genotype in this dataset

# Create scaled elevation variable
traits_corn %>% 
  mutate(elevation_sc = scale(elevation)[,1])-> traits_corn

## Aboveground biomass model ####

# Most complex model: all fixed effects interactions and up to 2-way
# interactions for random slopes
agb_mod_corn <- lmer(sqrt(agb_scam) ~ weight_init + date_cloned_grp + origin_lab +
                       (age + comp + co2 + salinity + elevation_sc)^5 + I(elevation_sc^2) +
                       (1+co2*elevation_sc + co2*comp + co2*salinity + elevation_sc*comp +
                          elevation_sc*salinity + comp*salinity|genotype) + (1|site_frame), data = traits_corn)

# Do stepwise model selection
agb_corn_evo_model <- get_model(lmerTest::step(agb_mod_corn)) 

# Check for convergence issues
summary(agb_corn_evo_model)
# All good

# Check assumptions
plot_model(agb_corn_evo_model, type = "diag")
# Looks good

# Get significance of model terms
anova(agb_corn_evo_model)

## Stem height model ####
# Need to remove pot 1830 because most of the stems were eaten by caterpillars
traits_corn %>% filter(pot_no != 1830)->traits_corn_height

# Fit most complex model
height_corn_mod <- lmer(height_scam_tot ~ weight_init + date_cloned_grp + origin_lab +
                                (age + comp + co2 + salinity + elevation_sc)^5 + I(elevation_sc^2) +
                                (1+co2*elevation_sc + co2*comp + co2*salinity +
                                   elevation_sc*comp + elevation_sc*salinity + comp*salinity|genotype) + (1|site_frame),
                        data = traits_corn_height)

# Do stepwise model selection
height_corn_evo_model <- get_model(step(height_corn_mod))

# Check for convergence errors
summary(height_corn_evo_model)
# All good

# Check assumptions
plot_model(height_corn_evo_model, type = "diag")

# Significance of model terms
anova(height_corn_evo_model)

## Stem width model ####

width_mod_corn <- lmer(width_scam_mid ~ weight_init + date_cloned_grp + origin_lab +
                         (age + comp + co2 + salinity + elevation_sc)^5 + I(elevation_sc^2) +
                         (1|site_frame) + (1+co2*elevation_sc + co2*comp + co2*salinity +
                                             elevation_sc*comp + elevation_sc*salinity + comp*salinity|genotype), data = traits_corn)

# Do stepwise model selection
width_corn_evo_model <- get_model(lmerTest::step(width_mod_corn))

# Check for convergence errors
summary(width_corn_evo_model) # All good

# Check assumptions.
plot_model(width_corn_evo_model, type = "diag") # All looks good!

# Get significance of model terms
anova(width_corn_evo_model)

## Stem density model ####

# Most complex model
density_mod_corn <- lmer(dens_scam_live ~ weight_init + date_cloned_grp + origin_lab +
                           (age + comp + co2 + salinity + elevation_sc)^5 + I(elevation_sc^2) +
                           (1|site_frame) + (1+co2*elevation_sc + co2*comp + co2*salinity +
                                               elevation_sc*comp + elevation_sc*salinity + comp*salinity|genotype), data = traits_corn)

# Get stepwise model
density_corn_evo_model <- get_model(lmerTest::step(density_mod_corn))

# Check for convergence warnings
summary(density_corn_evo_model) # All good

# Check assumptions.
plot_model(density_corn_evo_model, type = "diag") # All looks good!

# Get significance of model terms
anova(density_corn_evo_model)

