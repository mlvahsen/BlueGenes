# Fit linear models to subset of data to assess the role of eco-evo factors in
# explaining S. americanus trait variation

# Load libraries 
library(here); library(tidyverse); library(emmeans);
library(lme4); library(lmerTest); library(blme);
library(brglm); library(lmtest); library(sjPlot);
library(kableExtra); library(geomtextpath); library(patchwork);
library(effectsize)

# Read in compiled trait data
source(here("supp_code", "CompileTraitData.R"))
##################
## AGB analysis ##
##################

# For the trait analysis, we are going to look at just those that survived in
# the top 4 levels. I think if we keep in the 25% that survived in level 5 and
# the 2.5% that survived in level 6, this could skew the parabolas of the
# genotypes.

# We will also take the same approach as above by first taking out blackwater
# and doing a no competition analysis. And then doing a competition analysis on
# just the corn data.

##
# No competition
##

## Data manipulation ####

# We are analyzing non-competition pots, from levels 1-4, from Corn or
# Kirkpatrick, that were extant at harvest.
bg_full %>% 
  filter(comp == 0 & location != "blackwater" & level < 5 &
            agb_scam > 0) -> traits_nocomp

# There are 230 observations in this data set
nrow(traits_nocomp)

# Count to see how many we have for each genotype in this data set
traits_nocomp %>% 
  group_by(genotype) %>% 
  summarize(n = n()) %>% 
  pull(n) %>% range()
# There are between 2 and 12 reps per genotype in this dataset

# Center and scale elevation values
traits_nocomp %>% 
  mutate(elevation_sc = scale(elevation)[,1],
         elevation_sc2 = scale(elevation^2)[,1]) -> traits_nocomp

## Fit model and model selection (random effects first, then fixed effects) ####

# Most complex fixed effects model with up to two-way interactions for random
# slopes
agb_mod <- lmer(sqrt(agb_scam) ~ weight_init + date_cloned_grp + origin_lab +
                  (age + location + co2 + salinity + elevation_sc)^5 +
                  I(elevation_sc^2) +
                  (1+co2*elevation_sc + salinity*co2 + salinity*elevation_sc|genotype) +
                  (1|site_frame), data = traits_nocomp)

# Use lmerTest::step for backwards model selection
step_agb <- get_model(lmerTest::step(agb_mod))

# Check to see if there are any convergence warnings
summary(step_agb) # All good

# Rename model
agb_evo_model <- step_agb

# Check model assumptions
plot_model(agb_evo_model, type = "diag") # All looks good!

plot_model(step_agb, type = "emm", terms = c("elevation_sc[all]", "co2", "salinity", "age"))
## Comparison to null model ####

# Now construct a null model that does not have age, location, or genotype in it
agb_mod_null <- lmer(sqrt(agb_scam) ~ weight_init + date_cloned_grp + origin_lab +
                  (co2 + salinity + elevation_sc)^3 + I(elevation_sc^2) +
                  (1|site_frame), data = traits_nocomp)

# Conduct stepwise model selection
agb_null_model <- get_model(lmerTest::step(agb_mod_null))

# Check for convergence warnings
summary(agb_null_model) # All good

# Compare conditional R2 from models
MuMIn::r.squaredGLMM(agb_null_model)
# Conditional R2 = 0.35
MuMIn::r.squaredGLMM(agb_evo_model)
# Conditional R2 = 0.52

# Extract type III Anova table from evolution model
anova(agb_evo_model)

## Calculate effect sizes of fixed effects ####

# Create effect size table. Here, we choose epsilon because it is an unbiased
# estimator of the population value eta: "Epsilon is analogous to adjusted R2
# (Allen, 2017, p. 382), and has been found to be less biased (Carroll &
# Nordholm, 1975)." -- from package documentation
epsilon_squared(agb_evo_model)

# Pull significant terms from anova table and do not include covariates
tibble(anova(agb_evo_model)) %>% 
  mutate(term = rownames(anova(agb_evo_model))) %>% 
  filter(`Pr(>F)` < 0.05 & term != "weight_init" & term!= "date_cloned_grp") %>% 
  pull(term) -> significant_agb_effects

tibble(epsilon_squared((agb_evo_model))) %>% 
  filter(Parameter %in% significant_agb_effects) %>% 
  mutate(category = c("eco", "eco", "eco", "eco-eco", "eco-evo", "eco-evo")) %>% 
  group_by(category) %>% 
  summarize(total_epsilon2 = sum(Epsilon2_partial)) -> agb_effectsize

#########################
## Root:Shoot analysis ##
#########################

# For the root:shoot analysis, we can only really do the no competition pots in
# the straightforward way (b/c species were not separated for the entire core;
# only separated in the top 10cm)

## Data manipulation ####

# Start with the same data set as the aboveground analysis. Calculate
# root-to-shoot ratio
traits_nocomp %>% 
  mutate(rs = total_bg / agb_scam) -> traits_nocomp

# Plot and see if there are outliers
traits_nocomp %>% 
  ggplot(aes(x = rs)) +
  geom_histogram() +
  geom_rug() +
  xlab("root-to-shoot ratio")

# There are two pots that have really high r:s and one that seems a bit high.
# Take a look at these further.
traits_nocomp %>% 
  filter(rs > 4.5)
# It seems like pot_no 98 is real -- just a lot of bg compared to ag (on level
# 1)

# The other two seem to be more of experimental artifacts (i.e. initial
# propagule weight is driving total belowground weight at the end), so we drop
# those two.

# There also appear to be some that are NAs
traits_nocomp %>% 
  filter(!complete.cases(rs))
# There were issues with data collection for pots 205 and 1700 (i.e. mislabeled
# or missing bags for sections of bgb) so these two are dropped from the
# analysis because we could not reconcile the issues.

# We also need to drop two pots #165 and #176 where we needed to harvest rhizome
# to "save" the genotype for archiving.

traits_nocomp %>% 
  filter(rs < 6 & pot_no !=165 & pot_no !=176) -> traits_nocomp_rs

# Calculate the total N at this point
nrow(traits_nocomp_rs)
# 224

## Fit model and model selection (fixed effects) ####

# Most complex fixed and random effects model
rs_mod <- lmer(log(rs) ~ weight_init + date_cloned_grp + origin_lab +
                 (age + location + co2 + salinity + elevation_sc)^5 + elevation_sc2 +
                 (1+salinity*elevation_sc + co2*elevation_sc + salinity*co2|genotype) +
                 (1|site_frame), data = traits_nocomp_rs)

# Get the best stepwise model
step_rs <- get_model(lmerTest::step(rs_mod))

# Check to see if there are any convergence errors
summary(step_rs) # All good

# Check assumptions
plot_model(step_rs, type = "diag") # all good

# Rename model object
rs_evo_model <- step_rs

## Comparison to null model ####

# Now construct a null model that does not have age, location, or genotype in it
rs_modelnull <- lmer(log(rs) ~ weight_init + date_cloned_grp + origin_lab +
                       (co2 + salinity + elevation_sc)^3 + I(elevation_sc^2) +
                       (1|site_frame),
                     data = traits_nocomp_rs)

# Get stepwise model
rs_null_model <- get_model(lmerTest::step(rs_modelnull))

# Compare R2 from models
MuMIn::r.squaredGLMM(rs_null_model)
# Conditional R2 = 0.61
MuMIn::r.squaredGLMM(rs_evo_model)
# Conditional R2 = 0.70

# Get overall significance of model terms
anova(rs_evo_model)

## Calculate effect sizes of fixed effects ####
epsilon_squared(rs_evo_model)

# Pull significant terms from anova table 
tibble(anova(rs_evo_model)) %>% 
  mutate(term = rownames(anova(rs_evo_model))) %>% 
  filter(`Pr(>F)` < 0.05 & term != "weight_init" & term!= "date_cloned_grp") %>% 
  pull(term) -> significant_rs_effects

tibble(epsilon_squared((rs_evo_model))) %>% 
  filter(Parameter %in% significant_rs_effects) %>% 
  mutate(category = c("eco", "eco", "eco-eco", "eco-evo")) %>% 
  group_by(category) %>% 
  summarize(total_epsilon2 = sum(Epsilon2_partial)) -> rs_effectsize

##########################################
## Root distribution parameter analysis ##
##########################################
## Data manipulation ####

# These data can be the same as the root-to-shoot data

## Fit model and model selection (fixed effects) ####

# Most complex model

# Fit the most complex model
beta_mod <- lmer(beta ~ weight_init + date_cloned_grp + origin_lab +
                   (age + location + co2 + salinity + elevation_sc)^5 + I(elevation_sc^2) +
                   (1|site_frame) + (1 + salinity*elevation_sc + co2*elevation_sc + salinity*co2|genotype), data = traits_nocomp_rs)

# Check for convergence errors
summary(beta_mod) # Issues here!

# Fit the most complex model where we can keep the fixed part of the random
# effects in the model without convergence errors. We don't want to drop out the
# fixed part of a random slope because that means that we would be estimating it
# to be zero. We tested these individually and this was the most complex that fit
beta_mod_SE <- lmer(beta ~ weight_init + date_cloned_grp + origin_lab +
                   (age + location + co2 + salinity + elevation_sc)^5 + I(elevation_sc^2) +
                   (1|site_frame) + (1 + salinity*elevation_sc|genotype), data = traits_nocomp_rs)

# Try lmerTest::step while keeping the fixed effect part of the random slope
step_beta <- get_model(lmerTest::step(beta_mod_SE, keep = "salinity:elevation_sc"))

# Check for convergence errors on this model
summary(step_beta) # All good

# Rename model object
beta_evo_model <- step_beta

## Comparison to null model ####

# Now construct a null model that does not have age, location, or genotype in it
beta_modelnull <- lmer(beta ~ weight_init + date_cloned_grp + origin_lab +
                         (co2 + salinity + elevation)^3 + I(elevation^2) +
                         (1|site_frame), data = traits_nocomp_rs)

# Do stepwise selection for the null model
beta_null_model <- get_model(lmerTest::step(beta_modelnull))

# Compare R2 from models
MuMIn::r.squaredGLMM(beta_null_model)
# Conditional R2 = 0.48
MuMIn::r.squaredGLMM(beta_evo_model)
# Conditional R2 = 0.77

## Calculate effect sizes of fixed effects ####
epsilon_squared(beta_evo_model)

# Pull significant terms from anova table 
tibble(anova(beta_evo_model)) %>% 
  mutate(term = rownames(anova(beta_evo_model))) %>% 
  filter(`Pr(>F)` < 0.05 & term != "weight_init" & term!= "date_cloned_grp") %>% 
  pull(term) -> significant_beta_effects

tibble(epsilon_squared((beta_evo_model))) %>% 
  filter(Parameter %in% significant_beta_effects) %>% 
  mutate(category = c("evo", "evo", "eco", "eco-evo")) %>% 
  group_by(category) %>% 
  summarize(total_epsilon2 = sum(Epsilon2_partial)) -> beta_effectsize

##################################
## Belowground biomass analysis ##
##################################

# For the belowground analysis, we can only really do the no competition pots in
# the straightforward way (b/c species were not separated for the entire core;
# only separated in the top 10cm). Use the same subset of data as for the
# root:shoot analysis and beta analysis.

## Fit model and model selection (fixed effects) ####

# Again, use the same data set as for root-to-shoot ratio and beta

# Fit most complex model 
bg_mod <- lmer(sqrt(total_bg) ~ weight_init + date_cloned_grp + origin_lab +
                 (age + location + co2 + salinity + elevation_sc)^5 + I(elevation_sc^2) +
                 (1+co2*elevation_sc + salinity*elevation_sc + co2*salinity|genotype) + (1|site_frame), data = traits_nocomp_rs)

# Do stepwise model selection
step_bg <- get_model(lmerTest::step(bg_mod))

# Check for convergence errors
summary(step_bg) # All good

# Check model assumptions
plot_model(step_bg, type = "diag") # All good

# Rename model object
bg_evo_model <- step_bg

## Comparison to null model ####

# Now construct a null model that does not have age, location, or genotype in it
bg_modelnull <- lmer(sqrt(total_bg) ~ weight_init + date_cloned_grp + origin_lab +
                       (co2 + salinity + elevation_sc)^3 + I(elevation_sc^2) +
                       (1|site_frame),
                     data = traits_nocomp_rs)

bg_null_model <- get_model(step(bg_modelnull))
bg_null_model <- lm(sqrt(total_bg) ~ weight_init + date_cloned_grp +
                        elevation_sc + I(elevation_sc^2) + co2, data = traits_nocomp_rs)
# Need to keep elevation term in the model with elevation2

# Compare R2 from models
MuMIn::r.squaredGLMM(bg_null_model)
# Conditional R2 = 0.22
MuMIn::r.squaredGLMM(bg_evo_model)
# Conditional R2 = 0.39

# Get anova table of full model
anova(bg_evo_model)

## Calculate effect sizes of fixed effects ####
epsilon_squared(bg_evo_model)

# Pull significant terms from anova table 
tibble(anova(bg_evo_model)) %>% 
  mutate(term = rownames(anova(bg_evo_model))) %>% 
  filter(`Pr(>F)` < 0.05 & term != "weight_init" & term!= "date_cloned_grp") %>% 
  pull(term) -> significant_bg_effects

tibble(epsilon_squared((bg_evo_model))) %>% 
  filter(Parameter %in% significant_bg_effects) %>% 
  mutate(category = c("eco", "eco", "eco-eco")) %>% 
  group_by(category) %>% 
  summarize(total_epsilon2 = sum(Epsilon2_partial)) -> bg_effectsize

##########################
## Stem height analysis ##
##########################

# Use the same data as for the aboveground biomass analysis

##
# No competition
##

## Fit model and model selection (fixed effects) ####

# Start with same data set as for aboveground biomass

# Most complex model
height_mod <- lmer(height_scam_tot ~ weight_init + date_cloned_grp + origin_lab +
                     (age + location + co2 + salinity + elevation_sc)^5 + I(elevation_sc^2) +
                     (1|site_frame) + (1+co2*elevation_sc + salinity*elevation_sc + co2*salinity|genotype), data = traits_nocomp)

# Do stepwise model selection
step_height <- get_model(lmerTest::step(height_mod))

# Check model assumptions
plot_model(step_height, type = "diag") # There appears to be a huge outlier. 

# Figure out which data point this is
predict(step_height) < 30
traits_nocomp[186,]
# There's a note that there was a bunch of caterpillar herbivory for pot 1830
# (top-down) so we should drop it for the height analysis
traits_nocomp %>% 
  filter(pot_no != 1830) -> traits_nocomp_height

# Refit the full model with the updated dataset
height_mod_noOut <- lmer(height_scam_tot ~ weight_init + date_cloned_grp + origin_lab +
                           (age + location + co2 + salinity + elevation_sc)^5 + I(elevation_sc^2) +
                           (1+co2*elevation_sc + salinity*elevation_sc + salinity*co2|genotype) + (1|site_frame), data = traits_nocomp_height)

# Stepwise model selection
step_height_noOut <- get_model(lmerTest::step(height_mod_noOut))

# Check for convergence issues
summary(step_height_noOut) # Convergence error

# Need to fit most complex model that still has the fixed effect part of the
# random slope. We went through multiple models and this was the most complex
# that fit without convergence errors and also fit the above criteria

height_mod_noOut_SE <- lmer(height_scam_tot ~ weight_init + date_cloned_grp + origin_lab +
                                (age + location + co2 + salinity + elevation_sc)^5 + I(elevation_sc^2) +
                                (1+salinity*elevation_sc|genotype) + (1|site_frame), data = traits_nocomp_height)

# Do stepwise model selection while making sure fixed effect terms are not
# dropped out if they are in random slopes
step_height_noOut_SE <- get_model(lmerTest::step(height_mod_noOut_SE,
                                                 keep = c("salinity:elevation_sc", "salinity", "elevation_sc")))

# Check for convergence errors
summary(step_height_noOut_SE) # All good

# Check model assumptions
plot_model(step_height_noOut_SE, type = "diag") # All good

# Rename model object
height_evo_model <- step_height_noOut_SE

## Comparison to null model ####

# Now construct a null model that does not have age, location, or genotype in it
height_model_noOut_null <- lmer(height_scam_tot ~ weight_init + date_cloned_grp + origin_lab +
                                (co2 + salinity + elevation_sc)^3 + I(elevation_sc^2) +
                                (1|site_frame),
                              data = traits_nocomp_height)

# Do stepwise model selection
height_null_model <- get_model(step(height_model_noOut_null))

# Compare R2 from models
MuMIn::r.squaredGLMM(height_null_model)
# Adjusted R2 = 0.38
MuMIn::r.squaredGLMM(height_evo_model)
# Conditional R2 = 0.69

# Get anova table for full evo model
anova(height_evo_model)

## Calculate effect sizes of fixed effects ####
epsilon_squared(height_evo_model)

# Pull significant terms from anova table 
tibble(anova(height_evo_model)) %>% 
  mutate(term = rownames(anova(height_evo_model))) %>% 
  filter(`Pr(>F)` < 0.05 & term != "weight_init" & term!= "date_cloned_grp") %>% 
  pull(term) -> significant_height_effects

tibble(epsilon_squared((height_evo_model))) %>% 
  filter(Parameter %in% significant_height_effects) %>% 
  mutate(category = c("eco", "eco", "eco-evo")) %>% 
  group_by(category) %>% 
  summarize(total_epsilon2 = sum(Epsilon2_partial)) -> height_effectsize

##########################
## Stem width analysis ##
##########################

# Use the same as for the aboveground biomass analysis

## Fit model and model selection (fixed effects) ####

# Most complex model
width_mod <- lmer(width_scam_mid ~ weight_init + date_cloned_grp + origin_lab +
                    (age + location + co2 + salinity + elevation_sc)^5 + I(elevation_sc^2) +
                    (1+co2*elevation_sc+salinity*elevation_sc + co2*salinity|genotype) + (1|site_frame), data = traits_nocomp)

# Stepwise selection procedure
step_width <- get_model(lmerTest::step(width_mod))

# Check for convergence errors
summary(step_width) # All good

# Check model assumptions
plot_model(step_width, type = "diag") # All looks good!

# Rename model object
width_evo_model <- step_width

## Comparison to null model ####

# Now construct a null model that does not have age, location, or genotype in it
width_model_null <- lmer(width_scam_mid ~ weight_init + date_cloned_grp + origin_lab +
                         (co2 + salinity + elevation_sc)^3 + I(elevation_sc^2) +
                         (1|site_frame), data = traits_nocomp)

# Do stepwise model selection
width_null_model <- get_model(step(width_model_null))

# Compare R2 from models
MuMIn::r.squaredGLMM(width_null_model)
# Conditional R2 = 0.42
MuMIn::r.squaredGLMM(width_evo_model)
# Conditional R2 = 0.63

# Get anova table
anova(width_evo_model)

## Calculate effect sizes of fixed effects ####
epsilon_squared(width_evo_model)

# Pull significant terms from anova table 
tibble(anova(width_evo_model)) %>% 
  mutate(term = rownames(anova(width_evo_model))) %>% 
  filter(`Pr(>F)` < 0.05 & term != "weight_init" & term!= "date_cloned_grp") %>% 
  pull(term) -> significant_width_effects

tibble(epsilon_squared((width_evo_model))) %>% 
  filter(Parameter %in% significant_width_effects) %>% 
  mutate(category = c("eco", "eco", "eco-evo", "eco-evo")) %>% 
  group_by(category) %>% 
  summarize(total_epsilon2 = sum(Epsilon2_partial)) -> width_effectsize

###########################
## Stem density analysis ##
###########################

## Fit model ####
# Stem density uses the same data as AGB
density_mod <- lmer(dens_scam_live ~ weight_init + date_cloned_grp + origin_lab +
                      (age + location + co2 + salinity + elevation_sc)^5 + I(elevation_sc^2) +
                      (1+salinity*elevation_sc + co2*elevation_sc + salinity*co2|genotype) + (1|site_frame), data = traits_nocomp)

# Do stepwise model selection
step_density <- get_model(lmerTest::step(density_mod))

# Check for convergence errors
summary(step_density) # not good

# After trying models with less random effects, no models with random slopes
# converge when making sure the environmental variable that is included in the
# random slope is also included in the fixed effects.
density_mod_int <- lmer(dens_scam_live ~ weight_init + date_cloned_grp + origin_lab +
                      (age + location + co2 + salinity + elevation_sc)^5 + I(elevation_sc^2) +
                      (1|genotype) + (1|site_frame), data = traits_nocomp)

# Do stepwise selection
step_density_int <- get_model(lmerTest::step(density_mod_int))

# Check assumptions
plot_model(step_density_int, type = "diag") # all good

# Rename model object
density_evo_model <- step_density_int

## Compare with null model ####
density_mod_null <- lmer(dens_scam_live ~ weight_init + date_cloned_grp + origin_lab +
                           (co2 + salinity + elevation_sc)^3 + I(elevation_sc^2) +
                           (1|site_frame), data = traits_nocomp)

density_null_model <- get_model(step(density_mod_null))

# Compare models
MuMIn::r.squaredGLMM(density_null_model)
# 0.19
MuMIn::r.squaredGLMM(density_evo_model)
# 0.44

## Calculate effect sizes of fixed effects ####
epsilon_squared(density_evo_model)

# Pull significant terms from anova table 
tibble(anova(density_evo_model)) %>% 
  mutate(term = rownames(anova(density_evo_model))) %>% 
  filter(`Pr(>F)` < 0.05 & term != "weight_init" & term!= "date_cloned_grp") %>% 
  pull(term) -> significant_density_effects

tibble(epsilon_squared((density_evo_model))) %>% 
  filter(Parameter %in% significant_density_effects) %>% 
  mutate(category = c("eco", "eco", "eco")) %>% 
  group_by(category) %>% 
  summarize(total_epsilon2 = sum(Epsilon2_partial)) -> density_effectsize

#############################


###########################
## Calcs for all models ###
###########################
## Get final models for each trait ####

list(
`AGB EVO` = agb_evo_model,
`AGB NULL` = agb_null_model,
`BG EVO` = bg_evo_model,
`BG NULL` = bg_null_model,
`RS EVO` = rs_evo_model,
`RS NULL` = rs_null_model,
`BETA EVO` = beta_evo_model,
`BETA NULL` = beta_null_model,
`HEIGHT EVO` = height_evo_model,
`HEIGHT NULL` = height_null_model,
`WIDTH EVO` = width_evo_model,
`WIDTH NULL` = width_null_model,
`DENSITY EVO` = density_evo_model,
`DENSITY NULL` = density_null_model
) -> all_models

# Create functional to calculate and extract marginal and conditional ICC values
calc_icc <- function(x){
  as.numeric(performance::icc(x)[1])
}

calc_icc_conditional <- function(x){
  as.numeric(performance::icc(x)[2])
}

# Calculate R2 values for each model
lapply(all_models, MuMIn::r.squaredGLMM)

# Calculate ICC for all evolution models
lapply(all_models, calc_icc)

# Get ANOVA table for each model
lapply(all_models, anova)

## Calculate how important genotype is in the model ####

# For models that do not have random slopes, we just need the ICC (density, width)
calc_icc_conditional(density_evo_model)
# 0.2642793
calc_icc_conditional(width_evo_model)
# 0.1920142

# For models with random slopes, create a model with just a random intercept and
# then calculate ICC. 
agb_for_genotype <- update(agb_evo_model, .~.-(co2|genotype) + (1|genotype))
calc_icc_conditional(agb_for_genotype)
# 0.1015913

rs_for_genotype <- update(rs_evo_model, .~.-(salinity|genotype) + (1|genotype))
calc_icc_conditional(rs_for_genotype)
# 0.09068409

beta_for_genotype <- update(beta_evo_model, .~.-(1+salinity*elevation_sc|genotype) -salinity:elevation_sc + (1|genotype))
calc_icc_conditional(beta_for_genotype)
# 0.09884053

bg_for_genotype <- update(bg_evo_model, .~.-(co2|genotype) + (1|genotype))
calc_icc_conditional(bg_for_genotype)
# 0.1010572

height_for_genotype <- update(height_evo_model, .~.-(1+salinity*elevation_sc|genotype) -salinity:elevation_sc + (1|genotype))
calc_icc_conditional(height_for_genotype)
# 0.2202527

## Look at effect sizes ####
agb_effectsize
rs_effectsize
beta_effectsize
bg_effectsize
height_effectsize
width_effectsize
density_effectsize

# Create a color scale for the table to fill in cells based on effect size by
# categories
png(here("supp_data/effect_size_colors.png"), height = 8, width = 8, res = 300, units = "in")
tibble(category = c(agb_effectsize$category,
                    rs_effectsize$category,
                    beta_effectsize$category,
                    bg_effectsize$category,
                    height_effectsize$category,
                    width_effectsize$category,
                    density_effectsize$category),
       ef = c(agb_effectsize$total_epsilon2,
              rs_effectsize$total_epsilon2,
              beta_effectsize$total_epsilon2,
              bg_effectsize$total_epsilon2,
              height_effectsize$total_epsilon2,
              width_effectsize$total_epsilon2,
              density_effectsize$total_epsilon2)) %>% 
  ggplot(aes(x = log(ef), y = log(ef), color = ef)) +
  geom_text(aes(label = round(ef,3))) +
  geom_point(size = 3, aes(x = log(ef)+0.2, y = log(ef))) +
  scale_color_gradient(low = "gray95", high = "gray15")
dev.off()
####################################
## Predicted means for text ########
####################################

## Aboveground biomass ####
# Calculate mean elevation at level 4
traits_nocomp %>% 
  filter(level == 4) %>% 
  pull(elevation_sc) %>% mean() -> level4_elevation

# Get predicted means at high level of sea-level rise for agb
summary(emmeans(agb_evo_model, ~co2:salinity:age:elevation_sc,
                at = list(elevation_sc = level4_elevation), type = "response",
                weights = "proportional"))-> agb_pred_means

# Collect predicted means for co2 treatments from fresh & modern
agb_pred_means %>% 
  filter(salinity == "fresh" & age == "modern") %>% 
  pull(response) -> fresh_modern

# Calculate percent increase due to elevated CO2
fresh_modern[2]/fresh_modern[1]
# 1.901671

# Collect predicted means for co2 treatments from fresh & modern
agb_pred_means %>% 
  filter(salinity == "salt" & age == "ancestral") %>% 
  pull(response) -> salt_ancestral

# Calculate percent increase
salt_ancestral[2]/salt_ancestral[1]
# 1.884215

## Root-to-shoot ratio ####

# Calculate the effect of CO2 for descendant root-to-shoot ratio at low
# elevation
summary(emmeans(rs_evo_model, ~age:co2,
                at = list(elevation_sc = level4_elevation), type = "response",
                weights = "proportional")) -> rs_pred_means1

rs_pred_means1$response[4]/rs_pred_means1$response[2]
# 1.239753

# Calculate mean elevation at level 1
traits_nocomp_rs %>% 
  filter(level == 1) %>% 
  pull(elevation_sc) %>% mean() -> level1_elevation

# Calculate the effect of age at high elevation
summary(emmeans(rs_evo_model, ~age:co2,
                at = list(elevation_sc = level1_elevation), type = "response",
                weights = "proportional")) -> rs_pred_means2

rs_pred_means2$response[3]/rs_pred_means2$response[1]
# 1.187871

# Calculate the effect of salinity at low elevation
summary(emmeans(rs_evo_model, ~salinity,
                at = list(elevation_sc = level4_elevation), type = "response",
                weights = "proportional")) -> rs_pred_means3

rs_pred_means3$response[2]/rs_pred_means3$response[1]
# 1.234271

## Stem height ####

summary(emmeans(height_evo_model, ~salinity|location,
                weights = "proportional"))$emmean -> height_pred_means

# Corn: fresh vs salt
height_pred_means[1]/height_pred_means[2]
# 1.087699

## Stem width ####
summary(emmeans(width_evo_model, ~salinity|location,
                weights = "proportional"))$emmean -> width_pred_means

# Corn: fresh vs salt
width_pred_means[1]/width_pred_means[2]
# 1.144501

## Belowground biomass ####

# Calculate effect of salinity at high elevation
summary(emmeans(bg_evo_model, ~salinity,
                at = list(elevation_sc = level1_elevation), type = "response",
                weights = "proportional")) -> bg_pred_means

bg_pred_means$response[1]/bg_pred_means$response[2]
#1.315405
