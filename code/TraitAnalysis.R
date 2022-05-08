# Read in libraries 
library(here); library(tidyverse); library(emmeans);
library(lme4); library(lmerTest); library(blme);
library(brglm); library(lmtest); library(sjPlot);
library(kableExtra); library(geomtextpath); library(patchwork)

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

bg_full %>% 
  filter(comp == 0 & location != "blackwater") %>% 
  filter(level < 5) %>% 
  filter(agb_scam > 0) -> traits_nocomp

# Count to see how many we have for each genotype in this data set
traits_nocomp %>% 
  group_by(genotype) %>% 
  summarize(n = n()) %>% 
  pull(n) %>% range()
# There are between 2 and 12 reps per genotype in this dataset

## Fit model and model selection (random effects first, then fixed effects) ####

# Most complex fixed effects model
agb_mod <- lmer(sqrt(agb_scam) ~ (age + location + co2 + salinity + elevation)^5 +
                  I(elevation^2) + (1 + co2*salinity + elevation*co2 + salinity*elevation|genotype) +
                  (1|site_frame), data = traits_nocomp)

# Use lmerTest::step for backwards model selection
step_agb <- get_model(lmerTest::step(agb_mod, alpha.random = 0.05, alpha.fixed = 0.05))

# Check to see if there are any convergence warnings
summary(step_agb)
# We have a convergence error that indicates one of the random effects is either
# correlated 1:1 or is estimated to be zero. In this case, it is the former.

# We could also fit this model using a different RE structure (interaction
# between genotype and co2 as intercept instead of random slope).

# This version of the model converges
step_agb_alt <- update(step_agb, .~.-(co2|genotype) + (1|co2:genotype) + (1|genotype))

# Now see if this model is better than an intercept only model at alpha = 0.05
step_agb_alt_drop <- update(step_agb_alt, .~.-(1|co2:genotype))
anova(step_agb_alt, step_agb_alt_drop, refit = F)
# This is not significantly different at the alpha = 0.10 level which is our
# criterion for dropping random effects above. Because (1) the first
# model didn't converge and (2) the alternative model was not significant, we
# will drop this random slope to be on the more conservative side of things.

# The final model for agb is just an intercept only model for genotype and the
# fixed effects as identified above.

agb_evo_model <- update(step_agb, .~.-(co2|genotype) + (1|genotype))

plot_model(agb_evo_model, type = "diag") # All looks good!

## Comparison to null model ####

# Now construct a null model that does not have age, location, or genotype in it
agb_mod_null <- lmer(sqrt(agb_scam) ~ weight_init + date_cloned_grp + origin_lab +
                  (co2 + salinity + elevation)^3 + I(elevation^2) +
                  (1|site_frame), data = traits_nocomp)

agb_null_model <- get_model(lmerTest::step(agb_mod_null))

# Compare R2 from models
MuMIn::r.squaredGLMM(agb_null_model)
# Conditional R2 = 0.35
MuMIn::r.squaredGLMM(agb_evo_model)
# Conditional R2 = 0.45

## Make plot of fixed effects ####

# Extract type III Anova table
anova(agb_evo_model)

## Calculate effect sizes for text ####
emmeans_agb <- summary(emmeans(agb_evo_model, ~co2|salinity:age,
                               at = list(elevation = 0.2), type = "response"))

# Calculate percent increase in aboveground biomass for elevated co2 vs ambient
# co2 for modern genotypes in freshwater
emmeans_agb %>% 
  filter(salinity == "fresh" & age == "modern") %>% 
  pull(response) -> means_fresh_modern
means_fresh_modern[2]/means_fresh_modern[1]

# Calculate percent increase in aboveground biomass for elevated co2 vs ambient
# co2 for ancestral genotypes in brackish
emmeans_agb %>% 
  filter(salinity == "salt" & age == "ancestral") %>% 
  pull(response) -> means_salt_ancestral
means_salt_ancestral[2]/means_salt_ancestral[1]



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

# We also need to drop two pots #165 and #176 where we needed to harvest rhizome
# to "save" the genotype for archiving.

traits_nocomp %>% 
  filter(rs < 6 & pot_no !=165 & pot_no !=176) -> traits_nocomp_rs

# There also appear to be some that are NAs
traits_nocomp %>% 
  filter(!complete.cases(rs))
# There were issues with data collection for pots 205 and 1700 (i.e. mislabeled
# or missing bags for sections of bgb) so these two are dropped from the
# analysis because we could not reconcile the issues.

# Calculate the total N at this point
nrow(traits_nocomp_rs)
# 224

## Fit model and model selection (fixed effects) ####

# Most complex model
rs_mod <- lmer(log(rs) ~ weight_init + date_cloned_grp + origin_lab +
                 (age + location + co2 + salinity + elevation)^5 + I(elevation^2) +
                 (1+salinity*elevation + co2*elevation + salinity*co2|genotype) +
                 (1|site_frame), data = traits_nocomp_rs)

# Get the best stepwise model
step_rs <- get_model(lmerTest::step(rs_mod))

# Check to see if there are any convergence errors
summary(step_rs) # all good

# Check assumptions
plot_model(step_rs, type = "diag") # all good

rs_evo_model <- step_rs

## Comparison to null model ####

# Now construct a null model that does not have age, location, or genotype in it
rs_modelnull <- lmer(log(rs) ~ weight_init + date_cloned_grp + origin_lab +
                       (co2 + salinity + elevation)^3 + I(elevation^2) +
                       (1|site_frame),
                     data = traits_nocomp_rs)

rs_null_model <- get_model(lmerTest::step(rs_modelnull))


# Compare R2 from models
MuMIn::r.squaredGLMM(rs_null_model)
# Conditional R2 = 0.61
MuMIn::r.squaredGLMM(rs_evo_model)
# Conditional R2 = 0.70

# Get overall significance of model terms
anova(rs_evo_model)

##########################################
## Root distribution parameter analysis ##
##########################################

## Data manipulation ####

# This can be the same as the root-to-shoot data

## Fit model and model selection (fixed effects) ####

# Most complex model
beta_mod <- lmer(beta ~ weight_init + date_cloned_grp + origin_lab +
                   (age + location + co2 + salinity + elevation)^5 + I(elevation^2) +
                   (1|site_frame) + (1 + salinity*elevation + co2*elevation + salinity*co2|genotype), data = traits_nocomp_rs)

# Try lmerTest::step
step_beta <- get_model(lmerTest::step(beta_mod))

# Check for convergence warnings
summary(step_beta)
# There is a convergence issue which it looks like is due to site_frame
# intercept being estimated at 0

# Drop that intercept and refit
step_beta_noFrame <- update(step_beta, .~.-(1|site_frame))
# Still convergence issues. So try adding in interaction terms that are used in
# the random slopes into the fixed effect.
step_beta_int <- update(step_beta, .~.+ co2:elevation +salinity:elevation)

# This fits! And this makes sense to keep because we want to preserve
# relationships that are in the random slopes and it follows with the order from
# random to fixed in model selection.

summary(step_beta_int)
# Now no convergence issues

# Check assumptions
plot_model(step_beta_int, type = "diag") # all good

beta_evo_model <- step_beta_int

## Comparison to null model ####

# Now construct a null model that does not have age, location, or genotype in it
beta_modelnull <- lmer(beta ~ weight_init + date_cloned_grp + origin_lab +
                         (co2 + salinity + elevation)^3 + I(elevation^2) +
                         (1|site_frame), data = traits_nocomp_rs)


beta_null_model <- get_model(lmerTest::step(beta_modelnull))

# Compare R2 from models
MuMIn::r.squaredGLMM(beta_null_model)
# Conditional R2 = 0.48
MuMIn::r.squaredGLMM(beta_evo_model)
# Conditional R2 = 0.75

# Create predicted vs observed plots for each model

## Plot of the fixed effects for supplement ####
plot_model(beta_mod1, terms = c("elevation[all]", "age"), type = "emm")

# Get emmeans values for plots
beta_EA_plot <- emmeans::emmeans(beta_mod1, ~elevation:age, at = list(elevation = seq(0.156, 0.544, length.out = 50)))

summary(beta_EA_plot) %>% 
  mutate(age = case_when(age == "ancestral" ~ "ancestral cohort (1900-1950)",
                         T ~ "descendant cohort (2000-2020)")) %>% 
  mutate(beta = emmean) %>% 
  ggplot(aes(x = elevation, y = beta, color = age)) +
  geom_point(data = traits_nocomp_plot, aes(x = elevation, y = beta), shape = 1, stroke = 0.8, alpha = 0.7, size = 0.8) +
  geom_textline(aes(label = age), fontface = 2, linewidth = 1.2, size = 3, textcolor = "black") +
  geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL, fill = age), alpha = 0.2, color = NA) +
  ylab(expression(paste('root distribution parameter (', beta, ")"))) +
  xlab('elevation (m NAVD88)') +
  scale_color_manual(values = c("#fb9a99","#e31a1c")) +
  scale_fill_manual(values = c("#fb9a99","#e31a1c")) +
  theme_bw() +
  theme(legend.position = "none") -> beta_figA

# Also make figure where we translate root distribution parameter to rooting
# profile

# Create vector of depths
depths <- seq(0, 70, 1)
# Extract mean for ancestral and descendant at mean elevation
pred_beta_age <- summary(emmeans(beta_mod1, ~elevation:age, at = list(elevation = c(mean(traits_nocomp$elevation), min(traits_nocomp$elevation), max(traits_nocomp$elevation)))))$emmean
# Combine these in a tibble
tibble(age = rep(c("ancestral", "modern"), each = length(depths)*3),
       cum_frac = c(1-pred_beta_age[1]^depths, 1-pred_beta_age[2]^depths, 1-pred_beta_age[3]^depths,
                    1-pred_beta_age[4]^depths, 1-pred_beta_age[5]^depths, 1-pred_beta_age[6]^depths),
       depth = rep(depths,6),
       group = rep(c("mean_anc", "low_anc", "high_anc", "mean_mod", "low_mod", "high_mod"), each = length(depths)),
       elevation = rep(c("mean", "low", "high","mean", "low", "high"), each = length(depths))) %>%
  mutate(elevation = factor(elevation, levels = c("low", "mean", "high"))) %>% 
  ggplot(aes(x = depth, y = cum_frac, color = age, group = group)) +
  geom_line(aes(linetype = elevation), size = 0.7) +
  coord_flip() +
  scale_y_reverse() +
  scale_x_reverse() +
  scale_linetype_manual(values = c("dotted", "solid", "dashed")) +
  scale_color_manual(values = c("#fb9a99","#e31a1c")) +
  xlab("depth below marsh surface (cm)") +
  ylab("cumulative proportion") +
  theme_bw() -> beta_figB

beta_figA + beta_figB + plot_annotation(tag_levels = "a")-> beta_fig

ggsave(here("figs", "SuppFig_beta.png"), beta_fig, height = 4, width = 9, units = "in")
##################################
## Belowground biomass analysis ##
##################################

# For the belowground analysis, we can only really do the no competition pots in
# the straightforward way (b/c species were not separated for the entire core;
# only separated in the top 10cm). Use the same subset of data as for the
# root:shoot analysis for now.

## Fit model and model selection (fixed effects) ####

# Most complex model
bg_mod <- lmer(sqrt(total_bg) ~ weight_init + date_cloned_grp + origin_lab +
                 (age + location + co2 + salinity + elevation)^5 + I(elevation^2) +
                 (1+co2*elevation + salinity*elevation + co2*salinity|genotype) + (1|site_frame), data = traits_nocomp_rs)

# Do stepwise model selection
step_bg <- get_model(lmerTest::step(bg_mod))

# Check for convergence errors
summary(step_bg) # all good

# Check model assumptions
plot_model(step_bg, type = "diag") # all good

bg_evo_model <- step_bg

## Comparison to null model ####

# Now construct a null model that does not have age, location, or genotype in it
bg_modelnull <- lmer(total_bg ~ weight_init + date_cloned_grp + origin_lab +
                       (co2 + salinity + elevation)^3 + I(elevation^2) +
                       (1|site_frame),
                     data = traits_nocomp_rs)

bg_null_model <- get_model(step(bg_modelnull))

# Compare R2 from models
MuMIn::r.squaredGLMM(bg_null_model)
# Conditional R2 = 0.23
MuMIn::r.squaredGLMM(bg_evo_model)
# Conditional R2 = 0.39

##########################
## Stem height analysis ##
##########################

# Use the same as for the AGB analysis

##
# No competition
##

## Fit model and model selection (fixed effects) ####

# Most complex model
height_mod <- lmer(height_scam_tot ~ weight_init + date_cloned_grp + origin_lab +
                     (age + location + co2 + salinity + elevation)^5 + I(elevation^2) +
                     (1|site_frame) + (1+co2*elevation + salinity*elevation + co2*salinity|genotype), data = traits_nocomp)

step_height <- get_model(lmerTest::step(height_mod, alpha.random = 0.05))

plot(step_height) # There appears to be a huge outlier. 

# Figure out which data point this is
traits_nocomp %>% 
  filter(height_scam_tot < 15)
# There's a note that there was a bunch of caterpillar herbivory for this pot
# (top-down) so we should drop it for the height analysis

traits_nocomp %>% 
  filter(pot_no != 1830) -> traits_nocomp_height

height_mod_noOut <- lmer(height_scam_tot ~ weight_init + date_cloned_grp + origin_lab +
                           (age + location + co2 + salinity + elevation)^5 + I(elevation^2) +
                           (1+co2*elevation + salinity*elevation + co2*salinity|genotype) + (1|site_frame), data = traits_nocomp_height)

step_height_noOut <- get_model(lmerTest::step(height_mod_noOut, alpha.random = 0.05))

# Check for convergence issues
summary(step_height_noOut) # all good

# Check for assumptions
plot_model(step_height_noOut, type = "diag") # all good

height_evo_model <- step_height_noOut

## Comparison to null model ####

# Now construct a null model that does not have age, location, or genotype in it
height_model_noOut_null <- lmer(height_scam_tot ~ weight_init + date_cloned_grp + origin_lab +
                                (co2 + salinity + elevation)^3 + I(elevation^2) +
                                (1|site_frame),
                              data = traits_nocomp_height)

height_null_model <- get_model(step(height_model_noOut_null))

# Compare R2 from models
MuMIn::r.squaredGLMM(height_null_model)
# Adjusted R2 = 0.38
MuMIn::r.squaredGLMM(height_evo_model)
# Conditional R2 = 0.69

##########################
## Stem width analysis ##
##########################

# Use the same as for the AGB analysis

##
# No competition
##

## Fit model and model selection (fixed effects) ####

# Most complex model
width_mod <- lmer(width_scam_mid ~ weight_init + date_cloned_grp + origin_lab +
                    (age + location + co2 + salinity + elevation)^5 + I(elevation^2) +
                    (1+co2*elevation+salinity*elevation + co2*salinity|genotype) + (1|site_frame), data = traits_nocomp)

# Stepwise selection procedure
step_width <- get_model(lmerTest::step(width_mod, alpha.random = 0.05))

# Check for convergence errors
summary(step_width) # all good

# Check model assumptions
plot_model(step_width, type = "diag") # All looks good!

width_evo_model <- step_width

## Comparison to null model ####

# Now construct a null model that does not have age, location, or genotype in it
width_model_null <- lmer(width_scam_mid ~ weight_init + date_cloned_grp + origin_lab +
                         (co2 + salinity + elevation)^3 + I(elevation^2) +
                         (1|site_frame), data = traits_nocomp)

width_null_model <- get_model(step(width_model_null))

# Compare R2 from models
MuMIn::r.squaredGLMM(width_null_model)
# Conditional R2 = 0.42
MuMIn::r.squaredGLMM(width_evo_model)
# Conditional R2 = 0.63

###########################
## Stem density analysis ##
###########################

## Fit model ####
# Stem density uses the same data as AGB
density_mod <- lmer(dens_scam_live ~ weight_init + date_cloned_grp + origin_lab +
                      (age + location + co2 + salinity + elevation)^5 + I(elevation^2) +
                      (1+co2*elevation + salinity*elevation + salinity*co2|genotype) + (1|site_frame), data = traits_nocomp)

step_density <- get_model(lmerTest::step(density_mod, alpha.random = 0.05))

# Check for convergence errors
summary(step_density) # all good

# Check assumptions
plot_model(step_density, type = "diag") # all good

density_evo_model <- step_density

## Compare with null model ####
density_mod_null <- lmer(dens_scam_live ~ weight_init + date_cloned_grp + origin_lab +
                           (co2 + salinity + elevation)^3 + I(elevation^2) +
                           (1|site_frame), data = traits_nocomp)

density_null_model <- get_model(step(density_mod_null))

# Compare models
MuMIn::r.squaredGLMM(density_null_model)
# 0.19
MuMIn::r.squaredGLMM(density_evo_model)
# 0.44

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

# Calculate ICC for all evolution models
lapply(all_models, calc_icc)

# Get ANOVA table for each model
lapply(all_models, anova)

## Calculate how important genotype is in the model ####
step(agb_evo_model)
step(width_evo_model)
step(density_evo_model)

rs_for_genotype <- update(rs_evo_model, .~.-(salinity|genotype) + (1|genotype))
step(rs_for_genotype)

beta_for_genotype <- update(beta_evo_model, .~.-(salinity + elevation + co2 + salinity:elevation + elevation:co2|genotype) + (1|genotype))
step(beta_for_genotype)

bg_for_genotype <- update(bg_evo_model, .~.-(co2|genotype) + (1|genotype))
step(bg_for_genotype)

height_for_genotype <- update(height_evo_model, .~.-(elevation + salinity + elevation:salinity |genotype) + (1|genotype))
step(height_for_genotype)

## Same for random intercepts ####
step(rs_evo_model)
step(beta_evo_model)
step(bg_evo_model)
step(height_evo_model)

traits_nocomp_rs %>% 
  ggplot(aes(x = salinity, y = total_bg, color = co2)) +
  geom_jitter() + facet_wrap(~age)

###########################
## Effect sizes ###########
###########################

## Calculate effect sizes for in-text ####

# Calculate mean elevation at level 4
bg_nocomp %>% 
  filter(level == 4) %>% 
  pull(elevation) %>% mean() -> level4_elevation

# Get predicted means at high level of sea-level rise for agb
summary(emmeans(agb_evo_model, ~co2:salinity:age:elevation,
                at = list(elevation = level4_elevation), type = "response"))-> agb_pred_means

# Collect predicted means for co2 treatments from fresh & modern
agb_pred_means %>% 
  filter(salinity == "fresh" & age == "modern") %>% 
  pull(response) -> fresh_modern

# Calculate percent increase
fresh_modern[2]/fresh_modern[1]
# 1.892717

# Collect predicted means for co2 treatments from fresh & modern
agb_pred_means %>% 
  filter(salinity == "salt" & age == "ancestral") %>% 
  pull(response) -> salt_ancestral

# Calculate percent increase
salt_ancestral[2]/salt_ancestral[1]
# 1.770453

# Calculate mean elevation at level 1
bg_nocomp %>% 
  filter(level == 1) %>% 
  pull(elevation) %>% mean() -> level1_elevation

# Do the same thing for salinity at high levels of SLR
summary(emmeans(rs_evo_model, ~salinity:elevation,
                at = list(elevation = level4_elevation), type = "response")) -> rs_pred_means2

rs_pred_means2$response[2]/rs_pred_means2$response[1]
# 1.236114

# Width and height for provenance x salinity
summary(emmeans(height_evo_model, ~salinity|location))$emmean -> height_pred_means

# Corn: fresh vs salt
height_pred_means[1]/height_pred_means[2]
# 1.087679

summary(emmeans(width_evo_model, ~salinity|location))$emmean -> width_pred_means

# Corn: fresh vs salt
width_pred_means[1]/width_pred_means[2]
# 1.149909

