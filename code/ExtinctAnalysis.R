# Conduct binomial regression to assess the role of eco-evo factors in
# predicting the survival of S. americanus plants in the mesocosm experiment

# Read in libraries 
library(here); library(tidyverse); library(emmeans);
library(lme4); library(lmerTest); library(blme);
library(brglm); library(lmtest); library(sjPlot)

# Read in compiled trait data
source(here("supp_code", "CompileTraitData.R"))

#####################
## Binary analysis ##
#####################

## Data manipulation ####

# We are going to have to analyze the competition pots separately, so first look
# at no competition only. Also drop blackwater for simplicity.

bg_full %>% 
  filter(comp == 0 & location != "blackwater") -> full_data_nocomp

# Recode extinction values to be survival values to make the graph easier to
# interpret

full_data_nocomp %>% 
  mutate(survive = ifelse(extinct == 0, 1, 0)) -> full_data_nocomp

## Fit model ####

# Fit a logistic GLMM
extinct_mod_nocomp <- glmer(survive ~ weight_init + date_cloned_grp + origin_lab +
                            (co2 + age + location + elevation + salinity)^5 +
                            I(elevation^2) + (1|genotype) + (1|site_frame),
                          data = full_data_nocomp, family = "binomial")

# Random effect estimated at zero
VarCorr(extinct_mod_nocomp)
# Looks like both random effects are estimated to be zero

# Now fit a glm without random effects
extinct_mod_nocomp_fixed <- glm(survive ~ weight_init + date_cloned_grp + origin_lab +
                                  (co2 + salinity + elevation + age + location)^5 +
                                  I(elevation^2), data = full_data_nocomp, family = "binomial")

# We will deal with this warning later...
#glm.fit: fitted probabilities numerically 0 or 1 occurred 

# Check significance of model terms
car::Anova(extinct_mod_nocomp_fixed)

# No 4-way or 5-way interactions are significant so drop to 3-way model
extinct_mod_nocomp_fixed3 <- glm(survive ~ weight_init + origin_lab + date_cloned_grp +
                                   (co2 + salinity + elevation + age + location)^3 +
                                   I(elevation^2), data = full_data_nocomp, family = "binomial")

# Check significance of terms
car::Anova(extinct_mod_nocomp_fixed3)

# There is some complete separation (i.e. all reps with the same covariate
# combinations have the same predicted value) likely at low elevations where
# everything went extinct. So we should fit this with a bias reduction method
# using the package 'brglm'

# Scale covariates for easier interpretation of model coefficients
full_data_nocomp %>% 
  mutate(elevation_sc = scale(elevation)[,1],
         elevation2_sc = scale(elevation^2)[,1],
         weight_init_sc = scale(weight_init)[,1]) -> full_data_nocomp

# Set contrasts to help with interpretation of main effect coefficients in the
# presence of interactions
options(contrasts = c("contr.sum", "contr.poly"))

# Fit model
extinct_mod_nocomp_fixed3_BR <- brglm(survive ~ weight_init_sc + origin_lab + date_cloned_grp +
                                        (co2 + salinity + elevation_sc + age + location)^3 +
                                        I(elevation_sc^2), data = full_data_nocomp,
                                      family = binomial(logit))

# Fit the following model for plotting Fig S1 (same model with unscaled covariates)
# extinct_mod_nocomp_fixed3_BR <- brglm(survive ~ weight_init + origin_lab + date_cloned_grp +
#                                         (co2 + salinity + elevation + age + location)^3 +
#                                         I(elevation^2), data = full_data_nocomp,
#                                       family = binomial(logit))

# Check to see if we are missing something major with excluding genotype
extinct_mod_BR_genotype <- brglm(survive ~ weight_init_sc + origin_lab + date_cloned_grp +
                                   (co2 + salinity + elevation_sc)^3 + genotype +
                                   I(elevation_sc^2), data = full_data_nocomp,
                                 family = binomial(logit))

emmeans::emmeans(extinct_mod_BR_genotype, ~ genotype, type = "response")
# Huge confidence intervals around genotype -- not strong patterns here

# Compare to a model without genotype
extinct_mod_BR_nogenotype <- brglm(survive ~ weight_init_sc + origin_lab + date_cloned_grp +
                                     (co2 + salinity + elevation_sc)^3 +
                                     I(elevation_sc^2), data = full_data_nocomp,
                                   family = binomial(logit))

lrtest(extinct_mod_BR_genotype, extinct_mod_BR_nogenotype)
# Genotype model is not better than the other model

## Does competition mediate this response? ####

# We can only look at this for Corn genotypes
bg_full %>% 
  filter(location == "corn") %>% 
  mutate(survive = case_when(agb_scam > 0 ~ 1,
                             T ~ 0)) %>% 
  mutate(age = case_when(grepl("ancestral", cohort) ~ "ancestral",
                         T ~ "modern"),
         elevation_sc = scale(elevation)[,1],
         weight_init_sc = scale(weight_init)) -> corn_only

# Fit a logistic GLMM
extinct_mod_corn <- glmer(survive ~ weight_init_sc + date_cloned_grp + origin_lab + 
                            (co2 + salinity + elevation_sc + age + comp)^5 + I(elevation_sc^2) +
                            (1|site_frame) + (1|genotype), data = corn_only, family = "binomial")

VarCorr(extinct_mod_corn)
# random effects estimated to be zero

# Fit up to 3 way interactions like the previous model. Also need to fit brglm
# to deal with complete separation

extinct_mod_corn_BR <- brglm(survive ~ weight_init_sc + date_cloned_grp + origin_lab + 
                               (co2 + salinity + elevation_sc + age + comp)^3 + I(elevation_sc^2),
                             data = corn_only, family = "binomial")


# Competition did not seem to drive or mediate the likelihood of extinction in
# our experiment (no significant interactions with competition and no
# significant main effect).
