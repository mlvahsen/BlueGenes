# Read in libraries 
library(here); library(tidyverse); library(emmeans);
library(lme4); library(lmerTest); library(blme);
library(brglm); library(lmtest); library(sjPlot)

# Read in compiled trait data
source(here("supp_code", "CompileTraitData.R"))

## Data subsetting ####

# For most analyses we are going to need to subset the data into the competition
# data set (all Corn Island -- these were the only genotypes that had
# competition; and the no competition pots). For the no competition pots, we are
# also dropping out the Blackwater genotypes because they can't really be
# separated into age and location cohorts. We are also going to focus in on
# levels 1-4 where most clones survived.

bg_corn <- bg_full %>% filter(location == "corn" & level < 5)
bg_nocomp <- bg_full %>% filter(comp == 0 & location != "blackwater" & level < 5)

# See what our replication is given this subsetting
bg_corn %>% 
  group_by(age, location) %>% 
  summarize(n = n())

bg_corn %>% 
  group_by(genotype) %>% 
  summarize(n = n()) %>% 
  arrange(n)
# 20 genotypes; 6-26 reps per genotype

bg_nocomp %>% 
  group_by(age, location) %>% 
  summarize(n = n()) %>% 
  arrange(n)

bg_nocomp %>% 
  group_by(genotype) %>% 
  summarize(n = n()) %>% 
  arrange(n) %>% print(n = Inf)
# 32 genotypes; 1 - 16 reps per genotype

#####################
## Binary analysis ##
#####################

## Data manipulation ####

# First, we want to see what drove plants to not survive (or fail to establish)
# in the experiment. First, we'll group together plants that just didn't survive
# and those that failed to establish in the first place.

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
extinct_mod_nocomp <- glm(survive ~ weight_init + date_planted_grp + origin_lab + (co2 + salinity + elevation + age + location)^3 +
                            I(elevation^2) + genotype, data = full_data_nocomp, family = "binomial")

# Random effect estimated at zero
summary(extinct_mod_nocomp)
# Looks like both are estimated to be zero

# Now fit a glm without random effects
extinct_mod_nocomp_fixed <- glm(survive ~ weight_init + date_planted_grp + origin_lab + (co2 + salinity + elevation + age + location)^5 +
                              I(elevation^2), data = full_data_nocomp, family = "binomial")

car::Anova(extinct_mod_nocomp_fixed)

# No 4-way or 5-way interactions are significant so drop to 3-way model
extinct_mod_nocomp_fixed3 <- glm(survive ~ weight_init + origin_lab + (co2 + salinity + elevation + age + location)^3 +
                                  I(elevation^2), data = full_data_nocomp, family = "binomial")

# Check significance of terms
car::Anova(extinct_mod_nocomp_fixed3)

# There is some complete separation (i.e. all reps with the same covariate
# combinations have the same predicted value) at low elevations where everything
# went extinct. So we should fit this with a bias reduction method using the
# package 'brglm'

extinct_mod_nocomp_fixed3_BR <- brglm(survive ~ weight_init + origin_lab +
                                        (co2 + salinity + elevation + age + location)^3 +
                                      I(elevation^2), data = full_data_nocomp,
                                    family = binomial(logit))

summary(extinct_mod_nocomp_fixed3_BR)

a <- plot_model(extinct_mod_nocomp_fixed3_BR, terms = c("elevation", "age", "salinity"), type = "emm")

plot_data <- tibble(survive = a$data$predicted,
                    elevation = a$data$x,
                    lower.ci = a$data$conf.low,
                    upper.ci = a$data$conf.high,
                    salinity = a$data$facet,
                    age = a$data$group)

plot_data %>% 
  ggplot(aes(x = elevation, y = survive, color = age)) +
  geom_line() +
  facet_wrap(~salinity) +
  geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci, fill = age), alpha = 0.2, color = NA) +
  geom_jitter(data = full_data_nocomp, aes(x = elevation, y = survive), height = 0.05, width = 0, alpha = 0.3)
  

b <- plot_model(extinct_mod_nocomp_fixed3_BR, terms = c("location", "age", "co2"), type = "emm")

plot_data_b <- tibble(survive = b$data$predicted,
                    location = c("kirkpatrick", "corn")[b$data$x],
                    lower.ci = b$data$conf.low,
                    upper.ci = b$data$conf.high,
                    co2 = b$data$facet,
                    age = b$data$group)

pd <- position_dodge(width = 0.2)

plot_data_b %>% 
  ggplot(aes(x = location, y = survive, color = age)) +
  geom_point(size = 2, position = pd) +
  geom_errorbar(aes(ymin = lower.ci, ymax = upper.ci), width = 0.2, position = pd) +
  facet_wrap(~co2) +
  geom_point(data = full_data_nocomp, aes(x = location, y = survive), alpha = 0.2,
             position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0.05))

c <- plot_model(extinct_mod_nocomp_fixed3_BR, terms = c("location", "salinity"), type = "emm")

plot_data_c <- tibble(survive = c$data$predicted,
                      location = c("kirkpatrick", "corn")[c$data$x],
                      lower.ci = c$data$conf.low,
                      upper.ci = c$data$conf.high,
                      salinity = c$data$group)

plot_data_c %>% 
  ggplot(aes(x = location, y = survive, color = salinity)) +
  geom_point(size = 2, position = pd) +
  geom_errorbar(aes(ymin = lower.ci, ymax = upper.ci), width = 0.2, position = pd) +
  geom_point(data = full_data_nocomp, aes(x = location, y = survive), alpha = 0.2,
             position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0.05))

# Do all significance tests as likelihood ratio tests following the principle of
# marginality

##
# 3-way interactions
##

# co2:salinity:age (ns)
extinct_mod_noCSA <- update(extinct_mod_nocomp_fixed3_BR, .~.-co2:salinity:age)
lrtest(extinct_mod_nocomp_fixed3_BR, extinct_mod_noCSA)

# co2:salinity:location (ns)
extinct_mod_noCSL <- update(extinct_mod_nocomp_fixed3_BR, .~.-co2:salinity:location)
lrtest(extinct_mod_nocomp_fixed3_BR, extinct_mod_noCSL)

# co2:salinity:elevation (ns)
extinct_mod_noCSE <- update(extinct_mod_nocomp_fixed3_BR, .~.-co2:salinity:elevation)
lrtest(extinct_mod_nocomp_fixed3_BR, extinct_mod_noCSE)

# co2:age:location (*)
extinct_mod_noCAL <- update(extinct_mod_nocomp_fixed3_BR, .~.-co2:age:location)
lrtest(extinct_mod_nocomp_fixed3_BR, extinct_mod_noCAL)

# co2:age:elevation (ns)
extinct_mod_noCAE <- update(extinct_mod_nocomp_fixed3_BR, .~.-co2:age:elevation)
lrtest(extinct_mod_nocomp_fixed3_BR, extinct_mod_noCAE)

# co2:location:elevation (ns)
extinct_mod_noCLE <- update(extinct_mod_nocomp_fixed3_BR, .~.-co2:location:elevation)
lrtest(extinct_mod_nocomp_fixed3_BR, extinct_mod_noCLE)

# salinity:age:location (ns)
extinct_mod_noSAL <- update(extinct_mod_nocomp_fixed3_BR, .~.-salinity:age:location)
lrtest(extinct_mod_nocomp_fixed3_BR, extinct_mod_noSAL)

# salinity:age:elevation (*)
extinct_mod_noSAE <- update(extinct_mod_nocomp_fixed3_BR, .~.-salinity:age:elevation)
lrtest(extinct_mod_nocomp_fixed3_BR, extinct_mod_noSAE)

# salinity:location:elevation (ns)
extinct_mod_noSLE <- update(extinct_mod_nocomp_fixed3_BR, .~.-salinity:location:elevation)
lrtest(extinct_mod_nocomp_fixed3_BR, extinct_mod_noSLE)

# age:location:elevation (ns)
extinct_mod_noALE <- update(extinct_mod_nocomp_fixed3_BR, .~.-age:location:elevation)
lrtest(extinct_mod_nocomp_fixed3_BR, extinct_mod_noALE)

##
# 2-way interactions
##

# co2:salinity (ns)
extinct_mod_forCS <- brglm(survive ~ weight_init + origin_lab +
                              (co2 + salinity + elevation + age + location)^2 +
                              (co2 + elevation + age + location)^3 +
                             (salinity + elevation + age + location)^3 +
                              I(elevation^2), data = full_data_nocomp,
                            family = binomial(logit))
extinct_mod_noCS <- update(extinct_mod_forCS, .~.-co2:salinity)
lrtest(extinct_mod_forCS, extinct_mod_noCS)

# co2:age (ns)
extinct_mod_forCA <- brglm(survive ~ weight_init + origin_lab +
                             (co2 + salinity + elevation + age + location)^2 +
                             (co2 + elevation + salinity + location)^3 +
                             (salinity + elevation + age + location)^3 +
                             I(elevation^2), data = full_data_nocomp,
                           family = binomial(logit))
extinct_mod_noCA <- update(extinct_mod_forCA, .~.-co2:age)
lrtest(extinct_mod_forCA, extinct_mod_noCA)

# co2:location (ns)
extinct_mod_forCL <- brglm(survive ~ weight_init + origin_lab +
                             (co2 + salinity + elevation + age + location)^2 +
                             (co2 + elevation + salinity + age)^3 +
                             (salinity + elevation + age + location)^3 +
                             I(elevation^2), data = full_data_nocomp,
                           family = binomial(logit))
extinct_mod_noCL <- update(extinct_mod_forCL, .~.-co2:location)
lrtest(extinct_mod_forCL, extinct_mod_noCL)

# co2:elevation (ns)
extinct_mod_forCE <- brglm(survive ~ weight_init + origin_lab +
                             (co2 + salinity + elevation + age + location)^2 +
                             (co2 + location + salinity + age)^3 +
                             (salinity + elevation + age + location)^3 +
                             I(elevation^2), data = full_data_nocomp,
                           family = binomial(logit))
extinct_mod_noCE <- update(extinct_mod_forCE, .~.-co2:elevation)
lrtest(extinct_mod_forCE, extinct_mod_noCE)

# salinity:age (*)
extinct_mod_forSA <- brglm(survive ~ weight_init + origin_lab +
                             (co2 + salinity + elevation + age + location)^2 +
                             (co2 + location + salinity + elevation)^3 +
                             (co2 + elevation + age + location)^3 +
                             I(elevation^2), data = full_data_nocomp,
                           family = binomial(logit))
extinct_mod_noSA <- update(extinct_mod_forSA, .~.-salinity:age)
lrtest(extinct_mod_forSA, extinct_mod_noSA)

# salinity:location (*)
extinct_mod_forSL <- brglm(survive ~ weight_init + origin_lab +
                             (co2 + salinity + elevation + age + location)^2 +
                             (co2 + location + age + elevation)^3 +
                             (co2 + elevation + age + salinity)^3 +
                             I(elevation^2), data = full_data_nocomp,
                           family = binomial(logit))
extinct_mod_noSL <- update(extinct_mod_forSL, .~.-salinity:location)
lrtest(extinct_mod_forSL, extinct_mod_noSL)

# salinity:elevation (***)
extinct_mod_forSE <- brglm(survive ~ weight_init + origin_lab +
                             (co2 + salinity + elevation + age + location)^2 +
                             (co2 + location + salinity + age)^3 +
                             (co2 + elevation + age + location)^3 +
                             I(elevation^2), data = full_data_nocomp,
                           family = binomial(logit))
extinct_mod_noSE <- update(extinct_mod_forSE, .~.-salinity:elevation)
lrtest(extinct_mod_forSE, extinct_mod_noSE)

# age:location (ns)
extinct_mod_forAL <- brglm(survive ~ weight_init + origin_lab +
                             (co2 + salinity + elevation + age + location)^2 +
                             (co2 + elevation + salinity + age)^3 +
                             (co2 + elevation + salinity + location)^3 +
                             I(elevation^2), data = full_data_nocomp,
                           family = binomial(logit))
extinct_mod_noAL <- update(extinct_mod_forAL, .~.-age:location)
lrtest(extinct_mod_forAL, extinct_mod_noAL)

# age:elevation (*)
extinct_mod_forAE <- brglm(survive ~ weight_init + origin_lab +
                             (co2 + salinity + elevation + age + location)^2 +
                             (co2 + location + salinity + elevation)^3 +
                             (co2 + age + salinity + location)^3 +
                             I(elevation^2), data = full_data_nocomp,
                           family = binomial(logit))
extinct_mod_noAE <- update(extinct_mod_forAE, .~.-age:elevation)
lrtest(extinct_mod_forAE, extinct_mod_noAE)

# location:elevation (*)
extinct_mod_forLE <- brglm(survive ~ weight_init + origin_lab +
                             (co2 + salinity + elevation + age + location)^2 +
                             (co2 + age + salinity + elevation)^3 +
                             (co2 + age + salinity + location)^3 +
                             I(elevation^2), data = full_data_nocomp,
                           family = binomial(logit))
extinct_mod_noLE <- update(extinct_mod_forLE, .~.-location:elevation)
lrtest(extinct_mod_forLE, extinct_mod_noLE)

##
# Main terms
##

# co2 (*)
extinct_mod_forC <- brglm(survive ~ weight_init + origin_lab +
                             (salinity + elevation + age + location)^3 + co2 +
                             I(elevation^2), data = full_data_nocomp,
                           family = binomial(logit))
extinct_mod_noC <- update(extinct_mod_forC, .~.-co2)
lrtest(extinct_mod_forC, extinct_mod_noC)

# salinity (ns)
extinct_mod_forS <- brglm(survive ~ weight_init + origin_lab +
                            (co2 + elevation + age + location)^3 + salinity +
                            I(elevation^2), data = full_data_nocomp,
                          family = binomial(logit))
extinct_mod_noS <- update(extinct_mod_forS, .~.-salinity)
lrtest(extinct_mod_forS, extinct_mod_noS)

# age (ns)
extinct_mod_forA <- brglm(survive ~ weight_init + origin_lab +
                            (co2 + elevation + salinity + location)^3 + age +
                            I(elevation^2), data = full_data_nocomp,
                          family = binomial(logit))
extinct_mod_noA <- update(extinct_mod_forA, .~.-age)
lrtest(extinct_mod_forA, extinct_mod_noA)

# location (ns)
extinct_mod_forL <- brglm(survive ~ weight_init + origin_lab +
                            (co2 + elevation + salinity + age)^3 + location +
                            I(elevation^2), data = full_data_nocomp,
                          family = binomial(logit))
extinct_mod_noL <- update(extinct_mod_forL, .~.-location)
lrtest(extinct_mod_forL, extinct_mod_noL)

# elevation (***)
extinct_mod_forE <- brglm(survive ~ weight_init + origin_lab +
                            (co2 + location + salinity + age)^3 + elevation +
                            I(elevation^2), data = full_data_nocomp,
                          family = binomial(logit))
extinct_mod_noE <- update(extinct_mod_forE, .~.-elevation)
lrtest(extinct_mod_forE, extinct_mod_noE)

##
# Covariate terms
##

# weight_init (***)
extinct_mod_noW <- update(extinct_mod_nocomp_fixed3_BR, .~.-weight_init)
lrtest(extinct_mod_nocomp_fixed3_BR, extinct_mod_noW)

# origin_lab (***)
extinct_mod_noO <- update(extinct_mod_nocomp_fixed3_BR, .~.-origin_lab)
lrtest(extinct_mod_nocomp_fixed3_BR, extinct_mod_noO)

# elevation^2 (***)
extinct_mod_noE2 <- update(extinct_mod_nocomp_fixed3_BR, .~.-I(elevation^2))
lrtest(extinct_mod_nocomp_fixed3_BR, extinct_mod_noE2)


##
# Make plots
##

plot_colors <- c("#fc8d62", "#8da0cb")
plot_shapes <- c(15, 17)

# A = elevation:age:salinity interaction
a <- plot_model(extinct_mod_nocomp_fixed3_BR, terms = c("elevation", "age", "salinity"), type = "emm")

plot_data_a <- tibble(survive = a$data$predicted,
                    elevation = a$data$x,
                    lower.ci = a$data$conf.low,
                    upper.ci = a$data$conf.high,
                    salinity = a$data$facet,
                    age = a$data$group)

plot_raw_data <- full_data_nocomp %>% 
  mutate(salinity = ifelse(salinity == "fresh", "freshwater site (6ppt)", "brackish site (8ppt)"))

plot_data_a %>% 
  mutate(salinity = ifelse(salinity == "fresh", "freshwater site (6ppt)", "brackish site (8ppt)")) %>% 
  ggplot(aes(x = elevation, y = survive, color = age)) +
  geom_line() +
  facet_wrap(~salinity) +
  geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci, fill = age), alpha = 0.2, color = NA) +
  geom_jitter(data = plot_raw_data, aes(x = elevation, y = survive), height = 0.05, width = 0, alpha = 0.3) +
  ylab("survival rate") + xlab("elevation (m NAVD88)") +
  scale_fill_manual(values = plot_colors, labels = c("ancestral (1900-1960)", "modern (2000-2020)")) +
  scale_color_manual(values = plot_colors, labels = c("ancestral (1900-1960)", "modern (2000-2020)")) +
  theme_bw() +
  theme(legend.position = "right", legend.background = element_rect(color = "gray47")) +
  guides(fill=guide_legend(title="age cohort"),
         color = guide_legend(title="age cohort"))-> plot_a

# B = elevation:location:salinity interaction
b <- plot_model(extinct_mod_nocomp_fixed3_BR, terms = c("elevation", "location", "salinity"), type = "emm")

plot_data_b <- tibble(survive = b$data$predicted,
                      elevation = b$data$x,
                      lower.ci = b$data$conf.low,
                      upper.ci = b$data$conf.high,
                      salinity = b$data$facet,
                      location = b$data$group)

plot_data_b %>% 
  mutate(location = factor(location, levels = c("corn", "kirkpatrick"))) %>% 
  mutate(salinity = ifelse(salinity == "fresh", "freshwater site (6ppt)", "brackish site (8ppt)")) %>% 
  ggplot(aes(x = elevation, y = survive, color = location, shape = location)) +
  geom_line() +
  facet_wrap(~salinity) +
  geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci, fill = location), alpha = 0.2, color = NA) +
  geom_jitter(data = plot_raw_data, aes(x = elevation, y = survive, shape = location),
              height = 0.05, width = 0, alpha = 0.7) +
  scale_color_manual(values = c("gray57", "black")) +
  scale_fill_manual(values = c("gray57", "black")) +
  scale_shape_manual(values = plot_shapes) +
  guides(fill=guide_legend(title="provenance"),
         color = guide_legend(title="provenance"),
         shape = guide_legend(title="provenance"))  +
  theme_bw() + theme(legend.background = element_rect(color = "gray47"))+
  labs(x = "elevation (m NAVD88)", y = "survival rate") -> plot_b
  

# C = location:age:co2 interaction
c <- plot_model(extinct_mod_nocomp_fixed3_BR, terms = c("location", "age", "co2"), type = "emm")

plot_data_c <- tibble(survive = c$data$predicted,
                      location = c("kirkpatrick", "corn")[c$data$x],
                      lower.ci = c$data$conf.low,
                      upper.ci = c$data$conf.high,
                      co2 = c$data$facet,
                      age = c$data$group)

pd <- position_dodge(width = 0.4)

my_labels <- c(ambient = "ambient~CO[2]", elevated = "elevated~CO[2]")
my_labeller <- as_labeller(my_labels,
                           default = label_parsed)

plot_data_c %>% 
  ggplot(aes(x = location, y = survive, color = age, shape = location)) +
  geom_point(size = 2, position = pd) +
  geom_errorbar(aes(ymin = lower.ci, ymax = upper.ci), width = 0.2, position = pd) +
  facet_wrap(~co2, labeller = my_labeller) +
  #geom_point(data = full_data_nocomp, aes(x = location, y = survive), alpha = 0.2,
             #position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0.05, dodge.width = 1)) +
  theme_bw() +
  theme(legend.position = "none") +
  guides(colour = "none",
         shape = "none") +
  labs(x = "provenance", y = "survival rate") +
  scale_shape_manual(values = plot_shapes) +
  scale_color_manual(values = plot_colors) -> plot_c

png("figs/Fig1_Extinction.png", res = 300, units = "in", height = 5, width = 8.2)
plot_a + plot_c + plot_b + guide_area() + plot_annotation(tag_levels = 'a')+
  plot_layout(guides = "collect", widths = c(4,3)) & theme(legend.justification = "left")
dev.off()

## Does competition mediate this response?

# We can only look at this for Corn genotypes
full_data_merged %>% 
  filter(location == "corn") %>% 
  mutate(extinct = case_when(agb_scam > 0 ~ 0,
                             T ~ 1)) %>% 
  mutate(age = case_when(grepl("ancestral", cohort) ~ "ancestral",
                         T ~ "modern")) -> corn_only

# Fit a logistic GLMM
extinct_mod_corn <- glm(extinct ~ weight_init + date_planted_grp + origin_lab + (salinity + elevation)^2 +
                          I(elevation^2) + co2 + age + comp, data = corn_only, family = "binomial")

# See much fewer significant interactions in this model

# Check to see if salinity by elevation interaction is in the same direction as
# the other model

plot_model(extinct_mod_corn, terms = c("elevation[all]", "salinity"), type = "emm")
# Yep!

# Competition did not seem to drive or mediate the likelihood of extinction in
# our experiment (no siginificant interactions with competition and no
# significant main effect). Careful with this because it's "confirming the
# null".


## Zeros that are just due to failing to establish ####

# We can see if there are different factors that predict whether or not plants
# will just fail to establish, not establish and then die.

# First, create an indicator for this in the no comp data set. This should be
# the pots that are NOT in the blue_genes dataset.

full_data_nocomp %>% 
  mutate(failed_establish = case_when(pot_no %in% blue_genes$pot_no ~ 0,
                                      T ~ 1)) -> full_data_nocomp
# Fit the same model as previously


establish_mod_nocomp <- glm(failed_establish ~ weight_init + origin_lab + (co2 + salinity + elevation + age + location)^2 +
                              I(elevation^2) + salinity:elevation:age + salinity:elevation:location, data = full_data_nocomp, family = "binomial")

car::Anova(establish_mod_nocomp)

# Looks like different interactions are significant so we'll need to think about
# what makes the most sense (and is the most justifiable) and go with that.

# Look at interactions
plot_model(establish_mod_nocomp, terms = c("elevation[all]", "salinity", "age"), type = "emm")
# Modern genotypes basically all failed to establish at low elevations at
# freshwater, but did well in brackish conditions. This pattern existed for
# ancestral genotypes but the differences between freshwater and brackish
# conditions were smaller.

plot_model(establish_mod_nocomp, terms = c("elevation[all]", "salinity", "location"), type = "emm")
# More deaths at freshwater overall, but for Kirk genotypes, there were
# basically no deaths at higher elevations at freshwater but there were some in
# more brackish conditions.

plot_model(establish_mod_nocomp, terms = c("co2"), type = "emm")
# Plants were more likely to fail to establish when in elevated vs ambient
# conditions

# Do the same for the corn dataset (i.e. with competition)
corn_only %>% 
  mutate(failed_establish = case_when(pot_no %in% blue_genes$pot_no ~ 0,
                                      T ~ 1)) -> corn_only 


establish_mod_corn <- glm(failed_establish ~ weight_init + origin_lab + (co2 + salinity + elevation + age + comp)^2 +
                            I(elevation^2), data = corn_only, family = "binomial")

car::Anova(establish_mod_corn)

# It looks like competition mediates it now
plot_model(establish_mod_corn, terms = c("age", "comp"), type = "emm")
# Only modern genotypes had an increased likelihood of failing to establish when
# in the presence of competition.

plot_model(establish_mod_corn, terms = c("elevation[all]", "salinity"), type = "emm")
# Basically all failed to establish at low levels at freshwater, but differences
# between the establishment rates between sites decreases with increasing
# elevation.

plot_model(establish_mod_corn, terms = c("co2", "comp"), type = "emm")
# Only in ambient conditions was there an increased likelihood of failing to
# establish when in the presence of competition.


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

full_data_nocomp %>% 
  filter(level < 5) %>% 
  filter(agb_scam > 0) %>% 
  # Also scale elevation to help with convergence
  mutate(elevation_sc = scale(elevation)[,1]) -> traits_nocomp

# Mean adjustment = 0.2464553
# Sd adjustment = 0.18788

# Count to see how many we have for each genotype in this data set
traits_nocomp %>% 
  group_by(genotype) %>% 
  summarize(n = n()) %>% 
  pull(n) %>% range()
# There are between 2 and 12 reps per genotype in this dataset

# Start with aboveground biomass analysis to just get a hang of how we are going
# to do this. We should write this model selection protocol into a function
# later...

## Fit model and model selection (fixed effects) ####

# Most complex model
agb_mod <- lmer(sqrt(agb_scam) ~ weight_init + date_planted_grp + origin_lab +
                  (age + location + co2 + salinity + elevation_sc)^5 + I(elevation_sc^2) +
                  (1|genotype), data = traits_nocomp)

# Try lmerTest::step for backwards model selection
step_agb <- lmerTest::step(agb_mod)

# This is the best model
agb_mod_step <- lmer(sqrt(agb_scam) ~ weight_init + I(elevation_sc^2) + 
                       age*co2*salinity*elevation_sc + (1 | genotype),
                     data = traits_nocomp)

emmeans::emmeans(agb_mod_step, ~elevation_sc:co2:age:salinity,
                 at = list(elevation_sc = seq(min(traits_nocomp$elevation_sc),
                                              max(traits_nocomp$elevation_sc),
                                              length.out = 10)),
                 type = "response")

plot_model(agb_mod_step, type = "pred", terms = c("elevation_sc",
                                                  "co2", "age", "salinity"))

summary(agb_mod)
# Site frame does not explain variation. This makes sense because we had the
# actual elevations of each chamber and the largest differences would be due to
# elevation.

# Take out (1|site_frame)
agb_mod2 <- update(agb_mod, .~.-(1|site_frame))
# Look for highest order significant interaction
car::Anova(agb_mod2)
# This would be age:co2:salinity:elevation (4-way interaction)

# Create a model that keeps all 3-way interactions as well as that 4-way
# interaction
agb_mod3 <- lmer(sqrt(agb_scam) ~ weight_init + date_planted_grp + origin_lab +
                   (age + location + co2 + salinity + elevation_sc)^3 + I(elevation_sc^2) +
                   co2:age:salinity:elevation_sc + (1|genotype), data = traits_nocomp)
# Look for highest order significant terms that are NOT contained in that
# interaction
car::Anova(agb_mod3)
# There are no other interactions that are not nested within that 4-way
# interaction. So make a model with main terms and that interaction.
agb_mod4 <- lmer(sqrt(agb_scam) ~ weight_init + date_planted_grp + origin_lab +
                   location + I(elevation_sc^2) + co2*age*salinity*elevation_sc +
                   (1|genotype), data = traits_nocomp)

# This is the final FIXED effects model. Check assumptions.
plot_model(agb_mod, type = "diag") # All looks good!

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

agb_mod5 <- update(agb_mod4, .~.-(1|genotype) + (1+co2*elevation|genotype))
agb_mod6 <- update(agb_mod4, .~.-(1|genotype) + (1+co2*salinity|genotype))
agb_mod7 <- update(agb_mod4, .~.-(1|genotype) + (1+elevation*salinity|genotype))
agb_mod8 <- update(agb_mod4, .~.-(1|genotype) + (1+co2+elevation|genotype))
agb_mod9 <- update(agb_mod4, .~.-(1|genotype) + (1+elevation + salinity|genotype))
agb_mod10 <- update(agb_mod4, .~.-(1|genotype) + (1+salinity + co2|genotype))
agb_mod11 <- update(agb_mod4, .~.-(1|genotype) + (1+co2|genotype))
agb_mod12<- update(agb_mod4, .~.-(1|genotype) + (1+elevation|genotype))
agb_mod13 <- update(agb_mod4, .~.-(1|genotype) + (1+salinity|genotype))

# The only one that did not have estimated zero was the co2*elevation model.
# Compare this to the base model.
MuMIn::r.squaredGLMM(agb_mod4)
MuMIn::r.squaredGLMM(agb_mod5)

plot_model(agb_mod5, type = "pred", pred.type = "re",
           terms = c("elevation_sc[all]", "co2", "genotype"))
# This doesn't seem meaningful so it seems like the best (and most parsimonious
# model) does not have a random slope.

## Comparison to null model ####

# Now construct a null model that does not have age, location, or genotype in it
agb_modelnull <- lm(sqrt(agb_scam) ~ weight_init + date_planted_grp + origin_lab +
                      (co2 + salinity + elevation_sc)^3 + I(elevation_sc^2),
                    data = traits_nocomp)

car::Anova(agb_modelnull)
# If no age, location, or genotype, we would say that elevation and co2
# influence biomass
plot_model(agb_modelnull, terms = "elevation_sc", type = "emm")
plot_model(agb_modelnull, terms = "co2", type = "emm")

# Compare R2 from models
summary(agb_modelnull)
# Adjusted R2 = 0.32
MuMIn::r.squaredGLMM(agb_mod4)
# Conditional R2 = 0.46

# Create predicted vs observed plots for each model
tibble(predicted = c(predict(agb_modelnull)^2, predict(agb_mod4)^2),
       observed = rep(traits_nocomp$agb_scam, 2),
       model = c(rep(c("null", "evolution"), each = nrow(traits_nocomp)))) %>% 
  ggplot(aes(x = predicted, y = observed)) +
  geom_point(alpha = 0.3) +
  facet_wrap(~model) +
  ylim(0,16) + xlim(0,16) +
  geom_abline(aes(slope = 1, intercept = 0))

# Repeat this process for the competition pots

##
# Corn only
##

## Plot of fixed effects ####
traits_nocomp %>% 
  mutate(salinity = case_when(salinity == "fresh" ~ "Freshwater Site (6ppt)", T ~ "Brackish Site (8ppt)")) %>% 
  mutate(age = case_when(age == "ancestral" ~ "Ancestral Cohort (1900-1960)", T ~ "Modern Cohort (2000-2020)")) %>% 
  mutate(co2 = case_when(co2 == "elevated" ~ "Elevated (700ppm)", T ~ "Ambient (400ppm)")) %>%
  ggplot(aes(x = elevation, y = agb_scam, color = co2)) +
  geom_point(size = 2, alpha = 0.5) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), se = F, size = 1.5) +
  facet_grid(salinity ~ age) +
  scale_color_manual(values = c("#78c679", "#006837"))+
  guides(color=guide_legend(title=expression(paste(CO[2], " level")))) +
  ylab("Aboveground Biomass (g)") +
  xlab("Elevation (m NAVD88)") +
  theme_bw() +
  theme(legend.position = "top")

## Plot of random effects ####
traits_nocomp %>% 
  ggplot(aes(x = reorder(genotype, agb_scam, FUN = median), y = agb_scam)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 1.5) +
  ylab("Aboveground Biomass (g)") +
  xlab("Genotype") +
  geom_hline(aes(yintercept = mean(agb_scam)), linetype = "dashed") +
  theme_bw() +
  theme(axis.text.x = element_blank())
#####

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
agb_mod_corn <- lmer(sqrt(agb_scam) ~ weight_init + date_planted_grp + origin_lab +
                       (age + comp + co2 + salinity + elevation_sc)^5 + I(elevation_sc^2) +
                       (1|genotype), data = traits_corn)

# Try lmerTest::step
lmerTest::step(agb_mod_corn) # Full model is the best

agb_corn_step_mod <- lmer(sqrt(agb_scam) ~ weight_init + I(elevation_sc^2) +
                            age*comp*salinity*co2*elevation_sc + (1|genotype), data = traits_corn)

summary(agb_mod)
# Site frame does not explain variation. This makes sense because we had the
# actual elevations of each chamber and the largest differences would be due to
# elevation.

# Take out (1|site_frame)
agb_mod_corn2 <- update(agb_mod_corn, .~.-(1|site_frame))
# Look for highest order significant interaction
car::Anova(agb_mod_corn2)
# This would be age:co2:salinity:elevation:comp (5-way interaction)

# The most complex model might be the best model in this case
agb_mod_corn3 <- lmer(sqrt(agb_scam) ~ weight_init + date_planted_grp + origin_lab +
                        (age + comp + co2 + salinity + elevation_sc)^4 + I(elevation_sc^2) +
                        age:co2:salinity:elevation_sc:comp + (1|genotype), data = traits_corn)

# This is the final FIXED effects model. Check assumptions.
plot_model(agb_mod_corn3, type = "diag") # All looks good!

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

agb_mod_corn4 <- update(agb_mod_corn3, .~.-(1|genotype) + (1+co2*elevation|genotype))
agb_mod_corn5 <- update(agb_mod_corn3, .~.-(1|genotype) + (1+co2*salinity|genotype))
agb_mod_corn6 <- update(agb_mod_corn3, .~.-(1|genotype) + (1+elevation*salinity|genotype))
agb_mod_corn7 <- update(agb_mod_corn3, .~.-(1|genotype) + (1+co2+elevation|genotype))
agb_mod_corn8 <- update(agb_mod_corn3, .~.-(1|genotype) + (1+elevation + salinity|genotype))
agb_mod_corn9 <- update(agb_mod_corn3, .~.-(1|genotype) + (1+salinity + co2|genotype))
agb_mod_corn10 <- update(agb_mod_corn3, .~.-(1|genotype) + (1+co2|genotype))
agb_mod_corn11<- update(agb_mod_corn3, .~.-(1|genotype) + (1+elevation|genotype))
agb_mod_corn12 <- update(agb_mod_corn3, .~.-(1|genotype) + (1+salinity|genotype))

# The only one that did not have convergence issues or 0 variance was the
# elevation random slope.
MuMIn::r.squaredGLMM(agb_mod_corn3)
MuMIn::r.squaredGLMM(agb_mod_corn11)

anova(agb_mod_corn3, agb_mod_corn11)

plot_model(agb_mod_corn11, type = "pred", pred.type = "re",
           terms = c("elevation_sc[all]", "genotype")) +
  scale_color_manual(values = rainbow(32))
# This doesn't seem meaningful so it seems like the best (and most parsimonious
# model) does not have a random slope.

## Comparison to null model ####

# Now construct a null model that does not have age, location, or genotype in it
agb_model_cornnull <- lm(sqrt(agb_scam) ~ weight_init + date_planted_grp + origin_lab +
                           (co2 + salinity + elevation_sc + comp)^4 + I(elevation_sc^2),
                         data = traits_corn)

car::Anova(agb_model_cornnull)
# If no age, location, or genotype, we would say that elevation x comp
# influences biomass
plot_model(agb_model_cornnull, terms = c("elevation_sc[all]", "comp"), type = "emm")

# Compare R2 from models
summary(agb_model_cornnull)
# Adjusted R2 = 0.43
MuMIn::r.squaredGLMM(agb_mod_corn3)
# Conditional R2 = 0.57

# Create predicted vs observed plots for each model
tibble(predicted = c(predict(agb_model_cornnull)^2, predict(agb_mod_corn3)^2),
       observed = rep(traits_corn$agb_scam, 2),
       model = c(rep(c("null", "evolution"), each = nrow(traits_corn)))) %>% 
  ggplot(aes(x = predicted, y = observed)) +
  geom_point(alpha = 0.3) +
  facet_wrap(~model) +
  ylim(0,16) + xlim(0,16) +
  geom_abline(aes(slope = 1, intercept = 0))


## Plot of fixed effects ####
traits_corn %>% 
  filter(comp == 0) %>% 
  mutate(salinity = case_when(salinity == "fresh" ~ "Freshwater Site (6ppt)", T ~ "Brackish Site (8ppt)")) %>% 
  mutate(age = case_when(age == "ancestral" ~ "Ancestral Cohort (1900-1960)", T ~ "Modern Cohort (2000-2020)")) %>% 
  mutate(co2 = case_when(co2 == "elevated" ~ "Elevated (700ppm)", T ~ "Ambient (400ppm)")) %>%
  ggplot(aes(x = elevation, y = agb_scam, color = co2)) +
  geom_point(size = 2, alpha = 0.5) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), se = F, size = 1.5) +
  facet_grid(salinity ~ age) +
  scale_color_manual(values = c("#78c679", "#006837"))+
  guides(color=guide_legend(title=expression(paste(CO[2], " level")))) +
  ylab("Aboveground Biomass (g)") +
  xlab("Elevation (m NAVD88)") +
  theme_bw() +
  theme(legend.position = "top") +
  ggtitle("WITHOUT competition")

traits_corn %>% 
  filter(comp == 1) %>% 
  mutate(salinity = case_when(salinity == "fresh" ~ "Freshwater Site (6ppt)", T ~ "Brackish Site (8ppt)")) %>% 
  mutate(age = case_when(age == "ancestral" ~ "Ancestral Cohort (1900-1960)", T ~ "Modern Cohort (2000-2020)")) %>% 
  mutate(co2 = case_when(co2 == "elevated" ~ "Elevated (700ppm)", T ~ "Ambient (400ppm)")) %>%
  ggplot(aes(x = elevation, y = agb_scam, color = co2)) +
  geom_point(size = 2, alpha = 0.5) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), se = F, size = 1.5) +
  facet_grid(salinity ~ age) +
  scale_color_manual(values = c("#78c679", "#006837"))+
  guides(color=guide_legend(title=expression(paste(CO[2], " level")))) +
  ylab("Aboveground Biomass (g)") +
  xlab("Elevation (m NAVD88)") +
  theme_bw() +
  theme(legend.position = "top") +
  ggtitle("WITH competition")

## Plot of random effects ####
traits_corn %>% 
  ggplot(aes(x = reorder(genotype, agb_scam, FUN = median), y = agb_scam)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 1.5) +
  ylab("Aboveground Biomass (g)") +
  xlab("Genotype") +
  geom_hline(aes(yintercept = mean(agb_scam)), linetype = "dashed") +
  theme_bw() +
  theme(axis.text.x = element_blank())


# So basically adding competition jacks everything up!!!

# Exploratory see how scam biomass influences SPPA height
full_data_merged %>% 
  filter(height_sppa > 0) %>% 
  ggplot(aes(x = agb_scam, y = height_sppa, color = elevation)) +
  geom_point(size = 3, alpha = 0.6) + 
  scale_color_gradientn(colours = rainbow(4)) +
  xlab(expression(paste(italic("S. americanus"), " aboveground biomass (g)"))) +
  ylab(expression(paste(italic("S. patens"), " maximum height (cm)")))

# Exploratory see how scam biomass influences SPPA biomass
full_data_merged %>% 
  filter(height_sppa > 0) %>% 
  ggplot(aes(x = agb_scam, y = agb_sppa, color = elevation)) +
  geom_point(size = 3, alpha = 0.6) + 
  scale_color_gradientn(colours = rainbow(4)) +
  xlab(expression(paste(italic("S. americanus"), " aboveground biomass (g)"))) +
  ylab(expression(paste(italic("S. patens"), " aboveground biomass (g)")))


#########################
## Root:Shoot analysis ##
#########################

# For the root:shoot analysis, we can only really do the no competition pots in
# the straightforward way (b/c species were not separated for the entire core;
# only separated in the top 10cm)

## Data manipulation ####
# Use the same data set as the aboveground analysis. Calculate root-to-shoot
# ratio
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
# those two

traits_nocomp %>% 
  filter(rs < 6) -> traits_nocomp_rs

# There also appear to be some that are NAs
traits_nocomp %>% 
  filter(!complete.cases(rs))
# Need to figure this out on data checking!!!

# Move forward with this for now, but I'll go back and fix these later

## Fit model and model selection (fixed effects) ####

# Most complex model
rs_mod <- lmer(log(rs) ~ weight_init + date_planted_grp + origin_lab +
                 (age + location + co2 + salinity + elevation_sc)^5 + I(elevation_sc^2) +
                 (1|genotype), data = traits_nocomp_rs)

# Try lmerTest::step
lmerTest::step(rs_mod)

# This is the best selected model
options(contrasts = c("contr.sum", "contr.poly"))
rs_step_mod <- lmer(log(rs) ~ salinity*elevation_sc + age*co2*elevation_sc +
                      (1 | genotype), data = traits_nocomp_rs)

traits_nocomp_rs %>% 
  ggplot(aes(x = elevation, y = beta, color = salinity)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x)

car::Anova(rs_step_mod, type = 3)
plot_model(rs_step_mod, terms = c("elevation_sc[all]", "co2", "age"), type = "emm")
plot_model(rs_step_mod, terms = c("elevation_sc[all]", "salinity"), type = "emm")


MuMIn::r.squaredGLMM(rs_step_mod)

summary(rs_mod)

# Look for highest order significant interaction
car::Anova(rs_mod)
# This would be age:co2:elevation (3-way interaction)

# Create a model that keeps all 2-way interactions as well as that 3-way
# interaction
rs_mod2 <- lmer(log(rs) ~ weight_init + date_planted_grp + origin_lab +
                  (age + location + co2 + salinity + elevation_sc)^2 + I(elevation_sc^2) +
                  co2:age:elevation_sc + (1|site_frame) + (1|genotype), data = traits_nocomp_rs)
# Look for highest order significant terms that are NOT contained in that
# interaction. 
car::Anova(rs_mod2)
# salinity:elevation_sc is significant. So keep that and drop out other
# non-significant interactions
rs_mod3 <- lmer(log(rs) ~ weight_init + date_planted_grp + origin_lab +
                  location + I(elevation_sc^2) + co2*age*elevation_sc +
                  salinity*elevation_sc + location +
                  (1|site_frame) + (1|genotype), data = traits_nocomp_rs)

car::Anova(rs_mod3)
# This is the final FIXED effects model. Check assumptions.
plot_model(rs_mod3, type = "diag") # All looks good!

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

rs_mod4 <- update(rs_mod3, .~.-(1|genotype) + (1+co2*elevation|genotype))
rs_mod5 <- update(rs_mod3, .~.-(1|genotype) + (1+co2*salinity|genotype))
rs_mod6 <- update(rs_mod3, .~.-(1|genotype) + (1+elevation*salinity|genotype))
rs_mod7 <- update(rs_mod3, .~.-(1|genotype) + (1+co2+elevation|genotype))
rs_mod8 <- update(rs_mod3, .~.-(1|genotype) + (1+elevation + salinity|genotype))
rs_mod9 <- update(rs_mod3, .~.-(1|genotype) + (1+salinity + co2|genotype))
rs_mod10 <- update(rs_mod3, .~.-(1|genotype) + (1+co2|genotype))
rs_mod11<- update(rs_mod3, .~.-(1|genotype) + (1+elevation|genotype))
rs_mod12 <- update(rs_mod3, .~.-(1|genotype) + (1+salinity|genotype))

# The only ones that did not have convergence issues were models 8, 11, and 12.
# Compare this to the base model.
MuMIn::r.squaredGLMM(rs_mod8)
MuMIn::r.squaredGLMM(rs_mod11)
MuMIn::r.squaredGLMM(rs_mod12)

# Compare to the original models
anova(rs_mod3, rs_mod8)# This has highest explained deviance
anova(rs_mod3, rs_mod11)
anova(rs_mod3, rs_mod12)

plot_model(rs_mod8, type = "pred", pred.type = "re",
           terms = c("elevation_sc[all]", "salinity", "genotype"))
# This seems meaningful, so let's keep it in.

# Final model is rs_mod8

## Comparison to null model ####

# Now construct a null model that does not have age, location, or genotype in it
rs_modelnull <- lmer(log(rs) ~ weight_init + date_planted_grp + origin_lab +
                       elevation_sc*salinity + co2*elevation_sc + I(elevation_sc^2) + (1|site_frame),
                     data = traits_nocomp_rs)

car::Anova(rs_modelnull)
# If no age, location, or genotype, we would say that elevation and co2
# influence biomass
plot_model(rs_modelnull, terms = c("elevation_sc[all]", "salinity"), type = "emm")

# Compare R2 from models
MuMIn::r.squaredGLMM(rs_modelnull)
# Conditional R2 = 0.43
MuMIn::r.squaredGLMM(rs_mod8)
# Conditional R2 = 0.74

# Create predicted vs observed plots for each model
tibble(predicted = c(exp(predict(rs_modelnull)), exp(predict(rs_mod8))),
       observed = rep(traits_nocomp_rs$rs, 2),
       model = c(rep(c("null", "evolution"), each = nrow(traits_nocomp_rs)))) %>% 
  ggplot(aes(x = predicted, y = observed)) +
  geom_point(alpha = 0.3) +
  facet_wrap(~model) +
  ylim(0,5) + xlim(0,5) +
  geom_abline(aes(slope = 1, intercept = 0))







## Plot of fixed effects ####

# Age:co2:elevation interaction
traits_nocomp_rs %>% 
  mutate(salinity = case_when(salinity == "fresh" ~ "Freshwater Site (6ppt)",
                              T ~ "Brackish Site (8ppt)")) %>% 
  mutate(age = case_when(age == "ancestral" ~ "Ancestral Cohort (1900-1960)",
                         T ~ "Modern Cohort (2000-2020)")) %>% 
  mutate(co2 = case_when(co2 == "elevated" ~ "Elevated (700ppm)", T ~ "Ambient (400ppm)")) %>%
  ggplot(aes(x = elevation, y = rs, color = co2)) +
  geom_point(size = 2, alpha = 0.5) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), se = F, size = 1.5) +
  facet_grid(~ age) +
  scale_color_manual(values = c("#78c679", "#006837"))+
  guides(color=guide_legend(title=expression(paste(CO[2], " level")))) +
  ylab("Root-to-Shoot Ratio") +
  xlab("Elevation (m NAVD88)") +
  theme_bw() +
  theme(legend.position = "top")

# Salinity:elevation interaction
traits_nocomp_rs %>% 
  mutate(salinity = case_when(salinity == "fresh" ~ "Freshwater Site (6ppt)",
                              T ~ "Brackish Site (8ppt)")) %>% 
  mutate(age = case_when(age == "ancestral" ~ "Ancestral Cohort (1900-1960)",
                         T ~ "Modern Cohort (2000-2020)")) %>% 
  mutate(co2 = case_when(co2 == "elevated" ~ "Elevated (700ppm)", T ~ "Ambient (400ppm)")) %>%
  ggplot(aes(x = elevation, y = rs, color = salinity)) +
  geom_point(size = 2, alpha = 0.5) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), se = F, size = 1.5) +
  scale_color_manual(values = c("#fb8072", "#80b1d3")) +
  ylab("Root-to-Shoot Ratio") +
  xlab("Elevation (m NAVD88)") +
  guides(color=guide_legend(title="Salinity")) +
  theme_bw() +
  theme(legend.position = "top")

# Random slopes plot
out <- plot_model(rs_mod8, type = "pred", pred.type = "re", 
                  terms = c("elevation_sc[all]", "salinity", "genotype")) +
  scale_color_manual(values = c("#78c679", "#006837"))

# Collect the predicted values
scale(traits_nocomp_rs$elevation)

rs_preds <- tibble(elevation = out$data$x * sd(traits_nocomp_rs$elevation) +
                     mean(traits_nocomp_rs$elevation),
                   rootshoot = out$data$predicted,
                   salinity = out$data$group,
                   genotype = out$data$facet)

# Example of different functional responses for 10 Kirkpatrick genotypes
rs_preds %>% 
  mutate(salinity = case_when(salinity == "fresh" ~ "Freshwater Site (6ppt)",
                              T ~ "Brackish Site (8ppt)")) %>%
  filter(genotype %in% c("km1", "km2", "km3", "km4", "km5", "km6",
                         "ka1", "ka2", "ka3", "ka4", "ka5", "ka6")) %>% 
  ggplot(aes(x = elevation, y = rootshoot, color = salinity)) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), se = F) +
  facet_wrap(~genotype, nrow = 2) +
  guides(color=guide_legend(title="Salinity")) +
  theme_bw() +
  ylab("Predicted Root-to-Shoot Ratio") +
  xlab("Elevation (m NAVD88)") +
  scale_color_manual(values = c("#fb8072", "#80b1d3"))+
  theme(legend.position = "top")


##########################################
## Root distribution parameter analysis ##
##########################################

## Data manipulation ####

# This is for the pots we are doing trait analysis for (no comp for right now
# for the same reasons as the root:shoot analysis)

merge(traits_nocomp_rs, betas_to_merge, by = "pot_no") -> traits_nocomp_rs

traits_nocomp_rs %>% 
  ggplot(aes(x = reorder(genotype, beta, FUN = median, na.rm = T), y = beta)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 1.5) +
  ylab("Root Distribution Parameter") +
  xlab("Genotype ID") +
  geom_hline(aes(yintercept = mean(beta, na.rm = T)), linetype = "dashed") +
  theme_bw() 




## Fit model and model selection (fixed effects) ####

# Most complex model
beta_mod <- lmer(beta ~ weight_init + date_planted_grp + origin_lab +
                   (age + location + co2 + salinity + elevation_sc)^5 + I(elevation_sc^2) +
                   (1|site_frame) + (1|genotype), data = traits_nocomp_rs)

# Try lmerTest::step
lmerTest::step(beta_mod, keep = c("co2", "salinity", "age", "elevation_sc", "location"))

beta_step_mod <- lmer(beta ~ age + location + co2 + salinity + elevation_sc +
                        (1 | genotype) + salinity:elevation_sc, data = traits_nocomp_rs)

MuMIn::r.squaredGLMM(beta_step_mod)
plot_model(beta_step_mod, terms =c("elevation_sc[all]", "salinity"), type = "emm")
plot_model(beta_step_mod, terms = c("age"), type = "emm")

# Look for highest order significant interaction
car::Anova(beta_mod)
# This would be a bunch of 2-way interactions 

# Create a model that keeps all significant 2-way interactions
beta_mod2 <- lmer(beta ~ weight_init + date_planted_grp + origin_lab +
                    salinity*elevation_sc + co2*elevation_sc + location*salinity + age + I(elevation_sc^2) +
                    (1|site_frame) + (1|genotype), data = traits_nocomp_rs)
# Look for highest order significant terms that are NOT contained in that
# interaction. 
car::Anova(beta_mod2)

# This is the final FIXED effects model. Check assumptions.
plot_model(beta_mod2, type = "diag") # All looks good!

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

beta_mod3 <- update(beta_mod2, .~.-(1|genotype) + (1+co2*elevation|genotype))
beta_mod4 <- update(beta_mod2, .~.-(1|genotype) + (1+co2*salinity|genotype))
beta_mod5 <- update(beta_mod2, .~.-(1|genotype) + (1+elevation*salinity|genotype)) # fit
beta_mod6 <- update(beta_mod2, .~.-(1|genotype) + (1+co2+elevation|genotype))
beta_mod7 <- update(beta_mod2, .~.-(1|genotype) + (1+elevation + salinity|genotype)) # fit
beta_mod8 <- update(beta_mod2, .~.-(1|genotype) + (1+salinity + co2|genotype))
beta_mod9 <- update(beta_mod2, .~.-(1|genotype) + (1+co2|genotype))
beta_mod10<- update(beta_mod2, .~.-(1|genotype) + (1+elevation|genotype)) # fit
beta_mod11 <- update(beta_mod2, .~.-(1|genotype) + (1+salinity|genotype)) # fit

# The only ones that did not have convergence issues were models 3, 7, 9, 10, 11.
# Compare this to the base model.
MuMIn::r.squaredGLMM(beta_mod3)
MuMIn::r.squaredGLMM(beta_mod7)
MuMIn::r.squaredGLMM(beta_mod9)
MuMIn::r.squaredGLMM(beta_mod10)
MuMIn::r.squaredGLMM(beta_mod11)

# Compare to the original models
anova(beta_mod2, beta_mod3)
anova(beta_mod2, beta_mod7)# This has highest explained deviance
anova(beta_mod2, beta_mod9)
anova(beta_mod2, beta_mod10)
anova(beta_mod2, beta_mod11)

plot_model(beta_mod5, type = "pred", pred.type = "re",
           terms = c("elevation_sc[all]", "salinity", "genotype"))
# This seems meaningful, so let's keep it in.

# Final model is beta_mod5

## Comparison to null model ####

# Now construct a null model that does not have age, location, or genotype in it
beta_modelnull <- lmer(beta ~ weight_init + date_planted_grp + origin_lab + I(elevation_sc^2) +
                         elevation_sc*salinity + co2*elevation_sc + I(elevation_sc^2) + (1|site_frame),
                       data = traits_nocomp_rs)

car::Anova(beta_modelnull)
# If no age, location, or genotype, we would say that elevation and co2
# influence biomass
plot_model(beta_modelnull, terms = c("elevation_sc[all]", "salinity"), type = "emm")
plot_model(beta_modelnull, terms = c("elevation_sc[all]", "co2"), type = "emm")

# Compare R2 from models
MuMIn::r.squaredGLMM(beta_modelnull)
# Conditional R2 = 0.49
MuMIn::r.squaredGLMM(beta_mod7)
# Conditional R2 = 0.72

# Create predicted vs observed plots for each model

# Subset out data set with no NAs for plotting
traits_nocomp_rs %>% 
  filter(complete.cases(beta)) -> traits_nocomp_rs_noNAs

tibble(predicted = c(predict(beta_modelnull), predict(beta_mod5)),
       observed = rep(traits_nocomp_rs_noNAs$beta, 2),
       model = c(rep(c("null", "evolution"), each = nrow(traits_nocomp_rs_noNAs)))) %>% 
  ggplot(aes(x = predicted, y = observed)) +
  geom_point(alpha = 0.3) +
  facet_wrap(~model) +
  ylim(0.75,0.97) + xlim(0.75,0.97) +
  geom_abline(aes(slope = 1, intercept = 0))


##################################
## Belowground biomass analysis ##
##################################

# For the belowground analysis, we can only really do the no competition pots in
# the straightforward way (b/c species were not separated for the entire core;
# only separated in the top 10cm). Use the same subset of data as for the
# root:shoot analysis for now.

## Data manipulation ####

## Fit model and model selection (fixed effects) ####

# Most complex model
bg_mod <- lmer(total_bg ~ weight_init + date_planted_grp + origin_lab +
                 (age + location + co2 + salinity + elevation_sc)^5 + I(elevation_sc^2) +
                 (1|genotype), data = traits_nocomp_rs)

summary(bg_mod)

# Try stepwise
lmerTest::step(bg_mod)

bg_step_mod <- lmer(total_bg ~ weight_init + date_planted_grp + co2 + I(elevation_sc^2) + 
                      (1 | genotype), data = traits_nocomp_rs)

# Look for highest order significant interaction. No significant interactions.
car::Anova(bg_mod)

# Create a model that keeps all main terms
bg_mod2 <- lmer(total_bg ~ weight_init + date_planted_grp + origin_lab +
                  age + location + co2 + salinity + elevation_sc + I(elevation_sc^2) +
                  (1|site_frame) + (1|genotype), data = traits_nocomp_rs)

# This is the final FIXED effects model. Check assumptions.
plot_model(bg_mod2, type = "diag") # All looks good!

# Overall significance of terms in final fixed effects model
car::Anova(bg_mod2)
# Of interest are just additive terms of co2 and elevation

## Model selection (random slopes) ####

# Now we are going to try and add up to 2-way interactions as random slopes.

# Here are the possibilities:
# co2 + elevation
# co2 + salinity
# salinity + elevation
# co2
# elevation
# salinity

bg_mod3 <- update(bg_mod2, .~.-(1|genotype) + (1+co2+elevation|genotype))
bg_mod4 <- update(bg_mod2, .~.-(1|genotype) + (1+elevation + salinity|genotype))
bg_mod5 <- update(bg_mod2, .~.-(1|genotype) + (1+salinity + co2|genotype))
bg_mod6 <- update(bg_mod2, .~.-(1|genotype) + (1+co2|genotype))
bg_mod7<- update(bg_mod2, .~.-(1|genotype) + (1+elevation|genotype))
bg_mod8 <- update(bg_mod2, .~.-(1|genotype) + (1+salinity|genotype))

# Nothing converged, so base model is the best

# Look at fixed effects
plot_model(bg_mod2, type = "emm", terms = "elevation_sc[all]")
plot_model(bg_mod2, type = "emm", terms = "co2")

## Comparison to null model ####

# Now construct a null model that does not have age, location, or genotype in it
bg_modelnull <- lmer(total_bg ~ weight_init + date_planted_grp + origin_lab +
                       co2 + salinity + elevation_sc + I(elevation_sc^2) +
                       (1|site_frame),
                     data = traits_nocomp_rs)

car::Anova(bg_modelnull)
# If no age, location, or genotype, we would say that elevation and co2
# influence biomass still

# Compare R2 from models
MuMIn::r.squaredGLMM(bg_modelnull)
# Conditional R2 = 0.23
MuMIn::r.squaredGLMM(bg_mod2)
# Conditional R2 = 0.31

# Create predicted vs observed plots for each model
tibble(predicted = c(predict(bg_modelnull), predict(bg_mod2)),
       observed = rep(traits_nocomp_rs$total_bg, 2),
       model = c(rep(c("null", "evolution"), each = nrow(traits_nocomp_rs)))) %>% 
  ggplot(aes(x = predicted, y = observed)) +
  geom_point(alpha = 0.3) +
  facet_wrap(~model) +
  geom_abline(aes(slope = 1, intercept = 0))


## Plot of fixed effects ####

# Effect of  interaction
traits_nocomp_rs %>% 
  mutate(salinity = case_when(salinity == "fresh" ~ "Freshwater Site (6ppt)",
                              T ~ "Brackish Site (8ppt)")) %>% 
  mutate(age = case_when(age == "ancestral" ~ "Ancestral Cohort (1900-1960)",
                         T ~ "Modern Cohort (2000-2020)")) %>% 
  mutate(co2 = case_when(co2 == "elevated" ~ "Elevated (700ppm)", T ~ "Ambient (400ppm)")) %>%
  ggplot(aes(x = elevation, y = total_bg, color = co2)) +
  geom_point(size = 2, alpha = 0.5) +
  scale_color_manual(values = c("#78c679", "#006837"))+
  guides(color=guide_legend(title=expression(paste(CO[2], " level")))) +
  ylab("Total belowground biomass (g)") +
  xlab("Elevation (m NAVD88)") +
  theme_bw() +
  theme(legend.position = "top")

predict_frame <- expand.grid(elevation = seq(min(traits_nocomp_rs$total_bg), max(traits_nocomp_rs$total_bg), length.out = 100),
                             co2 = c("ambient", "elevated"))

# Get predicted means of each co2 treatment
bg_means <- emmeans(bg_mod2, ~elevation_sc:co2, at = list(elevation_sc = seq(min(traits_nocomp_rs$elevation_sc),
                                                                             max(traits_nocomp_rs$elevation_sc), length.out = 113)))
# Get true elevations
elevation <- seq(min(traits_nocomp_rs$elevation),
                 max(traits_nocomp_rs$elevation), length.out = 113)

# Make a data frame of predictions
predict_frame <- data.frame(total_bg = summary(bg_means)$emmean,
                            elevation = rep(elevation,2),
                            co2 = rep(c("ambient", "elevated"), each = 113)) %>% 
  mutate(co2 = case_when(co2 == "elevated" ~ "Elevated (700ppm)", T ~ "Ambient (400ppm)"))

ggplot(dta, aes(math, num_awards, col = prog)) +
  geom_point() +
  geom_smooth(method = "glm", se = FALSE,
              method.args = list(family = "poisson"), linetype = "dashed")+
  geom_line(data=pframe)  ## use prediction data here
## (inherits aesthetics etc. from main ggplot call)

# Salinity:elevation interaction
traits_nocomp_rs %>% 
  mutate(salinity = case_when(salinity == "fresh" ~ "Freshwater Site (6ppt)",
                              T ~ "Brackish Site (8ppt)")) %>% 
  mutate(age = case_when(age == "ancestral" ~ "Ancestral Cohort (1900-1960)",
                         T ~ "Modern Cohort (2000-2020)")) %>% 
  mutate(co2 = case_when(co2 == "elevated" ~ "Elevated (700ppm)", T ~ "Ambient (400ppm)")) -> plot_data

ggplot(plot_data, aes(x = elevation, y = total_bg, color = co2)) +
  geom_point(size = 2, alpha = 0.5) +
  scale_color_manual(values = c("#78c679", "#006837")) +
  guides(color=guide_legend(title=expression(paste(CO[2], " level")))) +
  geom_line(aes(x = predict_frame$elevation, y = predict_frame$total_bg, color = predict_frame$co2), size = 1.5) +
  ylab("Belowground biomass (g)") +
  xlab("Elevation (m NAVD88)") +
  theme_bw() +
  theme(legend.position = "top") 


##########################
## Stem height analysis ##
##########################

# Use the same as for the AGB analysis

##
# No competition
##

## Data manipulation ####

## Fit model and model selection (fixed effects) ####

# Most complex model
height_mod <- lmer(height_scam_tot ~ weight_init + date_planted_grp + origin_lab +
                     (age + location + co2 + salinity + elevation_sc)^5 + I(elevation_sc^2) +
                     (1|site_frame) + (1|genotype), data = traits_nocomp)

summary(height_mod)
plot(height_mod) # There appears to be a huge outlier (obs # 186)
plot_model(height_mod, type = "diag")

# Look at overall differences and then drop and see if outlier mattered
car::Anova(height_mod)

traits_nocomp_noOut <- traits_nocomp[-186,]
height_mod_noOut <- lmer(height_scam_tot ~ weight_init + date_planted_grp + origin_lab +
                           (age + location + co2 + salinity + elevation_sc)^5 + I(elevation_sc^2) +
                           (1|genotype), data = traits_nocomp_noOut)
plot(height_mod_noOut)
car::Anova(height_mod_noOut) # That seemed to be driving somethings. So let's drop it!

# Try step
lmerTest::step(height_mod_noOut)

height_step_mod <- lmer(height_scam_tot ~ date_planted_grp + elevation_sc + (1 | genotype), data = traits_nocomp_noOut)

# Create a model that keeps location:salinity; also drop site_frame because it
# does not converge now
height_mod_noOut2 <- lmer(height_scam_tot ~ weight_init + date_planted_grp + origin_lab +
                            age + co2 + elevation_sc + I(elevation_sc^2) +
                            location*salinity + (1|genotype), data = traits_nocomp_noOut)

car::Anova(height_mod_noOut2)

# This is the final FIXED effects model. Check assumptions.
plot_model(height_mod_noOut2, type = "diag") # All looks good!

## Model selection (random slopes) ####

# Now we are going to try and add up to 2-way interactions as random slopes.

# Here are the possibilities:
# co2 + elevation
# co2 + salinity
# salinity + elevation
# co2
# elevation
# salinity

height_mod_noOut3 <- update(height_mod_noOut2, .~.-(1|genotype) + (1+co2+elevation|genotype))
height_mod_noOut4 <- update(height_mod_noOut2, .~.-(1|genotype) + (1+elevation + salinity|genotype))
height_mod_noOut5 <- update(height_mod_noOut2, .~.-(1|genotype) + (1+salinity + co2|genotype))
height_mod_noOut6 <- update(height_mod_noOut2, .~.-(1|genotype) + (1+co2|genotype))
height_mod_noOut7<- update(height_mod_noOut2, .~.-(1|genotype) + (1+elevation|genotype))
height_mod_noOut8 <- update(height_mod_noOut2, .~.-(1|genotype) + (1+salinity|genotype))

# The only ones that did not have estimated zeros were 3, 6, and 7
MuMIn::r.squaredGLMM(height_mod_noOut3)
MuMIn::r.squaredGLMM(height_mod_noOut6)
MuMIn::r.squaredGLMM(height_mod_noOut7)

# Compare to base model
anova(height_mod_noOut2, height_mod_noOut3)
anova(height_mod_noOut2, height_mod_noOut6)
anova(height_mod_noOut2, height_mod_noOut7) 

# None improve the model significantly

plot_model(height_mod_noOut3, type = "pred", pred.type = "re",
           terms = c("elevation_sc[all]", "co2", "genotype"))
# This doesn't seem meaningful so it seems like the best (and most parsimonious
# model) does not have a random slope.

# Look at random intercepts
plot_model(height_mod_noOut2, type = "re")

## Comparison to null model ####

# Now construct a null model that does not have age, location, or genotype in it
height_model_noOut_null <- lm(height_scam_tot ~ weight_init + date_planted_grp + origin_lab +
                                co2 + salinity + elevation_sc + I(elevation_sc^2),
                              data = traits_nocomp_noOut)

car::Anova(height_model_noOut_null)
# If no age, location, or genotype, we would say that elevation was the only
# thing that influenced height.
plot_model(height_model_noOut_null, terms = "elevation_sc[all]", type = "emm")

# Compare R2 from models
summary(height_model_noOut_null)
# Adjusted R2 = 0.38
MuMIn::r.squaredGLMM(height_mod_noOut2)
# Conditional R2 = 0.61

# Create predicted vs observed plots for each model
observed <- rep(traits_nocomp_noOut[!is.na(traits_nocomp_noOut$height_scam_tot), "height_scam_tot"], 2)

tibble(predicted = c(predict(height_model_noOut_null), predict(height_mod_noOut2)),
       observed = observed,
       model = c(rep(c("null", "evolution"), each =length(observed)/2))) %>% 
  ggplot(aes(x = predicted, y = observed)) +
  geom_point(alpha = 0.3) +
  facet_wrap(~model)  +
  geom_abline(aes(slope = 1, intercept = 0))

# Repeat this process for the competition pots

##
# Corn only
##

## Fit model and model selection (fixed effects) ####

# Most complex model
height_mod_corn <- lmer(height_scam_tot ~ weight_init + date_planted_grp + origin_lab +
                          (age + comp + co2 + salinity + elevation_sc)^5 + I(elevation_sc^2) +
                          (1|genotype), data = traits_corn)

# Try step
lmerTest::step(height_mod_corn)

traits_corn %>% filter(height_scam_tot > 20) ->traits_corn_sub
height_corn_step_mod <- lmer(height_scam_tot ~ weight_init + age*comp*co2*salinity + (1 | genotype)  + 
                               age*comp*co2*elevation_sc, data = traits_corn_sub)

plot_model(height_corn_step_mod, terms = c("salinity", "age", "comp", "co2"), type = "emm")

traits_corn %>% 
  ggplot(aes(x = salinity, y = height_scam_tot, color = age)) +
  geom_boxplot() +
  facet_wrap(~comp + co2) +
  geom_smooth(method = "lm")


car::Anova(height_corn_step_mod)
# This would be age:co2:elevation:comp (4-way interaction)

# Fit model with that 4-way interaction and all possible 3-way interactions
height_mod_corn3 <- lmer(height_scam_tot ~ weight_init + date_planted_grp + origin_lab +
                           (age + comp + co2 + salinity + elevation_sc)^3 + I(elevation_sc^2) +
                           age:co2:elevation_sc:comp + (1|genotype), data = traits_corn)

car::Anova(height_mod_corn3) # Should drop out 4-way interaction

height_mod_corn4 <- lmer(height_scam_tot ~ weight_init + date_planted_grp + origin_lab +
                           (age + comp + co2 + salinity + elevation_sc)^3 + I(elevation_sc^2) +
                           (1|genotype), data = traits_corn)

car::Anova(height_mod_corn4) # Drop all 3-way interactions but comp:co2:salinity

height_mod_corn5 <- lmer(height_scam_tot ~ weight_init + date_planted_grp + origin_lab +
                           (age + comp + co2 + salinity + elevation_sc)^2 + I(elevation_sc^2) +
                           comp:co2:salinity + (1|genotype), data = traits_corn)

car::Anova(height_mod_corn5) 
# Now drop any 2-way interactions that are not included and are not
# co2:elevation_sc

height_mod_corn6 <- lmer(height_scam_tot ~ weight_init + date_planted_grp +
                           origin_lab + age + I(elevation_sc^2) +
                           comp*co2*salinity + co2*elevation_sc +
                           (1|genotype), data = traits_corn)

car::Anova(height_mod_corn6)

# This is the final FIXED effects model. Check assumptions.
plot_model(height_mod_corn6, type = "diag") # All looks good!

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

height_mod_corn7 <- update(height_mod_corn6, .~.-(1|genotype) + (1+co2*elevation|genotype))
height_mod_corn8 <- update(height_mod_corn6, .~.-(1|genotype) + (1+co2*salinity|genotype))
height_mod_corn9 <- update(height_mod_corn6, .~.-(1|genotype) + (1+elevation*salinity|genotype))
height_mod_corn10 <- update(height_mod_corn6, .~.-(1|genotype) + (1+co2+elevation|genotype))
height_mod_corn11 <- update(height_mod_corn6, .~.-(1|genotype) + (1+elevation + salinity|genotype))
height_mod_corn12 <- update(height_mod_corn6, .~.-(1|genotype) + (1+salinity + co2|genotype))
height_mod_corn13 <- update(height_mod_corn6, .~.-(1|genotype) + (1+co2|genotype))
height_mod_corn14<- update(height_mod_corn6, .~.-(1|genotype) + (1+elevation|genotype))
height_mod_corn15 <- update(height_mod_corn6, .~.-(1|genotype) + (1+salinity|genotype))

# The only one that did not have convergence issues or 0 variance was the
# elevation random slope.
MuMIn::r.squaredGLMM(height_mod_corn6)
MuMIn::r.squaredGLMM(height_mod_corn14)

anova(height_mod_corn6, height_mod_corn14)

plot_model(height_mod_corn14, type = "pred", pred.type = "re",
           terms = c("elevation_sc[all]", "genotype")) +
  scale_color_manual(values = rainbow(32))
# This doesn't seem meaningful so it seems like the best (and most parsimonious
# model) does not have a random slope.

## Comparison to null model ####

# Now construct a null model that does not have age, location, or genotype in it
height_model_cornnull <- lm(height_scam_tot ~ weight_init + date_planted_grp +
                              origin_lab + I(elevation_sc^2) +
                              comp*co2*salinity + co2*elevation_sc, data = traits_corn)

car::Anova(height_model_cornnull)
# If no age, location, or genotype, we would say that co2 x elevation
# influences biomass
plot_model(height_model_cornnull, terms = c("elevation_sc[all]", "co2"), type = "emm")

# Compare R2 from models
summary(height_model_cornnull)
# Adjusted R2 = 0.40
MuMIn::r.squaredGLMM(height_mod_corn6)
# Conditional R2 = 0.60

# Create predicted vs observed plots for each model
observed_height <- rep(traits_corn[!is.na(traits_corn$height_scam_tot), "height_scam_tot"], 2)

tibble(predicted = c(predict(height_model_cornnull), predict(height_mod_corn6)),
       observed = observed_height,
       model = c(rep(c("null", "evolution"), each = length(observed)/2))) %>% 
  ggplot(aes(x = predicted, y = observed)) +
  geom_point(alpha = 0.3) +
  facet_wrap(~model) +
  ylim(10,80) + xlim(10,80) +
  geom_abline(aes(slope = 1, intercept = 0))


## Plot of fixed effects ####
traits_corn %>% 
  filter(comp == 0) %>% 
  mutate(salinity = case_when(salinity == "fresh" ~ "Freshwater Site (6ppt)", T ~ "Brackish Site (8ppt)")) %>% 
  mutate(age = case_when(age == "ancestral" ~ "Ancestral Cohort (1900-1960)", T ~ "Modern Cohort (2000-2020)")) %>% 
  mutate(co2 = case_when(co2 == "elevated" ~ "Elevated (700ppm)", T ~ "Ambient (400ppm)")) %>%
  ggplot(aes(x = elevation, y = agb_scam, color = co2)) +
  geom_point(size = 2, alpha = 0.5) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), se = F, size = 1.5) +
  facet_grid(salinity ~ age) +
  scale_color_manual(values = c("#78c679", "#006837"))+
  guides(color=guide_legend(title=expression(paste(CO[2], " level")))) +
  ylab("Aboveground Biomass (g)") +
  xlab("Elevation (m NAVD88)") +
  theme_bw() +
  theme(legend.position = "top") +
  ggtitle("WITHOUT competition")

traits_corn %>% 
  filter(comp == 1) %>% 
  mutate(salinity = case_when(salinity == "fresh" ~ "Freshwater Site (6ppt)", T ~ "Brackish Site (8ppt)")) %>% 
  mutate(age = case_when(age == "ancestral" ~ "Ancestral Cohort (1900-1960)", T ~ "Modern Cohort (2000-2020)")) %>% 
  mutate(co2 = case_when(co2 == "elevated" ~ "Elevated (700ppm)", T ~ "Ambient (400ppm)")) %>%
  ggplot(aes(x = elevation, y = agb_scam, color = co2)) +
  geom_point(size = 2, alpha = 0.5) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), se = F, size = 1.5) +
  facet_grid(salinity ~ age) +
  scale_color_manual(values = c("#78c679", "#006837"))+
  guides(color=guide_legend(title=expression(paste(CO[2], " level")))) +
  ylab("Aboveground Biomass (g)") +
  xlab("Elevation (m NAVD88)") +
  theme_bw() +
  theme(legend.position = "top") +
  ggtitle("WITH competition")

## Plot of random effects ####
traits_corn %>% 
  filter(complete.cases(height_scam_tot)) %>% 
  ggplot(aes(x = reorder(genotype, height_scam_tot, FUN = median), y = height_scam_tot)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 1.5) +
  ylab("Mean stem height (cm)") +
  xlab("Genotype") +
  geom_hline(aes(yintercept = mean(height_scam_tot, na.rm = T)), linetype = "dashed") +
  theme_bw() #+
#theme(axis.text.x = element_blank())




##########################
## Stem width analysis ##
##########################

# Use the same as for the AGB analysis

##
# No competition
##

## Data manipulation ####

## Fit model and model selection (fixed effects) ####

# Most complex model
width_mod <- lmer(width_scam_mid ~ weight_init + date_planted_grp + origin_lab +
                    (age + location + co2 + salinity + elevation_sc)^5 + I(elevation_sc^2) +
                    (1|genotype), data = traits_nocomp)

# Step mod

lmerTest::step(width_mod, keep = c("co2", "salinity", "elevation_sc", "location", "age"))

width_step_mod <- lmer(width_scam_mid ~ date_planted_grp + location + elevation_sc + 
                         salinity + (1 | genotype) + location*elevation_sc + location*salinity, data = traits_nocomp)

summary(width_mod) # Looks like we could probably drop site_frame
plot_model(width_mod, type = "diag") # Looks good

# Look at overall differences and then drop and see if outlier mattered
car::Anova(width_mod)

# Create a model that keeps location:salinity and location:elevation_sc; also drop site_frame because it
# does not converge now
width_mod2 <- lmer(width_scam_mid ~ weight_init + date_planted_grp + origin_lab +
                     age + co2 + I(elevation_sc^2) + location*salinity + location*elevation_sc + (1|genotype), data = traits_nocomp)

car::Anova(width_mod2)

# This is the final FIXED effects model. Check assumptions.
plot_model(width_mod2, type = "diag") # All looks good!

## Model selection (random slopes) ####

# Now we are going to try and add up to 2-way interactions as random slopes.

# Here are the possibilities:
# co2 + elevation
# co2 + salinity
# salinity + elevation
# co2
# elevation
# salinity

width_mod3 <- update(width_mod2, .~.-(1|genotype) + (1+co2+elevation|genotype))
width_mod4 <- update(width_mod2, .~.-(1|genotype) + (1+elevation + salinity|genotype))
width_mod5 <- update(width_mod2, .~.-(1|genotype) + (1+salinity + co2|genotype))
width_mod6 <- update(width_mod2, .~.-(1|genotype) + (1+co2|genotype))
width_mod7<- update(width_mod2, .~.-(1|genotype) + (1+elevation|genotype))
width_mod8 <- update(width_mod2, .~.-(1|genotype) + (1+salinity|genotype))

# The only ones that did not have estimated zeros were 3, 4, 7 and 8
MuMIn::r.squaredGLMM(width_mod3)
MuMIn::r.squaredGLMM(width_mod4) # This seems to take into account the most
MuMIn::r.squaredGLMM(width_mod7)
MuMIn::r.squaredGLMM(width_mod8)

# Compare to base model
anova(width_mod2, width_mod3)
anova(width_mod2, width_mod4)
anova(width_mod2, width_mod7) 
anova(width_mod2, width_mod8)

# None improve the model significantly

plot_model(width_mod3, type = "pred", pred.type = "re",
           terms = c("elevation_sc[all]", "co2", "genotype"))
# This doesn't seem meaningful so it seems like the best (and most parsimonious
# model) does not have a random slope.

# Look at random intercepts
plot_model(width_mod2, type = "re")

# Look at fixed effects
plot_model(width_mod2, type = "emm", terms = c("salinity", "location"))
# Kirkpatrick genotypes are sensitive to saline conditions
plot_model(width_mod2, type = "emm", terms = c("elevation_sc[all]", "location"))

## Comparison to null model ####

# Now construct a null model that does not have age, location, or genotype in it
width_model_null <- lm(width_scam_mid ~ weight_init + date_planted_grp + origin_lab +
                         co2 + salinity + elevation_sc + I(elevation_sc^2),
                       data = traits_nocomp)

car::Anova(width_model_null)
# If no age, location, or genotype, we would say that elevation was the only
# thing that influenced height.
plot_model(width_model_null, terms = "elevation_sc[all]", type = "emm")

# Compare R2 from models
summary(width_model_null)
# Adjusted R2 = 0.42
MuMIn::r.squaredGLMM(width_mod2)
# Conditional R2 = 0.63

# Create predicted vs observed plots for each model
observed <- rep(traits_nocomp[!is.na(traits_nocomp$width_scam_mid), "width_scam_mid"], 2)

tibble(predicted = c(predict(width_model_null), predict(width_mod2)),
       observed = observed,
       model = c(rep(c("null", "evolution"), each =length(observed)/2))) %>% 
  ggplot(aes(x = predicted, y = observed)) +
  geom_point(alpha = 0.3) +
  xlim(1.5,5.5) + ylim(1.5,5.5) +
  facet_wrap(~model)  +
  geom_abline(aes(slope = 1, intercept = 0))

# Repeat this process for the competition pots

##
# Corn only
##

## Fit model and model selection (fixed effects) ####

# Most complex model
width_mod_corn <- lmer(width_scam_mid ~ weight_init + date_planted_grp + origin_lab +
                         (age + comp + co2 + salinity + elevation_sc)^3 + I(elevation_sc^2) +
                         (1|site_frame) + (1|genotype), data = traits_corn)

# Try step mod
lmerTest::step(width_mod_corn)

# Step mod
width_corn_step_mod <- lmer(width_scam_mid ~ weight_init + comp + co2 + salinity + 
                              elevation_sc + (1 | genotype) + 
                              co2*elevation_sc + I(elevation_sc^2) +
                              age*co2*salinity + age*salinity*elevation_sc, data = traits_corn)

plot_model(width_corn_step_mod, terms = c("co2", "age", "salinity"), type = "emm")

pairs(emmeans::emmeans(width_corn_step_mod, ~ co2|salinity|age))

summary(width_mod_corn)

car::Anova(width_mod_corn)
# This would be age:salinity:elevation and age:co2:salinity (3-way interactions)

# Fit model with that 4-way interaction and all possible 3-way interactions
width_mod_corn2 <- lmer(width_scam_mid ~ weight_init + date_planted_grp + origin_lab +
                          (age + comp + co2 + salinity + elevation_sc)^2 + I(elevation_sc^2) +
                          age:salinity:elevation_sc + age:co2:salinity + (1|genotype), data = traits_corn)

car::Anova(width_mod_corn2) # Keep 3 way interactions (and dependents) as well as co2:elevation

width_mod_corn3 <- lmer(width_scam_mid ~ weight_init + date_planted_grp + origin_lab +
                          comp + salinity*co2*age + salinity*age*elevation_sc + elevation_sc*co2 + I(elevation_sc^2) +
                          (1|genotype), data = traits_corn)

car::Anova(width_mod_corn3) # Final fixed effects model


# This is the final FIXED effects model. Check assumptions.
plot_model(width_mod_corn3, type = "diag") # All looks good!

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

width_mod_corn4 <- update(width_mod_corn3, .~.-(1|genotype) + (1+co2*elevation|genotype))
width_mod_corn5 <- update(width_mod_corn3, .~.-(1|genotype) + (1+co2*salinity|genotype))
width_mod_corn6 <- update(width_mod_corn3, .~.-(1|genotype) + (1+elevation*salinity|genotype))
width_mod_corn7 <- update(width_mod_corn3, .~.-(1|genotype) + (1+co2+elevation|genotype))
width_mod_corn8 <- update(width_mod_corn3, .~.-(1|genotype) + (1+elevation + salinity|genotype))
width_mod_corn9 <- update(width_mod_corn3, .~.-(1|genotype) + (1+salinity + co2|genotype))
width_mod_corn10 <- update(width_mod_corn3, .~.-(1|genotype) + (1+co2|genotype))
width_mod_corn11<- update(width_mod_corn3, .~.-(1|genotype) + (1+elevation|genotype))
width_mod_corn12 <- update(width_mod_corn3, .~.-(1|genotype) + (1+salinity|genotype))

# The only ones that did not have convergence issues or 0 variance were 6,8,11,12.
MuMIn::r.squaredGLMM(width_mod_corn6)
MuMIn::r.squaredGLMM(width_mod_corn8)
MuMIn::r.squaredGLMM(width_mod_corn11)
MuMIn::r.squaredGLMM(width_mod_corn12)

anova(width_mod_corn6, width_mod_corn3)
anova(width_mod_corn8, width_mod_corn3)
anova(width_mod_corn11, width_mod_corn3)
anova(width_mod_corn12, width_mod_corn3)

plot_model(width_mod_corn6, type = "pred", pred.type = "re",
           terms = c("elevation_sc[all]", "salinity", "genotype")) +
  scale_color_manual(values = rainbow(32))
# This doesn't seem meaningful so it seems like the best (and most parsimonious
# model) does not have a random slope.

## Comparison to null model ####

# Now construct a null model that does not have age, location, or genotype in it
width_model_cornnull <- lm(width_scam_mid ~ weight_init + date_planted_grp +
                             origin_lab + I(elevation_sc^2) + comp + 
                             (salinity + elevation_sc + co2)^2, data = traits_corn)

car::Anova(width_model_cornnull)
# If no age, location, or genotype, we would say that co2 x elevation
# influences width as well as salinity x elevation and comp
plot_model(width_model_cornnull, terms = c("elevation_sc[all]", "co2"), type = "emm")

# Compare R2 from models
summary(width_model_cornnull)
# Adjusted R2 = 0.54
MuMIn::r.squaredGLMM(width_mod_corn3)
# Conditional R2 = 0.67

# Create predicted vs observed plots for each model
observed_width <- rep(traits_corn[!is.na(traits_corn$width_scam_mid), "width_scam_mid"], 2)

tibble(predicted = c(predict(width_model_cornnull), predict(width_mod_corn3)),
       observed = observed_width,
       model = c(rep(c("null", "evolution"), each = length(observed)/2))) %>% 
  ggplot(aes(x = predicted, y = observed)) +
  geom_point(alpha = 0.3) +
  facet_wrap(~model) + xlim(1.5,5.5) + ylim(1.5,5.5) +
  geom_abline(aes(slope = 1, intercept = 0))

## Plot of fixed effects ####

plot_model(width_mod_corn3, type = "emm", terms = c("salinity", "age", "co2"))

traits_corn %>% 
  filter(complete.cases(width_scam_mid)) -> test_data

traits_corn %>% 
  filter(comp == 0) %>% 
  mutate(salinity = case_when(salinity == "fresh" ~ "Freshwater Site (6ppt)", T ~ "Brackish Site (8ppt)")) %>% 
  mutate(age = case_when(age == "ancestral" ~ "Ancestral Cohort (1900-1960)", T ~ "Modern Cohort (2000-2020)")) %>% 
  mutate(co2 = case_when(co2 == "elevated" ~ "Elevated (700ppm)", T ~ "Ambient (400ppm)")) %>%
  ggplot(aes(x = elevation, y = agb_scam, color = co2)) +
  geom_point(size = 2, alpha = 0.5) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), se = F, size = 1.5) +
  facet_grid(salinity ~ age) +
  scale_color_manual(values = c("#78c679", "#006837"))+
  guides(color=guide_legend(title=expression(paste(CO[2], " level")))) +
  ylab("Aboveground Biomass (g)") +
  xlab("Elevation (m NAVD88)") +
  theme_bw() +
  theme(legend.position = "top") +
  ggtitle("WITHOUT competition")

traits_corn %>% 
  filter(comp == 1) %>% 
  mutate(salinity = case_when(salinity == "fresh" ~ "Freshwater Site (6ppt)", T ~ "Brackish Site (8ppt)")) %>% 
  mutate(age = case_when(age == "ancestral" ~ "Ancestral Cohort (1900-1960)", T ~ "Modern Cohort (2000-2020)")) %>% 
  mutate(co2 = case_when(co2 == "elevated" ~ "Elevated (700ppm)", T ~ "Ambient (400ppm)")) %>%
  ggplot(aes(x = elevation, y = agb_scam, color = co2)) +
  geom_point(size = 2, alpha = 0.5) +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), se = F, size = 1.5) +
  facet_grid(salinity ~ age) +
  scale_color_manual(values = c("#78c679", "#006837"))+
  guides(color=guide_legend(title=expression(paste(CO[2], " level")))) +
  ylab("Aboveground Biomass (g)") +
  xlab("Elevation (m NAVD88)") +
  theme_bw() +
  theme(legend.position = "top") +
  ggtitle("WITH competition")

## Plot of random effects ####
traits_corn %>% 
  filter(complete.cases(height_scam_tot)) %>% 
  ggplot(aes(x = reorder(genotype, height_scam_tot, FUN = median), y = height_scam_tot)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 1.5) +
  ylab("Mean stem height (cm)") +
  xlab("Genotype") +
  geom_hline(aes(yintercept = mean(height_scam_tot, na.rm = T)), linetype = "dashed") +
  theme_bw() #+
#theme(axis.text.x = element_blank())






## Other exploratory stuff ####
# Exploratory see how scam biomass influences SPPA height
full_data_merged %>% 
  filter(height_sppa > 0) %>% 
  ggplot(aes(x = agb_scam, y = height_sppa, color = elevation)) +
  geom_point(size = 3, alpha = 0.6) + 
  scale_color_gradientn(colours = rainbow(4)) +
  xlab(expression(paste(italic("S. americanus"), " aboveground biomass (g)"))) +
  ylab(expression(paste(italic("S. patens"), " maximum height (cm)")))

# Exploratory see how scam biomass influences SPPA biomass
full_data_merged %>% 
  filter(height_sppa > 0) %>% 
  ggplot(aes(x = agb_scam, y = agb_sppa, color = elevation)) +
  geom_smooth(method = "lm", formula = y ~ I(log(x+0.00001))) +
  geom_point(size = 3, alpha = 0.6) + 
  scale_color_gradientn(colours = rainbow(4)) +
  xlab(expression(paste(italic("S. americanus"), " aboveground biomass (g)"))) +
  ylab(expression(paste(italic("S. patens"), " aboveground biomass (g)")))



###########################
## Stem density analysis ##
###########################

# Stem density uses the same data as AGB
density_mod <- lmer(dens_scam_live ~ weight_init + date_planted_grp + origin_lab +
                      (age + location + co2 + salinity + elevation_sc)^5 + I(elevation_sc^2) +
                      (1|genotype), data = traits_nocomp)

lmerTest::step(density_mod, keep = c("co2", "salinity", "age", "location", "elevation_sc"))

density_step_mod <- lmer(dens_scam_live ~ age + location + co2 + salinity + elevation_sc + 
                           I(elevation_sc^2) + (1 | genotype), data = traits_nocomp)

car::Anova(density_step_mod)

plot_model(density_step_mod, terms = c("elevation_sc[all]", "co2"), type = "emm")



#############################
## Competition in top 10cm ##
#############################

# Filter out top 10cm that is SCAM or SPPA only (no mixed roots)
bg_clean %>% 
  filter(segment_top == 0 & species %in% c("scam", "sppa")) %>% 
  group_by(pot_no, species) %>% 
  summarize(bg_biomass10 = sum(weight)) %>% 
  spread(key = species, value = bg_biomass10) -> bg_10cm

# Merge in with full_data
merge(bg_10cm, full_data_merged, by.x = "pot_no") %>% 
  rename(bg_scam10 = scam, bg_sppa10 = sppa)-> full_data_top10

# Filter for just corn genotypes and where scam bg is greater than 0 and levels
# 1-4
full_data_top10 %>% 
  filter(location == "corn" & bg_scam10 > 0 & level < 5) -> data_top10_corn

# Calculate an aboveground to belowground ratio for the top 10cm and filter out
# those that don't make sense biologically
data_top10_corn %>% 
  mutate(rs10 = bg_scam10 / agb_scam) %>% 
  filter(rs10 < 3)-> data_top10_corn

# Look to see if there are any rs10 outliers
data_top10_corn %>% 
  ggplot(aes(x = rs10)) +
  geom_histogram()

# Make age predictor variable
data_top10_corn %>% 
  mutate(age = case_when(substr(genotype, 2, 2) == "a" ~ "ancestral",
                         T ~ "modern")) -> data_top10_corn

# Scale elevation
data_top10_corn %>% 
  mutate(elevation_sc = scale(elevation)[,1]) -> data_top10_corn

# Fit a linear mixed model with comp for the top 10cm of belowground biomass
bg10_mod_corn <- lmer(sqrt(bg_scam10) ~ weight_init + date_planted_grp + origin_lab +
                        (age + comp + co2 + salinity + elevation_sc)^5 + I(elevation_sc^2) +
                        (1|genotype), data = data_top10_corn)

plot(bg10_mod_corn)
qqnorm(resid(bg10_mod_corn))
qqline(resid(bg10_mod_corn))

# Do stepwise selection
lmerTest::step(bg10_mod_corn)

# Stepwise model
bg10_mod_corn_step <- lmer(sqrt(bg_scam10) ~ weight_init + comp + elevation_sc +
                             I(elevation_sc^2) + age*co2*salinity +
                             (1 | genotype), data = data_top10_corn)

# Look at fixed effects
plot_model(bg10_mod_corn_step, type = "emm", terms = c("co2", "age", "salinity"))
plot_model(bg10_mod_corn_step, type = "emm", terms = "elevation_sc[all]")
plot_model(bg10_mod_corn_step, type = "emm", terms = "comp")

# Random effects
MuMIn::r.squaredGLMM(bg10_mod_corn_step)
plot_model(bg10_mod_corn_step, type = "re")

# See if we can add any random slopes

bg10_mod_corn_step1 <- update(bg10_mod_corn_step, .~.-(1|genotype) + (1+elevation_sc|genotype))
bg10_mod_corn_step2 <- update(bg10_mod_corn_step, .~.-(1|genotype) + (1+co2|genotype))
bg10_mod_corn_step3 <- update(bg10_mod_corn_step, .~.-(1|genotype) + (1+salinity|genotype)) # *
bg10_mod_corn_step4 <- update(bg10_mod_corn_step, .~.-(1|genotype) + (1+comp|genotype))
bg10_mod_corn_step5 <- update(bg10_mod_corn_step, .~.-(1|genotype) + (1+co2*salinity|genotype)) # *
bg10_mod_corn_step6 <- update(bg10_mod_corn_step, .~.-(1|genotype) + (1+co2+salinity|genotype))

anova(bg10_mod_corn_step, bg10_mod_corn_step3, refit = F)
anova(bg10_mod_corn_step, bg10_mod_corn_step5, refit = F)

# Plot of fixed effects
data_top10_corn %>% 
  ggplot(aes(x = age, y = bg_scam10, color = age)) +
  geom_boxplot() +
  geom_jitter(height = 0, width = 0.1) +
  facet_grid(salinity ~ co2)

emmip(bg10_mod_corn_step, co2 ~ salinity | age, CIs = T, type = "response",
      style = "factor") +
  geom_jitter(aes(x = salinity, y = bg_scam10, colour = co2), 
              data = data_top10_corn, pch = 16, width = 0.1, alpha = 0.2) +
  facet_wrap(~ age) +
  theme_bw()

# Do the same with root:shoot
rs10_mod_corn <- lmer(log(rs10) ~ weight_init + date_planted_grp + origin_lab +
                        (age + comp + co2 + salinity + elevation_sc)^5 + I(elevation_sc^2) +
                        (1|genotype), data = data_top10_corn)


# Do stepwise selection
lmerTest::step(rs10_mod_corn)

# Stepwise model
rs10_mod_corn_step <- lmer(log(rs10) ~ weight_init + date_planted_grp + salinity +
                             elevation_sc + I(elevation_sc^2) + (1 | genotype), data = data_top10_corn)

emmip(rs10_mod_corn_step, ~ salinity, CIs = T, type = "response",
      style = "factor") +
  geom_jitter(aes(x = salinity, y = rs10), 
              data = data_top10_corn, pch = 16, width = 0.1, alpha = 0.2) +
  theme_bw()
