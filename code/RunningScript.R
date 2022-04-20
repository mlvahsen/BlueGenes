# Read in libraries 
library(here); library(tidyverse); library(emmeans);
library(lme4); library(lmerTest); library(blme);
library(brglm); library(lmtest); library(sjPlot);
library(kableExtra); library(geomtextpath); library(patchwork)

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
extinct_mod_nocomp <- glm(survive ~ weight_init + date_cloned_grp + origin_lab + (co2 + salinity + elevation + age + location)^3 +
                            I(elevation^2) + genotype, data = full_data_nocomp, family = "binomial")

# Random effect estimated at zero
summary(extinct_mod_nocomp)
# Looks like both are estimated to be zero

# Now fit a glm without random effects
extinct_mod_nocomp_fixed <- glm(survive ~ weight_init + date_cloned_grp + origin_lab + (co2 + salinity + elevation + age + location)^5 +
                              I(elevation^2), data = full_data_nocomp, family = "binomial")

car::Anova(extinct_mod_nocomp_fixed)

# No 4-way or 5-way interactions are significant so drop to 3-way model
extinct_mod_nocomp_fixed3 <- glm(survive ~ weight_init + origin_lab + date_cloned_grp + (co2 + salinity + elevation + age + location)^3 +
                                  I(elevation^2), data = full_data_nocomp, family = "binomial")

# Check significance of terms
car::Anova(extinct_mod_nocomp_fixed3)

# There is some complete separation (i.e. all reps with the same covariate
# combinations have the same predicted value) at low elevations where everything
# went extinct. So we should fit this with a bias reduction method using the
# package 'brglm'

extinct_mod_nocomp_fixed3_BR <- brglm(survive ~ weight_init + origin_lab + date_cloned_grp +
                                        (co2 + salinity + elevation + age + location)^3 +
                                      I(elevation^2), data = full_data_nocomp,
                                    family = binomial(logit))

# Do all significance tests as likelihood ratio tests following the principle of
# marginality

##
# 3-way interactions
##

extinct_tab <- matrix(NA, 29, 4)

get_lr_results <- function(large_model, small_model, rep, term){
  out <- lrtest(large_model, small_model)
  df <- abs(out$Df[2])
  chisq <- out$Chisq[2]
  p <- out$`Pr(>Chisq)`[2]
  return(c(term, df, chisq, p))
}

# co2:salinity:age (ns)
extinct_mod_noCSA <- update(extinct_mod_nocomp_fixed3_BR, .~.-co2:salinity:age)
extinct_tab[1,] <- get_lr_results(extinct_mod_nocomp_fixed3_BR, extinct_mod_noCSA, 1, "co2:salinity:age")

# co2:salinity:location (ns)
extinct_mod_noCSL <- update(extinct_mod_nocomp_fixed3_BR, .~.-co2:salinity:location)
extinct_tab[2,] <- get_lr_results(extinct_mod_nocomp_fixed3_BR, extinct_mod_noCSL, 2, "co2:salinity:location")

# co2:salinity:elevation (ns)
extinct_mod_noCSE <- update(extinct_mod_nocomp_fixed3_BR, .~.-co2:salinity:elevation)
extinct_tab[3,] <- get_lr_results(extinct_mod_nocomp_fixed3_BR, extinct_mod_noCSE, 3, "co2:salinity:elevation")

# co2:age:location (*)
extinct_mod_noCAL <- update(extinct_mod_nocomp_fixed3_BR, .~.-co2:age:location)
extinct_tab[4,] <- get_lr_results(extinct_mod_nocomp_fixed3_BR, extinct_mod_noCAL, 4, "co2:age:location")

# co2:age:elevation (ns)
extinct_mod_noCAE <- update(extinct_mod_nocomp_fixed3_BR, .~.-co2:age:elevation)
extinct_tab[5,] <- get_lr_results(extinct_mod_nocomp_fixed3_BR, extinct_mod_noCAE, 5, "co2:age:elevation")

# co2:location:elevation (ns)
extinct_mod_noCLE <- update(extinct_mod_nocomp_fixed3_BR, .~.-co2:location:elevation)
extinct_tab[6,] <- get_lr_results(extinct_mod_nocomp_fixed3_BR, extinct_mod_noCLE, 6, "co2:location:elevation")

# salinity:age:location (ns)
extinct_mod_noSAL <- update(extinct_mod_nocomp_fixed3_BR, .~.-salinity:age:location)
extinct_tab[7,] <- get_lr_results(extinct_mod_nocomp_fixed3_BR, extinct_mod_noSAL, 7, "salinity:age:location")

# salinity:age:elevation (*)
extinct_mod_noSAE <- update(extinct_mod_nocomp_fixed3_BR, .~.-salinity:age:elevation)
extinct_tab[8,] <- get_lr_results(extinct_mod_nocomp_fixed3_BR, extinct_mod_noSAE, 8, "salinity:age:elevation")

# salinity:location:elevation (ns)
extinct_mod_noSLE <- update(extinct_mod_nocomp_fixed3_BR, .~.-salinity:location:elevation)
extinct_tab[9,] <- get_lr_results(extinct_mod_nocomp_fixed3_BR, extinct_mod_noSLE, 9, "salinity:location:elevation")

# age:location:elevation (ns)
extinct_mod_noALE <- update(extinct_mod_nocomp_fixed3_BR, .~.-age:location:elevation)
extinct_tab[10,] <- get_lr_results(extinct_mod_nocomp_fixed3_BR, extinct_mod_noALE, 10, "age:location:elevation")

##
# 2-way interactions
##

# co2:salinity (ns)
extinct_mod_forCS <- brglm(survive ~ weight_init + origin_lab + date_cloned_grp +
                              (co2 + salinity + elevation + age + location)^2 +
                              (co2 + elevation + age + location)^3 +
                             (salinity + elevation + age + location)^3 +
                              I(elevation^2), data = full_data_nocomp,
                            family = binomial(logit))
extinct_mod_noCS <- update(extinct_mod_forCS, .~.-co2:salinity)
extinct_tab[11,] <- get_lr_results(extinct_mod_forCS, extinct_mod_noCS, 11, "co2:salinity")

# co2:age (ns)
extinct_mod_forCA <- brglm(survive ~ weight_init + origin_lab +date_cloned_grp +
                             (co2 + salinity + elevation + age + location)^2 +
                             (co2 + elevation + salinity + location)^3 +
                             (salinity + elevation + age + location)^3 +
                             I(elevation^2), data = full_data_nocomp,
                           family = binomial(logit))
extinct_mod_noCA <- update(extinct_mod_forCA, .~.-co2:age)
extinct_tab[12,] <- get_lr_results(extinct_mod_forCA, extinct_mod_noCA, 12, "co2:age")

# co2:location (ns)
extinct_mod_forCL <- brglm(survive ~ weight_init + origin_lab + date_cloned_grp +
                             (co2 + salinity + elevation + age + location)^2 +
                             (co2 + elevation + salinity + age)^3 +
                             (salinity + elevation + age + location)^3 +
                             I(elevation^2), data = full_data_nocomp,
                           family = binomial(logit))
extinct_mod_noCL <- update(extinct_mod_forCL, .~.-co2:location)
extinct_tab[13,] <- get_lr_results(extinct_mod_forCL, extinct_mod_noCL, 13, "co2:location")

# co2:elevation (ns)
extinct_mod_forCE <- brglm(survive ~ weight_init + origin_lab + date_cloned_grp +
                             (co2 + salinity + elevation + age + location)^2 +
                             (co2 + location + salinity + age)^3 +
                             (salinity + elevation + age + location)^3 +
                             I(elevation^2), data = full_data_nocomp,
                           family = binomial(logit))
extinct_mod_noCE <- update(extinct_mod_forCE, .~.-co2:elevation)
extinct_tab[14,] <- get_lr_results(extinct_mod_forCE, extinct_mod_noCE, 14, "co2:elevation")

# salinity:age (*)
extinct_mod_forSA <- brglm(survive ~ weight_init + origin_lab + date_cloned_grp +
                             (co2 + salinity + elevation + age + location)^2 +
                             (co2 + location + salinity + elevation)^3 +
                             (co2 + elevation + age + location)^3 +
                             I(elevation^2), data = full_data_nocomp,
                           family = binomial(logit))
extinct_mod_noSA <- update(extinct_mod_forSA, .~.-salinity:age)
extinct_tab[15,] <- get_lr_results(extinct_mod_forSA, extinct_mod_noSA, 15, "salinity:age")

# salinity:location (*)
extinct_mod_forSL <- brglm(survive ~ weight_init + origin_lab + date_cloned_grp +
                             (co2 + salinity + elevation + age + location)^2 +
                             (co2 + location + age + elevation)^3 +
                             (co2 + elevation + age + salinity)^3 +
                             I(elevation^2), data = full_data_nocomp,
                           family = binomial(logit))
extinct_mod_noSL <- update(extinct_mod_forSL, .~.-salinity:location)
extinct_tab[16,] <- get_lr_results(extinct_mod_forSL, extinct_mod_noSL, 16, "salinity:location")

# salinity:elevation (***)
extinct_mod_forSE <- brglm(survive ~ weight_init + origin_lab + date_cloned_grp +
                             (co2 + salinity + elevation + age + location)^2 +
                             (co2 + location + salinity + age)^3 +
                             (co2 + elevation + age + location)^3 +
                             I(elevation^2), data = full_data_nocomp,
                           family = binomial(logit))
extinct_mod_noSE <- update(extinct_mod_forSE, .~.-salinity:elevation)
extinct_tab[17,] <- get_lr_results(extinct_mod_forSE, extinct_mod_noSE, 17, "salinity:elevation")

# age:location (ns)
extinct_mod_forAL <- brglm(survive ~ weight_init + origin_lab + date_cloned_grp +
                             (co2 + salinity + elevation + age + location)^2 +
                             (co2 + elevation + salinity + age)^3 +
                             (co2 + elevation + salinity + location)^3 +
                             I(elevation^2), data = full_data_nocomp,
                           family = binomial(logit))
extinct_mod_noAL <- update(extinct_mod_forAL, .~.-age:location)
extinct_tab[18,] <- get_lr_results(extinct_mod_forAL, extinct_mod_noAL, 18, "age:location")

# age:elevation (*)
extinct_mod_forAE <- brglm(survive ~ weight_init + origin_lab + date_cloned_grp +
                             (co2 + salinity + elevation + age + location)^2 +
                             (co2 + location + salinity + elevation)^3 +
                             (co2 + age + salinity + location)^3 +
                             I(elevation^2), data = full_data_nocomp,
                           family = binomial(logit))
extinct_mod_noAE <- update(extinct_mod_forAE, .~.-age:elevation)
extinct_tab[19,] <- get_lr_results(extinct_mod_forAE, extinct_mod_noAE, 19, "age:elevation")

# location:elevation (*)
extinct_mod_forLE <- brglm(survive ~ weight_init + origin_lab + date_cloned_grp +
                             (co2 + salinity + elevation + age + location)^2 +
                             (co2 + age + salinity + elevation)^3 +
                             (co2 + age + salinity + location)^3 +
                             I(elevation^2), data = full_data_nocomp,
                           family = binomial(logit))
extinct_mod_noLE <- update(extinct_mod_forLE, .~.-location:elevation)
extinct_tab[20,] <- get_lr_results(extinct_mod_forLE, extinct_mod_noLE, 20, "location:elevation")

##
# Main terms
##

# co2 (*)
extinct_mod_forC <- brglm(survive ~ weight_init + origin_lab + date_cloned_grp +
                             (salinity + elevation + age + location)^3 + co2 +
                             I(elevation^2), data = full_data_nocomp,
                           family = binomial(logit))
extinct_mod_noC <- update(extinct_mod_forC, .~.-co2)
extinct_tab[21,] <- get_lr_results(extinct_mod_forC, extinct_mod_noC, 21, "co2")

# salinity (ns)
extinct_mod_forS <- brglm(survive ~ weight_init + origin_lab + date_cloned_grp +
                            (co2 + elevation + age + location)^3 + salinity +
                            I(elevation^2), data = full_data_nocomp,
                          family = binomial(logit))
extinct_mod_noS <- update(extinct_mod_forS, .~.-salinity)
extinct_tab[22,] <- get_lr_results(extinct_mod_forS, extinct_mod_noS, 22, "salinity")

# age (ns)
extinct_mod_forA <- brglm(survive ~ weight_init + origin_lab + date_cloned_grp +
                            (co2 + elevation + salinity + location)^3 + age +
                            I(elevation^2), data = full_data_nocomp,
                          family = binomial(logit))
extinct_mod_noA <- update(extinct_mod_forA, .~.-age)
extinct_tab[23,] <- get_lr_results(extinct_mod_forA, extinct_mod_noA, 23, "age")

# location (ns)
extinct_mod_forL <- brglm(survive ~ weight_init + origin_lab + date_cloned_grp +
                            (co2 + elevation + salinity + age)^3 + location +
                            I(elevation^2), data = full_data_nocomp,
                          family = binomial(logit))
extinct_mod_noL <- update(extinct_mod_forL, .~.-location)
extinct_tab[24,] <- get_lr_results(extinct_mod_forL, extinct_mod_noL, 24, "location")

# elevation (***)
extinct_mod_forE <- brglm(survive ~ weight_init + origin_lab + date_cloned_grp +
                            (co2 + location + salinity + age)^3 + elevation +
                            I(elevation^2), data = full_data_nocomp,
                          family = binomial(logit))
extinct_mod_noE <- update(extinct_mod_forE, .~.-elevation)
extinct_tab[25,] <- get_lr_results(extinct_mod_forE, extinct_mod_noE, 25, "elevation")

##
# Covariate terms
##

# weight_init (***)
extinct_mod_noW <- update(extinct_mod_nocomp_fixed3_BR, .~.-weight_init)
extinct_tab[26,] <- get_lr_results(extinct_mod_nocomp_fixed3_BR, extinct_mod_noW, 26, "weight_init")

# origin_lab (***)
extinct_mod_noO <- update(extinct_mod_nocomp_fixed3_BR, .~.-origin_lab)
extinct_tab[27,] <- get_lr_results(extinct_mod_nocomp_fixed3_BR, extinct_mod_noO, 27, "origin_lab")

# elevation^2 (***)
extinct_mod_noE2 <- update(extinct_mod_nocomp_fixed3_BR, .~.-I(elevation^2))
extinct_tab[28,] <- get_lr_results(extinct_mod_nocomp_fixed3_BR, extinct_mod_noE2, 28, "elevation^2")

# date_cloned_grp (***)
extinct_mod_noD <- update(extinct_mod_nocomp_fixed3_BR, .~.-date_cloned_grp)
extinct_tab[29,] <- get_lr_results(extinct_mod_nocomp_fixed3_BR, extinct_mod_noD, 29, "date_cloned_grp")

# Make tabular results into a tibble
colnames(extinct_tab) <- c("term", "df", "chisq", "p")
as_tibble(extinct_tab) %>% 
  mutate(rank = 1:29) %>% 
  arrange(-rank) %>% 
  select(-rank) %>% 
  mutate(chisq = round(as.numeric(chisq), 2),
         p = round(as.numeric(p),3)) %>% 
  mutate(sig = case_when(0.01 < p & p < 0.05 ~ "*",
                         0.001 < p & p < 0.01 ~ "**",
                         p < 0.001 ~ "***",
                         T ~ "")) -> extinct_tab_out

knitr::kable(extinct_tab_out, "simple")

##
# Make plots
##

age_colors <- c("#fb9a99", "#e31a1c")
loc_colors <- c("#cab2d6", "#6a3d9a")
#plot_shapes <- c(15, 17)

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
  scale_fill_manual(values = age_colors, labels = c("ancestral (1900-1960)", "modern (2000-2020)")) +
  scale_color_manual(values = age_colors, labels = c("ancestral (1900-1960)", "modern (2000-2020)")) +
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
  ggplot(aes(x = elevation, y = survive, color = location)) +
  geom_line() +
  facet_wrap(~salinity) +
  geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci, fill = location), alpha = 0.2, color = NA) +
  geom_jitter(data = plot_raw_data, aes(x = elevation, y = survive),
              height = 0.05, width = 0, alpha = 0.7) +
  scale_color_manual(values = loc_colors) +
  scale_fill_manual(values = loc_colors) +
  guides(fill=guide_legend(title="provenance"),
         color = guide_legend(title="provenance"))  +
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


plot_data_c %>% 
  mutate(provenance = location) %>% 
  ggplot(aes(x = co2, y = survive, color = age)) +
  geom_point(size = 2, position = pd) +
  geom_errorbar(aes(ymin = lower.ci, ymax = upper.ci), width = 0.2, position = pd) +
  facet_wrap(~provenance) +
  #geom_point(data = full_data_nocomp, aes(x = location, y = survive), alpha = 0.2,
             #position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0.05, dodge.width = 1)) +
  theme_bw() +
  theme(legend.position = "none") +
  guides(colour = "none",
         shape = "none") +
  labs(x = expression(CO[2]), y = "survival rate") +
  scale_color_manual(values = age_colors) -> plot_c 

png("figs/Fig1_Extinction.png", res = 300, units = "in", height = 5, width = 8.2)
plot_a + plot_c + plot_b + guide_area() + plot_annotation(tag_levels = 'a')+
  plot_layout(guides = "collect", widths = c(4,3)) & theme(legend.justification = "left")
dev.off()

# Calculate predicted probability of survival for differing levels of elevation
# (averaged across all other treatments)

emmeans(extinct_mod_nocomp_fixed3_BR, ~elevation,
        at = list(elevation = c(0.1,0.2)), type = "response")

## Does competition mediate this response? ####

# We can only look at this for Corn genotypes
bg_full %>% 
  filter(location == "corn") %>% 
  mutate(survive = case_when(agb_scam > 0 ~ 1,
                             T ~ 0)) %>% 
  mutate(age = case_when(grepl("ancestral", cohort) ~ "ancestral",
                         T ~ "modern")) -> corn_only

# Fit a logistic GLMM
extinct_mod_corn <- glmer(survive ~ weight_init + date_cloned_grp + origin_lab + 
                          (co2 + salinity + elevation + age + comp)^5 + I(elevation^2) +
                            (1|site_frame) + (1|genotype), data = corn_only, family = "binomial")

summary(extinct_mod_corn)
# random effects estimated to be zero

extinct_mod_corn_fixed <- glm(survive ~ weight_init + date_cloned_grp + origin_lab + 
                                  (co2 + salinity + elevation + age + comp)^4 + I(elevation^2),
                              data = corn_only, family = "binomial")

car::Anova(extinct_mod_corn_fixed)

# Need to fit brglm to deal with complete separation

extinct_mod_corn_BR <- brglm(survive ~ weight_init + date_cloned_grp + origin_lab + 
                             (co2 + salinity + elevation + age + comp)^3 + I(elevation^2),
                           data = corn_only, family = "binomial")

# Competition doesn't seem to drastically influence things (no significant
# effect or significant interactions)

# Just run the significance tests and we can put that table in the supplement

# Do all significance tests as likelihood ratio tests following the principle of
# marginality

##
# 3-way interactions
##

extinct_tab_comp <- matrix(NA, 29, 4)

# co2:salinity:age (ns)
extinct_mod_corn_noCSA <- update(extinct_mod_corn_BR, .~.-co2:salinity:age)
extinct_tab_comp[1,] <- get_lr_results(extinct_mod_corn_BR, extinct_mod_corn_noCSA, 1, "co2:salinity:age")

# co2:salinity:comp (ns)
extinct_mod_corn_noCSL <- update(extinct_mod_corn_BR, .~.-co2:salinity:comp)
extinct_tab_comp[2,] <- get_lr_results(extinct_mod_corn_BR, extinct_mod_corn_noCSL, 2, "co2:salinity:comp")

# co2:salinity:elevation (ns)
extinct_mod_corn_noCSE <- update(extinct_mod_corn_BR, .~.-co2:salinity:elevation)
extinct_tab_comp[3,] <- get_lr_results(extinct_mod_corn_BR, extinct_mod_corn_noCSE, 3, "co2:salinity:elevation")

# co2:age:comp (*)
extinct_mod_corn_noCAL <- update(extinct_mod_corn_BR, .~.-co2:age:comp)
extinct_tab_comp[4,] <- get_lr_results(extinct_mod_corn_BR, extinct_mod_corn_noCAL, 4, "co2:age:comp")

# co2:age:elevation (ns)
extinct_mod_corn_noCAE <- update(extinct_mod_corn_BR, .~.-co2:age:elevation)
extinct_tab_comp[5,] <- get_lr_results(extinct_mod_corn_BR, extinct_mod_corn_noCAE, 5, "co2:age:elevation")

# co2:comp:elevation (ns)
extinct_mod_corn_noCLE <- update(extinct_mod_corn_BR, .~.-co2:comp:elevation)
extinct_tab_comp[6,] <- get_lr_results(extinct_mod_corn_BR, extinct_mod_corn_noCLE, 6, "co2:comp:elevation")

# salinity:age:comp (ns)
extinct_mod_corn_noSAL <- update(extinct_mod_corn_BR, .~.-salinity:age:comp)
extinct_tab_comp[7,] <- get_lr_results(extinct_mod_corn_BR, extinct_mod_corn_noSAL, 7, "salinity:age:comp")

# salinity:age:elevation (*)
extinct_mod_corn_noSAE <- update(extinct_mod_corn_BR, .~.-salinity:age:elevation)
extinct_tab_comp[8,] <- get_lr_results(extinct_mod_corn_BR, extinct_mod_corn_noSAE, 8, "salinity:age:elevation")

# salinity:comp:elevation (ns)
extinct_mod_corn_noSLE <- update(extinct_mod_corn_BR, .~.-salinity:comp:elevation)
extinct_tab_comp[9,] <- get_lr_results(extinct_mod_corn_BR, extinct_mod_corn_noSLE, 9, "salinity:comp:elevation")

# age:comp:elevation (ns)
extinct_mod_corn_noALE <- update(extinct_mod_corn_BR, .~.-age:comp:elevation)
extinct_tab_comp[10,] <- get_lr_results(extinct_mod_corn_BR, extinct_mod_corn_noALE, 10, "age:comp:elevation")

##
# 2-way interactions
##

# co2:salinity (ns)
extinct_mod_forCS <- brglm(survive ~ weight_init + origin_lab + date_cloned_grp +
                             (co2 + salinity + elevation + age + comp)^2 +
                             (co2 + elevation + age + comp)^3 +
                             (salinity + elevation + age + comp)^3 +
                             I(elevation^2), data = corn_only,
                           family = binomial(logit))
extinct_mod_corn_noCS <- update(extinct_mod_forCS, .~.-co2:salinity)
extinct_tab_comp[11,] <- get_lr_results(extinct_mod_forCS, extinct_mod_corn_noCS, 11, "co2:salinity")

# co2:age (ns)
extinct_mod_forCA <- brglm(survive ~ weight_init + origin_lab + date_cloned_grp +
                             (co2 + salinity + elevation + age + comp)^2 +
                             (co2 + elevation + salinity + comp)^3 +
                             (salinity + elevation + age + comp)^3 +
                             I(elevation^2), data = corn_only,
                           family = binomial(logit))
extinct_mod_corn_noCA <- update(extinct_mod_forCA, .~.-co2:age)
extinct_tab_comp[12,] <- get_lr_results(extinct_mod_forCA, extinct_mod_corn_noCA, 12, "co2:age")

# co2:comp (ns)
extinct_mod_forCL <- brglm(survive ~ weight_init + origin_lab + date_cloned_grp +
                             (co2 + salinity + elevation + age + comp)^2 +
                             (co2 + elevation + salinity + age)^3 +
                             (salinity + elevation + age + comp)^3 +
                             I(elevation^2), data = corn_only,
                           family = binomial(logit))
extinct_mod_corn_noCL <- update(extinct_mod_forCL, .~.-co2:comp)
extinct_tab_comp[13,] <- get_lr_results(extinct_mod_forCL, extinct_mod_corn_noCL, 13, "co2:comp")

# co2:elevation (ns)
extinct_mod_forCE <- brglm(survive ~ weight_init + origin_lab +date_cloned_grp+
                             (co2 + salinity + elevation + age + comp)^2 +
                             (co2 + comp + salinity + age)^3 +
                             (salinity + elevation + age + comp)^3 +
                             I(elevation^2), data = corn_only,
                           family = binomial(logit))
extinct_mod_corn_noCE <- update(extinct_mod_forCE, .~.-co2:elevation)
extinct_tab_comp[14,] <- get_lr_results(extinct_mod_forCE, extinct_mod_corn_noCE, 14, "co2:elevation")

# salinity:age (*)
extinct_mod_forSA <- brglm(survive ~ weight_init + origin_lab + date_cloned_grp +
                             (co2 + salinity + elevation + age + comp)^2 +
                             (co2 + comp + salinity + elevation)^3 +
                             (co2 + elevation + age + comp)^3 +
                             I(elevation^2), data = corn_only,
                           family = binomial(logit))
extinct_mod_corn_noSA <- update(extinct_mod_forSA, .~.-salinity:age)
extinct_tab_comp[15,] <- get_lr_results(extinct_mod_forSA, extinct_mod_corn_noSA, 15, "salinity:age")

# salinity:comp (*)
extinct_mod_forSL <- brglm(survive ~ weight_init + origin_lab +date_cloned_grp +
                             (co2 + salinity + elevation + age + comp)^2 +
                             (co2 + comp + age + elevation)^3 +
                             (co2 + elevation + age + salinity)^3 +
                             I(elevation^2), data = corn_only,
                           family = binomial(logit))
extinct_mod_corn_noSL <- update(extinct_mod_forSL, .~.-salinity:comp)
extinct_tab_comp[16,] <- get_lr_results(extinct_mod_forSL, extinct_mod_corn_noSL, 16, "salinity:comp")

# salinity:elevation (***)
extinct_mod_forSE <- brglm(survive ~ weight_init + origin_lab +date_cloned_grp+
                             (co2 + salinity + elevation + age + comp)^2 +
                             (co2 + comp + salinity + age)^3 +
                             (co2 + elevation + age + comp)^3 +
                             I(elevation^2), data = corn_only,
                           family = binomial(logit))
extinct_mod_corn_noSE <- update(extinct_mod_forSE, .~.-salinity:elevation)
extinct_tab_comp[17,] <- get_lr_results(extinct_mod_forSE, extinct_mod_corn_noSE, 17, "salinity:elevation")

# age:comp (ns)
extinct_mod_forAL <- brglm(survive ~ weight_init + origin_lab +date_cloned_grp+
                             (co2 + salinity + elevation + age + comp)^2 +
                             (co2 + elevation + salinity + age)^3 +
                             (co2 + elevation + salinity + comp)^3 +
                             I(elevation^2), data = corn_only,
                           family = binomial(logit))
extinct_mod_corn_noAL <- update(extinct_mod_forAL, .~.-age:comp)
extinct_tab_comp[18,] <- get_lr_results(extinct_mod_forAL, extinct_mod_corn_noAL, 18, "age:comp")

# age:elevation (*)
extinct_mod_forAE <- brglm(survive ~ weight_init + origin_lab +date_cloned_grp+
                             (co2 + salinity + elevation + age + comp)^2 +
                             (co2 + comp + salinity + elevation)^3 +
                             (co2 + age + salinity + comp)^3 +
                             I(elevation^2), data = corn_only,
                           family = binomial(logit))
extinct_mod_corn_noAE <- update(extinct_mod_forAE, .~.-age:elevation)
extinct_tab_comp[19,] <- get_lr_results(extinct_mod_forAE, extinct_mod_corn_noAE, 19, "age:elevation")

# comp:elevation (*)
extinct_mod_forLE <- brglm(survive ~ weight_init + origin_lab +date_cloned_grp+
                             (co2 + salinity + elevation + age + comp)^2 +
                             (co2 + age + salinity + elevation)^3 +
                             (co2 + age + salinity + comp)^3 +
                             I(elevation^2), data = corn_only,
                           family = binomial(logit))
extinct_mod_corn_noLE <- update(extinct_mod_forLE, .~.-comp:elevation)
extinct_tab_comp[20,] <- get_lr_results(extinct_mod_forLE, extinct_mod_corn_noLE, 20, "comp:elevation")

##
# Main terms
##

# co2 (*)
extinct_mod_forC <- brglm(survive ~ weight_init + origin_lab +date_cloned_grp+
                            (salinity + elevation + age + comp)^3 + co2 +
                            I(elevation^2), data = corn_only,
                          family = binomial(logit))
extinct_mod_corn_noC <- update(extinct_mod_forC, .~.-co2)
extinct_tab_comp[21,] <- get_lr_results(extinct_mod_forC, extinct_mod_corn_noC, 21, "co2")

# salinity (ns)
extinct_mod_forS <- brglm(survive ~ weight_init + origin_lab +date_cloned_grp+
                            (co2 + elevation + age + comp)^3 + salinity +
                            I(elevation^2), data = corn_only,
                          family = binomial(logit))
extinct_mod_corn_noS <- update(extinct_mod_forS, .~.-salinity)
extinct_tab_comp[22,] <- get_lr_results(extinct_mod_forS, extinct_mod_corn_noS, 22, "salinity")

# age (ns)
extinct_mod_forA <- brglm(survive ~ weight_init + origin_lab +date_cloned_grp+
                            (co2 + elevation + salinity + comp)^3 + age +
                            I(elevation^2), data = corn_only,
                          family = binomial(logit))
extinct_mod_corn_noA <- update(extinct_mod_forA, .~.-age)
extinct_tab_comp[23,] <- get_lr_results(extinct_mod_forA, extinct_mod_corn_noA, 23, "age")

# comp (ns)
extinct_mod_forL <- brglm(survive ~ weight_init + origin_lab +date_cloned_grp+
                            (co2 + elevation + salinity + age)^3 + comp +
                            I(elevation^2), data = corn_only,
                          family = binomial(logit))
extinct_mod_corn_noL <- update(extinct_mod_forL, .~.-comp)
extinct_tab_comp[24,] <- get_lr_results(extinct_mod_forL, extinct_mod_corn_noL, 24, "comp")

# elevation (***)
extinct_mod_forE <- brglm(survive ~ weight_init + origin_lab +date_cloned_grp+
                            (co2 + comp + salinity + age)^3 + elevation +
                            I(elevation^2), data = corn_only,
                          family = binomial(logit))
extinct_mod_corn_noE <- update(extinct_mod_forE, .~.-elevation)
extinct_tab_comp[25,] <- get_lr_results(extinct_mod_forE, extinct_mod_corn_noE, 25, "elevation")

##
# Covariate terms
##

# weight_init (***)
extinct_mod_corn_noW <- update(extinct_mod_corn_BR, .~.-weight_init)
extinct_tab_comp[26,] <- get_lr_results(extinct_mod_corn_BR, extinct_mod_corn_noW, 26, "weight_init")

# origin_lab (***)
extinct_mod_corn_noO <- update(extinct_mod_corn_BR, .~.-origin_lab)
extinct_tab_comp[27,] <- get_lr_results(extinct_mod_corn_BR, extinct_mod_corn_noO, 27, "origin_lab")

# elevation^2 (***)
extinct_mod_corn_noE2 <- update(extinct_mod_corn_BR, .~.-I(elevation^2))
extinct_tab_comp[28,] <- get_lr_results(extinct_mod_corn_BR, extinct_mod_corn_noE2, 28, "elevation^2")

# date_cloned_grp (***)
extinct_mod_corn_noD <- update(extinct_mod_corn_BR, .~.-date_cloned_grp)
extinct_tab_comp[29,] <- get_lr_results(extinct_mod_corn_BR, extinct_mod_corn_noD, 29, "date_cloned_grp")

# Make tabular results into a tibble
colnames(extinct_tab_comp) <- c("term", "df", "chisq", "p")
as_tibble(extinct_tab_comp) %>% 
  mutate(rank = 1:29) %>% 
  arrange(-rank) %>% 
  select(-rank) %>% 
  mutate(chisq = round(as.numeric(chisq), 2),
         p = round(as.numeric(p),3)) %>% 
  mutate(sig = case_when(0.01 < p & p < 0.05 ~ "*",
                         0.001 < p & p < 0.01 ~ "**",
                         p < 0.001 ~ "***",
                         T ~ "")) -> extinct_tab_comp_out

knitr::kable(extinct_tab_comp_out, "simple")


# Competition did not seem to drive or mediate the likelihood of extinction in
# our experiment (no significant interactions with competition and no
# significant main effect). Careful with this because it's "confirming the
# null".


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

keep_model_terms_nocomp <- c("co2", "salinity", "elevation", "location", "age")
keep_model_terms_nocomp_null <- c("co2", "salinity", "elevation")

keep_model_terms_corn <- c("co2", "salinity", "elevation", "comp", "age")
keep_model_terms_corn_null <- c("co2", "salinity", "elevation", "comp")

# Most complex fixed effects model
agb_mod <- lmer(sqrt(agb_scam) ~ 
                  (age + location + co2 + salinity + elevation)^5 + I(elevation^2) +
                  (1|genotype) + (1|site_frame), data = traits_nocomp)

lmerTest::step(agb_mod, keep = keep_model_terms_nocomp)
# Try lmerTest::step for backwards model selection
step_agb <- get_model(lmerTest::step(agb_mod))

# This is the best model
agb_mod_step <- lmer(sqrt(agb_scam) ~ weight_init + date_cloned_grp + age +  
                       age*co2*salinity*elevation + I(elevation^2) + (1|genotype),
                      data = traits_nocomp)

# This is the final FIXED effects model. Check assumptions.
plot_model(agb_mod_step, type = "diag") # All looks good!

## Model selection (random slopes) ####

# Now we are going to try and add up to 2-way interactions as random slopes.

# If "boundary (singular) fit: see help('isSingular')" error than the model is
# too complicated 
agb_mod1 <- update(agb_mod_step, .~.+ (1+elevation|genotype:co2))#*
agb_mod2 <- update(agb_mod_step, .~.+ (1+elevation|genotype:salinity))
agb_mod3 <- update(agb_mod_step, .~.+(1|genotype:co2) + (1|genotype:salinity) + (1|genotype:co2:salinity))
agb_mod4 <- update(agb_mod_step, .~.+(0+elevation|genotype) + (1|genotype:co2))
agb_mod5 <- update(agb_mod_step, .~.+(0+elevation|genotype) + (1|genotype:salinity))
agb_mod6 <- update(agb_mod_step, .~.+(1|genotype:co2) + (1|genotype:salinity))
agb_mod7 <- update(agb_mod_step, .~.+(1|genotype:co2))#*
agb_mod8<- update(agb_mod_step, .~.+(0+elevation|genotype))
agb_mod9 <- update(agb_mod_step, .~.+(1|genotype:salinity))

# Compare random effects models that fit to the base model with just genotype as
# a random intercept 
anova(agb_mod_step, agb_mod1, refit = FALSE)
anova(agb_mod_step, agb_mod7, refit = FALSE)
# agb_mod1 is best

## Comparison to null model ####

# Now construct a null model that does not have age, location, or genotype in it
agb_mod_null <- lmer(sqrt(agb_scam) ~ weight_init + date_cloned_grp + origin_lab +
                  (co2 + salinity + elevation)^3 + I(elevation^2) +
                  (1|site_frame), data = traits_nocomp)

get_model(lmerTest::step(agb_mod_null))

# Fit reduced null model
agb_null_step <- lm(sqrt(agb_scam) ~ weight_init + co2 + salinity + 
                      elevation + I(elevation^2), data = traits_nocomp)

# Compare R2 from models
MuMIn::r.squaredGLMM(agb_null_step)
# Conditional R2 = 0.35
MuMIn::r.squaredGLMM(agb_mod1)
# Conditional R2 = 0.55

## Make plot of fixed effects ####

# The highest order interaction that we want to show is age*co2*salinity*elevation
agb_plot <- plot_model(agb_mod1, terms = c("elevation[all]", "co2", "age", "salinity"), type = "emm")

# Collect data from plot_model object
plot_agb_data <- tibble(elevation = c(agb_plot[[1]]$data$x, agb_plot[[2]]$data$x),
                        co2 = c(agb_plot[[1]]$data$group, agb_plot[[2]]$data$group),
                        age = c(agb_plot[[1]]$data$facet, agb_plot[[2]]$data$facet),
                        salinity = c(agb_plot[[1]]$data$panel, agb_plot[[2]]$data$panel),
                        agb_scam = c(agb_plot[[1]]$data$predicted, agb_plot[[2]]$data$predicted),
                        lower_ci = c(agb_plot[[1]]$data$conf.low, agb_plot[[2]]$data$conf.low),
                        upper_ci = c(agb_plot[[1]]$data$conf.high, agb_plot[[2]]$data$conf.high))

# Change the labels in raw data for plots
traits_nocomp %>% 
  mutate(age = case_when(age == "modern" ~ "descendant cohort (2000-2020)",
                         T ~ "ancestral cohort (1900-1950)")) %>% 
  mutate(salinity = case_when(salinity == "fresh" ~ "freshwater site (4ppt)",
                              T ~ "brackish site (6ppt)")) -> traits_nocomp_plot

# Create plot
plot_agb_data %>% 
  mutate(age = case_when(age == "modern" ~ "descendant cohort (2000-2020)",
                         T ~ "ancestral cohort (1900-1950)")) %>% 
  mutate(salinity = case_when(salinity == "fresh" ~ "freshwater site (4ppt)",
                              T ~ "brackish site (6ppt)")) %>% 
  ggplot(aes(x = elevation, y = agb_scam, color = co2)) +
  geom_point(data = traits_nocomp_plot, shape = 1, alpha = 0.6, stroke = 0.7) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = co2), alpha = 0.2, color = NA) +
  geom_textline(aes(label = co2), linewidth = 1.2, fontface = 2, hjust = 0.2) +
  facet_grid(salinity ~ age) +
  scale_color_manual(values = c("#b2df8a","#33a02c")) +
  scale_fill_manual(values = c("#b2df8a","#33a02c")) +
  ylab("aboveground biomass (g)") +
  xlab("elevation (m NAVD88)") +
  theme_bw() +
  theme(legend.position = "none")-> agb_plot_4way
  
ggsave(here("figs", "Fig2.png"), agb_plot_4way, height = 4.7, width = 5.6, units = "in")

# Extract type III Anova table
anova(agb_mod4)

## Calculate effect sizes for text ####
emmeans_agb <- summary(emmeans(agb_mod1, ~co2|salinity:age, at = list(elevation = 0.2), type = "response"))

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

## Random slopes figure ####

# Refit with summed contrasts for plotting (this doesn't change the model)
options(contrasts = c("contr.sum", "contr.poly"))
agb_mod1 <- lmer(sqrt(agb_scam) ~ weight_init + date_cloned_grp + age +  
       age*co2*salinity*elevation + I(elevation^2) + (1+elevation|genotype:co2) + (1|genotype),
     data = traits_nocomp)
options(contrasts = c("contr.treatment", "contr.poly"))
agb_slopes <- (coef(agb_mod1)$`genotype:co2`[,1] + coef(agb_mod1)$`genotype:co2`[,8]%*%t(seq(0.2,0.5, length.out = 20)) +
  coef(agb_mod1)$`genotype:co2`[,9]%*%t(seq(0.2,0.5, length.out = 20))^2)^2 

tibble(as.data.frame(agb_slopes)) %>% 
  mutate(genotype_co2 = rownames(coef(agb_mod1)$`genotype:co2`[,]),
         genotype = sub(":.*", "", genotype_co2),
         co2 = sub(".*:", "", genotype_co2)) %>% 
  gather(key = index, value = agb, V1:V20) %>% 
  mutate(elevation = rep(seq(0.2,0.5, length.out = 20), each = 61)) %>% 
  ggplot(aes(x = elevation, y = agb, group = genotype_co2, color = co2)) +
  geom_line() 

expand.grid(salinity = unique(traits_nocomp$salinity),
            date_cloned_grp = unique(traits_nocomp$date_cloned_grp),
            genotype = unique(row.names(ranef(agb_mod1)$genotype)),
            elevation = seq(0.2, 0.5, length.out = 20),
            co2 = unique(traits_nocomp$co2)) %>% 
  mutate(age = case_when(substr(genotype, 2, 2) == "a" ~ "ancestral",
                         T ~ "modern"),
         weight_init = mean(traits_nocomp$weight_init)) -> newData

newData %>% filter(genotype != "ca6" | co2 != "elevated") -> newData_ca6

predict(agb_mod1, newData_ca6) -> predictions

newData_ca6 %>% 
  mutate(predict_agb = predictions^2) %>% 
  group_by(elevation, co2, genotype, age) %>% 
  summarize(mean_agb = mean(predict_agb)) %>% 
  ggplot(aes(x = elevation, y = mean_agb, group = genotype)) +
  geom_line(alpha = 0.4) +
  facet_wrap(~co2) +
  ylab("predicted aboveground biomass (g)") +
  xlab("elevation (m NAVD88)")

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
# There were issues with data collection for pots 205 and 1700 (i.e. mislabeled
# or missing bags for sections of bgb) so these two are dropped from the
# analysis because we could not reconcile the issues.

## Fit model and model selection (fixed effects) ####

# Most complex model
rs_mod <- lmer(log(rs) ~ weight_init + date_cloned_grp + origin_lab +
                 (age + location + co2 + salinity + elevation)^5 + I(elevation^2) +
                 (1|genotype) + (1|site_frame), data = traits_nocomp_rs)

# Get the best stepwise model
get_model(lmerTest::step(rs_mod))


# This is the best selected model
rs_step_mod <- lmer(log(rs) ~ date_cloned_grp + co2 + location * elevation +
                      salinity * elevation + (1 | genotype), data = traits_nocomp_rs)

# Check assumptions
plot_model(rs_step_mod, type = "diag")

## Model selection (random slopes) ####

# Now we are going to try and add up to 2-way interactions as random slopes.

rs_mod1 <- update(rs_step_mod, .~.+ (1+elevation|genotype:co2))
rs_mod2 <- update(rs_step_mod, .~.+ (1+elevation|genotype:salinity))#*
rs_mod3 <- update(rs_step_mod, .~.+(1|genotype:co2) + (1|genotype:salinity) + (1|genotype:co2:salinity))#*
rs_mod4 <- update(rs_step_mod, .~.+(0+elevation|genotype) + (1|genotype:co2))
rs_mod5 <- update(rs_step_mod, .~.+(0+elevation|genotype) + (1|genotype:salinity))#*
rs_mod6 <- update(rs_step_mod, .~.+(1|genotype:co2) + (1|genotype:salinity))
rs_mod7 <- update(rs_step_mod, .~.+(1|genotype:co2))
rs_mod8<- update(rs_step_mod, .~.+(0+elevation|genotype))#*
rs_mod9 <- update(rs_step_mod, .~.+(1|genotype:salinity))#*

# Compare to the original models (all ns)
anova(rs_step_mod, rs_mod2, refit = F) # This is the best model, but not significant
anova(rs_step_mod, rs_mod3, refit = F)
anova(rs_step_mod, rs_mod5, refit = F)
anova(rs_step_mod, rs_mod8, refit = F)
anova(rs_step_mod, rs_mod9, refit = F)
# None are significant

## Comparison to null model ####

# Now construct a null model that does not have age, location, or genotype in it
rs_modelnull <- lmer(log(rs) ~ weight_init + date_cloned_grp + origin_lab +
                       (co2 + salinity + elevation)^3 + I(elevation^2) +
                       (1|site_frame),
                     data = traits_nocomp_rs)

get_model(lmerTest::step(rs_modelnull, keep = keep_model_terms_nocomp_null))


rs_mod_null_step <- lm(log(rs) ~ date_cloned_grp + co2 + salinity * elevation,
                  data = traits_nocomp_rs)

# Compare R2 from models
MuMIn::r.squaredGLMM(rs_mod_null_step)
# Conditional R2 = 0.61
MuMIn::r.squaredGLMM(rs_step_mod)
# Conditional R2 = 0.70

# Get overall significance of model terms
anova(rs_step_mod)


##########################################
## Root distribution parameter analysis ##
##########################################

## Data manipulation ####

# This is for the pots we are doing trait analysis for (no comp for right now
# for the same reasons as the root:shoot analysis)


## Fit model and model selection (fixed effects) ####

# Most complex model
beta_mod <- lmer(beta ~ weight_init + date_cloned_grp + origin_lab +
                   (age + location + co2 + salinity + elevation)^5 + I(elevation^2) +
                   (1|site_frame) + (1|genotype), data = traits_nocomp_rs)

# Try lmerTest::step
get_model(lmerTest::step(beta_mod))

beta_step_mod <- lmer(beta ~ location + salinity + elevation + (1 | genotype) +  
                        salinity:elevation + location:elevation, data = traits_nocomp_rs)

# Check assumptions
plot_model(beta_step_mod, type = "diag")

## Model selection (random slopes) ####

# Now we are going to try and add up to 2-way interactions as random slopes.

beta_mod1 <- update(beta_step_mod, .~.+(1+elevation|genotype:salinity)) # fit
beta_mod2 <- update(beta_step_mod, .~.+(1+elevation|genotype:co2))
beta_mod3 <- update(beta_step_mod, .~.+(0+elevation|genotype) + (1|genotype:salinity)) 
beta_mod4 <- update(beta_step_mod, .~.+ (1|genotype:salinity)+ (1|genotype:co2)) # fit
beta_mod5 <- update(beta_step_mod, .~.+(1|genotype:co2))
beta_mod6<- update(beta_step_mod, .~.+(0+elevation|genotype)) 
beta_mod7 <- update(beta_step_mod, .~.+(1|genotype:salinity)) # fit

# Compare to the original models
anova(beta_step_mod, beta_mod1, refit = F) # this is the best model
anova(beta_step_mod, beta_mod4, refit = F) 
anova(beta_step_mod, beta_mod7, refit = F)

# Final model is beta_mod1

## Comparison to null model ####

# Now construct a null model that does not have age, location, or genotype in it
beta_modelnull <- lmer(beta ~ weight_init + date_cloned_grp + origin_lab +
                         (co2 + salinity + elevation)^3 + I(elevation^2) +
                         (1|site_frame), data = traits_nocomp_rs)


get_model(lmerTest::step(beta_modelnull))

beta_modelnull_step <- lm(beta ~ salinity + elevation + co2*elevation + 
                            salinity*elevation, data = traits_nocomp_rs)

# Compare R2 from models
MuMIn::r.squaredGLMM(beta_modelnull_step)
# Conditional R2 = 0.47
MuMIn::r.squaredGLMM(beta_mod1)
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
bg_mod <- lmer(total_bg ~ weight_init + date_cloned_grp + origin_lab +
                 (age + location + co2 + salinity + elevation)^5 + I(elevation^2) +
                 (1|genotype) + (1|site_frame), data = traits_nocomp_rs)

# Try stepwise
get_model(lmerTest::step(bg_mod))

bg_step_mod <- lmer(total_bg ~ weight_init + date_cloned_grp + co2 +  
                      elevation + I(elevation^2) + (1 | genotype), data = traits_nocomp_rs)

## Model selection (random slopes) ####

bg_mod1 <- update(bg_step_mod, .~.+(1|genotype:co2))
bg_mod2<- update(bg_step_mod, .~.+(0+elevation|genotype)) 

# Compare to base model (ns)
anova(bg_step_mod, bg_mod1, refit = F)

## Comparison to null model ####

# Now construct a null model that does not have age, location, or genotype in it
bg_modelnull <- lmer(total_bg ~ weight_init + date_cloned_grp + origin_lab +
                       (co2 + salinity + elevation)^3 + I(elevation^2) +
                       (1|site_frame),
                     data = traits_nocomp_rs)

get_model(step(bg_modelnull))

# Fit stepwise null model
bg_modelnull_step <- lm(total_bg ~ weight_init + date_cloned_grp + co2 + 
                          elevation + I(elevation^2), data = traits_nocomp_rs)

# Compare R2 from models
MuMIn::r.squaredGLMM(bg_modelnull_step)
# Conditional R2 = 0.23
MuMIn::r.squaredGLMM(bg_step_mod)
# Conditional R2 = 0.33

# Get emmeans values for plots
bgb_EC_plot <- summary(emmeans::emmeans(bg_step_mod, ~elevation:co2, at = list(elevation = seq(0.156, 0.544, length.out = 50))))

tibble(co2 = bgb_EC_plot$co2,
       total_bg = bgb_EC_plot$emmean,
       elevation = bgb_EC_plot$elevation,
       lower = bgb_EC_plot$lower.CL,
       upper = bgb_EC_plot$upper.CL) %>% 
  ggplot(aes(x = elevation, y = total_bg, color = co2)) +
  geom_textline(label = bgb_EC_plot$co2, fontface = 2, linewidth = 1.2, size = 3) +
  theme_bw() +
  theme(legend.position = "none") +
  geom_point(data = traits_nocomp_rs, aes(x = elevation, y = total_bg, color = co2), shape = 1, stroke = 0.8, alpha = 0.7, size = 0.8) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = co2), alpha = 0.2, color = NA) +
  scale_color_manual(values = c("#b2df8a", "#33a02c")) +
  scale_fill_manual(values = c("#b2df8a", "#33a02c")) +
  ylab("belowground biomass (g)") +
  xlab("elevation (m NAVD88)") -> bg_plot
  
ggsave(here("figs", "SuppFig_bg.png"), bg_plot, height = 3, width = 4, units = "in")

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
                     (age + location + co2 + salinity + elevation_sc)^5 + I(elevation_sc^2) +
                     (1|site_frame) + (1|genotype), data = traits_nocomp)

summary(height_mod)
plot(height_mod) # There appears to be a huge outlier (obs # 186)
plot_model(height_mod, type = "diag")

# Look at overall differences and then drop and see if outlier mattered
car::Anova(height_mod)

traits_nocomp_noOut <- traits_nocomp[-186,]
height_mod_noOut <- lmer(height_scam_tot ~ weight_init + date_cloned_grp + origin_lab +
                           (age + location + co2 + salinity + elevation)^5 + I(elevation^2) +
                           (1|genotype) + (1|site_frame), data = traits_nocomp_noOut)
plot(height_mod_noOut)
car::Anova(height_mod_noOut) # That seemed to be driving somethings. So let's drop it!

# Try step
get_model(lmerTest::step(height_mod_noOut))

height_step_mod <- lmer(height_scam_tot ~ location * salinity + elevation + (1 | genotype), data = traits_nocomp_noOut)

## Model selection (random slopes) ####

height_mod1 <- update(height_step_mod, .~.+(0+elevation|genotype) + (1|genotype:salinity))
height_mod2<- update(height_step_mod, .~.+(0+elevation|genotype))
height_mod3 <- update(height_step_mod, .~.+(1|genotype:salinity))

# None of these fit

## Comparison to null model ####

# Now construct a null model that does not have age, location, or genotype in it
height_model_noOut_null <- lmer(height_scam_tot ~ weight_init + date_cloned_grp + origin_lab +
                                (co2 + salinity + elevation)^3 + I(elevation^2) +
                                (1|site_frame),
                              data = traits_nocomp_noOut)

get_model(step(height_model_noOut_null))

# Fit step model
height_model_noOut_null_step <- lm(height_scam_tot ~ weight_init + salinity + 
                                     elevation, data = traits_nocomp_noOut)

# Compare R2 from models
MuMIn::r.squaredGLMM(height_model_noOut_null)
# Adjusted R2 = 0.39
MuMIn::r.squaredGLMM(height_step_mod)
# Conditional R2 = 0.59

## Fixed effects plot for supplement ####
# Create fixed effects model for supplement

# Get emmeans values for plots
height_ES_plot <- summary(emmeans::emmeans(height_step_mod, ~elevation:salinity, at = list(elevation = seq(0.156, 0.544, length.out = 50))))

tibble(height_scam_tot = height_ES_plot$emmean,
       elevation = height_ES_plot$elevation,
       lower = height_ES_plot$lower.CL,
       upper = height_ES_plot$upper.CL,
       salinity = height_ES_plot$salinity) %>% 
  ggplot(aes(x = elevation, y = height_scam_tot, color = salinity)) +
  geom_textline(label = rep(c("freshwater site", "brackish site"), each = 50), linewidth = 1.2, fontface = 2, size = 3) +
  theme_bw() +
  geom_point(data = traits_nocomp_noOut, aes(x = elevation, y = height_scam_tot, color = salinity), shape = 1, stroke = 0.8, alpha = 0.7, size = 0.8) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = salinity), alpha = 0.2, color = NA) +
  ylab("mean stem height (cm)") +
  xlab("elevation (m NAVD88)") +
  scale_color_manual(values = c("#fdbf6f", "#ff7f00")) +
  scale_fill_manual(values = c("#fdbf6f", "#ff7f00")) +
  theme(legend.position = "none") -> height_plot

ggsave(here("figs", "SuppFig_height.png"), height_plot, height = 3, width = 4, units = "in")

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
                    (1|genotype) + (1|site_frame), data = traits_nocomp)

# Step mod

get_model(lmerTest::step(width_mod))

width_step_mod <- lmer(width_scam_mid ~ date_cloned_grp + location + salinity +  
                         elevation + (1 | genotype) + location:salinity + location:elevation, data = traits_nocomp)

# This is the final FIXED effects model. Check assumptions.
plot_model(width_step_mod, type = "diag") # All looks good!

## Model selection (random slopes) ####
# Now we are going to try and add up to 2-way interactions as random slopes.

width_mod1<- update(width_step_mod, .~.+(0+elevation|genotype))
width_mod2 <- update(width_step_mod, .~.+(1|genotype:salinity))

# Compare to base model (ns)
anova(width_step_mod, width_mod1, refit = F)
anova(width_step_mod, width_mod2, refit = F)

## Comparison to null model ####

# Now construct a null model that does not have age, location, or genotype in it
width_model_null <- lmer(width_scam_mid ~ weight_init + date_cloned_grp + origin_lab +
                         (co2 + salinity + elevation)^5 + I(elevation^2) +
                         (1|site_frame),
                       data = traits_nocomp)

get_model(step(width_model_null))

width_model_null_step <- lm(width_scam_mid ~ weight_init + salinity + 
                              elevation, traits_nocomp)

# Compare R2 from models
MuMIn::r.squaredGLMM(width_model_null_step)
# Adjusted R2 = 0.42
MuMIn::r.squaredGLMM(width_step_mod)
# Conditional R2 = 0.63

## Fixed effects plot for the supplement #### 

# Need to show interactions of provenance by salinity and provenance by
# elevation

# Get emmeans values for plots
width_EL_plot <- summary(emmeans::emmeans(width_step_mod, ~elevation:location, at = list(elevation = seq(0.156, 0.544, length.out = 50))))
width_SL_plot <- summary(emmeans::emmeans(width_step_mod, ~salinity:location))

tibble(width_scam_mid = width_EL_plot$emmean,
       elevation = width_EL_plot$elevation,
       lower = width_EL_plot$lower.CL,
       upper = width_EL_plot$upper.CL,
       location = rep(c("corn", "kirkpatrick"), each = 50)) %>% 
  ggplot(aes(x = elevation, y = width_scam_mid, color = location)) +
  geom_line(size = 1.2) +
  theme_bw() +
  geom_point(data = traits_nocomp, aes(x = elevation, y = width_scam_mid, color = location), shape = 1,
             stroke = 0.8, alpha = 0.6, size = 0.8) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = location), alpha = 0.2, color = NA) +
  ylab("mean stem width (mm)") +
  xlab("elevation (m NAVD88)") + 
  scale_color_manual(values = c("#cab2d6", "#6a3d9a")) +
  scale_fill_manual(values = c("#cab2d6", "#6a3d9a")) +
  theme(legend.position = "none") -> width_plot_a

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
  ylab("") +
  scale_color_manual(values = c("#cab2d6", "#6a3d9a")) +
  labs(color = "provenance") -> width_plot_b

width_plot_a + width_plot_b + plot_layout(widths = c(3,2), guides = "collect") +
  plot_annotation(tag_levels = "a") -> width_plot

ggsave(here("figs", "SuppFig_width.png"), width_plot, height = 3, width = 7, units = "in")

###########################
## Stem density analysis ##
###########################

## Fit model ####
# Stem density uses the same data as AGB
density_mod <- lmer(dens_scam_live ~ weight_init + date_cloned_grp + origin_lab +
                      (age + location + co2 + salinity + elevation)^5 + I(elevation^2) +
                      (1|genotype) + (1|site_frame), data = traits_nocomp)

get_model(lmerTest::step(density_mod))

density_step_mod <- lmer(dens_scam_live ~ co2 + elevation + I(elevation^2) +
                           (1 | genotype), data = traits_nocomp)

# Check for random slopes

density_mod1 <- update(density_step_mod, .~.+(0+elevation|genotype))
density_mod2 <- update(density_step_mod, .~.+(1|genotype:co2))

# Compare with base genotype model
anova(density_step_mod, density_mod2, refit = F)

# Compare with null model

density_mod_null <- lmer(dens_scam_live ~ weight_init + date_cloned_grp + origin_lab +
                           (co2 + salinity + elevation)^3 + I(elevation^2) +
                           (1|site_frame), data = traits_nocomp)

get_model(step(density_mod_null))

# Fit stepwise null model
density_mod_null_step <- lm(dens_scam_live ~ co2 + elevation + I(elevation^2), 
   data = traits_nocomp)

# Compare models
MuMIn::r.squaredGLMM(density_mod_null_step)
# 0.19
MuMIn::r.squaredGLMM(density_step_mod)
# 0.44

## Fixed effects plot for supplement ####
density_EC_plot <- summary(emmeans(density_step_mod, ~co2:elevation, at = list(elevation = seq(0.156, 0.544, length.out = 50))))

tibble(dens_scam_live = density_EC_plot$emmean,
  elevation = density_EC_plot$elevation,
  co2 = density_EC_plot$co2,
  lower = density_EC_plot$lower.CL,
  upper = density_EC_plot$upper.CL) %>% 
  ggplot(aes(x = elevation, y = dens_scam_live, color = co2)) +
  geom_point(data = traits_nocomp, aes(x = elevation, y = dens_scam_live, color = co2), shape = 1,
             stroke = 0.8, alpha = 0.6, size = 0.8) +
  geom_textline(linewidth = 1.2, label = rep(c("ambient","elevated"), each = 50), hjust = "ymax", fontface = 2, size = 3) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = co2), color = NA, alpha = 0.2) +
  ylab("stem density") +
  xlab("elevation (m NAVD88)") +
  scale_color_manual(values = c("#b2df8a","#33a02c")) +
  scale_fill_manual(values = c("#b2df8a","#33a02c")) +
  theme_bw() +
  theme(legend.position = "none") -> density_plot

ggsave(here("figs", "SuppFig_density.png"), density_plot, height = 3, width = 4, units = "in")


#############################

## Get final models for each trait ####

list(
# AGB no comp
agb_mod_step = agb_mod1,
# AGB no comp null
agb_null_step = agb_null_step,
# RS no comp
rs_step_mod = rs_step_mod,
# RS no comp null
rs_mod_null_step = rs_mod_null_step,
# BGB no comp
bg_step_mod = bg_step_mod,
# BGB no comp null
bg_modelnull_step = bg_modelnull_step,
# Beta no comp 
beta_step_mod = beta_mod1,
# Beta no comp null
beta_modelnull_step = beta_modelnull_step,
# Height no comp
height_step_mod = height_step_mod,
# Height no comp null
height_model_noOut_null_step = height_model_noOut_null_step,
# Width no comp
width_step_mod = width_step_mod,
# Width no comp null
width_model_null_step = width_model_null_step,
# Density no comp
density_step_mod = density_step_mod,
# Density no comp null
density_mod_null_step = density_mod_null_step) -> all_models

# Calculate R2 values for each model
lapply(all_models, MuMIn::r.squaredGLMM)

# Create functional to calculate and extract marginal and conditional ICC values
calc_icc <- function(x){
  as.numeric(performance::icc(x)[1])
}

calc_icc_conditional <- function(x){
  as.numeric(performance::icc(x)[2])
}

# Collect all final competition models
models_nocomp <- list(agb_mod = all_models$agb_mod_step,
                      bg_mod = all_models$bg_step_mod,
                      rs_mod = all_models$rs_step_mod,
                      width_mod = all_models$width_step_mod,
                      height_mod = all_models$height_step_mod,
                      density_mod = all_models$density_step_mod,
                      beta_mod = all_models$beta_step_mod)


###########################
## Extra stuff ############
###########################

## Calculate effect sizes for in-text ####

# Calculate mean elevation at level 4
bg_nocomp %>% 
  filter(level == 4) %>% 
  pull(elevation) %>% mean() -> level4_elevation

# Get predicted means at high level of sea-level rise for agb
summary(emmeans(models_nocomp$agb_mod, ~co2:salinity:age:elevation,
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
# 1.87033

# Calculate mean elevation at level 1
bg_nocomp %>% 
  filter(level == 1) %>% 
  pull(elevation) %>% mean() -> level1_elevation

# Get predicted means at low level of sea-level rise for rs
summary(emmeans(models_nocomp$rs_mod, ~location:elevation,
                at = list(elevation = level1_elevation), type = "response")) -> rs_pred_means

rs_pred_means$response[1]/rs_pred_means$response[2]

# Do the same thing for salinity at high levels of SLR
summary(emmeans(models_nocomp$rs_mod, ~salinity:elevation,
                at = list(elevation = level4_elevation), type = "response")) -> rs_pred_means2

rs_pred_means2$response[2]/rs_pred_means2$response[1]

# Width and height for provenance x salinity
summary(emmeans(models_nocomp$height_mod, ~salinity|location))$emmean -> height_pred_means

# Corn: fresh vs salt
height_pred_means[1]/height_pred_means[2]

summary(emmeans(models_nocomp$width_mod, ~salinity|location))$emmean -> width_pred_means

# Corn: fresh vs salt
width_pred_means[1]/width_pred_means[2]


## Fig S7 effect of SPPA biomass on SCAM density ####
bg_corn %>% 
  filter(comp == 1) %>% 
  ggplot(aes(x = agb_sppa, y = dens_scam_live, color = elevation)) +
  geom_point(aes(fill = elevation), size = 3, alpha = 0.7) + 
  ylab(expression(paste(italic("S. americanus "), "stem density"))) +
  xlab(expression(paste(italic("S. patens "), "aboveground biomass (g)"))) +
  theme_bw() -> FigS7

ggsave(here("figs", "FigS7.png"), FigS7, height = 3.5, width = 4.6, units = "in")

## Fig S4: eco-eco for root:shoot and beta ####

# Get emmeans values for plots
rs_ES_plot <- summary(emmeans::emmeans(models_nocomp$rs_mod, ~elevation:salinity, at = list(elevation = seq(0.156, 0.544, length.out = 50))))

tibble(rs = exp(rs_ES_plot$emmean),
       elevation = rs_ES_plot$elevation,
       lower = exp(rs_ES_plot$lower.CL),
       upper = exp(rs_ES_plot$upper.CL),
       salinity = rep(c("fresh", "salt"), each = 50)) %>% 
  ggplot(aes(x = elevation, y = rs, color = salinity)) +
  geom_point(data = traits_nocomp_rs, aes(x = elevation, y = rs, color = salinity), shape = 1,
             stroke = 0.8, alpha = 0.6, size = 0.8) +
  geom_textline(size = 3, fontface = "bold", label= rep(c("freshwater site", "brackish site"), each = 50),
                linewidth = 1.2, hjust = 0.2, textcolor = "black") +
  theme_bw() +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = salinity), alpha = 0.2, color = NA) +
  ylab("root-to-shoot ratio") +
  xlab("") + 
  scale_color_manual(values = c("#fdbf6f", "#ff7f00")) +
  scale_fill_manual(values = c("#fdbf6f", "#ff7f00")) +
  theme(legend.position = "none") -> rs_interaction

beta_ES_plot <- summary(emmeans::emmeans(models_nocomp$beta_mod, ~elevation:salinity, at = list(elevation = seq(0.156, 0.544, length.out = 50))))

tibble(beta = beta_ES_plot$emmean,
       elevation = beta_ES_plot$elevation,
       lower = beta_ES_plot$lower.CL,
       upper = beta_ES_plot$upper.CL,
       salinity = rep(c("fresh", "salt"), each = 50)) %>% 
  ggplot(aes(x = elevation, y = beta, color = salinity)) +
  geom_point(data = traits_nocomp_rs, aes(x = elevation, y = beta, color = salinity), shape = 1,
             stroke = 0.8, alpha = 0.6, size = 0.8) +
  geom_textline(size = 3, fontface = "bold", label= rep(c("freshwater site", "brackish site"), each = 50),
                linewidth = 1.2, hjust = 0.05, textcolor = "black") +
  theme_bw() +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = salinity), alpha = 0.2, color = NA) +
  ylab(expression(paste("root distribution parameter (", beta, ")"))) +
  xlab("elevation (m NAVD88)") + 
  scale_color_manual(values = c("#fdbf6f", "#ff7f00")) +
  scale_fill_manual(values = c("#fdbf6f", "#ff7f00")) +
  theme(legend.position = "none") -> beta_interaction

FigS4 <- rs_interaction / beta_interaction + plot_annotation(tag_levels = "a") & theme(plot.margin = margin(t = 0,r = 5,b = 0,l = 5))

ggsave(here("figs", "FigS4_ecoeco.png"), FigS4, height = 6, width = 4.5, units = "in")



## Fig S5: effect of CO2 on root-to-shoot ratio ####
# Get emmeans values for plots
rs_C_plot <- summary(emmeans::emmeans(models_nocomp$rs_mod, ~co2, type = "response"))

tibble(rs = rs_C_plot$response,
       lower = rs_C_plot$lower.CL,
       upper = rs_C_plot$upper.CL,
       co2 = c("ambient", "elevated")) %>% 
  ggplot(aes(x = co2, y = rs, color = co2)) +
  geom_jitter(data = traits_nocomp_rs, aes(x = co2, y = rs, color = co2), shape = 1,
             stroke = 0.8, alpha = 0.6, size = 0.8)+
  theme_bw() +
  geom_errorbar(aes(ymin = lower, ymax = upper), position=position_dodge(width=0.3), width = 0.2, color = "black") +
  geom_point(size = 3, position=position_dodge(width=0.3)) +
  geom_point(size = 3, position=position_dodge(width=0.3), shape = 1, stroke = 0.8, color = "black") +
  ylab("root-to-shoot ratio") +
  xlab("") + 
  scale_color_manual(values = c("#b2df8a", "#33a02c")) +
  scale_fill_manual(values = c("#b2df8a", "#33a02c")) +
  theme(legend.position = "none") -> rs_co2

ggsave(here("figs", "FigS5_rs_co2.png"), rs_co2, height = 3.5, width = 3, units = "in")

## Fig S6: effect of elevation on stem height ####
height_E_plot <- summary(emmeans::emmeans(models_nocomp$height_mod, ~elevation, at = list(elevation = seq(0.156, 0.544, length.out = 50))))

tibble(height_scam_tot = height_E_plot$emmean,
       elevation = height_E_plot$elevation,
       lower = height_E_plot$lower.CL,
       upper = height_E_plot$upper.CL) %>% 
  ggplot(aes(x = elevation, y = height_scam_tot)) +
  geom_point(data = traits_nocomp_noOut, aes(x = elevation, y = height_scam_tot), shape = 1,
             stroke = 0.8, alpha = 0.6, size = 0.8) +
  geom_line() +
  theme_bw() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  ylab("mean stem height (cm)") +
  xlab("elevation (m NAVD88)") -> height_elevation

ggsave(here("figs", "FigS6_height_elevation.png"), height_elevation, height = 4, width = 4, units = "in")
