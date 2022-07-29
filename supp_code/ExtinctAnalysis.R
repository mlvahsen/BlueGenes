
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

test <- brglm(survive ~ weight_init + origin_lab + date_cloned_grp +
                (co2 + salinity + elevation + age + location)^3 + genotype +
                I(elevation^2), data = full_data_nocomp,
              family = binomial(logit))

test_drop <- brglm(survive ~ weight_init + origin_lab + date_cloned_grp +
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


