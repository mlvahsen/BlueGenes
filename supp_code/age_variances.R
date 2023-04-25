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

genotype_means_abgDat <- traits_nocomp %>% 
  group_by(genotype, age) %>% 
  summarize(agb_scam = mean(agb_scam),
            dens_scam_live = mean(dens_scam_live),
            width_scam_mid = mean(width_scam_mid)) %>% 
  gather(key = trait, value = value, agb_scam:width_scam_mid)

genotype_means_height <- traits_nocomp_height %>% 
  group_by(genotype, age) %>% 
  summarize(height_scam_tot = mean(height_scam_tot)) %>% 
  gather(key = trait, value = value, height_scam_tot)

genotype_means_bgDat <- traits_nocomp_rs %>% 
  group_by(genotype, age) %>% 
  summarize(total_bg = mean(total_bg),
            rs = mean(rs),
            beta = mean(beta)) %>% 
  gather(key = trait, value = value, total_bg:beta)

all_genotype_means <- rbind(genotype_means_abgDat, genotype_means_bgDat, genotype_means_height)

all_genotype_means %>% 
  mutate(`age cohort` = case_when(age == "modern" ~ "descendant",
                               T ~ age)) -> all_genotype_means

# AGB graph
size_set <- 4

agb_annotate <- round(car::leveneTest(value ~ factor(age), data = all_genotype_means %>% filter(trait == "agb_scam"))$`Pr(>F)`[1],3)

all_genotype_means %>% 
  filter(trait == "agb_scam") %>% 
  ggplot(aes(x = value)) +
  geom_density(aes(fill = `age cohort`, color = `age cohort`), alpha = 0.2) +
  geom_rug(aes(color = `age cohort`)) +
  scale_color_manual(values = c("#fb9a99", "#e31a1c")) +
  scale_fill_manual(values = c("#fb9a99", "#e31a1c")) +
  theme_classic(base_size = 14) +
  labs(y = "probability density",
       x = "aboveground biomass (g)") +
  annotate("text", x = Inf, y = Inf, label = paste("p = ", agb_annotate),
           vjust = 2, hjust = 1.1, size = size_set)-> agb_graph

# BG graph
bg_annotate <- round(car::leveneTest(value ~ factor(age), data = all_genotype_means %>% filter(trait == "total_bg"))$`Pr(>F)`[1],3)

all_genotype_means %>% 
  filter(trait == "total_bg") %>% 
  ggplot(aes(x = value)) +
  geom_density(aes(fill = `age cohort`, color = `age cohort`), alpha = 0.2) +
  geom_rug(aes(color = `age cohort`)) +
  scale_color_manual(values = c("#fb9a99", "#e31a1c")) +
  scale_fill_manual(values = c("#fb9a99", "#e31a1c")) +
  theme_classic(base_size = 14) +
  labs(y = "probability density",
       x = "belowground biomass (g)") +
  annotate("text", x = Inf, y = Inf, label = paste("p = ", bg_annotate),
           vjust = 2, hjust = 1.1, size = size_set) -> bgb_graph

# RS graph
rs_annotate <- round(car::leveneTest(value ~ factor(age), data = all_genotype_means %>% filter(trait == "rs"))$`Pr(>F)`[1],3)

all_genotype_means %>% 
  filter(trait == "rs") %>% 
  ggplot(aes(x = value)) +
  geom_density(aes(fill = `age cohort`, color = `age cohort`), alpha = 0.2) +
  geom_rug(aes(color = `age cohort`)) +
  scale_color_manual(values = c("#fb9a99", "#e31a1c")) +
  scale_fill_manual(values = c("#fb9a99", "#e31a1c")) +
  theme_classic(base_size = 14) +
  labs(y = "probability density",
       x = "root-to-shoot ratio") +
  annotate("text", x = Inf, y = Inf, label = paste("p = ", rs_annotate),
           vjust = 2, hjust = 1.1, size = size_set) -> rs_graph

# Stem density graph
dens_annotate <- round(car::leveneTest(value ~ factor(age), data = all_genotype_means %>% filter(trait == "dens_scam_live"))$`Pr(>F)`[1],3)

all_genotype_means %>% 
  filter(trait == "dens_scam_live") %>% 
  ggplot(aes(x = value)) +
  geom_density(aes(fill = `age cohort`, color = `age cohort`), alpha = 0.2) +
  geom_rug(aes(color = `age cohort`)) +
  scale_color_manual(values = c("#fb9a99", "#e31a1c")) +
  scale_fill_manual(values = c("#fb9a99", "#e31a1c")) +
  theme_classic(base_size = 14) +
  labs(y = "probability density",
       x = "stem density") +
  annotate("text", x = Inf, y = Inf, label = paste("p = ", dens_annotate),
           vjust = 2, hjust = 1.1, size = size_set) -> dens_graph

# Stem height graph
height_annotate <- round(car::leveneTest(value ~ factor(age), data = all_genotype_means %>% filter(trait == "height_scam_tot"))$`Pr(>F)`[1],3)

all_genotype_means %>% 
  filter(trait == "height_scam_tot") %>% 
  ggplot(aes(x = value)) +
  geom_density(aes(fill = `age cohort`, color = `age cohort`), alpha = 0.2) +
  geom_rug(aes(color = `age cohort`)) +
  scale_color_manual(values = c("#fb9a99", "#e31a1c")) +
  scale_fill_manual(values = c("#fb9a99", "#e31a1c")) +
  theme_classic(base_size = 14) +
  labs(y = "probability density",
       x = "mean stem height (cm)") +
  annotate("text", x = Inf, y = Inf, label = paste("p = ", paste(height_annotate, "0", sep = "")),
           vjust = 2, hjust = 1.1, size = size_set) -> height_graph

# Stem width graph
width_annotate <- round(car::leveneTest(value ~ factor(age), data = all_genotype_means %>% filter(trait == "width_scam_mid"))$`Pr(>F)`[1],3)

all_genotype_means %>% 
  filter(trait == "width_scam_mid") %>% 
  ggplot(aes(x = value)) +
  geom_density(aes(fill = `age cohort`, color = `age cohort`), alpha = 0.2) +
  geom_rug(aes(color = `age cohort`)) +
  scale_color_manual(values = c("#fb9a99", "#e31a1c")) +
  scale_fill_manual(values = c("#fb9a99", "#e31a1c")) +
  theme_classic(base_size = 14) +
  labs(y = "probability density",
       x = "mean stem width (mm)") +
  annotate("text", x = Inf, y = Inf, label = paste("p = ", width_annotate),
           vjust = 2, hjust = 3.3, size = size_set) -> width_graph

# Beta graph
beta_annotate <- round(car::leveneTest(value ~ factor(age), data = all_genotype_means %>% filter(trait == "beta"))$`Pr(>F)`[1],3)

all_genotype_means %>% 
  filter(trait == "beta") %>% 
  ggplot(aes(x = value)) +
  geom_density(aes(fill = `age cohort`, color = `age cohort`), alpha = 0.2) +
  geom_rug(aes(color = `age cohort`)) +
  scale_color_manual(values = c("#fb9a99", "#e31a1c")) +
  scale_fill_manual(values = c("#fb9a99", "#e31a1c")) +
  theme_classic(base_size = 14) +
  labs(y = "probability density",
       x = "root distribution parameter") +
  annotate("text", x = Inf, y = Inf, label = paste("p = ", beta_annotate),
           vjust = 2, hjust = 3.3, size = size_set) -> beta_graph

# Bring plots together
png("figs/FigSXX_AgeVariances.png", height = 8.8, width = 12.3, units = "in", res = 300)
agb_graph + bgb_graph + rs_graph + rs_graph +
  beta_graph + height_graph + width_graph + dens_graph + 
  plot_layout(guides = "collect") + plot_annotation(tag_levels = "a")
dev.off()
