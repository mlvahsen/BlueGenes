# This code is for generating Supplemental Figures S2-S7. Note: Fig S1 can be
# generated using the code "supp_code/SurfaceWater_Nutrients.R"

library(tidyverse); library(patchwork); library(sjPlot); library(emmeans)

## Read in trait data ####
source("supp_code/CompileTraitData.R")
## Get tidal data to plot flooding on the x-axis ####

# Read in tidal data
tidal_1 <- read_csv(here("supp_data", "tidal_1.csv"))
tidal_2 <- read_csv(here("supp_data", "tidal_2.csv"))
tidal_3 <- read_csv(here("supp_data", "tidal_3.csv"))
tidal_4 <- read_csv(here("supp_data", "tidal_4.csv"))

# Bind all together
tidal_all <- rbind(tidal_1, tidal_2, tidal_3, tidal_4)

# Change column names
mean_tide <- mean(tidal_all$`Verified (m)`)

## Fig S2: trait distributions for extant experimental units ####

plot_multi_histogram <- function(df, feature, label_column) {
  plt <- ggplot(df, aes(x=eval(parse(text=feature)), fill=eval(parse(text=label_column)))) +
    #geom_histogram(alpha=0.7, position="identity", aes(y = ..density..), color="black") +
    geom_density(alpha=0.7) +
    #geom_vline(aes(xintercept=mean(eval(parse(text=feature)))), color="black", linetype="dashed", size=1) +
    labs(x=feature, y = "Density") + scale_fill_manual(values = c("#f1eef6","#d0d1e6","#a6bddb","#74a9cf","#08519c","#08306b"))
  plt + guides(fill=guide_legend(title="Flooding Level"))
}

bg_alive <- bg_full %>% filter(agb_scam > 0 & location != "blackwater") %>%
  mutate(rs = total_bg/agb_scam) %>% filter(rs < 15) %>% 
  mutate(level = case_when(level == 1 ~ "1 (least flooded)",
                           level == 6 ~ "6 (most flooded)",
                           T ~ level))

plot_multi_histogram(bg_alive, 'agb_scam', 'level') +
  xlab("aboveground biomass (g)") +
  theme_bw(base_size = 14) -> a

plot_multi_histogram(bg_alive, 'total_bg', 'level') +
  xlab("belowground biomass (g)") +
  theme_bw(base_size = 14) -> b

plot_multi_histogram(bg_alive, 'dens_scam_live', 'level') +
  xlab("stem density") +
  theme_bw(base_size = 14) -> c

plot_multi_histogram(bg_alive, 'width_scam_mid', 'level') +
  xlab("mean stem width (mm)") +
  theme_bw(base_size = 14) -> d

plot_multi_histogram(bg_alive, 'height_scam_tot', 'level') +
  xlab("mean stem height (cm)") +
  theme_bw(base_size = 14) -> e

plot_multi_histogram(bg_alive, 'rs', 'level') +
  xlab("root-to-shoot ratio") +
  theme_bw(base_size = 14) -> f

plot_multi_histogram(bg_alive, 'beta', 'level') +
  xlab("root distribution parameter") +
  theme_bw(base_size = 14) -> g

png("figs/FigS2_TraitDists.png", height = 8.8, width = 11.2, res = 300, units = "in")
a+b+f+g+e+d+c+ plot_layout(guides = "collect", nrow = 3) +
  plot_annotation(tag_levels = "a")
dev.off()
## Fig S3: eco-evo of extinction/survival ####
extinct_mod_nocomp_fixed3_BR <- readRDS("derived_data/extinct_mod.rda")

age_colors <- c("#fb9a99", "#e31a1c")
loc_colors <- c("#cab2d6", "#6a3d9a")

# A = elevation:age:salinity interaction
a <- plot_model(extinct_mod_nocomp_fixed3_BR, terms = c("elevation", "age", "salinity"), type = "emm")

plot_data_a <- tibble(survive = a$data$predicted,
                      elevation = a$data$x,
                      lower.ci = a$data$conf.low,
                      upper.ci = a$data$conf.high,
                      salinity = a$data$facet,
                      age = a$data$group)

# Create data frame and recode
bg_full %>% 
  filter(comp == 0 & location != "blackwater") -> full_data_nocomp

full_data_nocomp %>% 
  mutate(survive = ifelse(extinct == 0, 1, 0)) -> full_data_nocomp

plot_raw_data <- full_data_nocomp %>% 
  mutate(salinity = ifelse(salinity == "fresh", "freshwater site (4 ppt)", "brackish site (6 ppt)"),
         flooding = mean_tide - elevation)

plot_data_a %>% 
  mutate(salinity = ifelse(salinity == "fresh", "freshwater site (4 ppt)", "brackish site (6 ppt)"),
         flooding = mean_tide - elevation) %>% 
  ggplot(aes(x = flooding, y = survive, color = age)) +
  geom_line() +
  facet_wrap(~salinity) +
  geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci, fill = age), alpha = 0.2, color = NA) +
  geom_jitter(data = plot_raw_data, aes(x = flooding, y = survive), height = 0.05, width = 0, alpha = 0.3) +
  ylab("survival rate") + xlab("average inundation (m)") +
  scale_fill_manual(values = age_colors, labels = c("ancestral (1900-1950)", "descendant (2000-2020)")) +
  scale_color_manual(values = age_colors, labels = c("ancestral (1900-1950)", "descendant (2000-2020)")) +
  theme_bw(base_size = 14) +
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
  mutate(salinity = ifelse(salinity == "fresh", "freshwater site (4 ppt)", "brackish site (6 ppt)"),
         flooding = mean_tide - elevation) %>% 
  ggplot(aes(x = flooding, y = survive, color = location)) +
  geom_line() +
  facet_wrap(~salinity) +
  geom_ribbon(aes(ymin = lower.ci, ymax = upper.ci, fill = location), alpha = 0.2, color = NA) +
  geom_jitter(data = plot_raw_data, aes(x = flooding, y = survive),
              height = 0.05, width = 0, alpha = 0.7) +
  scale_color_manual(values = loc_colors) +
  scale_fill_manual(values = loc_colors) +
  guides(fill=guide_legend(title="provenance"),
         color = guide_legend(title="provenance"))  +
  theme_bw(base_size = 14) + theme(legend.background = element_rect(color = "gray47"))+
  labs(x = "average inundation (m)", y = "survival rate") -> plot_b

# C = location:age:co2 interaction
c <- plot_model(extinct_mod_nocomp_fixed3_BR, terms = c("location", "age", "co2"), type = "emm")

plot_data_c <- tibble(survive = c$data$predicted,
                      location = c("kirkpatrick", "corn")[c$data$x],
                      lower.ci = c$data$conf.low,
                      upper.ci = c$data$conf.high,
                      co2 = c$data$facet,
                      age = c("ancestral", "descendant")[c$data$group])

pd <- position_dodge(width = 0.4)

plot_data_c %>% 
  mutate(provenance = location) %>% 
  ggplot(aes(x = co2, y = survive, color = age)) +
  geom_point(size = 2, position = pd) +
  geom_errorbar(aes(ymin = lower.ci, ymax = upper.ci), width = 0.2, position = pd) +
  facet_wrap(~provenance) +
  #geom_point(data = full_data_nocomp, aes(x = location, y = survive), alpha = 0.2,
  #position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0.05, dodge.width = 1)) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none") +
  guides(colour = "none",
         shape = "none") +
  labs(x = expression(CO[2]), y = "survival rate") +
  scale_color_manual(values = age_colors) -> plot_c 

png("figs/FigS3.png", res = 300, units = "in", height = 5, width = 9)
plot_a + plot_c + plot_b + guide_area() + plot_annotation(tag_levels = 'a')+
  plot_layout(guides = "collect", widths = c(4,3)) & theme(legend.justification = "left")
dev.off()

# Calculate predicted probability of survival for differing levels of elevation
# (averaged across all other treatments)

emmeans(extinct_mod_nocomp_fixed3_BR, ~elevation,
        at = list(elevation = c(0.1,0.2)), type = "response")

## Fig S4: effect of salinity and elevation on r:s and bgb ####
rs_evo_model <- readRDS("derived_data/rs_model.rda")
bg_evo_model <- readRDS("derived_data/bg_model.rda")

# Get no comp data and create scaled elevation values
bg_full %>% 
  filter(comp == 0 & location != "blackwater" & level < 5 &
           agb_scam > 0) -> traits_nocomp

traits_nocomp %>% 
  mutate(elevation_sc = scale(elevation)[,1],
         elevation_sc2 = scale(elevation^2)[,1]) -> traits_nocomp

traits_nocomp %>% 
  mutate(rs = total_bg / agb_scam) %>% 
  filter(rs < 6 & pot_no !=165 & pot_no !=176) -> traits_nocomp_rs

# root-shoot-ratio ~ salinity x elevation 
rs_SE <- summary(emmeans(rs_evo_model, ~elevation_sc:salinity, weights = "proportional",
                         at = list(elevation_sc = seq(min(traits_nocomp$elevation_sc),
                                                      max(traits_nocomp$elevation_sc), length.out = 50)), type = "response"))

# Collect data from plot_model object
plot_rs_data <- tibble(elevation_sc = rs_SE$elevation_sc,
                       salinity = rs_SE$salinity,
                       rs = rs_SE$response,
                       lower_ci = rs_SE$lower.CL,
                       upper_ci = rs_SE$upper.CL) %>% 
  # unscale elevation data
  mutate(elevation = elevation_sc * sd(traits_nocomp_rs$elevation) + mean(traits_nocomp_rs$elevation),
         flooding = mean_tide - elevation,
         salinity = case_when(salinity == "fresh" ~ "freshwater site (4 ppt)",
                              T ~ "brackish site (6 ppt)"))

# Change the labels in raw data for plots
traits_nocomp_rs %>% 
  mutate(salinity = case_when(salinity == "fresh" ~ "freshwater site (4 ppt)",
                              T ~ "brackish site (6 ppt)"),
         flooding = mean_tide - elevation) -> traits_nocomp_rs_plot

plot_rs_data %>% 
  ggplot(aes(x = flooding, y = rs, color = salinity)) +
  geom_point(data = traits_nocomp_rs_plot, shape = 1, alpha = 0.6, stroke = 1) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = salinity), alpha = 0.2, color = NA) +
  geom_line(size = 1.2) +
  ylab("root-to-shoot ratio") +
  xlab("") +
  theme_bw(base_size = 14) +
  theme(legend.position = "top") +
  scale_color_manual(values = c("#ff7f00", "#fdbf6f")) +
  scale_fill_manual(values = c("#ff7f00", "#fdbf6f")) -> rs_SE_plot

# bgb ~ salinity x elevation 
bg_SE <- summary(emmeans(bg_evo_model, ~elevation_sc:salinity, weights = "proportional",
                         at = list(elevation_sc = seq(min(traits_nocomp$elevation_sc),
                                                      max(traits_nocomp$elevation_sc), length.out = 50)), type = "response"))

# Collect data from plot_model object
plot_bg_data <- tibble(elevation_sc = bg_SE$elevation_sc,
                       salinity = bg_SE$salinity,
                       total_bg = bg_SE$response,
                       lower_ci = bg_SE$lower.CL,
                       upper_ci = bg_SE$upper.CL) %>% 
  # unscale elevation data
  mutate(elevation = elevation_sc * sd(traits_nocomp_rs$elevation) + mean(traits_nocomp_rs$elevation),
         flooding = mean_tide - elevation,
         salinity = case_when(salinity == "fresh" ~ "freshwater site (4 ppt)",
                              T ~ "brackish site (6 ppt)"))

plot_bg_data %>% 
  ggplot(aes(x = flooding, y = total_bg, color = salinity)) +
  geom_point(data = traits_nocomp_rs_plot, shape = 1, alpha = 0.6, stroke = 1) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, fill = salinity), alpha = 0.2, color = NA) +
  geom_line(size = 1.2) +
  ylab("belowground biomass (g)") +
  xlab("average inundation (m)") +
  theme_bw(base_size = 14) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("#ff7f00", "#fdbf6f")) +
  scale_fill_manual(values = c("#ff7f00", "#fdbf6f")) -> bg_SE_plot

png("figs/FigS4.png", height = 8.7, width = 5.5, res = 300, units = "in")
rs_SE_plot / bg_SE_plot +
  plot_annotation(tag_levels = "a")
dev.off()

## Fig S5: effect of age and elevation on root depth distribution ####
beta_evo_model <- readRDS("derived_data/beta_model.rda")

# Get emmeans values for plots
beta_EA_plot <- emmeans::emmeans(beta_evo_model, ~age|elevation_sc,
                                 weights = "proportional",
                                 at = list(elevation_sc = seq(min(traits_nocomp_rs$elevation_sc),
                                                              max(traits_nocomp_rs$elevation_sc), length.out = 50)))

traits_nocomp_rs %>% 
  mutate(age = case_when(age == "ancestral" ~ "ancestral cohort (1900-1950)",
                         T ~ "descendant cohort (2000-2020)"),
         flooding = mean_tide - elevation) -> traits_nocomp_plot

summary(beta_EA_plot) %>% 
  mutate(age = case_when(age == "ancestral" ~ "ancestral cohort (1900-1950)",
                         T ~ "descendant cohort (2000-2020)")) %>% 
  mutate(beta = emmean,
         elevation = elevation_sc*sd(traits_nocomp_rs$elevation) + mean(traits_nocomp_rs$elevation),
         flooding = mean_tide - elevation) %>% 
  ggplot(aes(x = flooding, y = beta, color = age)) +
  geom_point(data = traits_nocomp_plot, aes(x = flooding, y = beta), shape = 1, stroke = 0.8, alpha = 0.7, size = 1.2) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL, fill = age), alpha = 0.2, color = NA) +
  ylab(expression(paste('root distribution parameter (', beta, ")"))) +
  xlab('average inundation (m)') +
  scale_color_manual(values = c("#fb9a99","#e31a1c")) +
  scale_fill_manual(values = c("#fb9a99","#e31a1c")) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none") -> beta_figA

# Also make figure where we translate root distribution parameter to rooting
# profile

# Get raw data
bgb_beta_analysis_zeros <- read_csv("derived_data/belowground_cumulative.csv")

bgb_beta_analysis_zeros %>% 
  select(pot_no, depth = depth_roots, cum_frac = cum_sum) -> raw_cumulative_fractions

merge(raw_cumulative_fractions, traits_nocomp_rs %>% select(pot_no, age)) %>% 
  mutate(age = case_when(age == "modern" ~ "descendant",
         T ~ age))-> belowground_rawdata

# Create vector of depths
depths <- seq(0, 70, 1)
# Extract mean for ancestral and descendant at mean elevation
pred_beta_age <- summary(emmeans(beta_evo_model, ~elevation_sc:age,
                                 at = list(elevation_sc = c(mean(traits_nocomp$elevation_sc),
                                                            min(traits_nocomp$elevation_sc),
                                                            max(traits_nocomp$elevation_sc)))))$emmean
# Combine these in a tibble
tibble(age = rep(c("ancestral", "descendant"), each = length(depths)*3),
       cum_frac = c(1-pred_beta_age[1]^depths, 1-pred_beta_age[2]^depths, 1-pred_beta_age[3]^depths,
                    1-pred_beta_age[4]^depths, 1-pred_beta_age[5]^depths, 1-pred_beta_age[6]^depths),
       depth = rep(depths,6),
       group = rep(c("mean_anc", "low_anc", "high_anc", "mean_mod", "low_mod", "high_mod"), each = length(depths)),
       elevation = rep(c("mean", "low", "high","mean", "low", "high"), each = length(depths))) %>%
  mutate(flooding = case_when(elevation == "low" ~ "high",
                              elevation == "high" ~ "low",
                              elevation == "mean" ~ "mean"),
         inundation = factor(flooding, levels = c("low", "mean", "high"))) %>% 
  ggplot(aes(x = depth, y = cum_frac, color = age, group = group)) +
  geom_line(aes(linetype = inundation), size = 1.2) +
  geom_line(data = belowground_rawdata, aes(x = depth, y = cum_frac, group = pot_no), alpha = 0.05) +
  coord_flip() +
  scale_y_reverse() +
  scale_x_reverse() +
  scale_linetype_manual(values = c("dotted", "solid", "dashed")) +
  scale_color_manual(values = c("#fb9a99","#e31a1c")) +
  xlab("depth below marsh surface (cm)") +
  ylab("cumulative proportion") +
  theme_bw(base_size = 14) -> beta_figB

beta_figA + beta_figB + plot_annotation(tag_levels = "a")-> beta_fig
ggsave(here("figs", "FigS5.png"), beta_fig, height = 4, width = 9, units = "in")

## Fig S6: effect of SPPA competition on stem height and width ####
width_corn_evo_model <- readRDS("derived_data/width_comp_model.rda")
height_corn_evo_model <- readRDS("derived_data/height_comp_model.rda")

# Filter for only Corn Island genotypes, within levels 1-4, for extant plants 
bg_full %>% 
  filter(location == "corn" & level < 5 & agb_scam > 0) -> traits_corn

# Get emmeans values for plots
width_comp_plot <- summary(emmeans::emmeans(width_corn_evo_model, ~comp,
                                            weights = "proportional"))

tibble(width_scam_mid = width_comp_plot$emmean,
       lower = width_comp_plot$lower.CL,
       upper = width_comp_plot$upper.CL,
       comp = c(0, 1)) %>% 
  ggplot(aes(x = as.factor(comp), y = width_scam_mid, color = as.factor(comp))) +
  geom_jitter(data = traits_corn, aes(x = as.factor(comp), y = width_scam_mid, color = as.factor(comp)), shape = 1,
              stroke = 0.8, alpha = 0.5, size = 0.8)+
  theme_bw(base_size = 14) +
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
height_comp_plot <- summary(emmeans::emmeans(height_corn_evo_model, ~comp),
                            weights = "proportional")

tibble(height_scam_tot = height_comp_plot$emmean,
       lower = height_comp_plot$lower.CL,
       upper = height_comp_plot$upper.CL,
       comp = c(0, 1)) %>% 
  ggplot(aes(x = as.factor(comp), y = height_scam_tot, color = as.factor(comp))) +
  geom_jitter(data = traits_corn, aes(x = as.factor(comp), y = height_scam_tot, color = as.factor(comp)), shape = 1,
              stroke = 0.8, alpha = 0.5, size = 0.8)+
  theme_bw(base_size = 14) +
  geom_errorbar(aes(ymin = lower, ymax = upper), position=position_dodge(width=0.3), width = 0.2, color = "black") +
  geom_point(size = 3, position=position_dodge(width=0.3)) +
  geom_point(size = 3, position=position_dodge(width=0.3), shape = 1, stroke = 0.8, color = "black") +
  ylab("mean stem height (cm)") +
  scale_color_manual(values = c("gray47", "black")) +
  scale_fill_manual(values = c("gray47", "black")) +
  theme(legend.position = "none") + 
  scale_x_discrete(labels=expression(paste("without ", italic("S. patens")), paste("with ", italic("S. patens")))) +
  xlab("") -> height_comp

FigS6 <- height_comp + width_comp + plot_annotation(tag_levels = "a")

ggsave(here("figs", "FigS6.png"), FigS6, height = 3.5, width = 7, units = "in")

## Fig S7: correlation between SPPA agb and SCAM stem density ####
traits_corn %>% 
  filter(comp == 1) %>% 
  mutate(sppa_code = ifelse(agb_sppa == 0, "dead", "alive"),
         flooding = mean_tide - elevation) %>% 
  ggplot(aes(x = agb_sppa, y = dens_scam_live, color = flooding, shape = sppa_code)) +
  geom_point(size = 4, alpha = 0.5, stroke = 2) +
  xlab(expression(paste(italic("S. patens"), " aboveground biomass (g)"))) +
  ylab(expression(paste(italic("S. americanus"), " stem density"))) +
  theme_bw(base_size = 14) +
  scale_shape_manual(values = c(19,4), name = "",
                     labels = c(expression(paste(italic("S. patens"), " alive")),
                                expression(paste(italic("S. patens"), " dead")))) +
  scale_color_gradient(high = "#810f7c", low = "#b3cde3") +
  labs(color = "inundation (m)")-> FigS7

png("figs/FigS7.png", height = 4.3, width = 6.3, res = 300, units = "in")
FigS7
dev.off()
