## Random intercept fig ####
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

# Create functional to calculate and extract adjusted and conditional ICC values
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

# Set number of bootstrap simulations
sims = 100

tibble(agb = as.numeric(bootMer(models_nocomp$agb_mod, FUN = calc_icc, nsim = sims)$t),
       bg = as.numeric(bootMer(models_nocomp$bg_mod, FUN = calc_icc, nsim = sims)$t),
       rs = as.numeric(bootMer(models_nocomp$rs_mod, FUN = calc_icc, nsim = sims)$t),
       width = as.numeric(bootMer(models_nocomp$width_mod, FUN = calc_icc, nsim = sims)$t),
       height = as.numeric(bootMer(models_nocomp$height_mod, FUN = calc_icc, nsim = sims)$t),
       density = as.numeric(bootMer(models_nocomp$density_mod, FUN = calc_icc, nsim = sims)$t),
       beta = as.numeric(bootMer(models_nocomp$beta_mod, FUN = calc_icc, nsim = sims)$t)) -> boot_samples

tibble(agb = as.numeric(bootMer(models_nocomp$agb_mod, FUN = calc_icc_conditional, nsim = sims)$t),
       bg = as.numeric(bootMer(models_nocomp$bg_mod, FUN = calc_icc_conditional, nsim = sims)$t),
       rs = as.numeric(bootMer(models_nocomp$rs_mod, FUN = calc_icc_conditional, nsim = sims)$t),
       width = as.numeric(bootMer(models_nocomp$width_mod, FUN = calc_icc_conditional, nsim = sims)$t),
       height = as.numeric(bootMer(models_nocomp$height_mod, FUN = calc_icc_conditional, nsim = sims)$t),
       density = as.numeric(bootMer(models_nocomp$density_mod, FUN = calc_icc, nsim = sims)$t),
       beta = as.numeric(bootMer(models_nocomp$beta_mod, FUN = calc_icc_conditional, nsim = sims)$t)) -> boot_samples_conditional

# Bind together two data sets and convert to long format
rbind(boot_samples, boot_samples_conditional) %>% 
  mutate(icc_type = rep(c("adjusted", "conditional"), each = nrow(boot_samples))) %>% 
  gather(key = trait, value = icc, agb:beta) %>% 
  filter(complete.cases(icc)) -> boot_samples_all

# Calculate mean iccs across all traits
mean_iccC <- boot_samples_all %>% filter(icc_type == "conditional") %>% pull(icc) %>% mean()
mean_iccM <- boot_samples_all %>% filter(icc_type == "adjusted") %>% pull(icc) %>% mean()

boot_samples_all %>% 
  ggplot(aes(x = reorder(trait, -icc, mean), y = icc)) +
  geom_violin(aes(fill = icc_type), position = position_dodge(width = 0.7), alpha = 0.4, color = NA) +
  stat_summary(fun = mean, fun.min = function(z) { quantile(z,0.025) },
               fun.max = function(z) { quantile(z,0.975) }, aes(colour = icc_type),
               position = position_dodge(width = 0.7)) +
  coord_flip() +
  geom_hline(aes(yintercept = mean_iccC), linetype = "dashed", color = "#1f78b4") +
  geom_hline(aes(yintercept = mean_iccM), linetype = "dashed", color = "#a6cee3") +
  ylab("intraclass correlation coefficient (ICC)") +
  xlab("trait") +
  labs(fill = "ICC type",
       color = "ICC type") +
  scale_color_manual(values = c("#a6cee3", "#1f78b4")) +
  scale_fill_manual(values = c("#a6cee3", "#1f78b4")) +
  theme_bw() +
  theme(legend.position = "bottom")-> random_intercept


## Aboveground Biomass G x E ####

# Refit with summed contrasts for plotting (this doesn't change the model)
expand.grid(salinity = unique(traits_nocomp$salinity),
            date_cloned_grp = unique(traits_nocomp$date_cloned_grp),
            genotype = unique(row.names(ranef(agb_mod1)$genotype)),
            elevation = seq(0.156, 0.544, length.out = 20),
            co2 = unique(traits_nocomp$co2)) %>% 
  mutate(age = case_when(substr(genotype, 2, 2) == "a" ~ "ancestral",
                         T ~ "modern"),
         weight_init = mean(traits_nocomp$weight_init)) -> newData

newData %>% filter(genotype != "ca6" | co2 != "elevated") -> newData_ca6

predict(agb_mod1, newData_ca6) -> predictions

my_labels <- c(ambient="ambient~CO[2]", elevated="elevated~CO[2]")
my_labeller <- as_labeller(my_labels,
                           default = label_parsed)

# Create color palette such that the genotype lines are colored by the order of
# the intercept in ambient CO2

newData_ca6 %>% 
  mutate(predict_agb = predictions^2) %>% 
  filter(elevation == min(newData_ca6$elevation) & co2 == "ambient") %>% 
  group_by(co2, genotype, age) %>% 
  summarize(mean_agb = mean(predict_agb)) %>% 
  arrange(mean_agb) %>% 
  ungroup() %>% 
  mutate(color_order = 1:31) %>% 
  select(genotype, color_order) -> color_labels_agb

left_join(newData_ca6, color_labels_agb) -> plot_data_agb

plot_data_agb %>% 
  mutate(predict_agb = predictions^2) %>% 
  group_by(elevation, co2, genotype, age, color_order) %>% 
  summarize(mean_agb = mean(predict_agb)) %>% 
  ggplot(aes(x = elevation, y = mean_agb, group = genotype, color = color_order)) +
  geom_line() +
  facet_wrap(~co2, labeller = my_labeller) +
  ylab("aboveground biomass (g)") +
  xlab("") +
  theme_bw() +
  theme(axis.text.x = element_blank()) +
  scale_color_gradientn(colours = rainbow(31)) +
  theme(legend.position = "none")-> co2_GxE

## Root distribution parameter G x E ####
# Refit with summed contrasts for plotting (this doesn't change the model)

beta_mod1 <- lmer(beta ~ location * elevation + elevation * salinity + (1+elevation|genotype:salinity) + (1|genotype),
                 data = traits_nocomp_rs)
expand.grid(salinity = unique(traits_nocomp$salinity),
            genotype = unique(row.names(ranef(beta_mod1)$genotype)),
            elevation = seq(0.156, 0.544, length.out = 20)) %>% 
  mutate(location = case_when(substr(genotype, 1, 1) == "c" ~ "corn",
                         T ~ "kirkpatrick")) -> newData

newData %>% filter(genotype != "ca8" | salinity != "salt") -> newData_ca8

predict(beta_mod1, newData_ca8) -> predictions_beta

newData_ca8 %>% 
  mutate(predict_beta = predictions_beta) %>% 
  mutate(salinity = case_when(salinity == "fresh" ~ "freshwater site (4 ppt)",
                              T ~ "brackish site (6 ppt)")) %>%
  filter(elevation == min(newData_ca8$elevation) & salinity == "brackish site (6 ppt)") %>% 
  arrange(predict_beta) %>% 
  ungroup() %>% 
  mutate(color_order = 1:30) %>% 
  select(genotype, color_order) -> color_labels_beta

newData_ca8 %>% 
  mutate(predict_beta = predictions_beta) %>% 
  mutate(salinity = case_when(salinity == "fresh" ~ "freshwater site (4 ppt)",
                              T ~ "brackish site (6 ppt)")) -> newData_ca8

left_join(newData_ca8, color_labels_beta) -> plot_data_beta
# Add last color for ca8 because it is currently NA
plot_data_beta %>% 
  mutate(color_order = ifelse(is.na(color_order), 31, color_order)) -> plot_data_beta

plot_data_beta %>% 
  ggplot(aes(x = elevation, y = predict_beta, group = genotype, color = color_order)) +
  geom_line() +
  facet_wrap(~salinity) +
  ylab("root distribution parameter") +
  xlab("elevation (m NAVD88)") +
  theme_bw() +
  scale_color_gradientn(colours = rainbow(31)) +
  theme(legend.position = "none") -> salinity_GxE

## Bring plots together ####
png("figs/FigX_GxE.png", height = 5.5, width = 9, res = 300, units = "in")
random_intercept + co2_GxE/salinity_GxE + plot_annotation(tag_levels = "a")
dev.off()
