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
rs_ES_plot <- summary(emmeans::emmeans(rs_evo_model, ~elevation:salinity, at = list(elevation = seq(0.156, 0.544, length.out = 50))))

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

bg_ES_plot <- summary(emmeans::emmeans(bg_evo_model, ~elevation:salinity, at = list(elevation = seq(0.156, 0.544, length.out = 50))))

tibble(total_bg = bg_ES_plot$emmean^2,
       elevation = bg_ES_plot$elevation,
       lower = bg_ES_plot$lower.CL^2,
       upper = bg_ES_plot$upper.CL^2,
       salinity = rep(c("fresh", "salt"), each = 50)) %>% 
  ggplot(aes(x = elevation, y = total_bg, color = salinity)) +
  geom_point(data = traits_nocomp_rs, aes(x = elevation, y = total_bg, color = salinity), shape = 1,
             stroke = 0.8, alpha = 0.6, size = 0.8) +
  geom_textline(size = 3, fontface = "bold", label= rep(c("freshwater site", "brackish site"), each = 50),
                linewidth = 1.2, hjust = 0.90, textcolor = "black") +
  theme_bw() +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = salinity), alpha = 0.2, color = NA) +
  ylab("belowground biomass") +
  xlab("elevation (m NAVD88)") + 
  scale_color_manual(values = c("#fdbf6f", "#ff7f00")) +
  scale_fill_manual(values = c("#fdbf6f", "#ff7f00")) +
  theme(legend.position = "none") -> bg_interaction

FigS4 <- rs_interaction / bg_interaction + plot_annotation(tag_levels = "a") & theme(plot.margin = margin(t = 0,r = 5,b = 0,l = 5))

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

## Old Stuff ####
## Root-to-shoot ratio G x E 

expand.grid(salinity = unique(traits_nocomp_rs$salinity),
            genotype = unique(row.names(ranef(rs_evo_model)$genotype)),
            elevation = seq(0.156, 0.544, length.out = 20),
            co2 = unique(traits_nocomp_rs$co2)) %>% 
  mutate(age = case_when(substr(genotype, 2, 2) == "a" ~ "ancestral",
                         T ~ "modern")) -> newData_rs

predict(rs_evo_model, newData_rs) -> predictions_rs_RE

newData_rs %>% 
  mutate(predict_rs = exp(predictions_rs_RE)) %>% 
  group_by(salinity, genotype) %>% 
  summarize(mean_rs = mean(predict_rs)) -> color_labels_rs

left_join(newData_rs, color_labels_rs) -> plot_data_rs

plot_data_rs %>% 
  mutate(salinity = case_when(salinity == "fresh" ~ "freshwater site",
                              T ~ "brackish site")) %>% 
  ggplot(aes(x = salinity, y = mean_rs, group = genotype, color = genotype)) +
  geom_line() +
  geom_point(size = 1.5) +
  ylab("root-to-shoot ratio") +
  xlab("salinity") +
  theme_bw() +
  scale_color_manual(values = rainbow(31)) +
  theme(legend.position = "none")

## Belowground biomass G x E 

expand.grid(salinity = unique(traits_nocomp_rs$salinity),
            genotype = unique(row.names(ranef(bg_evo_model)$genotype)),
            elevation = seq(0.156, 0.544, length.out = 20),
            co2 = unique(traits_nocomp_rs$co2),
            weight_init = mean(traits_nocomp$weight_init),
            date_cloned_grp = unique(traits_nocomp_rs$date_cloned_grp)) -> newData_bg

predict(bg_evo_model, newData_bg) -> predictions_bg_RE

newData_bg %>% 
  mutate(predict_bg = predictions_bg_RE^2) %>% 
  group_by(co2, genotype) %>% 
  summarize(mean_bg = mean(predict_bg)) -> color_labels_bg

left_join(newData_bg, color_labels_bg) -> plot_data_bg

plot_data_bg %>% 
  ggplot(aes(x = co2, y = mean_bg, group = genotype, color = genotype)) +
  geom_line() +
  geom_point(size = 1.5) +
  ylab("belowground biomass") +
  xlab(expression(paste(CO[2], " treatment"))) +
  theme_bw() +
  scale_color_manual(values = rainbow(31)) +
  theme(legend.position = "none")

traits_nocomp_rs %>% 
  group_by(co2, genotype) %>% 
  summarize(mean_bg = mean(total_bg)) %>% 
  ggplot(aes(x = co2, y = mean_bg)) +
  geom_point()





plot_data_bg %>% 
  spread(co2, mean_bg) %>% 
  mutate(diff = elevated - ambient) %>% 
  select(genotype, diff) %>% 
  unique() %>% 
  arrange(diff) %>% 
  ggplot(aes(diff)) +
  geom_density() + 
  geom_rug()


## Aboveground Biomass G x E 

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

## Root distribution parameter G x E 
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


plot_data_height %>% 
  mutate(salinity = case_when(salinity == "salt" ~ "brackish site",
                              T ~ "freshwater site")) %>% 
  mutate(predict_height = predictions_height_RE) %>% 
  group_by(elevation, salinity, genotype, color_order) %>% 
  summarize(mean_height = mean(predict_height)) %>% 
  ggplot(aes(x = elevation, y = mean_height, group = genotype, color = color_order)) +
  geom_line() +
  facet_wrap(~salinity) +
  ylab("mean stem height (cm)") +
  xlab("elevation (cm NAVD88)") +
  theme_bw() +
  scale_color_gradientn(colours = rainbow(31)) +
  theme(legend.position = "none")

