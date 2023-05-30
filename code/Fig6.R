# Figure 6: genotype and GxE plot

## Read in trait data ####
source("supp_code/CompileTraitData.R")

bg_full %>% 
  filter(comp == 0 & location != "blackwater" & level < 5 &
           agb_scam > 0) -> traits_nocomp

# Center and scale elevation values
traits_nocomp %>% 
  mutate(elevation_sc = scale(elevation)[,1],
         elevation_sc2 = scale(elevation^2)[,1]) -> traits_nocomp

# Drop outlier for height data
traits_nocomp %>% 
  filter(pot_no != 1830) -> traits_nocomp_height

traits_nocomp %>% 
  mutate(rs = total_bg / agb_scam) -> traits_nocomp

traits_nocomp %>% 

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

## Make ICC plot ####
# Collect all final no competition models
agb_evo_model <- readRDS("derived_data/agb_model.rda")
bg_evo_model <- readRDS("derived_data/bg_model.rda")
rs_evo_model <- readRDS("derived_data/rs_model.rda")
width_evo_model <- readRDS("derived_data/width_model.rda")
height_evo_model <- readRDS("derived_data/height_model.rda")
density_evo_model <- readRDS("derived_data/density_model.rda")
beta_evo_model <- readRDS("derived_data/beta_model.rda")


models_nocomp <- list(agb_mod = agb_evo_model,
                      bg_mod = bg_evo_model,
                      rs_mod = rs_evo_model,
                      width_mod = width_evo_model,
                      height_mod = height_evo_model,
                      density_mod = density_evo_model,
                      beta_mod = beta_evo_model)

# Create functional to calculate and extract marginal and conditional ICC values
calc_icc <- function(x){
  as.numeric(performance::icc(x)[1])
}

calc_icc_conditional <- function(x){
  as.numeric(performance::icc(x)[2])
}

# Set number of bootstrap simulations
sims = 1000

# Calculate ICC over bootstrapped simulations
tibble(agb = as.numeric(bootMer(models_nocomp$agb_mod, FUN = calc_icc, nsim = sims)$t),
       bgb = as.numeric(bootMer(models_nocomp$bg_mod, FUN = calc_icc, nsim = sims)$t),
       rs = as.numeric(bootMer(models_nocomp$rs_mod, FUN = calc_icc, nsim = sims)$t),
       width = as.numeric(bootMer(models_nocomp$width_mod, FUN = calc_icc, nsim = sims)$t),
       height = as.numeric(bootMer(models_nocomp$height_mod, FUN = calc_icc, nsim = sims)$t),
       density = as.numeric(bootMer(models_nocomp$density_mod, FUN = calc_icc, nsim = sims)$t),
       beta = as.numeric(bootMer(models_nocomp$beta_mod, FUN = calc_icc, nsim = sims)$t)) -> boot_samples

# Calculate conditional ICC over bootstrapped simulations
tibble(agb = as.numeric(bootMer(models_nocomp$agb_mod, FUN = calc_icc_conditional, nsim = sims)$t),
       bgb = as.numeric(bootMer(models_nocomp$bg_mod, FUN = calc_icc_conditional, nsim = sims)$t),
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

# Calculate average conditional ICC for each trait
boot_samples_all %>%
  filter(icc_type == "conditional") %>%
  group_by(trait) %>% summarize(m = mean(icc))

# Create plot of conditional ICC for each trait
boot_samples_all %>% 
  filter(icc_type == "conditional") %>% 
  ggplot(aes(x = reorder(trait, icc, mean), y = icc)) +
  stat_summary(fun = mean, fun.min = function(z) { quantile(z,0.025) },
               fun.max = function(z) { quantile(z,0.975) },
               position = position_dodge(width = 0.7), color = "black", size = 0.8) +
  ylab("intraclass correlation coefficient (ICC)") +
  xlab("trait") +
  theme_bw(base_size = 14) +
  theme(legend.position = "top") +
  ylim(0, 0.5) +
  scale_x_discrete(labels = c("r:s", "agb", "bgb", "width","\u03B2", "density", "height"))-> random_intercept

## GxE graphs ####

# Height graph

# Predict height for all possibe combinations of predictor variables
expand.grid(salinity = unique(traits_nocomp_height$salinity),
            genotype = unique(row.names(ranef(height_evo_model)$genotype)),
            elevation_sc = seq(min(traits_nocomp_height$elevation_sc),
                            max(traits_nocomp_height$elevation_sc), length.out = 20),
            weight_init = mean(traits_nocomp_height$weight_init)) %>% 
  mutate(location = case_when(substr(genotype, 1, 1) == "c" ~ "corn",
                              T ~ "kirkpatrick"),
         elevation = rep(seq(min(traits_nocomp_height$elevation),
                         max(traits_nocomp_height$elevation), length.out = 20), each = 62)) %>% 
  mutate(flooding = mean_tide - elevation)-> newData_height

predict(height_evo_model, newData_height) -> predictions_height_RE

# Group predictions by elevation, salinity, and genotype for graphing
newData_height %>% 
  mutate(predict_height = predictions_height_RE) %>% 
  mutate(salinity = case_when(salinity == "salt" ~ "brackish site",
                              T ~ "freshwater site")) %>% 
  group_by(elevation, salinity, genotype) %>% 
  summarize(mean_height = mean(predict_height)) -> height_dat

# Create height plot
height_dat %>% 
  spread(salinity, mean_height) %>% 
  mutate(diff = `freshwater site` - `brackish site`,
         flooding = mean_tide - elevation) %>% 
  ggplot(aes(flooding, diff, group = genotype, color = genotype)) +
  geom_line() + 
  geom_hline(aes(yintercept = 0), color = "black", size = 1, linetype = "dashed") +
  ylab(expression(paste(height[fresh.]-height[brack.]))) +
  xlab("average inundation (m)")+
  theme_bw(base_size = 14) +
  theme(legend.position = "none") -> height_GxE

# Beta graph - elevation*salinity
expand.grid(salinity = unique(traits_nocomp_rs$salinity),
            genotype = unique(row.names(ranef(beta_evo_model)$genotype)),
            date_cloned_grp = unique(traits_nocomp_rs$date_cloned_grp),
            elevation_sc = seq(min(traits_nocomp_height$elevation_sc),
                               max(traits_nocomp_height$elevation_sc), length.out = 20),
            weight_init = mean(traits_nocomp_rs$weight_init)) %>% 
  mutate(location = case_when(substr(genotype, 1, 1) == "c" ~ "corn",
                              T ~ "kirkpatrick"),
         age = case_when(substr(genotype, 2, 2) == "a" ~ "ancestral",
                                T ~ "modern"),
         elevation = rep(seq(min(traits_nocomp_height$elevation),
                             max(traits_nocomp_height$elevation), length.out = 20), each = 186),
         flooding = mean_tide - elevation)-> newData_beta

predict(beta_evo_model, newData_beta) -> predictions_beta_RE

# Group predictions by elevation, salinity, and genotype for graphing
newData_beta %>% 
  mutate(predict_beta = predictions_beta_RE) %>% 
  mutate(salinity = case_when(salinity == "salt" ~ "brackish site",
                              T ~ "freshwater site")) %>% 
  group_by(elevation, salinity, genotype) %>% 
  summarize(mean_beta = mean(predict_beta)) -> beta_dat

beta_dat %>% 
  spread(salinity, mean_beta) %>% 
  mutate(diff = `freshwater site` - `brackish site`,
         flooding = mean_tide - elevation) %>% 
  ggplot(aes(flooding, diff, group = genotype, color = genotype)) +
  geom_line() + 
  geom_hline(aes(yintercept = 0), color = "black", size = 1, linetype = "dashed") +
  ylab(expression(paste(beta[fresh.]-beta[brack.]))) +
  xlab("average inundation (m)")+
  theme_bw(base_size = 14) +
  theme(legend.position = "none") -> beta_GxE

# Aboveground biomass graph
expand.grid(salinity = unique(traits_nocomp_rs$salinity),
            co2 = unique(traits_nocomp_rs$co2),
            genotype = unique(row.names(ranef(rs_evo_model)$genotype)),
            date_cloned_grp = unique(traits_nocomp$date_cloned_grp),
            elevation_sc = mean(traits_nocomp_rs$elevation_sc),
            weight_init = mean(traits_nocomp$weight_init)) %>% 
  mutate(age = case_when(substr(genotype, 2, 2) == "a" ~ "ancestral",
                         T ~ "modern")) -> newData_agb

predict(agb_evo_model, newData_agb) -> predictions_agb_RE

# Group predictions by genotype and co2
newData_agb %>% 
  mutate(predict_agb = predictions_agb_RE^2) %>% 
  group_by(co2, genotype) %>% 
  summarize(mean_agb = mean(predict_agb)) -> color_labels_agb

left_join(newData_agb, color_labels_agb) -> plot_data_agb

plot_data_agb %>% 
  spread(co2, mean_agb) %>% 
  mutate(diff = elevated - ambient) %>% 
  select(genotype, diff) %>% 
  unique() %>% 
  arrange(diff) %>% 
  ggplot(aes(reorder(genotype, diff), diff)) +
  geom_point(aes(color = genotype), size = 3, alpha = 0.5) +
  geom_point(shape = 1, size = 3) +
  theme_bw(base_size = 14) +
  geom_hline(aes(yintercept = 0), color = "black", size = 1, linetype = "dashed") +
  theme(axis.text.x = element_blank(),
        legend.position = "none") +
  ylab(expression(paste(agb[elev.]-agb[amb.]))) +
  xlab("genotype") -> agb_GxE

# Root-to-shoot graph
expand.grid(salinity = unique(traits_nocomp_rs$salinity),
            genotype = unique(row.names(ranef(rs_evo_model)$genotype)),
            elevation_sc = mean(traits_nocomp_rs$elevation_sc),
            co2 = unique(traits_nocomp_rs$co2)) %>% 
  mutate(age = case_when(substr(genotype, 2, 2) == "a" ~ "ancestral",
                         T ~ "modern")) -> newData_rs

predict(rs_evo_model, newData_rs) -> predictions_rs_RE

# Group predictions by salinity and genotype
newData_rs %>% 
  mutate(predict_rs = exp(predictions_rs_RE)) %>% 
  group_by(salinity, genotype) %>% 
  summarize(mean_rs = mean(predict_rs)) -> color_labels_rs

left_join(newData_rs, color_labels_rs) -> plot_data_rs

plot_data_rs %>% 
  spread(salinity, mean_rs) %>% 
  mutate(diff = fresh - salt) %>% 
  select(genotype, diff) %>% 
  unique() %>% 
  arrange(diff) %>% 
  ggplot(aes(reorder(genotype, diff), diff)) +
  geom_point(aes(color = genotype), size = 3, alpha = 0.5) +
  geom_point(shape = 1, size = 3) +
  theme_bw(base_size = 14) +
  geom_hline(aes(yintercept = 0), color = "black", size = 1, linetype = "dashed") +
  theme(axis.text.x = element_blank(),
        legend.position = "none") +
  ylab(expression(paste(r:s[fresh.]-r:s[brack.]))) +
  xlab("genotype") -> rs_GxE

# Belowground biomass graph
expand.grid(salinity = unique(traits_nocomp_rs$salinity),
            genotype = unique(row.names(ranef(bg_evo_model)$genotype)),
            elevation_sc = mean(traits_nocomp_rs$elevation_sc),
            co2 = unique(traits_nocomp_rs$co2),
            weight_init = mean(traits_nocomp$weight_init),
            date_cloned_grp = unique(traits_nocomp_rs$date_cloned_grp)) -> newData_bg

predict(bg_evo_model, newData_bg) -> predictions_bg_RE

# Group new predictions by co2 and genotype
newData_bg %>% 
  mutate(predict_bg = predictions_bg_RE^2) %>% 
  group_by(co2, genotype) %>% 
  summarize(mean_bg = mean(predict_bg)) -> color_labels_bg

left_join(newData_bg, color_labels_bg) -> plot_data_bg

plot_data_bg %>% 
  spread(co2, mean_bg) %>% 
  mutate(diff = elevated - ambient) %>% 
  select(genotype, diff) %>% 
  unique() %>% 
  arrange(diff) %>% 
  ggplot(aes(reorder(genotype, diff), diff)) +
  geom_point(aes(color = genotype), size = 3, alpha = 0.5) +
  geom_point(shape = 1, size = 3) +
  theme_bw(base_size = 14) +
  geom_hline(aes(yintercept = 0), color = "black", size = 1, linetype = "dashed") +
  theme(axis.text.x = element_blank(),
        legend.position = "none") +
  ylab(expression(paste(bgb[elev.]-bgb[amb.]))) +
  xlab("genotype") -> bg_GxE

## Bring plots together ####

# Create custom design
design <- c(
  area(1, 1, 2, 2),
  area(1, 3, 1, 3),
  area(2, 3, 2, 3),
  area(3, 1, 3, 1),
  area(3, 2, 3, 2),
  area(3, 3, 3, 3)
)

# Preview design of plot
# plot(design)

png("figs/Fig6.png", height = 7, width = 9.5, units = "in", res = 300)

random_intercept+beta_GxE+height_GxE+
  agb_GxE+bg_GxE+rs_GxE+
  plot_layout(design = design) +
  plot_annotation(tag_levels = "a")

dev.off()
