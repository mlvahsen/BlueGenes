# Compile trait data for Blue Genes

# Blue Genes Analysis -- Clean code!!! #

####################
## Data Wrangling ##
####################

## Preliminaries ####

# Load libraries
library(tidyverse); library(lubridate); library(here)

# Set graphical theme
theme_set(theme_bw())

# Read in elevation data
source(here("supp_code", "Calculate_Pot_Elevation.R"))

# Read in initial conditions data
inits <- read_csv(here("supp_data","InitialConditions_FinalReps_NOEDITS.csv"))

# Read in the rest of the data
census <- read_csv(here("supp_data","Biweekly_toSept16.csv"))
abg <- read_csv(here("supp_data","AbovegroundBiomass_data_17Jun20.csv"))
reap_pot <- read_csv(here("supp_data","Harvest_FullData_PotLevel.csv"))
reap_stem <- read_csv(here("supp_data","Harvest_FullData_StemLevel.csv"))
pot_ids <- read_csv(here("supp_data","HarvestPot_ID.csv")) # for matching up pot_no with data sheets that don't have pot_no
rna_harvest <- read.csv(here("supp_data","RNA_harvest.csv"))

## Data cleaning and formatting ##### 

# Add pot numbers to initial conditions data frame
inits %>% 
  rename(level = flood_level) %>% 
  mutate(id = paste(site_frame, level, position, sep = "_")) -> inits

elevation_data$site_level <- paste(elevation_data$site_frame, elevation_data$level, sep = "_")
inits$site_level <- paste(inits$site_frame, inits$level, sep = "_")

# Add true elevation to initial conditions data
merge(elevation_data, inits, by = "site_level") -> inits_elev

# Add pot numbers
inits_elev %>% 
  mutate(id = paste(site_level, position, sep = "_")) %>% 
  dplyr::select(-site_frame.y, -level.y) %>% 
  rename(site_frame = site_frame.x,
         level = level.x) -> inits_elev

pot_ids %>% 
  mutate(id = paste(`site_ frame`, level, pos, sep = "_")) %>% 
  dplyr::select(id, pot_no) -> pot_ids_only

merge(inits_elev, pot_ids_only, by = "id") -> inits_ids

# Add aboveground biomass

# Get total weight of biomass that was weighed in bags
abg %>% 
  dplyr::select(-Date) %>% # Helps deal with the double entry problem
  dplyr::select(-Recorder) %>% 
  #filter(notes != "dh rhizome" & notes != "dh stem" | is.na(notes)) %>% 
  filter(notes != "dh rhizome" & notes != "dh stem" | is.na(notes)) %>% 
  filter(!grepl("dh rhizome",notes)) %>% 
  dplyr::select(-notes) %>% 
  mutate(allom = ifelse(allom == "N/A", NA, allom)) %>% 
  # Deal with the double entry problem
  distinct() %>% 
  group_by(pot_no, species) %>% 
  summarize(weight = sum(weight)) %>% 
  spread(key = species, value = weight) %>% 
  mutate(scam = ifelse(is.na(scam), 0, scam),
         sppa = ifelse(is.na(sppa), 0, sppa)) -> abg_no_zeros

# Figure out which pots have zero biomass (i.e. no biomass was weighed)
pot_ids %>% 
  filter(!pot_no %in% abg_no_zeros$pot_no) %>% 
  mutate(scam = 0,
         sppa = 0) %>% 
  dplyr::select(pot_no, scam, sppa) -> abg_zeros

# Bind together
abg_all <- rbind(abg_no_zeros, abg_zeros)

# Merge abg with initial conditions data frame
full_data <- merge(inits_ids, abg_all, by = "pot_no")

# Match up initial conditions with harvest census
reap_pot %>% 
  dplyr::select(pot_no, scam_live_density, scam_dead_density, max_sppa_height,
         soil_to_pot, dh_count) %>% 
  mutate(max_sppa_height = ifelse(max_sppa_height == "N/A", NA, max_sppa_height)) -> reap_sum

# Merge reap_sum with full_data
full_data <- merge(full_data, reap_sum, by = "pot_no")

# Now deal with pot level harvest data
# There are some non-numerics in the width data
reap_stem[is.na(as.numeric(reap_stem$mid_wd)),"mid_wd"]

reap_stem %>% 
  mutate(mid_wd = as.numeric(ifelse(mid_wd == "N/A", NA, mid_wd)),
         bas_wd = as.numeric(ifelse(bas_wd == "N/A", NA, bas_wd))) %>% 
  # Filter out drainage hole stems and those that are cut when calculating
  # average height and widths
  filter(dh == 0 & cut == 0) %>% 
  group_by(pot_no) %>% 
  summarize(mean_tot_ht = mean(tot_ht, na.rm = T),
            mean_grn_ht = mean(grn_ht, na.rm = T),
            mean_mid_wd = mean(mid_wd, na.rm = T),
            mean_bas_wd = mean(bas_wd, na.rm = T),
            # Also calculate proportion twisted
            prop_twist = length(which(twist == 1)) / length(twist)) -> reap_height_sum

# Create rows for data that does not have measurements for scam
pot_ids %>% 
  filter(!pot_no %in% reap_height_sum$pot_no) %>% 
  dplyr::select(pot_no) %>% 
  mutate(mean_tot_ht = NA,
         mean_grn_ht = NA,
         mean_mid_wd = NA,
         mean_bas_wd = NA,
         prop_twist = NA) -> reap_height_zeros

# Bind together data with heights and data with no heights
reap_all <- rbind(reap_height_sum, reap_height_zeros)

# Merge with full data
full_data <- merge(full_data, reap_all, by = "pot_no")

# Check to see if we have any instances where we don't have height/width data
# but do have biomass data
full_data %>% filter(is.na(mean_tot_ht)) %>% 
  dplyr::select(scam_live_density, pot_no) %>% 
  filter(scam_live_density > 0)
# Harvest date on this is 13 October so this must be the pot that MLV pulled
# late that actually had biomass in it

# Check to see where there is biomass but there were no live stems
full_data[which(full_data$scam != 0 & full_data$scam_live_density == 0),] %>% 
  dplyr::select(pot_no, scam, scam_live_density, mean_tot_ht)
# These are all just dead biomass so this makes sense

# Check to see where there was no biomass but there were live stems
full_data[which(full_data$scam == 0 & full_data$scam_live_density != 0),] %>% 
  dplyr::select(pot_no, scam, scam_live_density, mean_tot_ht)
# Good!

# Group date_cloned and date_planted data information so it can be used in the
# models in a meaningful way. Clean up data and dplyr::select what we want for models.
full_data %>% 
  mutate(date_cloned_grp = factor(ifelse(date_cloned %in% c("5/22/19", "5/28/19", "5/27/19"), 1,
                                         ifelse(date_cloned %in% c("6/9/19", "6/10/19"), 2, 3))),
         date_planted_grp = factor(ifelse(date_planted %in% c("6/13/19", "6/14/19", "6/17/19", "6/12/19"), 1, 2))) %>% 
  dplyr::select(pot_no,site_frame,level,position,
         co2,
         salinity,
         comp,
         elevation,
         genotype,
         location,
         depth_seed = depth,
         cohort,
         origin_lab = origin,
         date_cloned_grp,
         date_planted_grp,
         weight_init = weight,
         width_init = width,
         nodes_init = nodes,
         agb_scam = scam,
         agb_sppa = sppa,
         dens_scam_live = scam_live_density,
         dens_scam_dead = scam_dead_density,
         height_sppa = max_sppa_height,
         height_scam_tot = mean_tot_ht,
         height_scam_grn = mean_grn_ht,
         width_scam_mid = mean_mid_wd,
         width_scam_bas = mean_bas_wd,
         prop_twist,
         soil_to_pot
  ) -> full_data

## Dealing with RNA stems #### 
# One stem was harvested for all pots in the rna_harvest data frame. We add one
# stem to the density count and then use the allometric data to calculate the
# biomass and then add that to the measured biomass for the rest of the stems.

# Use allometric stems to create a calibration to predict stems that were
# harvested for RNA
abg %>% 
  filter(allom %in% c(1,2)) %>% 
  mutate(stem_id = paste(pot_no, allom, sep = "_")) %>% 
  dplyr::select(stem_id, weight) %>% 
  # deals with double entry of some biomass data
  distinct() -> allom_biomass

# Get height and width from census data for allom stems
# Fix one reap_stem data entry issue
reap_stem %>% 
  filter(allom %in% c(1,2)) %>% 
  mutate(stem_id = paste(pot_no, allom, sep = "_")) %>% 
  dplyr::select(stem_id, tot_ht, grn_ht, mid_wd, bas_wd) -> allom_hw

dim(allom_biomass)
dim(allom_hw) # doesn't quite match up so filter out ids that match up

merge(allom_biomass, allom_hw, by.x = "stem_id") -> allom_data

# Fit log-transformed linear model accounting for heights, widths and interaction
# between height and width
hw_mod <- lm(log(weight) ~ tot_ht * mid_wd, data = allom_data)
#plot(hw_mod)
# Fit log-transformed linear model accounting for heights only
h_mod <- lm(log(weight) ~ tot_ht, data = allom_data)
#plot(h_mod)

# Predict weights of RNA stems
# Filter data that has heights and widths (all but 2)
rna_harvest %>% 
  filter(complete.cases(mid_wd)) -> rna_hw
# Filter data that only has heights
rna_harvest %>% 
  filter(!complete.cases(mid_wd)) -> rna_h

# Predict mass for stems we have height and widths
rna_hw$biomass <- exp(predict(hw_mod, newdata = rna_hw %>% dplyr::select(tot_ht, mid_wd)))
rna_h$biomass <- exp(predict(h_mod, newdata = rna_h %>% dplyr::select(tot_ht)))
rna_pred <- rbind(rna_hw, rna_h) 

# Sum biomass and dplyr::select pot_id and additional biomass column
rna_pred %>% 
  group_by(pot_no) %>% 
  summarize(biomass_add = sum(biomass)) -> rna_pred

# Add additional biomass
full_data[full_data$pot_no %in% rna_pred$pot_no,"agb_scam"] <- full_data[full_data$pot_no %in% rna_pred$pot_no,"agb_scam"] + rna_pred$biomass_add

# Add in extra density if RNA harvested stem (+1 stem)
full_data %>% 
  mutate(dens_scam_live = ifelse(pot_no %in% rna_harvest$pot_no, dens_scam_live + 1, dens_scam_live)) -> full_data

## Belowground ####
# Read in belowground data
source(here("supp_code", "bgb_data_fix.R"))

# Group data by pot to get total belowground
bg_clean %>% 
  group_by(pot_no) %>% 
  summarize(total_bg = sum(weight)) -> bg_total

full_data_merged <- merge(bg_total, full_data, by = "pot_no", all.y = T)

# Also calculate root distribution parameters values for non-competition pots
# Get pot numbers for no competition pots
full_data %>% 
  filter(comp == 0) %>% 
  pull(pot_no)-> no_comp_pots

bg_clean %>% 
  filter(pot_no %in% no_comp_pots) %>% 
  # Group everything by layer
  group_by(pot_no, segment_top) %>% 
  summarize(biomass = sum(weight))-> bg_nocomp

# Make summary of total biomass and number of layers for each pot
bg_nocomp %>%
  group_by(pot_no) %>%  
  summarize(total = sum(biomass),
            n = length(biomass)) %>% 
  # There are a couple pots that had sectioning issues (with just 1-2 sections)
  # so we drop these
  filter(n > 2) -> totals

# Figure out which pots we dropped
bg_nocomp %>%
  group_by(pot_no) %>%  
  summarize(total = sum(biomass),
            n = length(biomass)) %>% 
  filter(n < 3) %>% pull(pot_no) -> dropped_pots

# Make %notin%
`%notin%` <- Negate(`%in%`)

# Make sure data is sorted by pot and rooting depth
bg_nocomp %>%
  ungroup() %>% 
  filter(pot_no %notin% dropped_pots) %>% 
  arrange(pot_no, segment_top) %>% 
  # Create repeated entries of the total belowground biomass in order to
  # calculate proportions
  mutate(total = rep(totals$total, totals$n)) %>% 
  # Calculate proportions
  mutate(prop = biomass / total) %>% 
  # Group by pot so that we can calculate cumulative proportion across depth
  group_by(pot_no) %>% 
  mutate(cum_sum = cumsum(prop)) %>% 
  ungroup() -> bgb_beta_analysis

# Create an additional set of data such that at depth 0, there is no biomass
bgb_beta_analysis %>% 
  dplyr::select(pot_no) %>% 
  unique() -> for_adding_zeros 

# Here depth is -5, but will get changed to 0 in the next step
for_adding_zeros$segment_top <- -10
for_adding_zeros$biomass <- 0
for_adding_zeros$cum_sum <- 0

# Order columns the same way as the for_adding_zeros dataframe
bgb_beta_analysis %>% 
  dplyr::select(pot_no, segment_top, biomass, cum_sum)->bgb_beta_analysis

bgb_beta_analysis_zeros <- rbind(for_adding_zeros, bgb_beta_analysis) %>% 
  mutate(depth_roots = case_when(segment_top == -10 ~ 0,
                                 segment_top == 0 ~ 10,
                                 segment_top == 10 ~ 20,
                                 T ~ 70)) %>% 
  arrange(pot_no, depth_roots) %>% 
  dplyr::select(pot_no, depth_roots, biomass, cum_sum)

## Run non-linear regression on cumulative probability for each pot

# Create output vector to hold information
beta <- NULL

for (i in 1:length(unique(for_adding_zeros$pot_no))){
  sub_data <- bgb_beta_analysis_zeros %>% filter(pot_no == unique(for_adding_zeros$pot_no)[i])
  mod_temp <- nls(cum_sum ~ 1 - beta ^ (depth_roots), data = sub_data, start = list(beta = 1))
  beta[i] <- coef(mod_temp)[1]
}

for_adding_zeros$beta <- beta

# Merge together with the rest of the data
for_adding_zeros %>% 
  dplyr::select(pot_no, beta) %>% 
  # Add NAs for dropped pots
  add_row(pot_no = dropped_pots[1], beta = NA) %>%
  add_row(pot_no = dropped_pots[2], beta = NA) %>%
  add_row(pot_no = dropped_pots[3], beta = NA) %>%
  arrange(pot_no) -> betas_to_merge

bg_full <- merge(betas_to_merge, full_data_merged, all.y = T, by = "pot_no")

# Create an extinct vs extant binary variable based on SCAM aboveground biomass
bg_full %>% 
  mutate(extinct = case_when(agb_scam > 0 ~ 0,
                             T ~ 1)) %>% 
  # Also make a variable for age cohort
  mutate(age = case_when(grepl("ancestral", cohort) ~ "ancestral",
                         T ~ "modern"))-> bg_full

# Fix spartina competition pots -- if there is no SPPA, then SPPA max height is
# zero (not NA)

bg_full %>% 
  mutate(height_sppa = ifelse(height_sppa == 0 & comp == 1, NA, height_sppa)) -> bg_full

## Remove all ancillary datasets and keep blue_genes data set and full_data.
## This is just housekeeping so I don't get confused by all of the data frames
rm(list=setdiff(ls(), c("bg_full")))

