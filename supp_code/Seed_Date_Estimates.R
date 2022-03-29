# Estimate age of seeds as predicted from 3 Kirkpatrick cores collected in 2002.
# Data from Gentile 2015 dissertation documents

library(tidyverse); library(lubridate); library(strex); library(here)

# Read in soil dating data
soil <- read_csv(here("supp_data","SoilDateEst.csv"))
soil <- filter(soil, complete.cases(Depth))
# Read in genotype data
genotype <- read_csv(here("supp_data","BlueGenesGenotypes.csv"))

# Pull out unique seed depths
genotype %>% 
  mutate(depth_top = case_when(depth_top == "modern" ~ "0",
                               T ~ depth_top)) %>% 
  dplyr::select(depth_top, date_collected) %>% 
  dplyr::mutate(year_collected = year(mdy(date_collected))) %>% 
  # For now we are going to assume that the ones that are unknown were collected
  # in 2002
  dplyr::mutate(year_collected = case_when(is.na(year_collected) ~ 2002,
                                    T ~ year_collected)) %>% 
  filter(complete.cases(depth_top)) %>% 
  dplyr::select(depth_top, year_collected) %>% 
  unique() %>% 
  mutate(depth_top = as.numeric(depth_top)) %>% 
  arrange(depth_top) -> unique_depths

soil_depth_params <- readRDS(here("supp_data","calibration_priors.rds"))
# Get random draws for regression coefficients relating seed depth to seed age
beta_draws <- mvtnorm::rmvnorm(1000, soil_depth_params$beta_prior_mean, soil_depth_params$beta_prior_covar)
sigma_draws <- rgamma(1000, soil_depth_params$sigma_prior$alpha, soil_depth_params$sigma_prior$beta)

predicted_age <- matrix(NA, nrow = 1000, ncol = nrow(unique_depths))

for(i in 1:1000){
  for (j in 1:nrow(unique_depths)){
    predicted_age_mean <- beta_draws[i,1] + beta_draws[i,2]*unique_depths$depth_top[j] + beta_draws[i,3]*unique_depths$depth_top[j]^2
    predicted_age[i,j] <- rnorm(1, predicted_age_mean, sigma_draws[i])
  }
}

colnames(predicted_age) <- unique_depths$depth_top

predage_medians <- apply(predicted_age, 2, median)
unique_depths$pred_age <- predage_medians

# Calculate predicted decade 
unique_depths %>% 
  mutate(pred_year = year_collected - pred_age) -> unique_depths

# Merge together
genotype %>% 
  mutate(depth_top = case_when(depth_top == "modern" ~ "0",
                               T ~ depth_top)) %>% 
  mutate(year_collected = year(mdy(date_collected))) %>% 
  # For now we are going to assume that the ones that are unknown were collected
  # in 2002
  mutate(year_collected = case_when(is.na(year_collected) ~ 2002,
                                    T ~ year_collected)) %>% 
  mutate(unique_id = paste(depth_top, year_collected, sep = "_")) -> genotype_year

unique_depths$unique_id <- paste(unique_depths$depth_top, unique_depths$year_collected, sep = "_")

merge(genotype_year, unique_depths, by = c("unique_id")) -> genotype_with_ages

# Now group by cohort and get age ranges
genotype_with_ages %>% 
  mutate(cohort = as.character(str_extract_non_numerics(bluegenes_code))) %>% 
  group_by(cohort) %>% 
  # round to nearest decade
  summarize(min_year = round(min(pred_year),-1),
            max_year = round(max(pred_year), -1))

# Get genotype info for supplementary table
genotype_with_ages %>% 
  filter(substr(genotype_with_ages$bluegenes_code, 1, 1) != "b") %>% 
  dplyr::select(bluegenes_code, lab, lat, long, pred_year) %>% 
  mutate(pred_year = round(pred_year)) %>% 
  arrange(bluegenes_code)







