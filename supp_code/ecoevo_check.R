# Load libraries
library(tidyverse); library(lme4)

set.seed(123)

## AGB ####
# Test how often there is an eco-evo interaction in AGB for 100 iterations using
# 80% of the data

# Create storage
store_r2 <- NULL
store_EcoEvo <- NULL

# Create %notin% operator
`%notin%` <- Negate(`%in%`)

# Set seed
set.seed(123)

# Iterate for 130 because some of these will be dropped due to NAs that I don't
# think are biasing the results (see below)

for (i in 1:130){
  # Calculate what 80% of sample size is
  perc80 <- 0.8 * nrow(traits_nocomp)
  # Create random sample of that length
  sample80id <- sample(1:nrow(traits_nocomp), perc80, replace = F)
  # Create training data
  traits_nocomp[sample80id,] -> train_data
  
  # Create test data
  traits_nocomp[1:nrow(traits_nocomp) %notin% sample80id,] -> test_data

  # Fit model to training data
  agb_mod <- lmer(sqrt(agb_scam) ~ weight_init + date_cloned_grp + origin_lab +
                    (age + location + co2 + salinity + elevation_sc)^5 +
                    I(elevation_sc^2) +
                    (1+co2+elevation_sc+salinity|genotype) +
                    (1|site_frame), data = train_data)
  
  # Use lmerTest::step for backwards model selection
  step_agb <- get_model(lmerTest::step(agb_mod))
  
  # Rename model
  agb_evo_model <- step_agb
  
  # Predict test data using trained model
  skip_to_next <- FALSE

  # "Error in levelfun(r, n, allow.new.levels = allow.new.levels) : new levels
  # detected in newdata" -- This gets rid of this error and I don't think biases
  # the datasets that are being subsetted
  tryCatch(predict(agb_evo_model, test_data), error = function(e) { skip_to_next <<- TRUE})
  if(skip_to_next) {next}

  # Predict AGB for 20% held out data
  predict(agb_evo_model, test_data) -> pred_agb

  # Calculate r2 for linear model between predicted and observed for held out
  # data
  summary(lm(pred_agb^2 ~ test_data$agb_scam))$r.squared -> r2
 
  # Pull out terms that were in the final model
  names(attr(anova(agb_evo_model), "hypotheses")) -> names
  
  # Create vector that tells us if there is a significant age cohort or
  # provenance interaction; ; if there was at least one interaction in the final
  # model there is at least one significant interaction
  grepl("age:", names, fixed = TRUE) | grepl("location:", names, fixed = TRUE) -> logical_vector
  
  # Count up the number of eco-evo interactions in the final model
  length(which(logical_vector)) -> count
  
  # Store the r2 in a vector
  store_r2[i] <- r2
  # Store if there was at least one eco-evo interaction in the final model
  store_EcoEvo[i] <- ifelse(count > 0, 1, 0)
  print(i)
}

# Pull the r2 for the train vs test regression
r2_100 <- as.numeric(na.omit(store_r2))[1:100]

# Pull the binary vector of whether or not there was a significant interaction
EcoEvo_100 <- as.numeric(na.omit(store_EcoEvo))[1:100]

# Get mean r2
mean(r2_100) # 0.4595134

# Get proportion of models that had significant eco-evo interaction
sum(EcoEvo_100)/length(EcoEvo_100) # 0.86

# Repeat process but randomize data set each time

# Create storage
store_r2 <- NULL
store_EcoEvo <- NULL

# Set seed
set.seed(123)

# Iterate for 130 because some of these will be dropped due to NAs that I don't
# think are biasing the results (see below)

for (i in 1:130){
  # Calculate what 80% of sample size is
  perc80 <- 0.8 * nrow(traits_nocomp)
  # Create random sample of that length
  sample80id <- sample(1:nrow(traits_nocomp), perc80, replace = F)
  # Create training data
  traits_nocomp[sample80id,] -> train_data
  
  # Re-sample the training data
  train_data$agb_scam <- sample(train_data$agb_scam, length(train_data$agb_scam), replace = F)

  
  # Create test data
  traits_nocomp[1:nrow(traits_nocomp) %notin% sample80id,] -> test_data
  
  # Fit model to training data
  agb_mod <- lmer(sqrt(agb_scam) ~ weight_init + date_cloned_grp + origin_lab +
                    (age + location + co2 + salinity + elevation_sc)^5 +
                    I(elevation_sc^2) +
                    (1+co2+elevation_sc+salinity|genotype) +
                    (1|site_frame), data = train_data)
  
  # Use lmerTest::step for backwards model selection
  step_agb <- get_model(lmerTest::step(agb_mod))
  
  # Rename model
  agb_evo_model <- step_agb
  
  # Predict test data using trained model
  skip_to_next <- FALSE
  
  # "Error in levelfun(r, n, allow.new.levels = allow.new.levels) : new levels
  # detected in newdata" -- This gets rid of this error and I don't think biases
  # the datasets that are being subsetted
  tryCatch(predict(agb_evo_model, test_data), error = function(e) { skip_to_next <<- TRUE})
  if(skip_to_next) {next}
  
  # Predict AGB for 20% held out data
  predict(agb_evo_model, test_data) -> pred_agb
  
  # Calculate r2 for linear model between predicted and observed for held out
  # data
  summary(lm(pred_agb^2 ~ test_data$agb_scam))$r.squared -> r2
  
  # Pull out terms that were in the final model
  names(attr(anova(agb_evo_model), "hypotheses")) -> names
  
  # Create vector that tells us if there is a significant age cohort or
  # provenance interaction; ; if there was at least one interaction in the final
  # model there is at least one significant interaction
  grepl("age:", names, fixed = TRUE) | grepl("location:", names, fixed = TRUE) -> logical_vector
  
  # Count up the number of eco-evo interactions in the final model
  length(which(logical_vector)) -> count
  
  # Store the r2 in a vector
  store_r2[i] <- r2
  # Store if there was at least one eco-evo interaction in the final model
  store_EcoEvo[i] <- ifelse(count > 0, 1, 0)
  print(i)
}

# Pull the r2 for the train vs test regression
r2_100_random <- as.numeric(na.omit(store_r2))[1:100]

# Pull the binary vector of whether or not there was a significant interaction
EcoEvo_100_random <- as.numeric(na.omit(store_EcoEvo))[1:100]

# Get mean r2
mean(r2_100_random) # 0.0867908

# Get proportion of models that had significant eco-evo interaction
sum(EcoEvo_100_random)/length(EcoEvo_100_random) # 0.14


## R:S ####
# Test how often there is an eco-evo interaction in R:S for 100 iterations using
# 80% of the data

# Create storage
store_r2 <- NULL
store_EcoEvo <- NULL

# Create %notin% operator
`%notin%` <- Negate(`%in%`)

# Set seed
set.seed(123)

# Iterate for 130 because some of these will be dropped due to NAs that I don't
# think are biasing the results (see below)

for (i in 1:130){
  # Calculate what 80% of sample size is
  perc80 <- 0.8 * nrow(traits_nocomp_rs)
  # Create random sample of that length
  sample80id <- sample(1:nrow(traits_nocomp_rs), perc80, replace = F)
  # Create training data
  traits_nocomp_rs[sample80id,] -> train_data
  
  # Create test data
  traits_nocomp_rs[1:nrow(traits_nocomp_rs) %notin% sample80id,] -> test_data
  
  # Fit model to training data
  rs_mod <- lmer(log(rs) ~ weight_init + date_cloned_grp + origin_lab +
                    (age + location + co2 + salinity + elevation_sc)^5 +
                    I(elevation_sc^2) +
                    (1+co2+elevation_sc+salinity|genotype) +
                    (1|site_frame), data = train_data)
  
  # Use lmerTest::step for backwards model selection
  step_rs <- get_model(lmerTest::step(rs_mod))
  
  # Rename model
  rs_evo_model <- step_rs
  
  # Predict test data using trained model
  skip_to_next <- FALSE
  
  # "Error in levelfun(r, n, allow.new.levels = allow.new.levels) : new levels
  # detected in newdata" -- This gets rid of this error and I don't think biases
  # the datasets that are being subsetted
  tryCatch(predict(rs_evo_model, test_data), error = function(e) { skip_to_next <<- TRUE})
  if(skip_to_next) {next}
  
  # Predict R:S for 20% held out data
  predict(rs_evo_model, test_data) -> pred_rs
  
  # Calculate r2 for linear model between predicted and observed for held out
  # data
  summary(lm(exp(pred_rs) ~ test_data$rs))$r.squared -> r2
  
  # Pull out terms that were in the final model
  names(attr(anova(rs_evo_model), "hypotheses")) -> names
  
  # Create vector that tells us if there is a significant age cohort or
  # provenance interaction; ; if there was at least one interaction in the final
  # model there is at least one significant interaction
  grepl("age:", names, fixed = TRUE) | grepl("location:", names, fixed = TRUE) -> logical_vector
  
  # Count up the number of eco-evo interactions in the final model
  length(which(logical_vector)) -> count
  
  # Store the r2 in a vector
  store_r2[i] <- r2
  # Store if there was at least one eco-evo interaction in the final model
  store_EcoEvo[i] <- ifelse(count > 0, 1, 0)
  print(i)
}

# Pull the r2 for the train vs test regression
r2_100 <- as.numeric(na.omit(store_r2))[1:100]

# Pull the binary vector of whether or not there was a significant interaction
EcoEvo_100 <- as.numeric(na.omit(store_EcoEvo))[1:100]

# Get mean r2
mean(r2_100) # 0.5141297

# Get proportion of models that had significant eco-evo interaction
sum(EcoEvo_100)/length(EcoEvo_100) # 0.72

# Repeat process but randomize data set each time

# Create storage
store_r2 <- NULL
store_EcoEvo <- NULL

# Set seed
set.seed(123)

# Iterate for 130 because some of these will be dropped due to NAs that I don't
# think are biasing the results (see below)

for (i in 1:130){
  # Calculate what 80% of sample size is
  perc80 <- 0.8 * nrow(traits_nocomp_rs)
  # Create random sample of that length
  sample80id <- sample(1:nrow(traits_nocomp_rs), perc80, replace = F)
  # Create training data
  traits_nocomp_rs[sample80id,] -> train_data
  
  # Re-sample the training data
  train_data$rs<- sample(train_data$rs, length(train_data$rs), replace = F)
  
  # Create test data
  traits_nocomp_rs[1:nrow(traits_nocomp_rs) %notin% sample80id,] -> test_data
  
  # Fit model to training data
  rs_mod <- lmer(log(rs) ~ weight_init + date_cloned_grp + origin_lab +
                    (age + location + co2 + salinity + elevation_sc)^5 +
                    I(elevation_sc^2) +
                    (1+co2+elevation_sc+salinity|genotype) +
                    (1|site_frame), data = train_data)
  
  # Use lmerTest::step for backwards model selection
  step_rs <- get_model(lmerTest::step(rs_mod))
  
  # Rename model
  rs_evo_model <- step_rs
  
  # Predict test data using trained model
  skip_to_next <- FALSE
  
  # "Error in levelfun(r, n, allow.new.levels = allow.new.levels) : new levels
  # detected in newdata" -- This gets rid of this error and I don't think biases
  # the datasets that are being subsetted
  tryCatch(predict(rs_evo_model, test_data), error = function(e) { skip_to_next <<- TRUE})
  if(skip_to_next) {next}
  
  # Predict AGB for 20% held out data
  predict(rs_evo_model, test_data) -> pred_rs
  
  # Calculate r2 for linear model between predicted and observed for held out
  # data
  summary(lm(exp(pred_rs) ~ test_data$rs))$r.squared -> r2
  
  # Pull out terms that were in the final model
  names(attr(anova(rs_evo_model), "hypotheses")) -> names
  
  # Create vector that tells us if there is a significant age cohort or
  # provenance interaction; ; if there was at least one interaction in the final
  # model there is at least one significant interaction
  grepl("age:", names, fixed = TRUE) | grepl("location:", names, fixed = TRUE) -> logical_vector
  
  # Count up the number of eco-evo interactions in the final model
  length(which(logical_vector)) -> count
  
  # Store the r2 in a vector
  store_r2[i] <- r2
  # Store if there was at least one eco-evo interaction in the final model
  store_EcoEvo[i] <- ifelse(count > 0, 1, 0)
  print(i)
}

# Pull the r2 for the train vs test regression
r2_100_random <- as.numeric(na.omit(store_r2))[1:100]

# Pull the binary vector of whether or not there was a significant interaction
EcoEvo_100_random <- as.numeric(na.omit(store_EcoEvo))[1:100]

# Get mean r2
mean(r2_100_random) # 0.07852436

# Get proportion of models that had significant eco-evo interaction
sum(EcoEvo_100_random)/length(EcoEvo_100_random) # 0.13

## beta ####
# Test how often there is an eco-evo interaction in R:S for 100 iterations using
# 80% of the data

# Create storage
store_r2 <- NULL
store_EcoEvo <- NULL

# Create %notin% operator
`%notin%` <- Negate(`%in%`)

# Set seed
set.seed(123)

# Iterate for 130 because some of these will be dropped due to NAs that I don't
# think are biasing the results (see below)

for (i in 1:130){
  # Calculate what 80% of sample size is
  perc80 <- 0.8 * nrow(traits_nocomp_rs)
  # Create random sample of that length
  sample80id <- sample(1:nrow(traits_nocomp_rs), perc80, replace = F)
  # Create training data
  traits_nocomp_rs[sample80id,] -> train_data
  
  # Create test data
  traits_nocomp_rs[1:nrow(traits_nocomp_rs) %notin% sample80id,] -> test_data
  
  # Fit model to training data
  beta_mod <- lmer(beta ~ weight_init + date_cloned_grp + origin_lab +
                   (age + location + co2 + salinity + elevation_sc)^5 +
                   I(elevation_sc^2) +
                   (1+co2+elevation_sc+salinity|genotype) +
                   (1|site_frame), data = train_data)
  
  # Use lmerTest::step for backwards model selection
  step_beta <- get_model(lmerTest::step(beta_mod))
  
  # Rename model
  beta_evo_model <- step_beta
  
  # Predict test data using trained model
  skip_to_next <- FALSE
  
  # "Error in levelfun(r, n, allow.new.levels = allow.new.levels) : new levels
  # detected in newdata" -- This gets rid of this error and I don't think biases
  # the datasets that are being subsetted
  tryCatch(predict(beta_evo_model, test_data), error = function(e) { skip_to_next <<- TRUE})
  if(skip_to_next) {next}
  
  # Predict R:S for 20% held out data
  predict(beta_evo_model, test_data) -> pred_beta
  
  # Calculate r2 for linear model between predicted and observed for held out
  # data
  summary(lm(pred_beta ~ test_data$beta))$r.squared -> r2
  
  # Pull out terms that were in the final model
  names(attr(anova(beta_evo_model), "hypotheses")) -> names
  
  # Create vector that tells us if there is a significant age cohort or
  # provenance interaction; ; if there was at least one interaction in the final
  # model there is at least one significant interaction
  grepl("age:", names, fixed = TRUE) | grepl("location:", names, fixed = TRUE) -> logical_vector
  
  # Count up the number of eco-evo interactions in the final model
  length(which(logical_vector)) -> count
  
  # Store the r2 in a vector
  store_r2[i] <- r2
  # Store if there was at least one eco-evo interaction in the final model
  store_EcoEvo[i] <- ifelse(count > 0, 1, 0)
  print(i)
}

# Pull the r2 for the train vs test regression
r2_100 <- as.numeric(na.omit(store_r2))[1:100]

# Pull the binary vector of whether or not there was a significant interaction
EcoEvo_100 <- as.numeric(na.omit(store_EcoEvo))[1:100]

# Get mean r2
mean(r2_100) # 0.5141297

# Get proportion of models that had significant eco-evo interaction
sum(EcoEvo_100)/length(EcoEvo_100) # 0.72

# Repeat process but randomize data set each time

# Create storage
store_r2 <- NULL
store_EcoEvo <- NULL

# Set seed
set.seed(123)

# Iterate for 130 because some of these will be dropped due to NAs that I don't
# think are biasing the results (see below)

for (i in 1:130){
  # Calculate what 80% of sample size is
  perc80 <- 0.8 * nrow(traits_nocomp_rs)
  # Create random sample of that length
  sample80id <- sample(1:nrow(traits_nocomp_rs), perc80, replace = F)
  # Create training data
  traits_nocomp_rs[sample80id,] -> train_data
  
  # Re-sample the training data
  train_data$beta<- sample(train_data$beta, length(train_data$beta), replace = F)
  
  # Create test data
  traits_nocomp_rs[1:nrow(traits_nocomp_rs) %notin% sample80id,] -> test_data
  
  # Fit model to training data
  beta_mod <- lmer(beta ~ weight_init + date_cloned_grp + origin_lab +
                   (age + location + co2 + salinity + elevation_sc)^5 +
                   I(elevation_sc^2) +
                   (1+co2+elevation_sc+salinity|genotype) +
                   (1|site_frame), data = train_data)
  
  # Use lmerTest::step for backwards model selection
  step_beta <- get_model(lmerTest::step(beta_mod))
  
  # Rename model
  beta_evo_model <- step_beta
  
  # Predict test data using trained model
  skip_to_next <- FALSE
  
  # "Error in levelfun(r, n, allow.new.levels = allow.new.levels) : new levels
  # detected in newdata" -- This gets rid of this error and I don't think biases
  # the datasets that are being subsetted
  tryCatch(predict(beta_evo_model, test_data), error = function(e) { skip_to_next <<- TRUE})
  if(skip_to_next) {next}
  
  # Predict AGB for 20% held out data
  predict(beta_evo_model, test_data) -> pred_beta
  
  # Calculate r2 for linear model between predicted and observed for held out
  # data
  summary(lm(pred_beta ~ test_data$beta))$r.squared -> r2
  
  # Pull out terms that were in the final model
  names(attr(anova(beta_evo_model), "hypotheses")) -> names
  
  # Create vector that tells us if there is a significant age cohort or
  # provenance interaction; ; if there was at least one interaction in the final
  # model there is at least one significant interaction
  grepl("age:", names, fixed = TRUE) | grepl("location:", names, fixed = TRUE) -> logical_vector
  
  # Count up the number of eco-evo interactions in the final model
  length(which(logical_vector)) -> count
  
  # Store the r2 in a vector
  store_r2[i] <- r2
  # Store if there was at least one eco-evo interaction in the final model
  store_EcoEvo[i] <- ifelse(count > 0, 1, 0)
  print(i)
}

# Pull the r2 for the train vs test regression
r2_100_random <- as.numeric(na.omit(store_r2))[1:100]

# Pull the binary vector of whether or not there was a significant interaction
EcoEvo_100_random <- as.numeric(na.omit(store_EcoEvo))[1:100]

# Get mean r2
mean(r2_100_random) # 0.1139169

# Get proportion of models that had significant eco-evo interaction
sum(EcoEvo_100_random)/length(EcoEvo_100_random) # 0.12

## BGB ####
# Test how often there is an eco-evo interaction in BGB for 100 iterations using
# 80% of the data

# Create storage
store_r2 <- NULL
store_EcoEvo <- NULL

# Create %notin% operator
`%notin%` <- Negate(`%in%`)

# Set seed
set.seed(123)

# Iterate for 130 because some of these will be dropped due to NAs that I don't
# think are biasing the results (see below)

for (i in 1:130){
  # Calculate what 80% of sample size is
  perc80 <- 0.8 * nrow(traits_nocomp_rs)
  # Create random sample of that length
  sample80id <- sample(1:nrow(traits_nocomp_rs), perc80, replace = F)
  # Create training data
  traits_nocomp_rs[sample80id,] -> train_data
  
  # Create test data
  traits_nocomp_rs[1:nrow(traits_nocomp_rs) %notin% sample80id,] -> test_data
  
  # Fit model to training data
  bgb_mod <- lmer(sqrt(total_bg) ~ weight_init + date_cloned_grp + origin_lab +
                     (age + location + co2 + salinity + elevation_sc)^5 +
                     I(elevation_sc^2) +
                     (1+co2+elevation_sc+salinity|genotype) +
                     (1|site_frame), data = train_data)
  
  # Use lmerTest::step for backwards model selection
  step_bgb <- get_model(lmerTest::step(bgb_mod))
  
  # Rename model
  bgb_evo_model <- step_bgb
  
  # Predict test data using trained model
  skip_to_next <- FALSE
  
  # "Error in levelfun(r, n, allow.new.levels = allow.new.levels) : new levels
  # detected in newdata" -- This gets rid of this error and I don't think biases
  # the datasets that are being subsetted
  tryCatch(predict(bgb_evo_model, test_data), error = function(e) { skip_to_next <<- TRUE})
  if(skip_to_next) {next}
  
  # Predict BGB for 20% held out data
  predict(bgb_evo_model, test_data) -> pred_bgb
  
  # Calculate r2 for linear model between predicted and observed for held out
  # data
  summary(lm(pred_bgb^2 ~ test_data$total_bg))$r.squared -> r2
  
  # Pull out terms that were in the final model
  names(attr(anova(bgb_evo_model), "hypotheses")) -> names
  
  # Create vector that tells us if there is a significant age cohort or
  # provenance interaction; ; if there was at least one interaction in the final
  # model there is at least one significant interaction
  grepl("age:", names, fixed = TRUE) | grepl("location:", names, fixed = TRUE) -> logical_vector
  
  # Count up the number of eco-evo interactions in the final model
  length(which(logical_vector)) -> count
  
  # Store the r2 in a vector
  store_r2[i] <- r2
  # Store if there was at least one eco-evo interaction in the final model
  store_EcoEvo[i] <- ifelse(count > 0, 1, 0)
  print(i)
}

# Pull the r2 for the train vs test regression
r2_100 <- as.numeric(na.omit(store_r2))[1:100]

# Pull the binary vector of whether or not there was a significant interaction
EcoEvo_100 <- as.numeric(na.omit(store_EcoEvo))[1:100]

# Get mean r2
mean(r2_100) # 0.1457327

# Get proportion of models that had significant eco-evo interaction
sum(EcoEvo_100)/length(EcoEvo_100) # 0.54

# Repeat process but randomize data set each time

# Create storage
store_r2 <- NULL
store_EcoEvo <- NULL

# Set seed
set.seed(123)

# Iterate for 130 because some of these will be dropped due to NAs that I don't
# think are biasing the results (see below)

for (i in 1:130){
  # Calculate what 80% of sample size is
  perc80 <- 0.8 * nrow(traits_nocomp_rs)
  # Create random sample of that length
  sample80id <- sample(1:nrow(traits_nocomp_rs), perc80, replace = F)
  # Create training data
  traits_nocomp_rs[sample80id,] -> train_data
  
  # Re-sample the training data
  train_data$total_bg<- sample(train_data$total_bg, length(train_data$total_bg), replace = F)
  
  # Create test data
  traits_nocomp_rs[1:nrow(traits_nocomp_rs) %notin% sample80id,] -> test_data
  
  # Fit model to training data
  bgb_mod <- lmer(sqrt(total_bg) ~ weight_init + date_cloned_grp + origin_lab +
                     (age + location + co2 + salinity + elevation_sc)^5 +
                     I(elevation_sc^2) +
                     (1+co2+elevation_sc+salinity|genotype) +
                     (1|site_frame), data = train_data)
  
  # Use lmerTest::step for backwards model selection
  step_bgb <- get_model(lmerTest::step(bgb_mod))
  
  # Rename model
  bgb_evo_model <- step_bgb
  
  # Predict test data using trained model
  skip_to_next <- FALSE
  
  # "Error in levelfun(r, n, allow.new.levels = allow.new.levels) : new levels
  # detected in newdata" -- This gets rid of this error and I don't think biases
  # the datasets that are being subsetted
  tryCatch(predict(bgb_evo_model, test_data), error = function(e) { skip_to_next <<- TRUE})
  if(skip_to_next) {next}
  
  # Predict BGB for 20% held out data
  predict(bgb_evo_model, test_data) -> pred_bgb
  
  # Calculate r2 for linear model between predicted and observed for held out
  # data
  summary(lm(pred_bgb^2 ~ test_data$total_bg))$r.squared -> r2
  
  # Pull out terms that were in the final model
  names(attr(anova(bgb_evo_model), "hypotheses")) -> names
  
  # Create vector that tells us if there is a significant age cohort or
  # provenance interaction; ; if there was at least one interaction in the final
  # model there is at least one significant interaction
  grepl("age:", names, fixed = TRUE) | grepl("location:", names, fixed = TRUE) -> logical_vector
  
  # Count up the number of eco-evo interactions in the final model
  length(which(logical_vector)) -> count
  
  # Store the r2 in a vector
  store_r2[i] <- r2
  # Store if there was at least one eco-evo interaction in the final model
  store_EcoEvo[i] <- ifelse(count > 0, 1, 0)
  print(i)
}

# Pull the r2 for the train vs test regression
r2_100_random <- as.numeric(na.omit(store_r2))[1:100]

# Pull the binary vector of whether or not there was a significant interaction
EcoEvo_100_random <- as.numeric(na.omit(store_EcoEvo))[1:100]

# Get mean r2
mean(r2_100_random) # 0.07715906

# Get proportion of models that had significant eco-evo interaction
sum(EcoEvo_100_random)/length(EcoEvo_100_random) # 0.12

## density ####
# Test how often there is an eco-evo interaction in AGB for 100 iterations using
# 80% of the data

# Create storage
store_r2 <- NULL
store_EcoEvo <- NULL

# Create %notin% operator
`%notin%` <- Negate(`%in%`)

# Set seed
set.seed(123)

# Iterate for 130 because some of these will be dropped due to NAs that I don't
# think are biasing the results (see below)

for (i in 1:130){
  # Calculate what 80% of sample size is
  perc80 <- 0.8 * nrow(traits_nocomp)
  # Create random sample of that length
  sample80id <- sample(1:nrow(traits_nocomp), perc80, replace = F)
  # Create training data
  traits_nocomp[sample80id,] -> train_data
  
  # Create test data
  traits_nocomp[1:nrow(traits_nocomp) %notin% sample80id,] -> test_data
  
  # Fit model to training data
  density_mod <- lmer(dens_scam_live ~ weight_init + date_cloned_grp + origin_lab +
                    (age + location + co2 + salinity + elevation_sc)^5 +
                    I(elevation_sc^2) +
                    (1+co2+elevation_sc+salinity|genotype) +
                    (1|site_frame), data = train_data)
  
  # Use lmerTest::step for backwards model selection
  step_density <- get_model(lmerTest::step(density_mod))
  
  # Rename model
  density_evo_model <- step_density
  
  # Predict test data using trained model
  skip_to_next <- FALSE
  
  # "Error in levelfun(r, n, allow.new.levels = allow.new.levels) : new levels
  # detected in newdata" -- This gets rid of this error and I don't think biases
  # the datasets that are being subsetted
  tryCatch(predict(density_evo_model, test_data), error = function(e) { skip_to_next <<- TRUE})
  if(skip_to_next) {next}
  
  # Predict Density for 20% held out data
  predict(density_evo_model, test_data) -> pred_density
  
  # Calculate r2 for linear model between predicted and observed for held out
  # data
  summary(lm(pred_density ~ test_data$dens_scam_live))$r.squared -> r2
  
  # Pull out terms that were in the final model
  names(attr(anova(density_evo_model), "hypotheses")) -> names
  
  # Create vector that tells us if there is a significant age cohort or
  # provenance interaction; ; if there was at least one interaction in the final
  # model there is at least one significant interaction
  grepl("age:", names, fixed = TRUE) | grepl("location:", names, fixed = TRUE) -> logical_vector
  
  # Count up the number of eco-evo interactions in the final model
  length(which(logical_vector)) -> count
  
  # Store the r2 in a vector
  store_r2[i] <- r2
  # Store if there was at least one eco-evo interaction in the final model
  store_EcoEvo[i] <- ifelse(count > 0, 1, 0)
  print(i)
}

# Pull the r2 for the train vs test regression
r2_100 <- as.numeric(na.omit(store_r2))[1:100]

# Pull the binary vector of whether or not there was a significant interaction
EcoEvo_100 <- as.numeric(na.omit(store_EcoEvo))[1:100]

# Get mean r2
mean(r2_100) # 0.2889353

# Get proportion of models that had significant eco-evo interaction
sum(EcoEvo_100)/length(EcoEvo_100) # 0.48

# Repeat process but randomize data set each time

# Create storage
store_r2 <- NULL
store_EcoEvo <- NULL

# Set seed
set.seed(123)

# Iterate for 130 because some of these will be dropped due to NAs that I don't
# think are biasing the results (see below)

for (i in 1:130){
  # Calculate what 80% of sample size is
  perc80 <- 0.8 * nrow(traits_nocomp)
  # Create random sample of that length
  sample80id <- sample(1:nrow(traits_nocomp), perc80, replace = F)
  # Create training data
  traits_nocomp[sample80id,] -> train_data
  
  # Re-sample the training data
  train_data$dens_scam_live <- sample(train_data$dens_scam_live, length(train_data$dens_scam_live), replace = F)
  
  # Create test data
  traits_nocomp[1:nrow(traits_nocomp) %notin% sample80id,] -> test_data
  
  # Fit model to training data
  density_mod <- lmer(dens_scam_live ~ weight_init + date_cloned_grp + origin_lab +
                    (age + location + co2 + salinity + elevation_sc)^5 +
                    I(elevation_sc^2) +
                    (1+co2+elevation_sc+salinity|genotype) +
                    (1|site_frame), data = train_data)
  
  # Use lmerTest::step for backwards model selection
  step_density <- get_model(lmerTest::step(density_mod))
  
  # Rename model
  density_evo_model <- step_density
  
  # Predict test data using trained model
  skip_to_next <- FALSE
  
  # "Error in levelfun(r, n, allow.new.levels = allow.new.levels) : new levels
  # detected in newdata" -- This gets rid of this error and I don't think biases
  # the datasets that are being subsetted
  tryCatch(predict(density_evo_model, test_data), error = function(e) { skip_to_next <<- TRUE})
  if(skip_to_next) {next}
  
  # Predict Density for 20% held out data
  predict(density_evo_model, test_data) -> pred_density
  
  # Calculate r2 for linear model between predicted and observed for held out
  # data
  summary(lm(pred_density ~ test_data$dens_scam_live))$r.squared -> r2
  
  # Pull out terms that were in the final model
  names(attr(anova(density_evo_model), "hypotheses")) -> names
  
  # Create vector that tells us if there is a significant age cohort or
  # provenance interaction; ; if there was at least one interaction in the final
  # model there is at least one significant interaction
  grepl("age:", names, fixed = TRUE) | grepl("location:", names, fixed = TRUE) -> logical_vector
  
  # Count up the number of eco-evo interactions in the final model
  length(which(logical_vector)) -> count
  
  # Store the r2 in a vector
  store_r2[i] <- r2
  # Store if there was at least one eco-evo interaction in the final model
  store_EcoEvo[i] <- ifelse(count > 0, 1, 0)
  print(i)
}

# Pull the r2 for the train vs test regression
r2_100_random <- as.numeric(na.omit(store_r2))[1:100]

# Pull the binary vector of whether or not there was a significant interaction
EcoEvo_100_random <- as.numeric(na.omit(store_EcoEvo))[1:100]

# Get mean r2
mean(r2_100_random) # 0.0867908

# Get proportion of models that had significant eco-evo interaction
sum(EcoEvo_100_random)/length(EcoEvo_100_random) # 0.14


## stem width ####
# Test how often there is an eco-evo interaction in width for 100 iterations using
# 80% of the data

# Create storage
store_r2 <- NULL
store_EcoEvo <- NULL

# Create %notin% operator
`%notin%` <- Negate(`%in%`)

# Set seed
set.seed(123)

# Iterate for 130 because some of these will be dropped due to NAs that I don't
# think are biasing the results (see below)

for (i in 1:130){
  # Calculate what 80% of sample size is
  perc80 <- 0.8 * nrow(traits_nocomp)
  # Create random sample of that length
  sample80id <- sample(1:nrow(traits_nocomp), perc80, replace = F)
  # Create training data
  traits_nocomp[sample80id,] -> train_data
  
  # Create test data
  traits_nocomp[1:nrow(traits_nocomp) %notin% sample80id,] -> test_data
  
  # Fit model to training data
  width_mod <- lmer(width_scam_mid ~ weight_init + date_cloned_grp + origin_lab +
                        (age + location + co2 + salinity + elevation_sc)^5 +
                        I(elevation_sc^2) +
                        (1+co2+elevation_sc+salinity|genotype) +
                        (1|site_frame), data = train_data)
  
  # Use lmerTest::step for backwards model selection
  step_width <- get_model(lmerTest::step(width_mod))
  
  # Rename model
  width_evo_model <- step_width
  
  # Predict test data using trained model
  skip_to_next <- FALSE
  
  # "Error in levelfun(r, n, allow.new.levels = allow.new.levels) : new levels
  # detected in newdata" -- This gets rid of this error and I don't think biases
  # the datasets that are being subsetted
  tryCatch(predict(width_evo_model, test_data), error = function(e) { skip_to_next <<- TRUE})
  if(skip_to_next) {next}
  
  # Predict Density for 20% held out data
  predict(width_evo_model, test_data) -> pred_width
  
  # Calculate r2 for linear model between predicted and observed for held out
  # data
  summary(lm(pred_width ~ test_data$width_scam_mid))$r.squared -> r2
  
  # Pull out terms that were in the final model
  names(attr(anova(width_evo_model), "hypotheses")) -> names
  
  # Create vector that tells us if there is a significant age cohort or
  # provenance interaction; ; if there was at least one interaction in the final
  # model there is at least one significant interaction
  grepl("age:", names, fixed = TRUE) | grepl("location:", names, fixed = TRUE) -> logical_vector
  
  # Count up the number of eco-evo interactions in the final model
  length(which(logical_vector)) -> count
  
  # Store the r2 in a vector
  store_r2[i] <- r2
  # Store if there was at least one eco-evo interaction in the final model
  store_EcoEvo[i] <- ifelse(count > 0, 1, 0)
  print(i)
}

# Pull the r2 for the train vs test regression
r2_100 <- as.numeric(na.omit(store_r2))[1:100]

# Pull the binary vector of whether or not there was a significant interaction
EcoEvo_100 <- as.numeric(na.omit(store_EcoEvo))[1:100]

# Get mean r2
mean(r2_100) # 0.2889353

# Get proportion of models that had significant eco-evo interaction
sum(EcoEvo_100)/length(EcoEvo_100) # 0.48

# Repeat process but randomize data set each time

# Create storage
store_r2 <- NULL
store_EcoEvo <- NULL

# Set seed
set.seed(123)

# Iterate for 130 because some of these will be dropped due to NAs that I don't
# think are biasing the results (see below)

for (i in 1:130){
  # Calculate what 80% of sample size is
  perc80 <- 0.8 * nrow(traits_nocomp)
  # Create random sample of that length
  sample80id <- sample(1:nrow(traits_nocomp), perc80, replace = F)
  # Create training data
  traits_nocomp[sample80id,] -> train_data
  
  # Re-sample the training data
  train_data$width_scam_mid <- sample(train_data$width_scam_mid, length(train_data$width_scam_mid), replace = F)
  
  # Create test data
  traits_nocomp[1:nrow(traits_nocomp) %notin% sample80id,] -> test_data
  
  # Fit model to training data
  width_mod <- lmer(width_scam_mid ~ weight_init + date_cloned_grp + origin_lab +
                        (age + location + co2 + salinity + elevation_sc)^5 +
                        I(elevation_sc^2) +
                        (1+co2+elevation_sc+salinity|genotype) +
                        (1|site_frame), data = train_data)
  
  # Use lmerTest::step for backwards model selection
  step_width <- get_model(lmerTest::step(width_mod))
  
  # Rename model
  width_evo_model <- step_width
  
  # Predict test data using trained model
  skip_to_next <- FALSE
  
  # "Error in levelfun(r, n, allow.new.levels = allow.new.levels) : new levels
  # detected in newdata" -- This gets rid of this error and I don't think biases
  # the datasets that are being subsetted
  tryCatch(predict(width_evo_model, test_data), error = function(e) { skip_to_next <<- TRUE})
  if(skip_to_next) {next}
  
  # Predict Density for 20% held out data
  predict(width_evo_model, test_data) -> pred_width
  
  # Calculate r2 for linear model between predicted and observed for held out
  # data
  summary(lm(pred_width ~ test_data$width_scam_mid))$r.squared -> r2
  
  # Pull out terms that were in the final model
  names(attr(anova(width_evo_model), "hypotheses")) -> names
  
  # Create vector that tells us if there is a significant age cohort or
  # provenance interaction; ; if there was at least one interaction in the final
  # model there is at least one significant interaction
  grepl("age:", names, fixed = TRUE) | grepl("location:", names, fixed = TRUE) -> logical_vector
  
  # Count up the number of eco-evo interactions in the final model
  length(which(logical_vector)) -> count
  
  # Store the r2 in a vector
  store_r2[i] <- r2
  # Store if there was at least one eco-evo interaction in the final model
  store_EcoEvo[i] <- ifelse(count > 0, 1, 0)
  print(i)
}

# Pull the r2 for the train vs test regression
r2_100_random <- as.numeric(na.omit(store_r2))[1:100]

# Pull the binary vector of whether or not there was a significant interaction
EcoEvo_100_random <- as.numeric(na.omit(store_EcoEvo))[1:100]

# Get mean r2
mean(r2_100_random) # 0.09430775

# Get proportion of models that had significant eco-evo interaction
sum(EcoEvo_100_random)/length(EcoEvo_100_random) # 0.16

## stem height ####
# Test how often there is an eco-evo interaction in height for 100 iterations using
# 80% of the data

# Create storage
store_r2 <- NULL
store_EcoEvo <- NULL

# Create %notin% operator
`%notin%` <- Negate(`%in%`)

# Set seed
set.seed(123)

# Iterate for 130 because some of these will be dropped due to NAs that I don't
# think are biasing the results (see below)

for (i in 1:130){
  # Calculate what 80% of sample size is
  perc80 <- 0.8 * nrow(traits_nocomp_height)
  # Create random sample of that length
  sample80id <- sample(1:nrow(traits_nocomp_height), perc80, replace = F)
  # Create training data
  traits_nocomp_height[sample80id,] -> train_data
  
  # Create test data
  traits_nocomp_height[1:nrow(traits_nocomp_height) %notin% sample80id,] -> test_data
  
  # Fit model to training data
  height_mod <- lmer(height_scam_tot ~ weight_init + date_cloned_grp + origin_lab +
                        (age + location + co2 + salinity + elevation_sc)^5 +
                        I(elevation_sc^2) +
                        (1+co2+elevation_sc+salinity|genotype) +
                        (1|site_frame), data = train_data)
  
  # Use lmerTest::step for backwards model selection
  step_height <- get_model(lmerTest::step(height_mod))
  
  # Rename model
  height_evo_model <- step_height
  
  # Predict test data using trained model
  skip_to_next <- FALSE
  
  # "Error in levelfun(r, n, allow.new.levels = allow.new.levels) : new levels
  # detected in newdata" -- This gets rid of this error and I don't think biases
  # the datasets that are being subsetted
  tryCatch(predict(height_evo_model, test_data), error = function(e) { skip_to_next <<- TRUE})
  if(skip_to_next) {next}
  
  # Predict Density for 20% held out data
  predict(height_evo_model, test_data) -> pred_height
  
  # Calculate r2 for linear model between predicted and observed for held out
  # data
  summary(lm(pred_height ~ test_data$height_scam_tot))$r.squared -> r2
  
  # Pull out terms that were in the final model
  names(attr(anova(height_evo_model), "hypotheses")) -> names
  
  # Create vector that tells us if there is a significant age cohort or
  # provenance interaction; ; if there was at least one interaction in the final
  # model there is at least one significant interaction
  grepl("age:", names, fixed = TRUE) | grepl("location:", names, fixed = TRUE) -> logical_vector
  
  # Count up the number of eco-evo interactions in the final model
  length(which(logical_vector)) -> count
  
  # Store the r2 in a vector
  store_r2[i] <- r2
  # Store if there was at least one eco-evo interaction in the final model
  store_EcoEvo[i] <- ifelse(count > 0, 1, 0)
  print(i)
}

# Pull the r2 for the train vs test regression
r2_100 <- as.numeric(na.omit(store_r2))[1:100]

# Pull the binary vector of whether or not there was a significant interaction
EcoEvo_100 <- as.numeric(na.omit(store_EcoEvo))[1:100]

# Get mean r2
mean(r2_100) # 0.4976839

# Get proportion of models that had significant eco-evo interaction
sum(EcoEvo_100)/length(EcoEvo_100) # 0.71

# Repeat process but randomize data set each time

# Create storage
store_r2 <- NULL
store_EcoEvo <- NULL

# Set seed
set.seed(123)

# Iterate for 130 because some of these will be dropped due to NAs that I don't
# think are biasing the results (see below)

for (i in 1:130){
  # Calculate what 80% of sample size is
  perc80 <- 0.8 * nrow(traits_nocomp_height)
  # Create random sample of that length
  sample80id <- sample(1:nrow(traits_nocomp_height), perc80, replace = F)
  # Create training data
  traits_nocomp_height[sample80id,] -> train_data
  
  # Re-sample the training data
  train_data$height_scam_tot <- sample(train_data$height_scam_tot, length(train_data$height_scam_tot), replace = F)
  
  # Create test data
  traits_nocomp_height[1:nrow(traits_nocomp_height) %notin% sample80id,] -> test_data
  
  # Fit model to training data
  height_mod <- lmer(height_scam_tot ~ weight_init + date_cloned_grp + origin_lab +
                        (age + location + co2 + salinity + elevation_sc)^5 +
                        I(elevation_sc^2) +
                        (1+co2+elevation_sc+salinity|genotype) +
                        (1|site_frame), data = train_data)
  
  # Use lmerTest::step for backwards model selection
  step_height <- get_model(lmerTest::step(height_mod))
  
  # Rename model
  height_evo_model <- step_height
  
  # Predict test data using trained model
  skip_to_next <- FALSE
  
  # "Error in levelfun(r, n, allow.new.levels = allow.new.levels) : new levels
  # detected in newdata" -- This gets rid of this error and I don't think biases
  # the datasets that are being subsetted
  tryCatch(predict(height_evo_model, test_data), error = function(e) { skip_to_next <<- TRUE})
  if(skip_to_next) {next}
  
  # Predict Density for 20% held out data
  predict(height_evo_model, test_data) -> pred_height
  
  # Calculate r2 for linear model between predicted and observed for held out
  # data
  summary(lm(pred_height ~ test_data$height_scam_tot))$r.squared -> r2
  
  # Pull out terms that were in the final model
  names(attr(anova(height_evo_model), "hypotheses")) -> names
  
  # Create vector that tells us if there is a significant age cohort or
  # provenance interaction; ; if there was at least one interaction in the final
  # model there is at least one significant interaction
  grepl("age:", names, fixed = TRUE) | grepl("location:", names, fixed = TRUE) -> logical_vector
  
  # Count up the number of eco-evo interactions in the final model
  length(which(logical_vector)) -> count
  
  # Store the r2 in a vector
  store_r2[i] <- r2
  # Store if there was at least one eco-evo interaction in the final model
  store_EcoEvo[i] <- ifelse(count > 0, 1, 0)
  print(i)
}

# Pull the r2 for the train vs test regression
r2_100_random <- as.numeric(na.omit(store_r2))[1:100]

# Pull the binary vector of whether or not there was a significant interaction
EcoEvo_100_random <- as.numeric(na.omit(store_EcoEvo))[1:100]

# Get mean r2
mean(r2_100_random) # 0.09539823

# Get proportion of models that had significant eco-evo interaction
sum(EcoEvo_100_random)/length(EcoEvo_100_random) # 0.14
