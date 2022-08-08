store_r2 <- NULL
store_EcoEvo <- NULL

for (i in 1:100){
  # Calculate what 80% of sample size is
  perc80 <- 0.8 * nrow(traits_nocomp_rs)
  # Create random sample of that length
  sample80id <- sample(1:nrow(traits_nocomp_rs), perc80, replace = F)
  # Create training data
  traits_nocomp_rs[sample80id,] -> train_data
  # Create %notin% operator
  `%notin%` <- Negate(`%in%`)
  # Create test data
  traits_nocomp_rs[1:nrow(traits_nocomp_rs) %notin% sample80id,] -> test_data
  
  # Fit model to training data
  agb_mod <- lmer(log(rs) ~ weight_init + date_cloned_grp + origin_lab +
                    (age + location + co2 + salinity + elevation_sc)^5 +
                    I(elevation_sc^2) +
                    (1+co2+salinity+elevation_sc|genotype) +
                    (1|site_frame), data = train_data)
  
  # Use lmerTest::step for backwards model selection
  step_agb <- get_model(lmerTest::step(agb_mod))
  
  anova(step_agb)
  
  # Check to see if there are any convergence warnings
  #plot(step_agb) # All good
  
  # Rename model
  agb_evo_model <- step_agb
  
  # Predict test data using trained model
  skip_to_next <- FALSE
  
  tryCatch(predict(agb_evo_model, test_data), error = function(e) { skip_to_next <<- TRUE})
  
  if(skip_to_next) { next }     
  
  predict(agb_evo_model, test_data) -> pred_agb
  
  plot(pred_agb, log(test_data$rs))
  summary(lm(pred_agb ~ log(test_data$rs)))$r.squared -> r2
  
  names(attr(anova(agb_evo_model), "hypotheses")) -> names
  
  grepl("age:", names, fixed = TRUE) | grepl("location:", names, fixed = TRUE) |
    grepl("age", names, fixed = TRUE) -> logical_vector
  
  length(which(logical_vector)) -> count
  
  store_r2[i] <- r2
  store_EcoEvo[i] <- ifelse(count > 0, 1, 0)
  print(i)
}


