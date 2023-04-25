diff_store <- NULL

for (i in 1:1000){
  
  traits_nocomp$agb_scam_random <- sample(traits_nocomp$agb_scam, nrow(traits_nocomp),
                                          replace = F)
  
  agb_test <- lmer(sqrt(agb_scam_random) ~ weight_init + date_cloned_grp + age + co2 +  
                     salinity + elevation_sc + I(elevation_sc^2) + (co2| genotype) +  
                     age:co2 + age:salinity + age:elevation_sc + co2:salinity +  
                     co2:elevation_sc + salinity:elevation_sc + age:co2:salinity +  
                     age:co2:elevation_sc + age:salinity:elevation_sc + co2:salinity:elevation_sc +  
                     age:co2:salinity:elevation_sc, data = traits_nocomp)
  
  marginal <- MuMIn::r.squaredGLMM(agb_test)[1]
  conditional <- MuMIn::r.squaredGLMM(agb_test)[2]
  
  diff_store[i] <- conditional - marginal
  
}

hist(diff_store)
quantile(diff_store, 0.95) -> upper

agb_test <- lmer(sqrt(agb_scam) ~ weight_init + date_cloned_grp + age + co2 +  
                   salinity + elevation_sc + I(elevation_sc^2) + (co2| genotype) +  
                   age:co2 + age:salinity + age:elevation_sc + co2:salinity +  
                   co2:elevation_sc + salinity:elevation_sc + age:co2:salinity +  
                   age:co2:elevation_sc + age:salinity:elevation_sc + co2:salinity:elevation_sc +  
                   age:co2:salinity:elevation_sc, data = traits_nocomp)

MuMIn::r.squaredGLMM(agb_test)

tibble(diff = diff_store) %>% 
  mutate(color_code = ifelse(diff > upper, "upper", "lower")) %>% 
  ggplot(aes(x = diff)) +
  geom_histogram(aes(fill = color_code), bins = 21, color = "gray28") +
  labs(x = "prop. variation explained by genotype or GxE") +
  theme_classic(base_size = 14) +
  geomtextpath::geom_textvline(aes(xintercept = 0.1532535), color = "black",
                               linetype = "solid", size = 4,
                               label = "true value", linewidth = 1) +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("gray87", "orange")) +
  ggtitle("Permutation results for aboveground biomass model") -> perm_test

ggsave("~/Desktop/permutation_test.png", perm_test, height = 5.18, width = 6.3, units = "in")

