traits_nocomp %>% 
  #filter(genotype == "cm10") %>% 
  group_by(genotype, salinity) %>% 
  summarize(n = n()) %>% 
  print(n = Inf)

bg_full %>% 
  filter(level == 6 & agb_scam > 0) %>% 
  group_by(genotype) %>% 
  summarize(n = n())

bg_full %>% 
  filter(level == 6 & agb_scam == 0) %>% 
  group_by(genotype) %>% 
  summarize(n = n()) %>% 
  print(n = Inf)


bg_full %>% 
  filter(level %in% c(5,6) & agb_scam > 0 ) %>% 
  select(level, agb_scam, dens_scam_live, total_bg,
         beta, width_scam_mid, height_scam_tot) %>% 
  mutate(rs = total_bg/agb_scam) %>%
  gather(key = trait, value = value, agb_scam:rs) %>% 
  ggplot(aes(x = value)) +
  geom_histogram() +
  facet_wrap(~trait, scales = "free")

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

library(patchwork)

png("figs/FigSXX_TraitDists.png", height = 8.8, width = 11.2, res = 300, units = "in")
a+b+f+g+e+d+c+ plot_layout(guides = "collect", nrow = 3) +
  plot_annotation(tag_levels = "a")
dev.off()
  
