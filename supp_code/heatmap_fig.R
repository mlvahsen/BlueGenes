
tibble(value = c("p < 0.001", "not in final model", "not in final model", "p < 0.01", "did not evaluate", "p < 0.05",
                 "p < 0.001", "p < 0.001", "not in final model", "p < 0.05", "p < 0.001", "not in final model",
                 "p < 0.001", "not in final model", "not in final model", "not in final model", "p < 0.01", "not in final model",
                 "p < 0.001", "p > 0.05 in final model", "not in final model", "p > 0.05 in final model", "did not evaluate", "p < 0.001",
                 "p < 0.001", "not in final model", "not in final model", "p < 0.05", "p < 0.001", "not in final model",
                 "p < 0.001", "not in final model", "p < 0.05", "p < 0.01", "p < 0.001", "not in final model",
                 "p < 0.001", "not in final model", "not in final model", "not in final model", "p < 0.001", "not in final model"),
       trait = rep(c("agb", "rs", "bg", "beta", "height", "width", "density"), each = 6),
       effect = rep(c("eco", "eco-eco", "evo", "eco-evo", "gen var", "g x e"), 7)) %>% 
  mutate(trait = factor(trait, levels = rev(c("agb", "bg", "rs", "beta", "height", "width", "density")))) %>% 
  mutate(effect = factor(effect, levels = c("eco", "eco-eco", "evo", "eco-evo", "gen var", "g x e"))) %>% 
  mutate(value = factor(value, levels = c("p < 0.001", "p < 0.01", "p < 0.05", "p > 0.05 in final model", "not in final model", "did not evaluate"))) %>% 
  mutate(value = as.factor(value)) %>% 
  ggplot(aes(x = effect, y = trait, fill = value)) +
  geom_tile() +
  scale_fill_manual(values = c("#08519c","#3182bd", "#6baed6","#bdd7e7","gray95", "gray80")) +
  scale_x_discrete(position = "top") +
  theme(legend.background = element_rect(color = "black")) +
  ylab("") + xlab("")-> heatmap

png("figs/FigX_heatmap.png", height = 5, width = 6, res = 300, units = "in")
heatmap
dev.off()
  

