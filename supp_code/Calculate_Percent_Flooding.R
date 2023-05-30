# This code calculate percent flooding for pots in Blue Genes experiment based on
# Annapolis tidal gauge data for that time period (used in Figure 2)

# Load libraries
library(tidyverse); library(here)

# Read in trait data to get actual elevations
source(here("supp_code", "CompileTraitData.R"))

# Read in tidal data
tidal_1 <- read_csv(here("supp_data", "tidal_1.csv"))
tidal_2 <- read_csv(here("supp_data", "tidal_2.csv"))
tidal_3 <- read_csv(here("supp_data", "tidal_3.csv"))
tidal_4 <- read_csv(here("supp_data", "tidal_4.csv"))

# Bind all together
tidal_all <- rbind(tidal_1, tidal_2, tidal_3, tidal_4)

# Change column names
tidal_all %>% 
  mutate(date = Date,
         time = `Time (GMT)`,
         elevation = `Verified (m)`,
         index = 1:nrow(tidal_all)) %>% 
  select(index, date, time, elevation) -> tidal_all

# Get average elevation of each level
bg_full %>% 
  group_by(level) %>% 
  summarize(mean_elev = mean(elevation)) -> elev_tab

# Calculate percent flooded for each level at average elevation
perc_elev_store <- NULL

for(i in 1:nrow(elev_tab)){
  perc_elev_store[i] <- length(which(tidal_all$elevation > elev_tab$mean_elev[i])) / nrow(tidal_all)
}

perc_elev_label <- paste(round(perc_elev_store*100, 1), "%", sep = "")

# Create vertical lines for each week
months_index <- tibble(x = c(1, 4561, 12001, 19441, 26641, 29041))

# Create plot of tides over time
tidal_all %>% 
  ggplot(aes(x = index, y = elevation)) +
  geom_line(color = "gray37") +
  # Add horizontal lines for each level
  geom_hline(data = elev_tab, aes(yintercept = mean_elev),
             color = c("#d0d1e6", "#a6bddb", "#74a9cf", "#3690c0", "#0570b0", "#034e7b"),
             size = 0.8) +
  ylab("tidal elevation (m NAVD88)") +
  xlab("") +
  geom_vline(data = months_index, aes(xintercept = x), color = "gray57", linetype = "dashed") +
  scale_x_continuous(breaks = months_index$x,
                     labels = c("Jun 12", "Jul 1", "Aug 1", "Sep 1", "Oct 1", "Oct 11")) +
  annotate("text", x = rep(30500, 6), y = elev_tab$mean_elev + 0.05, label = perc_elev_label,
           fontface = "bold", color = c("#d0d1e6", "#a6bddb", "#74a9cf", "#3690c0", "#0570b0", "#034e7b")) +
  theme(axis.text.x = element_text(size=12),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(size = 12)) -> tidal_plot

# Save plot
ggsave(here("figs", "Fig2_tidal_chart.png"), tidal_plot, height = 2.8, width = 9, units = "in")

# Calculate percent flooded for a variety of elevations NADV88 in-text
perc_elev_store_text <- NULL

for(i in 1:2){
  perc_elev_store_text[i] <- length(which(tidal_all$elevation > c(0.1, 0.2)[i])) / nrow(tidal_all)
}

perc_elev_store_text
# 0.8364413 0.6446721
