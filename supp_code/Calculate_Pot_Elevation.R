# This script uses both RTK-GPS and hand measurements to get actual heights of
# pots within frames. We can then use this data as a covariate in models rather
# than levels. This should overall reduce frame effects as well.

# All measurements taken by Helena Kleiner in Spring 2020

# Load libraries
library(tidyverse); library(here)

# Set plotting theme
theme_set(theme_classic())

# Read in data
# First data set is RTK-GPS measures of the top rack
top_level <- read_csv(here("supp_data", "organ_rack_height.csv"))
# Second data is the differences between racks (e.g. distance between top level
# and second level from the top)
level_diffs <- read_csv(here("supp_data", "shelf_heights_measured.csv"))

# Start with top_level and tidy data
top_level %>% 
  dplyr::select(id = Point.ID,
         elev_1 = Elevation_1,
         elev_2 = Elevation_2,
         elev_3 = Elevation_3) %>% 
  # Add 25 cm to all elevation values (screws for all pots but G6 were 30 cm
  # down and peat was 5 cm down from the top of the pot). Need to add 9 cm to
  # the measurement of G6 because screws had to be placed lower because we
  # couldn't move that rack
  mutate(elev_1 = ifelse(id == "G6", elev_1 + 0.09 + 0.25, elev_1 + 0.25),
         elev_2 = ifelse(id == "G6", elev_2 + 0.09 + 0.25, elev_2 + 0.25),
         elev_3 = ifelse(id == "G6", elev_3 + 0.09 + 0.25, elev_3 + 0.25)) %>% 
  group_by(id) %>% 
  mutate(mean_elev = mean(c(elev_1, elev_2, elev_3)),
         sd_elev = sd(c(elev_1, elev_2, elev_3)),
         se_elev = sd_elev / sqrt(3)) -> top_level_summary

# Plot top levels with standard errors
top_level_summary %>% 
  ggplot(aes(x = id, y = mean_elev)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_elev - se_elev, ymax = mean_elev + se_elev), width = 0.5) +
  ylab("mean elevation (m NAVD88)")

# We could use the information about variability as an errors in variables
# problem to be complete about how variation in elevation drives differences in
# responses, but I think for now we can just take the mean elevation for each
# rack.

top_level_summary %>% 
  dplyr::select(id, mean_elev) -> top_level_mean

range(top_level_mean$mean_elev)
# We have up to ~5cm difference across frames so it's probably good we are
# accounting for this!

# Now format the level differences
level_diffs %>% 
  # Take out measurements that had "DISREGARD" in the notes
  filter(!grepl("DISREGARD",level_diffs$notes)) %>% 
  group_by(site_frame) %>% 
  # There should be 5 for each frame (1-2, 2-3, 3-4, 4-5, 5-6)
  dplyr::count() %>% 
  filter(n != 5)

# Check what's up with gcrew_3 and gcrew_4
level_diffs %>% 
  filter(!grepl("DISREGARD",level_diffs$notes)) %>% 
  filter(site_frame %in% c("gcrew_3", "gcrew_4")) %>% 
  pull(notes)
# Ok, it looks like we should average measurements taken with rows that have the
# same level based on the notes

level_diffs %>% 
  filter(!grepl("DISREGARD",level_diffs$notes)) %>% 
  group_by(site_frame, level) %>% 
  # Convert to meters
  summarize(measurement_R = mean(measurement_R) * 0.01,
            measurement_C = mean(measurement_C * 0.01),
            measurement_L = mean(measurement_L) * 0.01) %>% 
  group_by(site_frame, level) %>% 
  # Now take the global mean across right, center and left
  summarize(mean_diff = mean(c(measurement_R, measurement_C, measurement_L))) %>% 
  ungroup() %>% 
  # Rename level to be numeric 2-6, where 2 represents the difference between 1 to 2
  mutate(level = ifelse(level == "1 to 2", 2,
                        ifelse(level == "2 to 3", 3,
                               ifelse(level == "3 to 4", 4,
                                      ifelse(level == "4 to 5", 5, 6))))) -> level_diffs_sum

# Create storage for loop that will go through and subtract elevations for each
# step
storage <- matrix(NA, nrow = 14, ncol = 6)

# Make sure both data frames are sorted before we loop through
top_level_summary <- top_level_summary[order(top_level_summary$id),]
level_diffs_sum <- level_diffs_sum[order(level_diffs_sum$site_frame, level_diffs_sum$level),]

# Set vector of frame names for level_diffs
frame_names <- unique(level_diffs_sum$site_frame)

# Loop through to get elevations for each level for each frame
for (i in 1:14){
  temp_top <- top_level_summary$mean_elev[i]
  temp_rest <- level_diffs_sum %>% filter(site_frame == frame_names[i]) %>% pull(mean_diff)
  storage[i,1] <- temp_top
  for (j in 2:6){
    storage[i,j] <- storage[i,j-1] - temp_rest[j-1]
  }
}

# Add back ids and levels and put in a tibble
all_elevations <- data.frame(storage)
colnames(all_elevations) <- paste("level", 1:6, sep = "_")
tibble(all_elevations) %>% 
  mutate(site_frame = frame_names) %>% 
  gather(key = level, value = elevation, level_1:level_6) %>% 
  mutate(level = substr(level, nchar(level), nchar(level))) -> elevation_data

# Remove all previous data frames and keep the last one
rm(list=setdiff(ls(), "elevation_data"))
