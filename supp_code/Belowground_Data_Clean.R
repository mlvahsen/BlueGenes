# Script for cleaning up belowground biomass data

library(tidyverse)

bg_check <- read_csv(here("supp_data","BelowgroundBiomass_Check.csv"))
bg_old <- read_csv(here("supp_data","BelowgroundBiomass_OLD.csv"))
ic <- read_csv(here("supp_data","HarvestPot_ID.csv"))

# create not in function
"%notin%" <- Negate("%in%")

# identify unique pots for each
bg_old_pots <- unique(bg_old$pot_no)
bg_check_pots <- unique(bg_check$pot_no)

# which pots are in bg_old but are not in bg_check
bg_old_pots[which(bg_old_pots %notin% bg_check_pots)] # none
# which pots are in bg_check but are not in bg_old
bg_check_pots[which(bg_check_pots %notin% bg_old_pots)]# pot 205

bg_check %>% 
  filter(pot_no == 205) 
# something seems wrong here... 
ic %>% 
  filter(`site_ frame` == "fresh_7") %>% 
  print(nrow = Inf)# no competition

# create unique layer x pot_no index to see what is missing in each
bg_old_pots_layer <- unique(paste(bg_old$pot_no, bg_old$segment_top, sep = "_"))
bg_check_pots_layer <- unique(paste(bg_check$pot_no, bg_check$segment_top, sep = "_"))

# which pot_layers are in bg_old but are not in bg_check
bg_old_pots_layer[which(bg_old_pots_layer %notin% bg_check_pots_layer)] 
# "165_10"  "1700_10" "1659_0"  "340_10"  "354_10"  "214_10"  "311_10"  "394_20" 

# which pot_layers are in bg_check but are not in bg_old
bg_check_pots_layer[which(bg_check_pots_layer %notin% bg_old_pots_layer)]# pot 205
# "122_10"  "198_10"  "205_0"   "205_10"  "205_20"  "212_20"  "217_10"  "219_10"  "219_20"  "1624_0"  "1624_10"
# "1642_0"  "1642_10" "1652_10" "1660_20" "1687_0"  "1697_10" "1698_10" "1882_20" "1884_10" "1884_20" "1885_10"
# "1886_0"  "1886_10" "1890_10"

# Start with check and then go through by hand
bg_check%>% 
  arrange(pot_no, segment_top, species) %>% 
  filter(pot_no %in% 314) 

# fresh_1
bg_check %>% 
  filter(pot_no %in% 1:48) %>% 
  mutate(species = case_when(species == "mixed fine roots" ~ "mixed roots",
                             species == "fine roots" ~ "mixed roots",
                             species == "scam" ~ "scam",
                             species == "sppa" ~ "sppa",
                             species == "no visible sppa" ~ "scam",
                             T ~ "scam")) -> bg_fresh_1

bg_check %>% 
  filter(pot_no %in% 49:96) %>% 
  mutate(species = case_when(species == "mixed fine roots" ~ "mixed roots",
                             species == "mixed rotts" ~ "mixed roots",
                             species == "scam" ~ "scam",
                             species == "sppa" ~ "sppa",
                             species == "SPPA" ~ "sppa",
                             species == "mixed roots" ~ "mixed roots",
                             T ~ "scam")) -> bg_fresh_2

bg_check %>% 
  filter(pot_no %in% 97:144) %>% 
  mutate(species = case_when(species == "mixed fine roots" ~ "mixed roots",
                             species == "mixed roots" ~ "mixed roots",
                             species == "scam" ~ "scam",
                             species == "sppa" ~ "sppa",
                             species == "fine roots" ~ "mixed roots",
                             notes == "mixed scam/sppa" ~ "mixed roots",
                             notes == "sppa & scam fine roots tangled together in these groups" ~ "mixed roots",
                             T ~ "scam")) -> bg_fresh_3

bg_check %>% 
  filter(pot_no %in% 145:192) %>% 
  mutate(species = case_when(species == "mixed fine roots" ~ "mixed roots",
                             species == "fine roots" ~ "mixed roots",
                             species == "mixed roots" ~ "mixed roots",
                             species == "mixed" ~ "mixed roots",
                             species == "scam" ~ "scam",
                             species == "sppa" ~ "sppa",
                             notes == "fine roots" ~ "mixed roots",
                             notes == "c3" ~ "scam",
                             notes == "c4" ~ "sppa",
                             T ~ "scam")) %>% 
  filter(pot_no != 154 | weight != 1.152) %>% 
  add_row(date = "11/05/19", pot_no = 165, segment_top = "10",
          segment_bottom = "20", weight = "1.984", species = NA, notes = NA) -> bg_fresh_4

bg_check %>% 
  filter(pot_no %in% 193:240) %>% 
  mutate(species = case_when(species == "mixed fine roots" ~ "mixed roots",
                             species == "mixed roots" ~ "mixed roots",
                             species == "scam" ~ "scam",
                             species == "sppa" ~ "sppa",
                             notes == "quite a lot of sppa/fine roots in scam rhizome" ~ "mixed roots",
                             T ~ "scam")) %>% 
  add_row(date = "01/16/20", pot_no = 214, segment_top = "10",
          segment_bottom = "20", weight = "0.691", species = NA, notes = NA) -> bg_fresh_5

bg_check %>% 
  filter(pot_no %in% 241:288) %>% 
  mutate(species = case_when(species == "mixed fine roots" ~ "mixed roots",
                             species == "mixed roots" ~ "mixed roots",
                             species == "scam" ~ "scam",
                             species == "sppa" ~ "sppa",
                             T ~ "scam")) -> bg_fresh_6

bg_check %>% 
  filter(pot_no %in% 289:336) %>% 
  mutate(species = case_when(species == "mixed fine roots" ~ "mixed roots",
                             species == "mixed roots" ~ "mixed roots",
                             species == "scam" ~ "scam",
                             species == "sppa" ~ "sppa",
                             T ~ "scam")) %>% 
  filter(pot_no != 304 | segment_bottom != "end") %>% 
  add_row(date = "01/15/20", pot_no = 311, segment_top = "10", segment_bottom = "20",
          weight = "2.281", species = NA, notes = NA) %>% 
  mutate(segment_bottom = ifelse(segment_top == "0" & weight == "0.084", NA, segment_bottom)) %>% 
  mutate(weight = ifelse(weight == "7/457", "7.457", weight)) -> bg_fresh_7

bg_tibble_fresh <- rbind(bg_fresh_1,
                         bg_fresh_2,
                         bg_fresh_3,
                         bg_fresh_4,
                         bg_fresh_5,
                         bg_fresh_6,
                         bg_fresh_7)

bg_tibble_fresh %>% 
  mutate(weight = as.numeric(weight),
         segment_top = as.numeric(segment_top),
         segment_bottom = as.numeric(segment_bottom)) -> bg_tibble_fresh


bg_tibble_fresh %>%
  group_by(pot_no, segment_top) %>%
  summarize(total_biomass = sum(weight)) %>%
  ggplot(aes(x = segment_top, y = total_biomass, group = pot_no)) + 
  geom_line(alpha = 0.2)

bg_tibble_fresh %>% 
  group_by(pot_no) %>% 
  summarize(total_bgb = sum(weight)) -> bg_tibble_fresh_sum

# Now GCREW
ic %>% 
  filter(`site_ frame` == "gcrew_7") %>% 
  print(nrow = Inf)# no competition

bg_check %>% 
  filter(pot_no %in% 337:384) %>% 
  mutate(species = case_when(species == "mixed fine roots" ~ "mixed roots",
                             species == "mixed fine root" ~ "mixed roots",
                             species == "mixed roots" ~ "mixed roots",
                             species == "scam" ~ "scam",
                             species == "sppa" ~ "sppa",
                             notes == "most sppa fine roots, too mixed with scam to get" ~ "mixed roots",
                             T ~ "scam")) %>% 
  mutate(weight = ifelse(pot_no == 339 & segment_top == "20", 2.906, weight)) %>% 
  add_row(date = "12/04/19", pot_no = 340, segment_top = "10", segment_bottom = "20",
          weight = "1.875", species = NA, notes = NA) %>% 
  mutate(pot_no = ifelse(pot_no == 348 & date == "12/17/19", 358, pot_no)) %>% 
  mutate(weight = ifelse(pot_no == 337 & segment_top == "0", 2.763, weight)) -> bg_gcrew_1

bg_check %>% 
  filter(pot_no %in% 385:1632 & pot_no != 810) %>% 
  mutate(species = case_when(species == "mixed fine roots" ~ "mixed roots",
                             species == "mixed roots" ~ "mixed roots",
                             species == "scam" ~ "scam",
                             species == "sppa" ~ "sppa",
                             T ~ "scam")) %>% 
  add_row(date = "11/25/19", pot_no = 394, segment_top = "20",
          segment_bottom = "end", weight = "2.092", species = NA, notes = NA) -> bg_gcrew_2

bg_check %>% 
  filter(pot_no %in% 1633:1680) %>% 
  mutate(species = case_when(species == "mixed fine roots" ~ "mixed roots",
                             species == "scam" ~ "scam",
                             species == "sppa" ~ "sppa",
                             T ~ "scam")) %>% 
  filter(notes != "stringy material - NOT ROOTS" | is.na(notes)) %>% 
  add_row(date = "12/03/19", pot_no = 1659, segment_top = "0",
          segment_bottom = "10", weight = "3.199", species = NA, notes = NA) -> bg_gcrew_3

bg_check %>% 
  filter(pot_no %in% c(1681:1828, 810)) %>% 
  mutate(species = case_when(species == "mixed fine roots" ~ "mixed roots",
                             species == "mixed roots" ~ "mixed roots",
                             species == "fine roots" ~ "mixed roots",
                             species == "scam" ~ "scam",
                             species == "sppa" ~ "sppa",
                             T ~ "scam")) %>% 
  mutate(segment_top = ifelse(segment_top == "1", "0", segment_top)) %>% 
  add_row(date = "11/07/19", pot_no = 1700, segment_top = "10", 
          segment_bottom = "20", weight = "2.600", species = NA, notes = "NA") %>% 
  mutate(pot_no = ifelse(pot_no == 810, 1810, pot_no))-> bg_gcrew_4

bg_check %>% 
  filter(pot_no %in% 1829:1876) %>% 
  mutate(species = case_when(species == "mixed fine roots" ~ "mixed roots",
                             species == "mixed roots" ~ "mixed roots",
                             species == "sppa?" ~ "sppa",
                             species == "scam" ~ "scam",
                             species == "sppa" ~ "sppa",
                             T ~ "scam")) %>% 
  filter(notes != "fiberous material--roots?" | is.na(notes)) %>% 
  mutate(pot_no = ifelse(pot_no == 1849 & date == "02/05/20", 1847, pot_no))-> bg_gcrew_5

bg_check %>% 
  filter(pot_no %in% 1877:1952) %>% 
  mutate(species = case_when(species == "mixed fine roots" ~ "mixed roots",
                             species == "mixed roots" ~ "mixed roots",
                             species == "mixed" ~ "mixed roots",
                             species == "mixed scam/sppa" ~ "mixed roots",
                             species == "fine roots" ~ "mixed roots",
                             species == "scam" ~ "scam",
                             species == "sppa" ~ "sppa",
                             T ~ "scam")) %>% 
  mutate(pot_no = ifelse(pot_no == 1880 & date == "02/12/20", 1888, pot_no)) %>% 
  mutate(pot_no = ifelse(pot_no == 1892 & date == "02/04/20", 1894, pot_no)) %>% 
  mutate(weight = ifelse(pot_no == 1899 & segment_top == "20", 1.560, weight)) -> bg_gcrew_6

bg_check %>% 
  filter(pot_no %in% c(1953:2000, 1679)) %>% 
  mutate(species = case_when(species == "mixed fine roots" ~ "mixed roots",
                             species == "mixed roots" ~ "mixed roots",
                             species == "scam" ~ "scam",
                             species == "sppa" ~ "sppa",
                             T ~ "scam")) %>% 
  filter(notes != "empty" | is.na(notes)) %>% 
  mutate(pot_no = ifelse(pot_no == 1965 & date == "01/07/20", 1956, pot_no)) %>% 
  mutate(pot_no = ifelse(pot_no == 1679 & date == "12/04/19", 1979, pot_no))-> bg_gcrew_7

bg_tibble_gcrew <- rbind(bg_gcrew_1,
                         bg_gcrew_2,
                         bg_gcrew_3,
                         bg_gcrew_4,
                         bg_gcrew_5,
                         bg_gcrew_6,
                         bg_gcrew_7)


bg_tibble_gcrew %>% 
  mutate(weight = as.numeric(weight),
         segment_top = as.numeric(segment_top),
         segment_bottom = as.numeric(segment_bottom)) -> bg_tibble_gcrew


bg_tibble_gcrew %>%
  group_by(pot_no, segment_top) %>%
  summarize(total_biomass = sum(weight)) %>%
  ggplot(aes(x = segment_top, y = total_biomass, group = pot_no)) + 
  geom_line(alpha = 0.2)

bg_tibble_fresh %>% 
  group_by(pot_no) %>% 
  summarize(total_bgb = sum(weight)) -> bg_tibble_fresh_sum

bg_all <- rbind(bg_tibble_fresh, bg_tibble_gcrew)

# Pots that still need figuring out: 78, 205, 354, 1700, 1885, 1984

# 78: missing 10-20 and 20-end; no notes indicating why (pg. 14); AD
# 205: additional 0-66.5 bag for 4 total bags; no notes indicating why (pg. 10); AD
# 354: missing segment_top and segment_bottom for one later and missing another; no notes indicated why (pg. 5); NJI
# 1700: missing 0-10 layer; note on recording page says "on outside" (pg. 2); HSK
# 1885: missing 20-end layer; no notes indicating why (pg. 10); AMD
# 1984: duplicated 0-10 bags and no 20-end bag; no notes indicating why (pg. 5); NJI

# 11 June 2021 - update from checking dried biomass boxes in SB
bg_all %>% 
  mutate(species = case_when(pot_no == 215 ~ "sppa",
                             T ~ species))

# MLV fixing from May 2021 trip to SB 
bg_all %>% 
  mutate(species = case_when(pot_no == 215 ~ "sppa",
                             T ~ species)) %>% 
  mutate(pot_no = case_when(pot_no == 343 & weight == 0.146 ~ 345,
                            T ~ pot_no)) %>% 
  mutate(species = case_when(pot_no == 28 ~ "sppa",
                             T ~ species)) %>% 
  mutate(species = case_when(pot_no == 160 ~ "sppa",
                             T ~ species)) %>% 
  mutate(species = case_when(pot_no == 318 & weight %in% c(0.127, 0.323) ~ "sppa",
                             T ~ species)) %>%
  mutate(segment_top = case_when(pot_no == 399 & weight > 3 ~ 0,
                             T ~ segment_top)) %>%
  mutate(segment_bottom = case_when(pot_no == 399 & weight > 3 ~ 10,
                                 T ~ segment_bottom)) %>%
  mutate(species = case_when(pot_no == 399 & weight %in% c(0.961, 0.207) ~ "mixed roots",
                                    T ~ species)) %>%
  mutate(species = case_when(pot_no == 1847 & segment_top %in% c(10, 20) ~ "mixed roots",
                             T ~ species)) %>%
  mutate(species = case_when(pot_no == 1967 & segment_top %in% c(10, 20) ~ "mixed roots",
                             T ~ species)) %>%
  filter(pot_no != 1679) %>% 
  add_row(date = "07/01/20", pot_no = 78, segment_top = 10, segment_bottom = 20,
          weight = 0.885, species = "mixed roots", notes = NA) %>% 
  add_row(date = "07/01/20", pot_no = 78, segment_top = 20, segment_bottom = 65,
          weight = 0.768, species = "mixed roots", notes = NA) %>% 
  add_row(date = "12/16/19", pot_no = 354, segment_top = 10, segment_bottom = 20,
          weight = 1.152, species = "scam", notes = NA) %>% 
  mutate(segment_top = case_when(pot_no == 354 & weight < 1 ~ 20,
                                 T ~ segment_top)) %>%
  mutate(segment_bottom = case_when(pot_no == 354 & weight < 1 ~ 67,
                                    T ~ segment_bottom)) %>%
  filter(pot_no != 1866) %>% 
  mutate(species = case_when(pot_no %in% c(121, 161, 355, 401, 406, 1664, 1687,
                                           1816, 1894, 1966, 1982) ~ "scam",
                                           T ~ species)) -> bg_all

bg_clean <- bg_all %>% 
  filter(pot_no %notin% c(205, 1700, 1885))
