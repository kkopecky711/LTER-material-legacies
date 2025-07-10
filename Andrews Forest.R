#### H.J. Andrews Forest ####
# Analysis of tree growth, mortality, and ingrowth as a function of standing and downed dead wood

# Working directory
setwd('/Users/kako4300/Library/CloudStorage/OneDrive-UCB-O365/Projects/LTER material legacy synthesis')

# Packages
library(tidyverse)
library(janitor)
library(glmmTMB)
library(ggeffects)
library(DHARMa)
library(broom.mixed)

# Dead tree data
hja_dead <- read_csv("Datasets/Andrews Forest/OHJA_downed wood summary.csv") %>% 
  clean_names() %>% 
  # filter(stand %in% c("RS01", "RS02", "RS03", "RS04", "RS05", "RS07", "RS08", "RS10", "RS12", "RS14", "RS15", "RS16", "RS17", "RS18", "RS20", "RS22", "RS23", "RS24", "RS26", "RS27", "RS29", "RS31", "RS33")) %>% 
  mutate(obs_id = paste(stand, plot, year, sep = "_"))

# Live tree data
hja_live <- read_csv("Datasets/Andrews Forest/PSP_Plot_Change_20250610.csv") %>% 
  clean_names() %>% 
  filter(d_year != "NA")

# Group species within live tree data
hja_live.spp_grouped <- hja_live %>%
  group_by(stand, plot, year) %>% 
  summarize(mortality = sum(mort_ba_spp),
            ingrowth = sum(ingrowth_ba_spp),
            growth = sum(growth_ba_spp)) %>% 
  mutate(plot = as.character(plot),
         obs_id = paste(stand, plot, year, sep = "_"))

# Step 1: Rename columns for clarity
dead <- hja_dead %>%
  rename(year_dead = year) %>% 
  select(-obs_id)

live <- hja_live.spp_grouped %>%
  rename(year_live = year)%>% 
  select(-obs_id)

# Step 2: Join on stand and plot
temporal_joined <- inner_join(live, dead, by = c("stand", "plot")) %>%
  # Step 3: Keep only cases where live year is after dead year
  filter(year_live > year_dead) %>%
  # Step 4: Calculate time since dead measurement
  mutate(years_since_dead = year_live - year_dead)

# Filter for only first two records after dead measurement
temporal_joined.first_two <- temporal_joined %>%
  group_by(stand, plot, year_dead) %>%
  arrange(years_since_dead, .by_group = TRUE) %>%
  slice_head(n = 2) %>%
  ungroup()

# Exploratory graph
ggplot(temporal_joined, aes(x = total_volume, y = growth)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic()


## 
hja.summary <- temporal_joined %>% 
  group_by(stand, plot, years_since_dead) %>% 
  summarize(mean_growth = mean(growth),
            mean_dead = mean(total_volume))

