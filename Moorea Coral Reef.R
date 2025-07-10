#### Moorea Coral Reef LTER ####
# Analysis of change in live coral cover ~ dead coral cover from dead coral removal experiment

# Packages
library(tidyverse)
library(janitor)
library(glmmTMB)
library(ggeffects)
library(DHARMa)

# Working directory
setwd('/Users/kako4300/Library/CloudStorage/OneDrive-UCB-O365/Projects/LTER material legacy synthesis')

# Read in datafile for live and dead coral annotations
regions_master <- read_csv("Datasets/Moorea Coral Reef/Regions_master.csv") %>% 
  clean_names

# Create dataframe of live coral cover
live_coral <- regions_master %>% 
  filter(tag_lab_class_name == "Pocillopora" | tag_lab_class_name == "Acropora") %>% 
  group_by(plot, treatment, tag_lab_date) %>% 
  summarize(live_coral.m_sq = sum(tag_lab_surf_area*0.0001)) %>% 
  mutate(year = case_when(tag_lab_date == "8/1/19" ~ 2019,
                                tag_lab_date == "8/1/20" ~ 2020,
                                tag_lab_date == "8/1/21" ~ 2021,
                                tag_lab_date == "8/1/22" ~ 2022,
                                tag_lab_date == "8/1/23" ~ 2023),
         years_since_start = case_when(tag_lab_date == "8/1/19" ~ 0,
                                       tag_lab_date == "8/1/20" ~ 1,
                                       tag_lab_date == "8/1/21" ~ 2,
                                       tag_lab_date == "8/1/22" ~ 3,
                                       tag_lab_date == "8/1/23" ~ 4),
         year = as.factor(year),
         obs_id = paste(year, plot, treatment, sep = "_")) %>% 
  select(-tag_lab_date)

# Create dataframe of dead coral cover
dead_coral <- regions_master %>%
  filter(tag_lab_class_name == "Dead coral") %>% 
  group_by(plot, treatment, tag_lab_date) %>% 
  summarize(dead_coral.m_sq = sum(tag_lab_surf_area*0.0001)) %>% 
  mutate(year = case_when(tag_lab_date == "8/1/19" ~ 2019,
                          tag_lab_date == "8/1/20" ~ 2020,
                          tag_lab_date == "8/1/21" ~ 2021,
                          tag_lab_date == "8/1/22" ~ 2022,
                          tag_lab_date == "8/1/23" ~ 2023),
         years_since_start = case_when(tag_lab_date == "8/1/19" ~ 0,
                                       tag_lab_date == "8/1/20" ~ 1,
                                       tag_lab_date == "8/1/21" ~ 2,
                                       tag_lab_date == "8/1/22" ~ 3,
                                       tag_lab_date == "8/1/23" ~ 4),
         year = as.factor(year),
         obs_id = paste(year, plot, treatment, sep = "_")) %>% 
  ungroup() %>% 
  select(c(obs_id, dead_coral.m_sq))

# Merge live and dead dataframes
live_dead <- live_coral %>% 
  left_join(dead_coral, 
            by = "obs_id")


# Change in coral cover
coral_change <- live_dead %>%
  arrange(plot, year) %>%
  group_by(plot) %>%
  mutate(
    live_coral_start = lag(live_coral.m_sq),
    live_coral_change_pct = (live_coral.m_sq - live_coral_start) / live_coral_start * 100,
    dead_coral_start = lag(dead_coral.m_sq)
  ) %>%
  ungroup()

coral_change_filtered <- coral_change %>%
  filter(!is.na(live_coral_change_pct) & !is.na(dead_coral_start))

# Fit GLMM
coral_glmm <- glmmTMB(
  live_coral_change_pct ~ dead_coral_start + (1 | treatment/plot) + (1 | year),
  data = coral_change_filtered,
  family = gaussian()
)
summary(coral_glmm)
tidy(coral_glmm, effects = "fixed", conf.int = TRUE, conf.level = 0.95)

# Assumptions check
simulateResiduals(coral_glmm, plot = TRUE)

# Visualization
# Extract predictions from model
preds <- ggpredict(coral_glmm, terms = "dead_coral_start")

# Plot model predictions + raw values
ggplot() +
  # Model prediction line with CI ribbon
  geom_ribbon(data = preds, 
              aes(x = x, 
                  ymin = conf.low, 
                  ymax = conf.high), 
              alpha = 0.2) +
  geom_line(data = preds, 
            aes(x = x, 
                y = predicted), 
            size = 1.2) +
  geom_point(data = coral_change_filtered, 
             aes(x = dead_coral_start, 
                 y = live_coral_change_pct),
             size = 2, alpha = 0.6) +
  labs(
    x = Dead~coral~cover~(m^2),
    y = "Change in live coral (% surf. area/yr)") +
  theme_classic(base_size = 13)

## Z-score scaling
coral_scaled <- coral_change_filtered %>%
  mutate(
    live_coral_change_pct_z = scale(live_coral_change_pct)[, 1],
    dead_coral_start_z = scale(dead_coral_start)[, 1]
  )

coral_glmm_z <- glmmTMB(
  live_coral_change_pct_z ~ dead_coral_start_z + (1 | treatment/plot) + (1 | year),
  data = coral_scaled,
  family = gaussian()
)
summary(coral_glmm_z)

## Visualization
# Extract predictions from model
preds_z <- ggpredict(coral_glmm_z, terms = "dead_coral_start_z")

# Plot model predictions + raw values
ggplot() +
  # Model prediction line with CI ribbon
  geom_ribbon(data = preds_z, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  geom_line(data = preds_z, aes(x = x, y = predicted), size = 1.2) +
  geom_point(data = coral_scaled, 
             aes(x = dead_coral_start_z, y = live_coral_change_pct_z),
             size = 2, alpha = 0.8) +
  labs(
    x = "Dead coral cover (Z-score)",
    y = "% Change in Live Coral (Z-score)") +
  theme_classic(base_size = 14)

# Extract effect size
coral_effect <- tidy(coral_glmm_z, effects = "fixed", conf.int = TRUE) |>
  filter(term == "dead_coral_start_z")

print(coral_effect)

write_csv(coral_effect, "Datasets/Effect sizes/mcr.effect_size.csv")
