#### Moorea Coral Reef LTER ####
# Analysis of change in live coral cover ~ dead coral cover from dead coral removal experiment

# Packages
library(tidyverse)
library(janitor)
library(glmmTMB)
library(ggeffects)
library(DHARMa)
library(broom.mixed)

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
coral_glmm.raw <- glmmTMB(
  live_coral_change_pct ~ dead_coral_start + (1 | treatment/plot) + (1 | year),
  data = coral_change_filtered,
  family = gaussian()
)
summary(coral_glmm.raw)

# Diagnostics
res <- simulateResiduals(coral_glmm.raw)
plot(res) # Tests are non-significant

# Extract raw effect size
mcr_effect.raw <- tidy(coral_glmm.raw, effects = "fixed", conf.int = TRUE, conf.level = 0.95) %>% 
  filter(term == "dead_coral_start")

write_csv(mcr_effect.raw, "Datasets/Effect sizes/Raw/mcr_effect.raw.csv")
  
## Visualization
# Extract predictions from model
preds <- ggpredict(coral_glmm.raw, terms = "dead_coral_start")

# Plot model predictions + raw values
ggplot() +
  geom_point(data = coral_change_filtered, 
             aes(x = dead_coral_start, 
                 y = live_coral_change_pct),
             color = "#20618D", 
             alpha = 0.6) +
  geom_ribbon(data = preds, 
              aes(x = x, 
                  ymin = conf.low, 
                  ymax = conf.high), 
              alpha = 0.2,
              fill = "#20618D") +
  geom_line(data = preds, 
            aes(x = x, 
                y = predicted), 
            linewidth = 0.75) +
  labs(
    x = Dead~coral~cover~(m^2),
    y = "Change in live coral (% surf. area/yr)") +
  theme_classic(base_size = 14)

## Standardized model
# Scale both variables to z-scores
coral_change_filtered <- coral_change_filtered %>%
  mutate(
    live_coral_change_pct.z = scale(live_coral_change_pct)[, 1],
    dead_coral_start.z = scale(dead_coral_start)[, 1]
  )

coral_glmm.z <- glmmTMB(
  live_coral_change_pct.z ~ dead_coral_start.z + (1 | treatment/plot) + (1 | year),
  data = coral_change_filtered,
  family = gaussian()
)
summary(coral_glmm.z)

# Diagnostics
res.z <- simulateResiduals(coral_glmm.z)
plot(res.z) # Tests are non-significant

## Visualization
# Extract predictions from model
preds.z <- ggpredict(coral_glmm.z, terms = "dead_coral_start.z")

# Plot model predictions + raw values
ggplot() +
  geom_point(data = coral_change_filtered, 
             aes(x = dead_coral_start.z, 
                 y = live_coral_change_pct.z),
             color = "darkgrey", 
             alpha = 0.6) +
  geom_ribbon(data = preds.z, 
              aes(x = x, 
                  ymin = conf.low, 
                  ymax = conf.high), 
              alpha = 0.2) +
  geom_line(data = preds.z, 
            aes(x = x, 
                y = predicted), 
            linewidth = 1.2) +
  scale_y_continuous(limits = c(-2, 3.4)) +
  labs(x = "Dead coral cover (Z-score)",
    y = "Change in live coral (Z-score)") +
  theme_classic(base_size = 14)

## Extract effect size
mcr_effect.z <- tidy(coral_glmm.z, effects = "fixed", conf.int = TRUE) |>
  filter(term == "dead_coral_start_z")

write_csv(mcr_effect.z, "Datasets/Effect sizes/Standardized/mcr_effect.z.csv")
