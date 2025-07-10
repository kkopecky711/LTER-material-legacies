#### SONGS kelp forest monitoring data ####
# Analysis of juvenile kelp recruits ~ cover of dead macroocystis holdfasts

# Packages
library(tidyverse)
library(janitor)
library(glmmTMB)
library(ggeffects)
library(DHARMa)
library(broom.mixed)
library(performance)

# Working directory
setwd('/Users/kako4300/Library/CloudStorage/OneDrive-UCB-O365/Projects/LTER material legacy synthesis')

# Data
songs <- read_csv("Datasets/SONGS/dmaho_cover_mapy_density_all_reef_2000_2023_masp_2009_2023.csv") %>% 
  clean_names()

## Recruits < 10cm
# Remove quadrats with 0% cover oh dead holdfasts
songs.no_zeroes <- songs %>% 
  filter(dmaho_percent_cover > 0,
         mapy_count_per_unit_area > 0)

# Exploratory graphs
ggplot(songs.no_zeroes, aes(x = dmaho_percent_cover, y = mapy_count_per_unit_area)) +
  geom_jitter() +
  geom_smooth(method = "lm") +
  theme_classic()

## Analysis
# Fit Poisson model
kelp_nbinom.glmm <- glmmTMB(
  mapy_count_per_unit_area ~ dmaho_percent_cover + 
    (1 | reef_code/transect_option_code) + (1 | year),
  data = songs.no_zeroes,
  family = nbinom2(link = "log") # negative binomial because response is count data with overdispersion
)
summary(kelp_nbinom.glmm)
tidy(kelp_nbinom.glmm, effects = "fixed", conf.int = TRUE, conf.level = 0.95)

# Get model predictions
preds <- ggpredict(kelp_nbinom.glmm, terms = "dmaho_percent_cover")

# Plot predictions + raw data
ggplot() +
  geom_jitter(data = songs.no_zeroes, 
             aes(x = dmaho_percent_cover, y = mapy_count_per_unit_area), 
             alpha = 0.5, 
             shape = 16) +
  geom_line(data = preds, 
            aes(x = x, y = predicted), 
            linewidth = 1.2) +
  geom_ribbon(data = preds, 
              aes(x = x, ymin = conf.low, ymax = conf.high), 
              alpha = 0.3) +
  labs(x = "% cover of dead holdfasts",
       y = Juvenile~kelp~density~(count/m^2)) +
  theme_classic(base_size = 14)


# Z-score model
# Create z-score versions of response and predictor
songs_z <- songs.no_zeroes %>%
  mutate(
    mapy_count_z = scale(mapy_count_per_unit_area)[, 1],
    dmaho_cover_z = scale(dmaho_percent_cover)[, 1]
  )

kelp_z_model <- glmmTMB(
  mapy_count_z ~ dmaho_cover_z + (1 | reef_code/transect_option_code/quadrat_option_code) + (1 | year),
  data = songs_z,
  family = gaussian(link = "identity")  # z-scored response → normal distribution
)

summary(kelp_z_model)

# Get predicted relationship from z-score model
preds_z <- ggpredict(kelp_z_model, terms = "dmaho_cover_z")

# Plot predictions with raw values as background
ggplot() +
  geom_jitter(data = songs_z, 
              aes(x = dmaho_cover_z, y = mapy_count_z), 
              alpha = 0.4) +
  geom_line(data = preds_z, 
            aes(x = x, y = predicted), 
            linewidth = 1.2) +
  geom_ribbon(data = preds_z, 
              aes(x = x, ymin = conf.low, ymax = conf.high), 
              alpha = 0.3) +
  labs(
    x = "Dead Holdfast Cover (z-score)",
    y = "Juvenile Kelp Density (z-score)"
  ) +
  theme_classic(base_size = 14)



## Effect size
kelp_effect <- tidy(kelp_z_model, effects = "fixed", conf.int = TRUE) %>%
  mutate(
    effect = "fixed",
    component = "cond",
    statistic = estimate / std.error
  ) %>%
  select(effect, component, term, estimate, std.error, statistic, p.value, conf.low, conf.high)%>%
  filter(term == "dmaho_cover_z")

write_csv(kelp_effect, "Datasets/Effect sizes/songs.effect_size.csv")
