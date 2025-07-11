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
# Fit neagative binomial GLMM
kelp_nbinom.glmm <- glmmTMB(
  mapy_count_per_unit_area ~ dmaho_percent_cover + 
    (1 | reef_code/transect_option_code) + (1 | year),
  data = songs.no_zeroes,
  family = nbinom2(link = "log") # negative binomial because response is count data with overdispersion
)
summary(kelp_nbinom.glmm)

# Diagnostics
res <- simulateResiduals(kelp_nbinom.glmm)
plot(res) # Tests are reasonable

# Extract effect size and CI, write to .csv file
songs_effect.raw <- tidy(kelp_nbinom.glmm, effects = "fixed", conf.int = TRUE, conf.level = 0.95) %>% 
  filter(term == "dmaho_percent_cover")

write_csv(songs_effect.raw, "Datasets/Effect sizes/Raw/songs_effect.raw.csv")

# Get model predictions
preds <- ggpredict(kelp_nbinom.glmm, terms = "dmaho_percent_cover")

# Plot predictions + raw data
ggplot() +
  geom_jitter(data = songs.no_zeroes, 
             aes(x = dmaho_percent_cover, y = mapy_count_per_unit_area), 
             color = "darkgrey", 
             alpha = 0.6, 
             width = 0.15) +
  geom_line(data = preds, 
            aes(x = x, y = predicted), 
            linewidth = 1.2) +
  geom_ribbon(data = preds, 
              aes(x = x, ymin = conf.low, ymax = conf.high), 
              alpha = 0.3) +
  labs(x = "% cover of dead holdfasts",
       y = Juvenile~kelp~density~(count/m^2/yr)) +
  theme_classic(base_size = 14)


# Z-score model
# Create z-score versions of response and predictor
songs_z <- songs.no_zeroes %>%
  mutate(
    mapy_count_z = scale(mapy_count_per_unit_area)[, 1],
    dmaho_cover_z = scale(dmaho_percent_cover)[, 1]
  )

kelp_glmm.z <- glmmTMB(
  mapy_count_z ~ dmaho_cover_z + (1 | reef_code/transect_option_code) + (1 | year),
  dispformula = ~ dmaho_cover_z,
  data = songs_z,
  family = gaussian()
)
summary(kelp_glmm.z)

# Diagnostics
res.z <- simulateResiduals(kelp_glmm.z)
plot(res.z) # Tests are reasonable

# Get predicted relationship from z-score model
preds_z <- ggpredict(kelp_glmm.z, terms = "dmaho_cover_z")

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
songs_effect.z <- tidy(kelp_glmm.z, effects = "fixed", conf.int = TRUE, conf.level = 0.95) %>% 
  filter(term == "dmaho_cover_z")

write_csv(songs_effect.z, "Datasets/Effect sizes/Standardized/songs_effect.z.csv")

