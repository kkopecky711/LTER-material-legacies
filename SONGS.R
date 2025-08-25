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
# Remove quadrats with 0% cover of dead holdfasts, otherwise data are highly skewed
songs.no_zeroes <- songs %>% 
  filter(dmaho_percent_cover > 0)

# Exploratory graphs
ggplot(songs.no_zeroes, aes(x = dmaho_percent_cover, y = mapy_count_per_unit_area)) +
  geom_jitter() +
  geom_smooth(method = "lm") +
  theme_classic()

## Analysis
# Fit negative binomial GLMM
kelp_glmm.raw <- glmmTMB(
  mapy_count_per_unit_area ~ dmaho_percent_cover + 
    (1 | reef_code/transect_option_code) + (1 | year),
  data = songs.no_zeroes,
  family = nbinom2(link = "log") # negative binomial because response is count data with overdispersion
)
summary(kelp_glmm.raw)

# Diagnostics
res <- simulateResiduals(kelp_glmm.raw)
plot(res) # Tests are reasonable

# compare with model that has zero-inflation component
# Fit negative binomial GLMM
kelp_glmm.raw.zi <- glmmTMB(
  mapy_count_per_unit_area ~ dmaho_percent_cover + 
    (1 | reef_code/transect_option_code) + (1 | year),
  ziformula = ~1,
  data = songs.no_zeroes,
  family = nbinom2(link = "log") # negative binomial because response is count data with overdispersion
)
summary(kelp_glmm.raw.zi)

# Compare models
AIC(kelp_glmm.raw, kelp_glmm.raw.zi) # zero-inflation component does not enhance model fit; keep original model structure

# Extract effect size and CI, write to .csv file
songs_effect.raw <- tidy(kelp_glmm.raw, effects = "fixed", conf.int = TRUE, conf.level = 0.95) %>% 
  filter(term == "dmaho_percent_cover")

write_csv(songs_effect.raw, "Datasets/Effect sizes/Raw/songs_effect.raw.csv")


## Visualization
# Get model predictions
preds <- ggpredict(kelp_glmm.raw, 
                   terms = "dmaho_percent_cover")

# Plot predictions + raw data
ggplot() +
  geom_jitter(data = songs.no_zeroes, 
             aes(x = dmaho_percent_cover, y = mapy_count_per_unit_area), 
             color = "#20618D", 
             alpha = 0.6, 
             width = 0.5) +
  geom_line(data = preds, 
            aes(x = x, y = predicted), 
            linewidth = 0.75) +
  geom_ribbon(data = preds, 
              aes(x = x, ymin = conf.low, ymax = conf.high), 
              alpha = 0.3,
              fill = "#20618D") +
  labs(x = "% cover of dead holdfasts",
       y = Juvenile~kelp~density~(count/m^2/yr)) +
  theme_classic(base_size = 14)


## Standardized model
# Create z-score versions of response and predictor
songs.no_zeroes <- songs.no_zeroes %>%
  mutate(
    mapy_count.z = scale(mapy_count_per_unit_area)[, 1],
    dmaho_cover.z = scale(dmaho_percent_cover)[, 1]
  )

# kelp_glmm.z <- glmmTMB(
#   mapy_count.z ~ dmaho_cover.z + (1 | reef_code/transect_option_code) + (1 | year),
#   dispformula = ~ dmaho_cover.z,
#   data = songs.no_zeroes,
#   family = gaussian()
# )
# summary(kelp_glmm.z)
# 
# # Diagnostics
# res.z <- simulateResiduals(kelp_glmm.z)
# plot(res.z) # Tests failed

# Try non-scaled response, scaled predictor, negative binomial model family
kelp_glmm.z <- glmmTMB(
  mapy_count_per_unit_area ~ dmaho_cover.z + (1 | reef_code/transect_option_code) + (1 | year),
  #dispformula = ~ dmaho_cover.z,
  data = songs.no_zeroes,
  family = nbinom2(link = "log")
)
summary(kelp_glmm.z)

# Diagnostics
res.z <- simulateResiduals(kelp_glmm.z)
plot(res.z) # Tests failed


## Visualization
# Get predicted relationship from z-score model
preds.z <- ggpredict(kelp_glmm.z, terms = "dmaho_cover.z")

# Plot predictions + raw data
ggplot() +
  geom_jitter(data = songs.no_zeroes, 
              aes(x = dmaho_cover.z, y = mapy_count.z), 
              color = "darkgrey", 
              alpha = 0.6, 
              width = 0.15) +
  geom_line(data = preds.z, 
            aes(x = x, y = predicted), 
            linewidth = 1.2) +
  geom_ribbon(data = preds.z, 
              aes(x = x, ymin = conf.low, ymax = conf.high), 
              alpha = 0.3) +
  labs(x = "Cover of dead holdfasts (Z-score)",
       y = Juvenile~kelp~density~(count/m^2/yr)) +
  theme_classic(base_size = 14)


## Effect size
songs_effect.z <- tidy(kelp_glmm.z, effects = "fixed", conf.int = TRUE, conf.level = 0.95) %>% 
  filter(term == "dmaho_cover.z")

write_csv(songs_effect.z, "Datasets/Effect sizes/Standardized/songs_effect.z.csv")

