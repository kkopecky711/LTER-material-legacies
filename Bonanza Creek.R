#### Bonanza Creek LTER ####
# Analysis of seed propagation from burned spruce trees following a large fire

# Packages
library(tidyverse)
library(janitor)
library(glmmTMB)
library(ggeffects)
library(DHARMa)
library(broom.mixed)

# Working directory
setwd('/Users/kako4300/Library/CloudStorage/OneDrive-UCB-O365/Projects/LTER material legacy synthesis')

# Load seed data
seeds <- read_csv("Datasets/Bonanza Creek/AK2004 sites seeds.csv") %>% 
  clean_names() %>% 
  select(burn, site, bs_stg_ba, total_m2) %>% 
  mutate(burn = as.factor(burn),
         site = as.factor(site)) %>% 
  drop_na()

# Exploratory graph of total seeds ~ basal area
ggplot(seeds, aes(x = bs_stg_ba, y = total_m2)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic()

## Analyses
# GLMM of seed density ~ standing basal area of burned trees
seed_model <- glmmTMB(
  total_m2 ~ bs_stg_ba + (1 | burn/site) + , # random effect of site nested within burn complex
  family = tweedie(link = "log"),
  data = seeds)

# False convergence warning; tried many alternative options, none of which resolved this. Proceeding anyway, but with caution

summary(seed_model) # Model output looks reasonable
tidy(seed_model, effects = "fixed", conf.int = TRUE, conf.level = 0.95)

## Visualization
# Generate predicted values
preds <- ggpredict(seed_model, terms = "bs_stg_ba")

# Plot predictions with 95% CI
ggplot(preds, aes(x = x, y = predicted)) +
  geom_line(size = 1.2) + 
            #color = "#3366AA") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
              #fill = "#3366AA", 
              alpha = 0.3) +
  geom_point(data = seeds,
             aes(x = bs_stg_ba,
                 y = total_m2),
             alpha = 0.6) +
  labs(
    x = "Burned stem basal area (m²)",
    y = "Seed density (no./m², log10)") +
  scale_y_continuous(trans = "log10") + 
  theme_classic(base_size = 14)

## Z-score standardization
seeds_z <- seeds %>%
  mutate(
    bs_stg_ba_z = scale(bs_stg_ba)[, 1],
    total_m2_z = scale(total_m2)[, 1]
  )

seed_model_z <- glmmTMB(
  total_m2_z ~ bs_stg_ba_z + (1 | burn/site),
  family = gaussian,  # Tweedie not needed since response is now continuous and standardized
  data = seeds_z
)
summary(seed_model_z)
# Visualization
preds_z <- ggpredict(seed_model_z, terms = "bs_stg_ba_z")

ggplot() +
  geom_point(data = seeds_z, aes(x = bs_stg_ba_z, y = total_m2_z), alpha = 0.6, size = 2) +
  geom_line(data = preds_z, aes(x = x, y = predicted), color = "black", size = 1.2) +
  geom_ribbon(data = preds_z, aes(x = x, ymin = conf.low, ymax = conf.high),
              fill = "#3366AA", alpha = 0.3) +
  labs(
    x = "Burned stem basal area (Z-score)",
    y = "Seed density (Z-score)"
  ) +
  theme_classic(base_size = 14)

# Effect size
boreal_effect <- tidy(seed_model_z, effects = "fixed", conf.int = TRUE) %>%
  filter(term == "bs_stg_ba_z")

print(boreal_effect)

write_csv(boreal_effect, "Datasets/Effect sizes/bnz.effect_size.csv")

