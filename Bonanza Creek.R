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
seed_glmm.raw <- glmmTMB(
  total_m2 ~ bs_stg_ba + (1 | burn/site), 
  family = tweedie(link = "log"),
  data = seeds)

# False convergence warning; tried many alternative options, none of which resolved this. Proceeding anyway, but with caution; diagnostics dests won't run due to non-convergence
# res <- simulateResiduals(seed_glmm.raw)
# plot(res)

summary(seed_glmm.raw) # Model output looks reasonable

bnz_effect.raw <- tidy(seed_glmm.raw, effects = "fixed", conf.int = TRUE, conf.level = 0.95) %>%
  filter(term == "bs_stg_ba")

write_csv(bnz_effect.raw, "Datasets/Effect sizes/Raw/bnz_effect.raw.csv")

## Visualization
# Generate predicted values
preds <- ggpredict(seed_glmm.raw, terms = "bs_stg_ba")

# Plot predictions with 95% CI and raw values
ggplot(preds, aes(x = x, y = predicted)) +
  geom_point(data = seeds,
             aes(x = bs_stg_ba,
                 y = total_m2),
             alpha = 0.6,
             color = "#6B6C58") +
  geom_line(linewidth = 0.75) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
              alpha = 0.3,
              fill = "#6B6C58") +
  labs(
    x = "Burned stem basal area (m²)",
    y = expression("Seed density (no./m"^2*", log"[10]*")")
  ) + 
  scale_y_continuous(trans = "log10") + 
  theme_classic(base_size = 14)

## Standardized model ----
seeds <- seeds %>%
  mutate(
    bs_stg_ba.z = scale(bs_stg_ba)[, 1],
    total_m2.z = scale(total_m2)[, 1]
  )

seed_glmm.z <- glmmTMB(
  total_m2.z ~ bs_stg_ba.z + (1 | burn/site),
  family = gaussian(),
  data = seeds
)
summary(seed_glmm.z)

# Diagnostics
res.z <- simulateResiduals(seed_glmm.z)
plot(res.z) # Tests show moderate deviation; attempts to modify model structure did not improve diagnostics; proceeding with caution when interpreting model output

## Visualization
# Generate predicted values
preds.z <- ggpredict(seed_glmm.z, terms = "bs_stg_ba.z")

# Plot predictions with 95% CI and raw values
ggplot(preds.z, aes(x = x, y = predicted)) +
  geom_point(data = seeds,
             aes(x = bs_stg_ba.z,
                 y = total_m2.z),
             alpha = 0.6,
             color = "darkgrey") +
  geom_line(linewidth = 1.2) + 
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
              alpha = 0.3) +
  labs(x = "Burned stem basal area (Z-score)",
    y = "Seed density (Z-score)") +
  theme_classic(base_size = 14)


## Effect size
bnz_effect.z <- tidy(seed_glmm.z, effects = "fixed", conf.int = TRUE) %>%
  filter(term == "bs_stg_ba_z")

write_csv(bnz_effect.z, "Datasets/Effect sizes/Standardized/bnz_effect.z.csv")

