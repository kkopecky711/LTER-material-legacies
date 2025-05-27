# Code for anlaysis of Santa Barbara Coastal LTER data
# Analysis of the inlfuence of dead giant kelp holdfasts on demographic processes in live macroalgae



# Packages
library(tidyverse)
library(janitor)
library(glmmTMB)
library(ggeffects) 

# Working directory
setwd('/Users/kako4300/Library/CloudStorage/OneDrive-UCB-O365/Projects/LTER material legacy synthesis')

#### Annual kelp forest monitoring data ####
# Read in datafile and filter for only algae observations 
sbc.juv_macro <- read_csv("Datasets/Santa Barbara Coastal/Annual_Quad_Swath_All_Years_20240823.csv") %>% 
  clean_names() %>% 
  filter(sp_code == "MPJ",
         count > -1,
         year > 2007) %>% 
  mutate(obvs_id = paste(year, site, transect, quad, side, sep = "_")) %>% 
  select(c(obvs_id, year, site, transect, quad, side, count)) %>% 
  rename(juv_density = count)

# Create dataframe with only dead macrocystis holdfasts
sbc.dead_macro <- read_csv("Datasets/Santa Barbara Coastal/SBS_Cover_All_Years_20240823.csv") %>%
  clean_names() %>% 
  filter(sp_code == "DMH") %>% 
  mutate(obvs_id = paste(year, site, transect, quad, side, sep = "_")) %>% 
  select(c(obvs_id, percent_cover)) %>% 
  rename(dead_cover = percent_cover)

# Merge juvenile macro and dead holdfast dataframes
sbc.juv_dead <- merge(sbc.juv_macro, sbc.dead_macro, by = "obvs_id")
 
# Exploratory graph of juvenile macrocystis count ~ % cover of dead holdfasts
ggplot(sbc.juv_dead, aes(x = dead_cover, y = juv_density)) +
  geom_jitter() +
  geom_smooth(method = "lm") +
  theme_minimal()

## Mixed effects model of juv count ~ dead cover, site as random effect
# Fit a zero-inflated negative binomial mixed model
zinb_model <- glmmTMB(
  juv_density ~ dead_cover + (1|site/transect),
  ziformula = ~1,  # model zero-inflation with an intercept
  family = nbinom2,
  data = sbc.juv_dead
)

# Run a negative binomial model with no random effects and no zero-inflation to compare model fit
nb_model <- glmmTMB(juv_density ~ dead_cover, family = nbinom2, data = sbc.juv_dead)

# Compare model fits with AIC scores
 AIC(nb_model, zinb_model)
 
# Zero-inflated model fits better, examine model outputs and run diagnostics
 summary(zinb_model)
 
library(DHARMa)

sim_res <- simulateResiduals(zinb_model)
plot(sim_res)
testDispersion(sim_res)
testZeroInflation(sim_res)

## Visualizations
# Generate predicted values across the range of dead_cover
preds <- ggpredict(
  zinb_model,
  terms = "dead_cover")

# Plot predicted values with 95% confidence intervals
ggplot(preds, aes(x = x, y = predicted)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
              alpha = 0.3,
              fill = "#3B4F8E") +
  geom_jitter(data = sbc.juv_dead, aes(x = dead_cover, y = juv_density),
              alpha = 0.6,
              size = 2,
              width = 0.1) +
  labs(
    x = "Percent cover of dead holdfasts",
    y = "Predicted juvenile kelp density (no./m²)"
  ) +
  theme_classic(base_size = 14)

#### Normalized model for ecosystem comparison ----
# Create a copy of the data
sbc.juv_dead.z <- sbc.juv_dead

# Z-score standardize the response and predictor
sbc.juv_dead.z$juv_density_z <- scale(sbc.juv_dead.z$juv_density)[, 1]
sbc.juv_dead.z$dead_cover_z  <- scale(sbc.juv_dead.z$dead_cover)[, 1]

zinb_z_gaussian <- glmmTMB(
  juv_density_z ~ dead_cover_z + (1 | site/transect),
  ziformula = ~1,
  family = gaussian,
  data = sbc.juv_dead.z
)

library(broom.mixed)

effect_z <- tidy(zinb_z_gaussian, effects = "fixed", conf.int = TRUE) |>
  filter(term == "dead_cover_z")

print(effect_z)
