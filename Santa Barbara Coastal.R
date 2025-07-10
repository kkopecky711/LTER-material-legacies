# Code for analysis of Santa Barbara Coastal LTER data
# Analysis of the influence of dead giant kelp holdfasts on demographic processes in live macroalgae

## Chat with Dan Reed
# SBC data does not have a good match up dead holdfasts and juvenile kelp
# SONGS data might be a better fit for the question


# Packages
library(tidyverse)
library(janitor)
library(glmmTMB)
library(ggeffects)
library(DHARMa)

# Working directory
setwd('/Users/kako4300/Library/CloudStorage/OneDrive-UCB-O365/Projects/LTER material legacy synthesis')

#### Annual kelp forest monitoring data ####
# Read in datafile and filter for only juvenile macrocystis observations 
sbc.juv_macro <- read_csv("Datasets/Santa Barbara Coastal/Annual_Quad_Swath_All_Years_20240823.csv") %>% 
  clean_names() %>% 
  filter(sp_code == "MPJ",
         count > -1,
         year > 2007) %>% 
  select(c(year, site, transect,  quad, side, count)) %>% 
  rename(juv_density = count)

# Create summarized dataframe by swath
sbc.juv_macro.swath <- sbc.juv_macro %>% 
  mutate(swath = if_else(quad %in% c(0,8, 16), 20, 40)) %>% 
  group_by(year, site, transect, swath, side) %>% 
  summarize(juv_density.mean = mean(juv_density)) %>% 
  mutate(obs_id = paste(year, site, transect, swath, side, sep = "_")) %>% 
  ungroup() 

# Create dataframe with only dead macrocystis holdfasts
sbc.dead_macro <- read_csv("Datasets/Santa Barbara Coastal/SBS_Cover_All_Years_20240823.csv") %>%
  clean_names() %>% 
  filter(sp_code == "DMH",
         percent_cover < 20) %>% 
  rename(swath = quad,
         dead_cover = percent_cover) %>% 
  mutate(obs_id = paste(year, site, transect, swath, side, sep = "_")) %>% 
  select(c(obs_id, dead_cover))

sbc.annual <- read_csv("Datasets/Santa Barbara Coastal/Annual_Cover_All_Years_algae-only.csv") %>% 
  clean_names() %>% 
  filter(sp_code == "MPJ")

# Remove any duplicate rows in juvenile macro and dead holdfast dataframes, then merge
sbc.juv_macro.swath <- sbc.juv_macro.swath[!duplicated(sbc.juv_macro.swath$obs_id), ]
sbc.dead_macro <- sbc.dead_macro[!duplicated(sbc.dead_macro$obs_id), ]

sbc.juv_dead <- inner_join(sbc.juv_macro.swath, sbc.dead_macro, by = "obs_id")
 
# Exploratory graph of juvenile macrocystis count ~ % cover of dead holdfasts
ggplot(sbc.juv_dead, aes(x = dead_cover, y = juv_density.mean)) +
  geom_jitter() +
  geom_smooth(method = "lm") +
  theme_minimal()

## Mixed effects model of juv count ~ dead cover, site as random effect
# Fit a zero-inflated negative binomial mixed model
zinb_model <- glmmTMB(
  juv_density.mean ~ dead_cover + (1|site/transect/swath) + (1 | year),
  ziformula = ~1,  # model zero-inflation with an intercept
  family = tweedie(link = "log"), # for skewed distribution of continuous values
  data = sbc.juv_dead
)

# Run a negative binomial model with no random effects and no zero-inflation to compare model fit
nb_model <- glmmTMB(juv_density.mean ~ dead_cover, 
                    family = tweedie(link = "log"), 
                    data = sbc.juv_dead)

# Compare model fits with AIC scores
 AIC(nb_model, zinb_model)
 
# Zero-inflated model fits better (lower AIC score), examine model outputs and run diagnostics
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
  geom_jitter(data = sbc.juv_dead, aes(x = dead_cover, y = juv_density.mean),
              alpha = 0.6,
              size = 2,
              width = 0.1) +
  labs(
    x = "Percent cover of dead holdfasts",
    y = "Juvenile kelp density (no./m²)"
  ) +
  theme_classic(base_size = 14)

#### Normalized model for ecosystem comparison ----
# Create a copy of the data
sbc.juv_dead.z <- sbc.juv_dead

# Z-score standardize the response and predictor
sbc.juv_dead.z$juv_density_z <- scale(sbc.juv_dead.z$juv_density.mean)[, 1]
sbc.juv_dead.z$dead_cover_z  <- scale(sbc.juv_dead.z$dead_cover)[, 1]

zinb_z_gaussian <- glmmTMB(
  juv_density_z ~ dead_cover_z + (1 | site/transect/swath) + (1 | year),
  ziformula = ~1,
  family = gaussian,
  data = sbc.juv_dead.z
)

# Z-score visualization
## Visualizations
# Generate predicted values across the range of dead_cover
preds_z <- ggpredict(
  zinb_z_gaussian,
  terms = "dead_cover_z")

# Plot predicted values with 95% confidence intervals
ggplot(preds_z, aes(x = x, y = predicted)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
              alpha = 0.3,
              fill = "#3B4F8E") +
  geom_jitter(data = sbc.juv_dead.z, aes(x = dead_cover_z, y = juv_density_z),
              alpha = 0.6,
              size = 2,
              width = 0.1) +
  labs(
    x = "Cover of dead holdfasts (Z-score)",
    y = "Juvenile kelp density (Z-score)"
  ) +
  theme_classic(base_size = 14)


library(broom.mixed)

effect_z <- tidy(zinb_z_gaussian, effects = "fixed", conf.int = TRUE) |>
  filter(term == "dead_cover_z")

print(effect_z)
