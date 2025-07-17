## Virginia Coastal Reserve oyster reef data ##

# Packages
library(tidyverse)
library(janitor)
library(glmmTMB)
library(ggeffects)
library(DHARMa)
library(broom.mixed)

# Working directory
setwd('/Users/kako4300/Library/CloudStorage/OneDrive-UCB-O365/Projects/LTER material legacy synthesis')

vcr.raw <- read_csv("Datasets/Virginia Coastal Reserve/Oyster_Count_Data_2022_03_10.csv")

# Filter original data for just live and dead oysters, create column for season
vcr.juv_dead <- vcr.raw %>% 
  clean_names() %>%
  filter(species == "Box Adult Oyster" | species == "Spat Oyster",
         restoration == "Reference") %>% 
  mutate(year = substr(date, nchar(date) - 1, nchar(date)),
         year = paste0("20", year),
         species = case_when(species == "Spat Oyster" ~ "Juvenile",
                             TRUE ~ "Dead"),
         obvs_id = paste(year, site, date, sample, sep = "_")) %>%
  select(obvs_id, site, year, date, species, species_count) %>% 
  pivot_wider(names_from = "species", 
              values_from = "species_count", 
              values_fn = sum)
  
# Calculate means across quadrats for site and year
vcr.juv_dead.summary <- vcr.juv_dead %>% 
  group_by(year, site) %>% 
  summarize(dead.mean = mean(Dead), 
            juv.mean = mean(Juvenile)) %>% 
  filter(dead.mean < 250)

## Analysis
# GLMM of mean annual juvenile density ~ mean dead density
oyster_glmm.raw <- glmmTMB(
  juv.mean ~ dead.mean + (1 | site) + (1 | year),
  family = Gamma(link = "log"),
  data = vcr.juv_dead.summary 
)
summary(oyster_glmm.raw)

# Diagnostics
res <- simulateResiduals(oyster_glmm.raw)
plot(res) # KS test non-significant; residual vs predicted test signnificant, but attempts to improve were not successful; proceeding with caution

# To report effect size and CI in main text and table
vcr_effect.raw <- tidy(oyster_glmm.raw, effects = "fixed", conf.int = TRUE, conf.level = 0.95) %>% 
  filter(term == "dead.mean")

write_csv(vcr_effect.raw, "Datasets/Effect sizes/Raw/vcr_effect.raw.csv")

## Visualization
# Generate datatable of model predictions
preds <- ggpredict(
  oyster_glmm.raw,
  bias_correction = TRUE,
  terms = "dead.mean")

# Plot predicted values with 95% confidence intervals and raw values; multiply predictor and response by 4 to get #/m^2 
ggplot(preds, aes(x = x*4, y = predicted*4)) +
  geom_point(data = vcr.juv_dead.summary, aes(x = dead.mean*4, y = juv.mean*4),
             alpha = 0.6,
             color = "darkgrey") +
  geom_line(linewidth = 1.2) +
  geom_ribbon(aes(ymin = conf.low*4, ymax = conf.high*4), 
              alpha = 0.3) +
  labs(
    x = "Dead oyster density (no./m²)",
    y = "Juvenile oyster density (no./m²/yr)"
  ) +
  theme_classic(base_size = 14)

## Standardized model ----
# Scale response and predictor to z-scores
vcr.juv_dead.summary$juv_density_z <- scale(vcr.juv_dead.summary$juv.mean)[, 1]
vcr.juv_dead.summary$dead_density_z  <- scale(vcr.juv_dead.summary$dead.mean)[, 1]

# Fit model to z-scored predictor and response
oyster_glmm.z <- glmmTMB(
  juv_density_z ~ dead_density_z + (1 | site) + (1 | year),
  family = gaussian(),
  data = vcr.juv_dead.summary
)
summary(oyster_glmm.z)

# Diagnostics
res.z <- simulateResiduals(oyster_glmm.z)
plot(res.z) # Tests are non-significant

## Visualization
preds_z <- ggpredict(
  oyster_glmm.z,
  terms = "dead_density_z")

# Plot predicted values with 95% confidence intervals
ggplot(preds_z, aes(x = x, y = predicted)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
              alpha = 0.3,
              fill = "#3B4F8E") +
  geom_point(data = vcr.juv_dead, aes(x = dead_density_z, y = juv_density_z),
             alpha = 0.6,
             size = 2) +
  labs(
    x = "Dead oyster density (Z-score)",
    y = "Juvenile oyster density (Z-score)"
  ) +
  theme_classic(base_size = 14)

# Extract effect size
vcr_effect.z <- tidy(oyster_glmm.z, effects = "fixed", conf.int = TRUE) |>
  filter(term == "dead_density_z")

write_csv(vcr_effect.z, "Datasets/Effect sizes/Standardized/vcr_effect.z.csv")
