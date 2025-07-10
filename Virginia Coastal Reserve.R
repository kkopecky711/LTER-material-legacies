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

# Filter original data for just live and dead oysters
vcr.juv_dead <- vcr.raw %>% 
  clean_names() %>%
  filter(species == "Box Adult Oyster" | species == "Spat Oyster",
         restoration == "Reference") %>% 
  mutate(year = substr(date, nchar(date) - 1, nchar(date)),
         year = paste0("20", year),
         date = mdy(date), # Convert the date column to Date format
         month = month(date), # Extract the month
         species = case_when(species == "Spat Oyster" ~ "Juvenile",
                             TRUE ~ "Dead"),
         season = case_when( # Assign seasons based on month
           month %in% c(12, 1, 2) ~ "Winter",
           month %in% c(3, 4, 5) ~ "Spring",
           month %in% c(6, 7, 8) ~ "Summer",
           month %in% c(9, 10, 11) ~ "Fall",
           TRUE ~ NA_character_),
         obvs_id = paste(year, season, site, sample, sep = "_")) %>%
  select(obvs_id, site, year, season, species, species_count) %>% 
  pivot_wider(names_from = "species", 
              values_from = "species_count", 
              values_fn = sum)
  
# Exploratory graph of Juvenile oyster density ~ dead oyster density in previous year
ggplot(vcr.juv_dead, aes(x = Dead, y = Juvenile, color = season)) +
  geom_point(alpha = 0.5,
             aes(color = season)) +
  geom_smooth(method = "lm",
              aes(group = season)) +
  #scale_x_continuous(expand = c(0.01,0)) +
  #scale_y_continuous(expand = c(0.01,0)) +
  labs(x = "Dead oyster density",
       y = "Juvenile oyster density") +
  theme_classic()

# Remove summer observations because they are not representative and Smith 2006 observations
vcr.juv_dead <- vcr.juv_dead %>% 
  filter(season != "Summer")

# Remove observations from Smith 2006, which contains a few high outliers in dead oyster density (> 350 m^2)
vcr.juv_dead.filtered <- vcr.juv_dead %>%
  filter(!(year == "2006" & site == "Smith"))


## GLMM with negative binomial distribution (due to overdispersion) and zero-inflation
mod_zinb <- glmmTMB(
  Juvenile ~ Dead + (1 | site) + (1 | year),
  ziformula = ~1,
  data = vcr.juv_dead.filtered,
  family = nbinom2() # negative binomial because variables are count data with overdispersion
)
summary(mod_zinb)

# Standard negative binomial model (no zero inflation)
mod_nb <- glmmTMB(
  Juvenile ~ Dead + (1 | site) + (1 | year),
  data = vcr.juv_dead.filtered,
  family = nbinom2
)

# Compare AIC scores to ensure zero-inflation improves model fit
AIC(mod_nb, mod_zinb) # zero-inflated has better fit

# Diagnostics
simres <- simulateResiduals(mod_zinb)
plot(simres)

# To report effect size and CI in main text and table
oyster_effect.raw <- tidy(mod_zinb, effects = "fixed", conf.int = TRUE, conf.level = 0.95) %>% 
  filter(term == "Dead")

write_csv(oyster_effect.raw, "Datasets/Effect sizes/Raw/vcr.effect_raw.csv")

## Visualization
# Generate datatable of model predictions
preds <- ggpredict(
  mod_zinb,
  terms = "Dead")

# Plot predicted values with 95% confidence intervals
ggplot(preds, aes(x = x, y = predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
              alpha = 0.3) +
  geom_point(data = vcr.juv_dead.filtered, aes(x = Dead, y = Juvenile),
              alpha = 0.6,
              size = 2) +
  #scale_y_continuous(limits = c(0, 5000)) +
  labs(
    x = "Dead oyster density (no./0.25 m²)",
    y = "Juvenile oyster density (no./0.25 m²)"
  ) +
  #scale_y_continuous(limits = c(0, 4000)) +
  theme_classic(base_size = 14)

#### Normalized model for ecosystem comparison ----
# Create a copy of the data
vcr.juv_dead.z <- vcr.juv_dead

# Z-score standardize the response and predictor
vcr.juv_dead.z$juv_density_z <- scale(vcr.juv_dead.z$Juvenile)[, 1]
vcr.juv_dead.z$dead_density_z  <- scale(vcr.juv_dead.z$Dead)[, 1]

zinb_z_gaussian <- glmmTMB(
  juv_density_z ~ dead_density_z + (1 | site/year),
  family = gaussian,
  data = vcr.juv_dead.z
)
summary(zinb_z_gaussian)

## Visualization
preds_z <- ggpredict(
  zinb_z_gaussian,
  terms = "dead_density_z")

# Plot predicted values with 95% confidence intervals
ggplot(preds_z, aes(x = x, y = predicted)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
              alpha = 0.3,
              fill = "#3B4F8E") +
  geom_point(data = vcr.juv_dead.z, aes(x = dead_density_z, y = juv_density_z),
             alpha = 0.6,
             size = 2) +
  labs(
    x = "Dead oyster density (Z-score)",
    y = "Juvenile oyster density (Z-score)"
  ) +
  theme_classic(base_size = 14)

# Extract effect size
library(broom.mixed)

oyster_effect <- tidy(zinb_z_gaussian, effects = "fixed", conf.int = TRUE) |>
  filter(term == "dead_density_z")

print(oyster_effect)

write_csv(oyster_effect, "Datasets/Effect sizes/vcr.effect_size.csv")
