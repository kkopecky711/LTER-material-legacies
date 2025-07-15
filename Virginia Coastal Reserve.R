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
  select(obvs_id, site, year, season, month, species, species_count) %>% 
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

vcr.juv_dead <- vcr.juv_dead %>% 
  mutate(month = as.factor(month))

ggplot(vcr.juv_dead, aes(x = month, y = Juvenile)) +
  geom_boxplot() +
  theme_minimal()

# Remove summer observations because they are not representative and Smith 2006 observations
vcr.juv_dead <- vcr.juv_dead %>% 
  filter(season != "Summer")

# Remove observations from Smith 2006, which contains a few high outliers in dead oyster density (> 350 m^2)
vcr.juv_dead.filtered <- vcr.juv_dead %>%
  filter(!(year == "2006" & site == "Smith"))


## GLMM with negative binomial distribution (due to overdispersion) and zero-inflation
oyster_glmm.raw <- glmmTMB(
  Juvenile ~ Dead + (1 | site) + (1 | year),
  ziformula = ~1,
  data = vcr.juv_dead.filtered,
  family = nbinom2() # negative binomial because variables are count data with overdispersion
)
summary(oyster_glmm.raw)

# Standard negative binomial model (no zero inflation)
mod_nb <- glmmTMB(
  Juvenile ~ Dead + (1 | site) + (1 | year),
  data = vcr.juv_dead.filtered,
  family = nbinom2
)

# Compare AIC scores to ensure zero-inflation improves model fit
AIC(mod_nb, oyster_glmm.raw) # zero-inflated has better fit

# Diagnostics
res <- simulateResiduals(oyster_glmm.raw)
plot(res)

# To report effect size and CI in main text and table
vcr_effect.raw <- tidy(oyster_glmm.raw, effects = "fixed", conf.int = TRUE, conf.level = 0.95) %>% 
  filter(term == "Dead")

write_csv(vcr_effect.raw, "Datasets/Effect sizes/Raw/vcr_effect.raw.csv")

## Visualization
# Generate datatable of model predictions
preds <- ggpredict(
  oyster_glmm.raw,
  bias_correction = TRUE,
  terms = "Dead")

# Plot predicted values with 95% confidence intervals and raw values; multiply predictor and response by 4 to get #/m^2 
ggplot(preds, aes(x = x*4, y = predicted*4)) +
  geom_point(data = vcr.juv_dead.filtered, aes(x = Dead*4, y = Juvenile*4),
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
# Z-score standardize the response and predictor
vcr.juv_dead$juv_density_z <- scale(vcr.juv_dead$Juvenile)[, 1]
vcr.juv_dead$dead_density_z  <- scale(vcr.juv_dead$Dead)[, 1]

# Fit model to z-scored predictor and response
oyster_glmm.z <- glmmTMB(
  juv_density_z ~ dead_density_z + (1 | site) + (1 | year),
  family = gaussian(),
  dispformula = ~ dead_density_z,
  data = vcr.juv_dead
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
