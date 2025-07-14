#### Konza Prairie LTER 

library(tidyverse)
library(janitor)
library(glmmTMB)
library(DHARMa)
library(ggeffects)
library(broom.mixed)

setwd('/Users/kako4300/Library/CloudStorage/OneDrive-UCB-O365/Projects/LTER material legacy synthesis')

#### Comparisons of one and two year burns with burn treatment as categorical predictor

# Annual burn watershed
knz.annual.1d <- read_csv("Datasets/Konza Prairie/PAB011.csv") %>% 
  clean_names() %>%  
  filter(watershed == "001d")%>%
  rename(year = recyear)

knz.annual.1d.means <- knz.annual.1d %>% 
  group_by(year, soiltype, transect) %>% 
  summarize(lvgrass.mean = mean(lvgrass)) %>% 
  mutate(burn_cat = "Burned")
  
# Two year watershed
knz.two_year.2d <- read_csv("Datasets/Konza Prairie/PAB041.csv") %>% 
  clean_names() %>% 
  filter(watershed == "002d") %>% 
  rename(year = recyear)

# Burn years
burn_years <- read_csv("Datasets/Konza Prairie/KFH011.csv") %>% 
  clean_names() %>% 
  filter(watershed == "2D") %>%
  select(year) %>% 
  distinct(year, .keep_all = TRUE) %>% 
  mutate(burn_year = "yes")

# Create complete set of burned and unburned years
all_years <- tibble(year = seq(1978, max(knz.two_year.2d$year)))

# Full join with all years and replace NAs with "no"
all_burn_years <- all_years %>%
  left_join(burn_years, by = "year") %>%
  mutate(burn_year = replace_na(burn_year, "no"))

# Merge burn data with biomass data
knz.two_year.no_burns <- knz.two_year.2d %>% 
  left_join(all_burn_years, by = "year") %>% 
  filter(burn_year == "no") %>% 
  select(-burn_year)

knz.two_year.no_burns.means <- knz.two_year.no_burns %>% 
  group_by(year, soiltype, transect) %>% 
  summarize(lvgrass.mean = mean(lvgrass)) %>% 
  mutate(burn_cat = "Unburned")

# Extract unique years from the non-burned dataset
target_years <- unique(knz.two_year.no_burns.means$year)

# Filter annual dataset to only those years
knz.annual.1d.means.filtered <- knz.annual.1d.means %>%
  filter(year %in% target_years)

# Merge one and two year watersheds, filter for only Florentine soil type
knz.merged <- rbind(knz.annual.1d.means.filtered, knz.two_year.no_burns.means)
knz.merged <- knz.merged %>% 
  mutate(burn_cat = as.factor(burn_cat)) %>% 
  filter(soiltype == "fl")

## GLMM of live grass biomass ~ burn category
grass_glmm_log <- glmmTMB(
  lvgrass.mean ~ burn_cat + (1 | year) + (1 | transect),
  data = knz.merged,
  family = gaussian(link = "log")
)
summary(grass_glmm_log)

# Diagnostics
res_gauss <- simulateResiduals(grass_glmm_log)
plot(res_gauss) # Tests are non-significant

## Visualization
# Generate datatable of model predictions
preds <- ggpredict(grass_glmm_log, terms = "burn_cat")

# Multiply response by 10 to scale up to g/m^2
ggplot() +
  geom_jitter(
    data = knz.merged,
    aes(x = burn_cat, y = lvgrass.mean*10),
    color = "darkgrey", 
    alpha = 0.6, 
    width = 0.15) +
  geom_point(
    data = preds,
    aes(x = x, y = predicted*10),
    size = 3) +
  geom_errorbar(
    data = preds,
    aes(x = x, ymin = conf.low*10, ymax = conf.high*10),
    width = 0) +
  geom_line(data = preds,
            aes(x = x, y = predicted*10, group = group)) +
  labs(x = "Burn Category",
    y = "Live grass biomass (g/m²/yr)") +
  theme_classic(base_size = 14)

# To report effect size and CI in main text and table 
knz_effect.raw <- tidy(grass_glmm_tw, effects = "fixed", conf.int = TRUE, conf.level = 0.95) %>% 
  filter(term == "burn_catUnburned")

write_csv(knz_effect.raw, "Datasets/Effect sizes/Raw/knz_effect.raw.csv")

## Z-score standardized model

# Scale response only to its Z-score
knz.merged <- knz.merged %>%
  mutate(
    lvgrass_z = scale(lvgrass.mean)[, 1]  # extract numeric vector from scale object
  )

grass_glmm_z <- glmmTMB(
  lvgrass_z ~ burn_cat + (1 | year) + (1 | transect),
  data = knz.merged,
  family = gaussian()
)
summary(grass_glmm_z)

# Diagnostics
res.z <- simulateResiduals(grass_glmm_z)
plot(res.z) # Tests are non-significant

# Extract effect size and CI, write to .csv file
knz_effect.z <- tidy(grass_glmm_z, effects = "fixed", conf.int = TRUE) %>% 
  filter(term == "burn_catUnburned")

write_csv(knz_effect.z, "Datasets/Effect sizes/Standardized/knz_effect.z.csv")
