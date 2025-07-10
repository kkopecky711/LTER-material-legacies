#### Konza Prairie LTER ####

## Notes
# Explore relationship of previous year dead biomass on current year productivity
# Should not group forbs and grass; grass could be considered the foundation species, and forbs neither produce legacies in the same way, nor are they annuals

# Look at relationship between dead biomass and biomass of forbs and woody 

# potential relationships between dead biomass and abiotic conditions (soil temperature, soil moisture, N content)

# Look into % cover data as well

# Watershed = fire treatment: 004 means burned every four years, 020 = every 20 years, b = statistical replicate
# Soiltype: tu = lowland, thick, moist; fl = upland, thin, rocky, drier

# Consider trying a running mean approach

library(tidyverse)
library(janitor)
library(glmmTMB)
library(DHARMa)
library(ggeffects)

setwd('/Users/kako4300/Library/CloudStorage/OneDrive-UCB-O365/Projects/LTER material legacy synthesis')

### Burn treatment as categorical predictor ----

#### Comparisons of one and two year burns

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

# Merge one and two year watersheds
knz.merged <- rbind(knz.annual.1d.means.filtered, knz.two_year.no_burns.means)
knz.merged <- knz.merged %>% 
  mutate(burn_cat = as.factor(burn_cat)) %>% 
  filter(soiltype == "fl")

## GLMM of live grass biomass ~ burn category
# grass_glmm <- glmmTMB(
#   lvgrass.mean ~ burn_cat + (1 | year) + (1 | transect),
#   data = knz.merged,
#   family = Gamma(link = "log")  # Gamma distribution for continuous, skewed, positive biomass
# )
# summary(grass_glmm)
# 
# # Diagnostics
# res <- simulateResiduals(grass_glmm)
# plot(res)

# Residuals not uniform, try log-transofrmed gaussian
# knz.merged <- knz.merged %>%
#   mutate(lvgrass_log = log1p(lvgrass.mean))

# grass_glmm_log <- glmmTMB(
#   lvgrass_log ~ burn_cat + (1 | year) + (1 | soiltype/transect),
#   data = knz.merged,
#   family = gaussian()
# )
# summary(grass_glmm_log)
# 
# res_log <- simulateResiduals(grass_glmm_log)
# plot(res_log)

# Use Tweedie distribution with log-link function
grass_glmm_tw <- glmmTMB(
  lvgrass.mean ~ burn_cat + (1 | year) + (1 | transect),
  data = knz.merged,
  family = tweedie(link = "log")
)
summary(grass_glmm_tw)

res_tw <- simulateResiduals(grass_glmm_tw)
plot(res_tw)

AIC(grass_glmm, grass_glmm_log, grass_glmm_tw)

# Tweedie model passes diagnostics and has lower AIC score, will stick with this version
preds <- ggpredict(grass_glmm_tw, terms = "burn_cat")

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
    y = "Live grass biomass (g/m²)") +
  theme_classic(base_size = 14)

# To report effect size and CI in main text
tidy(grass_glmm_tw, effects = "fixed", conf.int = TRUE, conf.level = 0.95)%>% 
  filter(term == "burn_catUnburned")

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

# Extract effect size and CI, write to .csv file
prairie_effect <- tidy(grass_glmm_z, effects = "fixed", conf.int = TRUE) %>% 
  filter(term == "burn_catUnburned")

print(prairie_effect)

write_csv(prairie_effect, "Datasets/Effect sizes/knz.effect_size.csv")


### Dead grass biomass as continuous predictor ----
# Exploratory graph of live grass biomass ~ previous year dead grass biomass, separated by burn frequency treatment
ggplot(knz, aes(x = pryrdead, y = lvgrass, color = watershed)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", 
              aes(color = watershed)) +
  #scale_x_continuous(expand = c(0.01,1)) +
  # scale_y_continuous(expand = c(0.01,1)) +
  labs(x = "Previous year dead biomass (g/m^2)",
       y = "Current year live biomass (g/m^2)") +
  theme_classic()

knz.four <- knz %>% 
  filter(watershed == "004b",
         pryrdead > 0)

ggplot(knz.four, aes(x = pryrdead, y = lvgrass)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm") +
  labs(x = "Previous year dead biomass (g/m^2)",
       y = "Current year live biomass (g/m^2)") +
  theme_classic()

#### Annual burns only -----
knz.annual <- knz %>% 
  filter(watershed == "001d",
         pryrdead != "NA") %>% 
  select(recyear, soiltype, transect, plotnum, lvgrass, pryrdead, cuyrdead) %>% 
  mutate(recyear = as.numeric(recyear),
         soiltype = as.factor(soiltype)) %>% 
  rename(year = recyear)

#### 4 year burns only -----
# Remove watersheds with annual and 20 year burns; remove all observations from burn years within 4-year watershed; remove observations with > 100 g/m^2 of woody biomass
knz.four_year <- knz %>% 
  filter(watershed == "004b",
         # pryrdead > 0,
         # lvgrass > 0,
         woody < 100) %>% 
  select(recyear, soiltype, transect, plotnum, lvgrass, pryrdead) %>% 
  mutate(recyear = as.numeric(recyear),
         soiltype = as.factor(soiltype)) %>% 
  rename(year = recyear)

hist(knz.four_year$lvgrass)
hist(log(knz.four_year$lvgrass))

hist(knz.four_year$pryrdead)
hist(log(knz.four_year$pryrdead))


# Create dataframe for years when burning was conducted
# Load burn year data and filter for only four-year watershed
burns.four_year <- read_csv("Datasets/Konza Prairie/KFH011.csv") %>% 
  clean_names() %>% 
  filter(watershed == "4B",
         code != "004A") %>% 
  select(year) %>% 
  distinct(year, .keep_all = TRUE) %>% 
  mutate(burn_year = "yes")

# Create complete set of burned and unburned years
all_years <- tibble(year = seq(1991, max(knz.four_year$year)))

# Full join with all years and replace NAs with "no"
burn_data_filled <- all_years %>%
  left_join(burns.four_year, by = "year") %>%
  mutate(burn_year = replace_na(burn_year, "no"))

# Merge burn data with biomass data
knz.four_year <- knz.four_year %>% 
  left_join(burn_data_filled, by = "year") 

knz.four_year <- knz.four_year %>% 
  mutate(burn_year = as.factor(burn_year))

# Create numerical variable for years since burn
knz.four_year <- knz.four_year %>% 
  arrange(year) %>%
  mutate(last_burn = ifelse(burn_year == "yes", year, NA)) %>%
  fill(last_burn, .direction = "down") %>%  
  mutate(years_since_burn = year - last_burn,
         years_since_burn = as.factor(years_since_burn)) 

# Filter for observations before 2001 and years since burn > 1
knz.four_year.2000 <- knz.four_year %>% 
  filter(year < 2001,
         years_since_burn != "0")

# GLMM of live biomass ~ previous year's dead biomass
knz.four_year.glmm <- glmmTMB(lvgrass ~ pryrdead + years_since_burn + (1|transect/plotnum) + (1 | year),
                              family = "gaussian",
                              data = knz.four_year.2000)
summary(knz.four_year.glmm)

predictions.knz.four_year <- ggpredict(knz.four_year.glmm, terms = ~ pryrdead*years_since_burn)
plot(predictions.knz.four_year)
predictions.knz.four_year <- as.data.frame(predictions.knz.four_year)

# Model visualization
ggplot()+
  geom_point(data = knz.four_year.2000, 
             aes(x = pryrdead, 
                 y = lvgrass),
             alpha = 0.6,
             size = 2) +
  geom_ribbon(data = predictions.knz.four_year, 
              aes(x = x, 
                  y = predicted, 
                  ymin =conf.low, 
                  ymax = conf.high,
                  group = group,
                  fill = group),
              alpha =0.5)+
  geom_line(data = predictions.knz.four_year, 
            aes(x = x , 
                y = predicted,
                group = group,
                color = group)) +
  scale_x_continuous(expand = c(0.01,1)) +
  scale_y_continuous(#limits = c(0,160), 
    expand = c(0.01,1)) +
  labs(x = Previous~year~dead~biomass~(g/m^2),
       y = Current~year~live~biomass~(g/m^2)) +
  #facet_wrap(~group) +
  theme_classic()

# Filter for observations before 2000 and years since burn = 3
knz.four_year.2000.3 <- knz.four_year.2000 %>% 
  filter(year < 2001,
         years_since_burn == "3",
         pryrdead > 0)

# GLMM of live biomass ~ previous year's dead biomass
knz.four_year.3.glmm <- glmmTMB(lvgrass ~ pryrdead + (1|transect),
                              family = "gaussian",
                              data = knz.four_year.2000.3)
summary(knz.four_year.3.glmm)

predictions.knz.four_year <- ggpredict(knz.four_year.glmm, terms = ~ pryrdead)
plot(predictions.knz.four_year)
predictions.knz.four_year <- as.data.frame(predictions.knz.four_year)

# Model visualization
ggplot()+
  geom_point(data = knz.four_year.2000, 
             aes(x = pryrdead, 
                 y = lvgrass),
             alpha = 0.6,
             size = 2) +
  geom_ribbon(data = predictions.knz.four_year, 
              aes(x = x, 
                  y = predicted, 
                  ymin =conf.low, 
                  ymax = conf.high),
              alpha =0.5)+
  geom_line(data = predictions.knz.four_year, 
            aes(x = x , 
                y = predicted,
                group = group)) +
  scale_x_continuous(expand = c(0.01,1)) +
  scale_y_continuous(#limits = c(0,160), 
    expand = c(0.01,1)) +
  labs(x = Previous~year~dead~biomass~(g/m^2),
       y = Current~year~live~biomass~(g/m^2)) +
  #facet_wrap(~years_since_burn) +
  theme_classic()

#### 20 year burns only ----
# Remove watersheds with annual and 20 year burns; remove all observations from burn years within 4-year watershed; remove observations with > 100 g/m^2 of woody biomass
knz.20_year <- knz %>% 
  mutate(woody = replace(woody, is.na(woody), 0)) %>% 
  filter(watershed == "020b") %>% 
  select(recyear, soiltype, transect, plotnum, lvgrass, pryrdead, cuyrdead) %>% 
  mutate(recyear = as.numeric(recyear),
         soiltype = as.factor(soiltype)) %>% 
  rename(year = recyear)

# Create dataframe for years when burning was conducted
# Load burn year data and filter for only four-year watershed
burns.20_year <- read_csv("Datasets/Konza Prairie/KFH011.csv") %>% 
  clean_names() %>% 
  filter(watershed == "20B") %>%
  select(year) %>% 
  distinct(year, .keep_all = TRUE) %>% 
  mutate(burn_year = "yes")

# Create complete set of burned and unburned years
all_years <- tibble(year = seq(1975, max(knz.20_year$year)))

# Full join with all years and replace NAs with "no"
burn_data_filled <- all_years %>%
  left_join(burns.20_year, by = "year") %>%
  mutate(burn_year = replace_na(burn_year, "no"))

# Merge burn data with biomass data
knz.20_year <- knz.20_year %>% 
  left_join(burn_data_filled, by = "year") 

knz.20_year <- knz.20_year %>% 
  mutate(burn_year = as.factor(burn_year))

# Create numerical variable for years since burn
knz.20_year.yrs_since <- knz.20_year %>% 
  arrange(year) %>%
  mutate(last_burn = case_when(year < 1991 ~ 1975,
                               year >= 1991 & year < 2017 ~ 1991,
                               TRUE ~ 2017)) %>%
  mutate(years_since_burn = year - last_burn) %>% 
  filter(lvgrass > 0,
         cuyrdead > 0,
         pryrdead > 0) %>% 
  drop_na()

# Generate breaks in years since burn 
breaks_seq <- seq(0, max(knz.20_year.yrs_since$years_since_burn, na.rm = TRUE) + 5, by = 5)

# Create matching labels
labels_seq <- paste(
  head(breaks_seq, -1),
  tail(breaks_seq, -1) - 1,
  sep = "-"
)

# Apply cut with matching breaks and labels
knz.20_year.yrs_since <- knz.20_year.yrs_since %>%
  mutate(
    burn_interval = cut(
      years_since_burn,
      breaks = breaks_seq,
      labels = labels_seq,
      right = FALSE,
      include.lowest = TRUE)
  ) %>% 
  filter(years_since_burn > 0,
         year < 2001,
         !burn_interval %in% c("20-24", "25-29"))

# GLMM of live grass biomass ~ previous year dead biomass, interaction with burn interval
knz_20yr.glmm.gamma <- glmmTMB(
  lvgrass ~ pryrdead*burn_interval + (1 | transect/plotnum/year),
  data = knz.20_year.yrs_since,
  family = Gamma(link = "log")
)
summary(knz_20yr.glmm.gamma)

## Vizualization
# Step 1: Get model predictions from ggeffects
preds_gamma <- ggpredict(knz_20yr.glmm.gamma, terms = c("pryrdead", "burn_interval"))

# Step 2: Plot with ggplot
ggplot(preds_gamma, aes(x = x, y = predicted, color = group)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group),
              alpha = 0.2, color = NA) +
  labs(
    x = "Previous Year Dead Biomass (g/m²)",
    y = "Predicted Live Biomass (g/m²)",
    color = "Burn Interval (yrs)",
    fill = "Burn Interval (yrs)"
  ) +
  theme_minimal(base_size = 14)


## Z-score standardized model
knz.20_year.yrs_since.z <- knz.20_year.yrs_since %>%
  group_by(burn_interval) %>%
  mutate(
    lvgrass_z = scale(lvgrass)[,1],
    pryrdead_z = scale(pryrdead)[,1]
  ) %>%
  ungroup()

# Build separate z-score models for each burn category
# Create an empty list to store models
z_models <- list()

# Get unique burn interval levels
burn_levels <- unique(knz.20_year.yrs_since.z$burn_interval)

# Fit one model per interval
for (b in burn_levels) {
  dat <- knz.20_year.yrs_since.z %>% filter(burn_interval == b)
  z_models[[as.character(b)]] <- glmmTMB(
    lvgrass_z ~ pryrdead_z + (1 | transect/plotnum/year),
    data = dat,
    family = gaussian()
  )
}

# Extract coefficients for each model
results <- lapply(names(z_models), function(name) {
  model <- z_models[[name]]
  coef <- summary(model)$coefficients$cond["pryrdead_z", ]
  data.frame(
    burn_interval = name,
    estimate = coef["Estimate"],
    std_error = coef["Std. Error"],
    p_value = coef["Pr(>|z|)"]
  )
}) %>%
  bind_rows()
