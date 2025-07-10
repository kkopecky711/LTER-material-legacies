#### Luquillo LTER ####
# Analysis of the effect of litterfall on seedling/sapling survival using Canopy Trimming Experiment data

# Packages
library(tidyverse)
library(janitor)
library(glmmTMB)
library(ggeffects)
library(DHARMa)

# Working directory
setwd('/Users/kako4300/Library/CloudStorage/OneDrive-UCB-O365/Projects/LTER material legacy synthesis')

# Load and clean seedling data
luq.seedlings <- read_csv("Datasets/Luquillo/CTESeedlingMeasurementData2003-20212.csv") %>% 
  clean_names() %>% 
  mutate(start_date = as.character(start_date),
         year = as.numeric(substr(start_date, 1, 4)))

# Treatment metadata
luq.cte_ttt <- read_csv("Datasets/Luquillo/CTE_Treatments_0.csv") %>% 
  clean_names() %>% 
  select(-c(latitude, longitude, subplot)) %>% 
  distinct()

# Merge with treatment metadata
luq.seedlings <- merge(luq.seedlings, luq.cte_ttt)

luq.seedlings <- luq.seedlings %>% 
  relocate(treatment, .before = block) %>% 
  select(-c(comments, new))

#### Seedlings < 10cm -----
seedling_lessthan10 <- luq.seedlings %>%   
  filter(lessthan10cm != "NA",
         treatment %in% c("Trim&clear", "Trim+Debris"),
         year %in% c(2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013)) %>% 
  group_by(treatment, year, block, plot) %>% 
  summarize(seedling.count = sum(lessthan10cm)) %>% 
  mutate(treatment = case_when(treatment == "Trim&clear" ~"Litter removed", 
                               TRUE ~ "Litter added"),
         treatment = factor(treatment),
         block = factor(block),
         plot = factor(plot),
         year = factor(year))

seedling_lessthan10_glmm <- glmmTMB(
  seedling.count ~ treatment + (1 | block/plot) + (1 | year),
  family = nbinom2,
  data = seedling_lessthan10
)
summary(seedling_lessthan10_glmm)
tidy(seedling_lessthan10_glmm, effects = "fixed", conf.int = TRUE, conf.level = 0.95)

## Visualization
# Get predictions
lessthan10_preds <- ggpredict(seedling_lessthan10_glmm, terms = "treatment")

# Plot
ggplot(lessthan10_preds, aes(x = x, y = predicted)) +
  geom_jitter(data = seedling_lessthan10,
              aes(x = treatment, y = seedling.count),
              width = 0.1, 
              height = 0,
              alpha = 0.6) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
                width = 0) +
  geom_line(aes(group = group)) +
  labs(x = "Litterfall treatment",
       y = "Seedling (< 10cm) counts") +
  theme_classic(base_size = 14)

## Z-score < 10cm seedling model
# Standardize seedling count
seedling_lessthan10$seedlings.z <- scale(seedling_lessthan10$seedling.count)

# Refit the model using standardized response
seedling_lessthan10_glmm.z <- glmmTMB(
  seedlings.z ~ treatment + (1 | block/plot) + (1 | year),
  family = gaussian,
  data = seedling_lessthan10
)
summary(seedling_lessthan10_glmm.z)

# Get predictions for Z-score model
lessthan10_preds.z <- ggpredict(seedling_lessthan10_glmm.z, terms = "treatment")

# Plot
ggplot(lessthan10_preds.z, aes(x = x, y = predicted)) +
  geom_jitter(data = seedling_lessthan10,
              aes(x = treatment, y = seedlings.z),
              width = 0.1, 
              height = 0,
              alpha = 0.4) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
                width = 0) +
  labs(x = "Litterfall Treatment",
       y = "Seedlings < 10cm (Z-score)") +
  theme_classic(base_size = 14)

## Effect size
tropical_effect <- tidy(seedling_lessthan10_glmm.z, effects = "fixed", conf.int = TRUE) %>%
  filter(term == "treatmentLitter added")

print(tropical_effect)

write_csv(tropical_effect, "Datasets/Effect sizes/luq.effect_size.csv")



#### Survival and mortality of seedlings < 10cm ----
# Create a summary table that sums up mortality and survival counts to the plot level
# Add columns for seedling survival and mortality as a binary variable
luq.seedlings_mortality <- luq.seedlings %>% 
  filter(#year %in% c(2004, 2005, 2006, 2007, 2008),
    #year == 2004 | year == 2005,
    dead_alive_notfound != "NF") %>% 
  mutate(mortality = if_else(dead_alive_notfound == "D", 1, 0),
         survival = if_else(dead_alive_notfound == "A", 1, 0),
         year = as.numeric(year))

luq.mortality_counts <- luq.seedlings_mortality %>%   
  group_by(treatment, year, block, plot) %>% 
  summarize(mortality.count = sum(mortality),
            survival.count = sum(survival)) %>% 
  mutate(treatment = factor(treatment),
         block = factor(block),
         plot = factor(plot),
         year = factor(year))

# Remove first year of experiment (pre-manipulation) and filter for only Control and trim + debris treatments
luq.seeding_counts <- luq.mortality_counts %>% 
  filter(year %in% c(2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013),
         treatment == "Trim&clear" | treatment == "Trim+Debris")

## Analysis
# Survival
# Fit GLMM for seedling survival
seedling_survival_glmm <- glmmTMB(
  survival.count ~ treatment + (1 | block/plot) + (1 | year),
  family = nbinom2,
  data = luq.seeding_counts
)
summary(seedling_survival_glmm)

## Visualization
# Get predictions
surv_preds <- ggpredict(seedling_survival_glmm, terms = "treatment")

# Plot
ggplot(surv_preds, aes(x = x, y = predicted)) +
  geom_jitter(data = luq.seeding_counts,
              aes(x = treatment, y = survival.count),
              width = 0.1, 
              height = 0,
              alpha = 0.6, 
              color = "darkgrey") +
  geom_point(color = "black",
             size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
                width = 0) +
  labs(x = "Litterfall Treatment",
       y = "Seedling counts") +
  theme_classic(base_size = 14)

## Z-score model
# Standardize survival count
luq.seeding_counts$survival.z <- scale(luq.seeding_counts$survival.count)

# Refit the model using standardized response
seedling_survival_glmm.z <- glmmTMB(
  survival.z ~ treatment + (1 | block/plot/year),
  family = gaussian,
  data = luq.seeding_counts
)
summary(seedling_survival_glmm.z)

# Get predictions for Z-score model
surv_preds.z <- ggpredict(seedling_survival_glmm.z, terms = "treatment")

# Plot
ggplot(surv_preds.z, aes(x = x, y = predicted)) +
  geom_jitter(data = luq.seeding_counts,
              aes(x = treatment, y = survival.z),
              width = 0.1, 
              height = 0,
              alpha = 0.6, 
              color = "darkgrey") +
  geom_point(color = "black",
             size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
                width = 0) +
  labs(x = "Litterfall Treatment",
       y = "Seedling counts (Z-score)") +
  theme_classic(base_size = 14)

## Effect size
tropical_effect <- tidy(seedling_survival_glmm.z, effects = "fixed", conf.int = TRUE) %>%
  filter(term == "treatmentTrim+Debris")

print(tropical_effect)

write_csv(tropical_effect, "Datasets/Effect sizes/luq.effect_size.csv")

#### Litter as a continuous variable ----
# Litterfall data
luq.litter <- read_csv("Datasets/Luquillo/CTE-Litterfall-data.csv") %>% 
  clean_names() %>% 
  select(-c(8:10)) %>% 
  mutate(date = as.character(date),
         year = as.numeric(substr(date, 1, 4)),
         litter_mass = leaves_in_g + wood_in_g + misc_in_g + fruits_seeds_flower_g) 


luq.litter <- merge(luq.litter, luq.cte_ttt)

luq.litter.first_cte <- luq.litter %>% 
  filter(year %in% c(2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013),
         treatment == "Trim&clear" | treatment == "Trim+Debris") %>% 
  group_by(treatment, year, block, plot) %>% 
  summarize(total_biomass = sum(litter_mass))

seedlings_litter <- merge(luq.seeding_counts, luq.litter.first_cte)
seedlings_litter <- seedlings_litter %>% 
  drop_na()

# Range of litter across treatments
ggplot(seedlings_litter, aes(x = total_biomass, fill = treatment)) +
  geom_histogram(position = "dodge",
                 binwidth = 100, color = "black")

# Seedling count ~ litter biomass
seedling_litter_glmm <- glmmTMB(
  mortality.count ~ total_biomass,
  family = nbinom2,
  data = seedlings_litter
)
summary(seedling_litter_glmm)

summary(seedlings_litter$total_biomass)
sd(seedlings_litter$total_biomass)

## Visualization
# Get predictions
mort_preds <- ggpredict(seedling_litter_glmm, terms = "total_biomass")

# Plot
ggplot(mort_preds, aes(x = x, y = predicted)) +
  geom_point(data = seedlings_litter,
              aes(x = total_biomass, y = mortality.count),
              alpha = 0.6) +
  geom_line(color = "black") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high),
              alpha = 0.5) +
  labs(x = "Biomass of litterfall",
       y = "Seedling mortality (count)") +
  theme_classic(base_size = 14)


