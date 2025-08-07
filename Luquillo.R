## Luquillo LTER ##
# Analysis of the effect of litterfall on seedling survival using Canopy Trimming Experiment data

# Packages
library(tidyverse)
library(janitor)
library(glmmTMB)
library(ggeffects)
library(DHARMa)
library(broom.mixed)

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

## Seedlings < 10cm 
# Create dataframe for seedlings < 10cm in years postimanipulation and pre-hurricane
seedling_lessthan10 <- luq.seedlings %>%   
  filter(lessthan10cm != "NA",
         treatment %in% c("Trim&clear", "Trim+Debris"),
         year %in% c(2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013)) %>% 
  group_by(treatment, year, block, plot) %>% 
  summarize(seedling.count = sum(lessthan10cm)) %>% 
  mutate(treatment = case_when(treatment == "Trim&clear" ~"Removed", 
                               TRUE ~ "Added"),
         treatment = factor(treatment),
         block = factor(block),
         plot = factor(plot),
         year = factor(year))


# Analysis
seedling_glmm.raw <- glmmTMB(
  seedling.count ~ treatment + (1 | block/plot) + (1 | year),
  family = nbinom2,
  data = seedling_lessthan10
)
summary(seedling_glmm.raw)

# Diagnostics
res <- simulateResiduals(seedling_glmm.raw)
plot(res) # Tests are non-significant

# Extract effect size and CI for main text and table
luq_effect.raw <- tidy(seedling_glmm.raw, effects = "fixed", conf.int = TRUE, conf.level = 0.95) %>% 
  filter(term == "treatmentLitter added") 

write_csv(luq_effect.raw, "Datasets/Effect sizes/Raw/luq_effect.raw.csv")

## Visualization
# Get predictions
lessthan10_preds <- ggpredict(seedling_glmm.raw, terms = "treatment")

# Plot
ggplot() +
  geom_jitter(data = seedling_lessthan10,
              aes(x = treatment, y = seedling.count),
              color = "#6B6C58", 
              alpha = 0.6, 
              width = 0.15) +
  geom_line(data = lessthan10_preds, 
            aes(x = x, y = predicted,
                group = group)) +
  geom_errorbar(data = lessthan10_preds, 
                aes(x = x, y = predicted, 
                    ymin = conf.low, ymax = conf.high), 
                width = 0) +
  geom_point(data = lessthan10_preds, 
             aes(x = x, y = predicted),
             size = 4,
             color = "black") +
  geom_point(data = lessthan10_preds, 
             aes(x = x, y = predicted),
             size = 3,
             color = "#6B6C58") +
  labs(x = "Canopy debris",
       y = "Seedling (< 10cm) counts/yr") +
  theme_classic(base_size = 14)

## Z-score < 10cm seedling model
# Standardize seedling count
seedling_lessthan10$seedlings.z <- scale(seedling_lessthan10$seedling.count)

# Refit the model using standardized response
seedling_glmm.z <- glmmTMB(
  seedlings.z ~ treatment + (1 | block/plot) + (1 | year),
  family = gaussian(),
  data = seedling_lessthan10
)
summary(seedling_lessthan10_glmm.z)

# Diagnostics
res.z <- simulateResiduals(seedling_glmm.z)
plot(res.z) # Tests are non-significant

# Get predictions for Z-score model
lessthan10_preds.z <- ggpredict(seedling_glmm.z, terms = "treatment")

# Plot
ggplot() +
  geom_jitter(data = seedling_lessthan10,
              aes(x = treatment, y = seedlings.z),
              color = "darkgrey", 
              alpha = 0.6, 
              width = 0.15) +
  geom_point(data = lessthan10_preds.z, 
             aes(x = x, y = predicted),
             size = 3) +
  geom_errorbar(data = lessthan10_preds.z, 
                aes(x = x, y = predicted, 
                    ymin = conf.low, ymax = conf.high), 
                width = 0) +
  geom_line(data = lessthan10_preds.z, 
            aes(x = x, y = predicted,
                group = group)) +
  labs(x = "Litterfall treatment",
       y = "Seedling counts (< 10cm) (Z-score)") +
  theme_classic(base_size = 14)


## Effect size
luq_effect.z <- tidy(seedling_glmm.z, effects = "fixed", conf.int = TRUE) %>%
  filter(term == "treatmentLitter added")

print(luq_effect.z)

write_csv(luq_effect.z, "Datasets/Effect sizes/Standardized/luq_effect.z.csv")
