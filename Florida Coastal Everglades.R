## Florida Coastal Everglades ##
# Analysis of hurricane litterfall deposition on mangrove root production

# Packages
library(tidyverse)
library(janitor)
library(glmmTMB)
library(ggeffects)
library(DHARMa)
library(broom.mixed)

# Working directory
setwd('/Users/kako4300/Library/CloudStorage/OneDrive-UCB-O365/Projects/LTER material legacy synthesis')

## Litterfall data
litterfall <- read_csv("Datasets/Florida Coastal Everglades/LT_PP_Castaneda_001.csv") %>% 
  clean_names()

# Filter for post-Irma observations leading up to date of root production measurement
litterfall.post_irma <- litterfall %>%
  filter((date > as.Date("2017-09-17") & date < as.Date("2019-03-01"))) %>% 
  group_by(sitename, plot_id, basket_id) %>% 
  summarize(total_weight = sum(total_weight)) # Sum up weight of litter accumulated over time in each basket

litterfall.summary <- litterfall.post_irma %>% 
  group_by(sitename, plot_id) %>% 
  summarize(mean_litter = mean(total_weight), # Average the amount of litter accumulated within each plot
            se_litter = sd(total_weight)/sqrt(n()))

## Root production data
root_production.irma <- read_csv("Datasets/Florida Coastal Everglades/FCE_1278_Root_Production_post-Irma.csv") %>%
  clean_names() %>%
  select(c("sitename", "month", "year", "point", "position", "location", "core_id", "root_size", "root_size_class", "root_production"))

# Summarize by plot
root_prod.summary <- root_production.irma %>%
  filter(root_size_class == "Fine",
         root_production > 0) %>% 
  mutate(plot_id = if_else(point %in% c("A", "B"), 1, 2)) %>% 
  group_by(sitename, plot_id) %>% 
  summarize(root_prod.mean = mean(root_production),
            root_prod.se = sd(root_production)/sqrt(n())) %>% 
  ungroup()

## Merge litterfall and root production dataframes
root_prod.litter <- root_prod.summary %>% 
  left_join(litterfall.summary, by = c("sitename", "plot_id"))

# Exploratory graph of root production ~ litterfall
ggplot(root_prod.litter, aes(x = mean_litter, y = root_prod.mean)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "Mean litter biomass",
       y = "Mean fine-root production") +
  theme_classic() 

## Analysis
# GLMM of mean fine root production ~ mean litterfall
root_glmm.raw <- glmmTMB(
  root_prod.mean ~ mean_litter + (1 | sitename),
  data = root_prod.litter,
  family = gaussian()
)
summary(root_glmm.raw)

# Diagnostics
res <- simulateResiduals(root_glmm.raw)
plot(res) # Tests are non-significant

# To report effect size and CI in main text
fce_effect.raw <- tidy(root_glmm.raw, effects = "fixed", conf.int = TRUE, conf.level = 0.95)%>% 
  filter(term == "mean_litter")

write_csv(fce_effect.raw, "Datasets/Effect sizes/Raw/fce_effect.raw.csv")

## Visualization
# Get model predictions and confidence intervals
preds <- ggpredict(root_glmm.raw, terms = "mean_litter")

# Plot
ggplot() +
  geom_point(data = root_prod.litter, 
             aes(x = mean_litter, y = root_prod.mean), 
             color = "darkgrey",
             alpha = 0.6) +
  geom_ribbon(data = preds, 
              aes(x = x, ymin = conf.low, ymax = conf.high), 
              alpha = 0.2) +
  geom_line(data = preds, 
            aes(x = x, y = predicted),
            linewidth = 1.2) +
  labs(
    x = "Mean litterfall (g)",
    y = "Mean root production (g/m²/yr)") +
  theme_classic(base_size = 14)

## Standardized model
# Scale both predictor and response to their z-scores
root_prod.litter <- root_prod.litter %>%
  mutate(
    root.z = scale(root_prod.mean)[,1],
    litter.z = scale(mean_litter)[,1]
  )

# Model relationship with same structure
root_glmm.z <- glmmTMB(
  root.z ~ litter.z + (1 | sitename),
  data = root_prod.litter,
  family = gaussian()
)
summary(root_glmm.z)

# Diagnostics
res.z <- simulateResiduals(root_glmm.z)
plot(res.z) # Tests are non-significant


## Visualization
# Get model predictions and confidence intervals
preds.z <- ggpredict(root_glmm.z, terms = "litter.z")

# Plot
ggplot() +
  geom_point(data = root_prod.litter, 
             aes(x = litter.z, y = root.z), 
             color = "darkgrey",
             alpha = 0.6) +
  geom_ribbon(data = preds.z, 
              aes(x = x, ymin = conf.low, ymax = conf.high), 
              alpha = 0.2) +
  geom_line(data = preds.z, 
            aes(x = x, y = predicted),
            linewidth = 1.2) +
  labs(
    x = "Mean litterfall (Z-score)",
    y = "Mean root production (Z-score)") +
  theme_classic(base_size = 14)

# Extract effect size and CI, write to .csv file
fce_effect.z <- tidy(root_glmm.z, effects = "fixed", conf.int = TRUE) %>% 
  filter(term == "litter_z")

write_csv(fce_effect.z, "Datasets/Effect sizes/Standardized/fce_effect.z.csv")

