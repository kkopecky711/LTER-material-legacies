#### Florida Coastal Everglades ####
# Analysis of hurricane sediment deposition on mangrove seedling development

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

# Filter for Wilma and Irma periods
litterfall.post_hurr <- litterfall %>%
  filter(
    (date > as.Date("2005-10-27") & date < as.Date("2006-02-01")) |
      (date > as.Date("2017-09-17") & date < as.Date("2019-03-01"))) %>% 
  mutate(year = substr(date, 1, 4),
         hurricane = if_else(year > 2006, "Irma", "Wilma")) %>% 
  group_by(sitename, hurricane, plot_id, basket_id) %>% 
  summarize(total_weight = sum(total_weight)) # Sum up weight of litter accumulated over time in each basket

litterfall.summary <- litterfall.post_hurr %>% 
  group_by(sitename, hurricane, plot_id) %>% 
  summarize(mean_litter = mean(total_weight), # Average the amount of litter accumulated within each plot
            se_litter = sd(total_weight)/sqrt(n())) %>% 
  filter(sitename != "TS/Ph8")

#### Root production ----
# Load in datasets for root production after hurricanes, then bind
root_production.wilma <- read_csv("Datasets/Florida Coastal Everglades/FCE1277_Root_Production_pre-hurricane.csv") %>% 
  clean_names() %>% 
  filter(year == "2006") %>% # filter for year after hurricane and only top later 
  mutate(hurricane = "Wilma") %>% 
  select(c("sitename", "hurricane", "month", "year", "point", "position", "location", "core_id", "root_size", "root_size_class", "root_production"))

root_production.irma <- read_csv("Datasets/Florida Coastal Everglades/FCE_1278_Root_Production_post-Irma.csv") %>%
  clean_names() %>% 
  mutate(hurricane = "Irma") %>% 
  select(c("sitename", "hurricane", "month", "year", "point", "position", "location", "core_id", "root_size", "root_size_class", "root_production"))

root_production <- rbind(root_production.wilma, root_production.irma)

# Summarize by plot
root_prod.summary <- root_production %>%
  filter(root_size_class == "Fine",
         sitename %in% c("SRS4", "SRS5", "SRS6")) %>% 
  mutate(hurricane = if_else(year > 2006, "Irma", "Wilma"),
         plot_id = if_else(point %in% c("A", "B"), 1, 2)) %>% 
  group_by(sitename, hurricane, plot_id) %>% 
  summarize(root_prod.mean = mean(root_production),
            root_prod.se = sd(root_production)/sqrt(n()))

## Merge litterfall and root production dataframes
root_prod.litter <- root_prod.summary %>% 
  left_join(litterfall.summary, by = c("sitename", "hurricane", "plot_id"))

# Exploratory graph of root production ~ litterfall
ggplot(root_prod.litter, aes(x = mean_litter, y = root_prod.mean)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "Mean litter biomass",
       y = "Mean fine-root production") +
  theme_classic() +
  facet_wrap(~hurricane)

## Analysis
# GLMM of mean fine root production ~ mean litterfall

root_prod.glmm <- glmmTMB(
  root_prod.mean ~ mean_litter + (1 | plot_id) + (1 | sitename),
  data = root_prod.litter,
  family = gaussian()
)
summary(root_prod.glmm)

# To report effect size and CI in main text
tidy(root_prod.glmm, effects = "fixed", conf.int = TRUE, conf.level = 0.95)%>% 
  filter(term == "mean_litter")

# Get model predictions and confidence intervals
preds <- ggpredict(root_prod.glmm, terms = "mean_litter")

# Plot
ggplot() +
  geom_ribbon(data = preds, aes(x = x, ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  geom_line(data = preds, aes(x = x, y = predicted)) +
  geom_point(data = root_prod.litter, 
             aes(x = mean_litter, y = root_prod.mean), 
             size = 2,
             alpha = 0.6) +
  labs(
    x = "Mean litterfall (g)",
    y = "Mean root production (g/m²/yr)") +
  theme_classic(base_size = 14)

# Extract coefficient and standard error
coef_est <- fixef(root_prod.glmm)$cond["mean_litter"]
coef_se  <- summary(root_prod.glmm)$coefficients$cond["mean_litter", "Std. Error"]

# Standard deviations of x and y
sd_x <- sd(root_prod.litter$mean_litter)
sd_y <- sd(root_prod.litter$root_prod.mean)

# Standardized effect and SE
beta_std <- coef_est * (sd_x / sd_y)
beta_std_se <- coef_se * (sd_x / sd_y)

# 95% CI (normal approximation)
lower_ci <- beta_std - 1.96 * beta_std_se
upper_ci <- beta_std + 1.96 * beta_std_se

library(tibble)

# Calculate z and p
z_stat <- beta_std / beta_std_se
p_val <- 2 * (1 - pnorm(abs(z_stat)))

# Format in same structure
mangrove_effect <- tibble(
  effect = "fixed",
  component = "cond",
  term = "mean_litter_z",  # mimic naming convention for z-scored terms
  estimate = beta_std,
  std.error = beta_std_se,
  statistic = z_stat,
  p.value = p_val,
  conf.low = lower_ci,
  conf.high = upper_ci
)

write_csv(mangrove_effect, "Datasets/Effect sizes/fce.effect_size.csv")



# Z-score standardized model
root_prod.scaled <- root_prod.litter %>%
  mutate(
    root_z = scale(root_prod.mean)[,1],
    litter_z = scale(mean_litter)[,1]
  )

root_prod.glmm.z <- glmmTMB(
  root_z ~ litter_z + (1 | plot_id) + (1 | sitename),
  data = root_prod.scaled,
  family = gaussian()
) # Model does not converge with any iteration of random effects structure (including no random effects)

# Try simpe linear model
root_prod.lm.z <- lm(root_z ~ litter_z, 
                     data = root_prod.scaled)
summary(root_prod.lm.z)

cor(root_prod.scaled$root_z, root_prod.scaled$litter_z, use = "complete.obs")


#### Root Biomass ----
# Wilma
root_biomass.wilma <- read_csv("Datasets/Florida Coastal Everglades/FCE1277_Root_Biomass_pre-hurricane.csv") %>% clean_names()

# Read in root biomass data, post-Irma
root_biomass.irma <- read_csv("Datasets/Florida Coastal Everglades/FCE_1278_Root_Biomass_post-Irma.csv") %>%
  clean_names()

# Summarize by plot
root_biomass.summary <- root_biomass.irma %>%
  mutate(hurricane = "Irma") %>% 
  filter(root_size_class == "Fine",
         sitename %in% c("SRS4", "SRS5", "SRS6")) %>% 
  group_by(sitename, plot_id) %>% 
  summarize(mean_biomass = mean(root_biomass),
            se_biomass = sd(root_biomass)/sqrt(n()))

# Filter litter dataframe for Irma only
litterfall.irma <- litterfall.summary %>% 
  ungroup() %>% 
  filter(hurricane == "Irma") %>% 
  select(-hurricane)

# Merge litter and root biomass dataframes
roots_biomass.litter.irma <- root_biomass.summary %>% 
  left_join(litterfall.irma, by = c("sitename", "plot_id"))

ggplot(roots_biomass.litter.irma, aes(x = mean_litter, y = mean_biomass)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "Mean litter biomass",
       y = "Mean root biomass",
       title = "Post-Irma") +
  theme_classic()





# Transects run from coast inland; 10m point should line up with the corners of the plots



# Look at average plot-level red mangrove root biomass ~ plot-level average litterfall (total or just red mangrove)
# Could also try max litterfall
# Fine root biomass most sensitive 
# Be careful not to use litterfall data that go past when root data were collected


