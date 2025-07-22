#### Georgia Coastal Ecosystems LTER ####
# Analysis of marshgrass productivity in the presence and absence of marshwrack disturbance

# Packages
library(tidyverse)
library(janitor)
library(glmmTMB)
library(ggeffects)
library(DHARMa)
library(broom.mixed)

# Working directory
setwd('/Users/kako4300/Library/CloudStorage/OneDrive-UCB-O365/Projects/LTER material legacy synthesis')

# Load and clean marshgrass biomass data
gce.biomass <- read_csv("Datasets/Georgia Coastal Ecosystems/GCE_Biomass_2000-2024.csv") %>% 
  clean_names()

gce.biomass_dist <- gce.biomass %>% 
  filter(quadrat_area == 0.25,
         total_plant_biomass_m2 != "NA") %>% 
  mutate(plot_disturbance = if_else(plot_disturbance == 1, "Yes", "No"),
         plot_disturbance = as.factor(plot_disturbance),
         plot_disturbance = ordered(plot_disturbance, levels = c("No", "Yes"))) %>% 
  select(c(1,2,4,7,16)) 

# Summarize data and visualize with exploratory graph
gce.biomass_dist.summary <- gce.biomass_dist %>% 
  group_by(site, year, plot_disturbance) %>% 
  summarize(biomass.mean = mean(total_plant_biomass_m2)) %>% #,
            #biomass.se = sd(total_plant_biomass_m2)/sqrt(n()))
  ungroup()

# Exploratory graph
ggplot(gce.biomass_dist.summary, aes(x = plot_disturbance, y = biomass.mean)) +
  geom_col() +
  geom_errorbar(aes(ymax = biomass.mean + biomass.se,
                    ymin = biomass.mean - biomass.se),
                width = 0) +
  theme_classic()

## Analysis
marsh_glmm.raw <- glmmTMB(
  biomass.mean ~ plot_disturbance + (1 | site) + (1 | year),
  data = gce.biomass_dist.summary,
  family = gaussian(link = "log")
)
summary(marsh_glmm.raw)

# Diagnostics
res <- simulateResiduals(marsh_glmm.raw)
plot(res) # Tests are non-significant

# To report effect size and CI in main text and table 
gce_effect.raw <- tidy(marsh_glmm.raw, effects = "fixed", conf.int = TRUE, conf.level = 0.95) %>% 
  filter(term == "plot_disturbance.L")

write_csv(gce_effect.raw, "Datasets/Effect sizes/Raw/gce_effect.raw.csv")

## Visualization
preds <- ggpredict(marsh_glmm.raw, terms = "plot_disturbance")

ggplot() +
  geom_jitter(data = gce.biomass_dist.summary,
              aes(x = plot_disturbance, y = biomass.mean),
              width = 0.15, 
              alpha = 0.6,
              color = "darkgrey") +
  geom_point(data = preds,
             aes(x = x, y = predicted),
             size = 3) +
  geom_errorbar(data = preds,
                aes(x = x, ymin = conf.low, ymax = conf.high),
                width = 0) +
  geom_line(data = preds,
            aes(x = x, y = predicted, group = group)) +
  labs(x = "Marshwrack disturbance",
    y = expression("Marshgrass biomass (g/m"^2*"/yr)")) +
  theme_classic(base_size = 14)

## Z-score standardization
# Standardize  response only (since predictor is binary and categorical)
gce.biomass_dist_z <- gce.biomass_dist.summary %>%
  mutate(biomass_z = scale(biomass.mean)[, 1])

marsh_glmm.z <- glmmTMB(
  biomass_z ~ plot_disturbance + (1 | site) + (1 | year),
  data = gce.biomass_dist_z,
  family = gaussian()
)
summary(marsh_glmm.z)

# Diagnostics
res.z <- simulateResiduals(marsh_glmm.z)
plot(res.z) # Tests are non-significant

## Visualization
preds.z <- ggpredict(marsh_glmm.z, terms = "plot_disturbance")

ggplot() +
  geom_jitter(data = gce.biomass_dist_z,
              aes(x = plot_disturbance, y = biomass_z),
              width = 0.15, 
              alpha = 0.6,
              color = "darkgrey") +
  geom_point(data = preds.z,
             aes(x = x, y = predicted),
             size = 3) +
  geom_errorbar(data = preds.z,
                aes(x = x, ymin = conf.low, ymax = conf.high),
                width = 0) +
  geom_line(data = preds.z,
            aes(x = x, y = predicted, group = group)) +
  labs(x = "Marshwrack disturbance",
       y = "Marshgrass biomass (Z-score)") +
  theme_classic(base_size = 14)

## Extract effect size
gce_effect.z <- tidy(marsh_model.z, effects = "fixed", conf.int = TRUE) %>% 
  filter(term == "plot_disturbance.L")

print(gce_effect.z)

write_csv(gce_effect.z, "Datasets/Effect sizes/Standardized/gce_effect.z.csv")
