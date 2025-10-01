#### H.J. Andrews Forest ####
# Analysis of tree growth, mortality, and ingrowth as a function of standing and downed dead wood

# Working directory
setwd('/Users/kako4300/Library/CloudStorage/OneDrive-UCB-O365/Projects/LTER material legacy synthesis')

# Packages
library(tidyverse)
library(janitor)
library(glmmTMB)
library(ggeffects)
library(DHARMa)
library(broom.mixed)

# Dead tree data
hja_dead <- read_csv("Datasets/Andrews Forest/OHJA_downed wood summary_v2.csv") %>% 
  clean_names() %>% 
  mutate(obs_id = paste(stand, plot, year, sep = "_")) %>% 
  rename(cwd_year = year)

## Douglas fir only ----
# Live tree data for Douglas fir (PSME)
hja_live <- read_csv("Datasets/Andrews Forest/PSP_Plot_Change_20year_PSME.csv") %>% 
  clean_names() %>% 
  mutate(obs_id = paste(stand, plot, cwd_year, sep = "_"))

# Now merge using stand, plot, and cwd_year
live_dead <- inner_join(hja_live, hja_dead,
                          by = c("stand", "plot", "cwd_year"))

# Calculate dead wood per hectare and individual tree growth rates, remove unneeded columns
live_dead.growth <- live_dead %>% 
  mutate(dw_mass_ha = total_mass/area_ha,
         dw_vol_ha = total_volume/area_ha,
         dw_cover_ha = total_cover/area_ha,
         dw_area_ha = total_area/area_ha,
         tree_growth_ind = growth_baph_spp / (tph0_spp * surv_prop_spp) / d_year) %>% 
  select(c(stand, plot, cwd_year, dw_mass_ha, dw_vol_ha, dw_cover_ha, dw_area_ha, tree_growth_ind))

## Analysis
tree_glmm.raw <- glmmTMB(
  tree_growth_ind ~ dw_mass_ha + (1 | stand),
  data = live_dead.growth,
  family = gaussian()
)
summary(tree_glmm.raw)

# Diagnostics
res <- simulateResiduals(tree_glmm.raw)
plot(res) # Tests are non-significant

# To report effect size and CI in main text and table
hja_effect.raw <- tidy(tree_glmm.raw, effects = "fixed", conf.int = TRUE, conf.level = 0.95) %>% 
  filter(term == "dw_mass_ha")

write_csv(hja_effect.raw, "Datasets/Effect sizes/Raw/hja_effect.raw.csv")

## Visualization
# Generate datatable of model predictions
preds <- ggpredict(tree_glmm.raw, terms = "dw_mass_ha")

# Plot predicted values with 95% confidence intervals and raw values
ggplot() +
  geom_point(data = live_dead.growth, 
             aes(x = dw_mass_ha, y = tree_growth_ind),
             alpha = 0.6,
             color = "#6B6C58") +
  geom_line(data = preds, 
            aes(x = x, y = predicted),
            linewidth = 0.75) +
  geom_ribbon(data = preds,
              aes(x = x, y = predicted, 
                  ymin = conf.low, ymax = conf.high), 
              alpha = 0.3,
              fill = "#6B6C58") +
  labs(x = "Dead wood mass (kg/ha)",
       y = "Tree growth (Douglas fir; m²/ind./yr)") +
  theme_classic(base_size = 14)

## Standardized model
# Scale predictor and response only to their Z-scores
live_dead.growth <- live_dead.growth %>%
  mutate(dw_mass_ha.z = scale(dw_mass_ha)[, 1],
         tree_growth_ind.z = scale(tree_growth_ind)[, 1])

tree_glmm.z <- glmmTMB(
  tree_growth_ind.z ~ dw_mass_ha.z + (1 | stand),
  data = live_dead.growth,
  family = gaussian()
)
summary(tree_glmm.z)

# Diagnostics
res.z <- simulateResiduals(tree_glmm.z)
plot(res.z) # Tests are non-significant

## Visualization
# Generate datatable of model predictions
preds.z <- ggpredict(tree_glmm.z, terms = "dw_mass_ha.z")

# Plot predicted values with 95% confidence intervals and raw values
ggplot() +
  geom_point(data = live_dead.growth, 
             aes(x = dw_mass_ha.z, y = tree_growth_ind.z),
             alpha = 0.6,
             color = "darkgrey") +
  geom_line(data = preds.z, 
            aes(x = x, y = predicted),
            linewidth = 1.2) +
  geom_ribbon(data = preds.z,
              aes(x = x, y = predicted, 
                  ymin = conf.low, ymax = conf.high), 
              alpha = 0.3) +
  scale_y_continuous(limits = c(-2, 3.4)) +
  labs(x = "Dead wood mass (Z-score)",
       y = "Individual tree growth (Z-score)") +
  theme_classic(base_size = 14)

# Extract effect size and CI, write to .csv file
hja_effect.z <- tidy(tree_glmm.z, effects = "fixed", conf.int = TRUE) %>% 
  filter(term == "dw_mass_ha.z")

write_csv(hja_effect.z, "Datasets/Effect sizes/Standardized/hja_effect.z.csv")

## Two-SD approach (predictor scaled by 2*SD; response scaled to 1 SD) ----
# Scale predictor by 2*SD and response to Z-score
live_dead.growth <- live_dead.growth %>%
  mutate(
    dw_mass_ha.2sd    = (dw_mass_ha - mean(dw_mass_ha, na.rm = TRUE)) / (2 * sd(dw_mass_ha, na.rm = TRUE)),
    tree_growth_ind.z = scale(tree_growth_ind)[, 1]
  )

tree_glmm.hybrid <- glmmTMB(
  tree_growth_ind.z ~ dw_mass_ha.2sd + (1 | stand),
  data = live_dead.growth,
  family = gaussian()
)
summary(tree_glmm.hybrid)

# Diagnostics
res.hybrid <- simulateResiduals(tree_glmm.hybrid)
plot(res.hybrid) # Tests should be comparable to your z–z model

## Visualization
# Generate datatable of model predictions
preds.hybrid <- ggpredict(tree_glmm.hybrid, terms = "dw_mass_ha.2sd")

# Plot predicted values with 95% confidence intervals and raw values
ggplot() +
  geom_point(data = live_dead.growth, 
             aes(x = dw_mass_ha.2sd, y = tree_growth_ind.z),
             alpha = 0.6,
             color = "darkgrey") +
  geom_line(data = preds.hybrid, 
            aes(x = x, y = predicted),
            linewidth = 1.2) +
  geom_ribbon(data = preds.hybrid,
              aes(x = x, y = predicted, 
                  ymin = conf.low, ymax = conf.high), 
              alpha = 0.3) +
  labs(x = "Dead wood mass (scaled by 2×SD; 1 unit = +2 SD)",
       y = "Individual tree growth (Z-score)") +
  theme_classic(base_size = 14)

# Extract effect size and CI, write to .csv file
# (Coefficient = SDs of Y per +2 SD change in X)
hja_effect.hybrid <- tidy(tree_glmm.hybrid, effects = "fixed", conf.int = TRUE) %>% 
  filter(term == "dw_mass_ha.2sd")

write_csv(hja_effect.hybrid, "Datasets/Effect sizes/Standardized/hja_effect.hybrid_2sdX_1sdY.csv")

