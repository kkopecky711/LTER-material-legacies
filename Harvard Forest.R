## Harvard Forest LTER Hemlock removal experiment ##

# Data include live and dead trees of different species and life stages
# Code co-written by Audrey Barker-Plotkin and Kai Kopecky 

# Packages
library(tidyverse)
library(janitor)
library(glmmTMB)
library(ggeffects)
library(DHARMa)
library(broom.mixed)

setwd('/Users/kako4300/Library/CloudStorage/OneDrive-UCB-O365/Projects/LTER material legacy synthesis')

# Dead wood sampled in 2005, 2007, 2009, 2011, 2013, 2015, 2017, 2021 (next 2025).
# Dead wood amounts inferred for 2004, 2014, and 2019 (to line up with tree census dates)
predictorsAll<-read.csv('Datasets/Harvard Forest/DeadByPlot.csv', header=TRUE)
str(predictorsAll)
predictorsAll<-predictorsAll %>% mutate_at(c(2,3,5), as.factor)

# Trees (>5 cm dbh) measured in 2004, 2009, 2014, 2019, 2024
# Sapling trees (>1.3 m tall but <5 cm dbh) measured 2004, 2007, 2009, 2011, 2013, 2015, 2017, 2019, 2021, 2023
responseAll<-read.csv('Datasets/Harvard Forest/TreeSapDensSpecies.csv', header=TRUE)
str(responseAll)
responseAll$species<-ifelse(is.na(responseAll$species), 'none', responseAll$species)
responseAll<-responseAll %>% mutate_at(c(2,3,5,8), as.factor)

#### Data Wrangling ####
levels(responseAll$species)
# group trees/saplings as hemlock (initial dominant), birch (post-disturbance dominant), and other
responseAll$spgroup<-ifelse(responseAll$species == 'TSCA', 'hemlock', ifelse(responseAll$species %in% c('BELE', 'BEAL'), 'birch', 'other'))

# Summarize for density (per hectare) and biomass (grams?)
responseL<-responseAll %>% 
  group_by (plot, trt, block, year, spgroup, stratum) %>%
  summarise (dens_ha = sum(dens.ha),
             biomass_gm2=sum(biomass.gm2, na.rm=TRUE))

# Create wide dataframe of repsonse variables where each species is a column
responseW<-responseL %>%
  pivot_wider(names_from = spgroup, values_from = c(dens_ha, biomass_gm2))
responseW[is.na(responseW)] <- 0 #not exactly since biomass is NA for saplings
responseW$dens_ha_total<-responseW$dens_ha_hemlock +responseW$dens_ha_birch + responseW$dens_ha_other 
responseW$biomass_gm2_total<-responseW$biomass_gm2_hemlock +responseW$biomass_gm2_birch + responseW$biomass_gm2_other 

# Create wide dataframe for predictor variable with two different metrics of deadwood (volume and mass) as columns
predictors<-predictorsAll[, c(1:7)] #the notes mess up the pivot
predictorsW<-predictors %>%
  pivot_wider(names_from = group, values_from = c(volm3ha, massgm2))

#use this data frame for looking at matched response & predictor at each time-point
pred.resp <- predictorsW %>% inner_join( responseW, 
                                         by=c('plot','year', 'trt', 'block'))

## Analyses
# Dataframe for just girdled and logged treatments and just saplings
gird_log.saps <- pred.resp %>% 
  filter(trt %in% c("girdled", "logged"),
         stratum == "Sapling",
         year > 2004) %>% 
  select(c("plot", "trt", "block", "year", "dens_ha_hemlock", "dens_ha_total"))

# time series of treatments, hemlock saplings only
ggplot(gird_log.saps, aes(x = year, y = dens_ha_hemlock, color = trt)) +
  geom_point() +
  geom_line() +
  labs(x = "Year",
       y = "Density of Hemlock saplings (no./ha)") +
  facet_wrap(~block) +
  theme_minimal()

# time series of treatments, all species of saplings
ggplot(gird_log.saps, aes(x = year, y = dens_ha_total, color = trt)) +
  geom_point() +
  geom_line() +
  labs(x = "Year",
       y = "Density of saplings (no./ha)") +
  facet_wrap(~block) +
  theme_minimal()

## GLMM of sapling density as a function of treatment for Hemlock only
# Set the logged treatment as the reference level
gird_log.saps$trt <- relevel(gird_log.saps$trt, ref = "logged") 

# Use Tweedie distribution with log-link
hemlock_glmm.raw <- glmmTMB(dens_ha_hemlock ~ trt + (1 | block/plot) + (1 | year),
                                    family = tweedie(link = "log"), 
                                    data = gird_log.saps)
summary(hemlock_glmm.raw)

# Zero-inflated model
hemlock_glmm.raw.zi <- glmmTMB(dens_ha_hemlock ~ trt + (1 | block/plot) + (1 | year),
                                    family = tweedie(link = "log"),
                                    ziformula = ~1,
                                    data = gird_log.saps)
summary(hemlock_glmm.raw.zi)

# Compare models to ensure zero-inflation improves model fit
AIC(hemlock_glmm.raw, hemlock_glmm.raw.zi) # Zero-inflated model has lower AIC score

# Diagnostics
sim <- simulateResiduals(hemlock_glmm.raw.zi)
plot(sim)
testZeroInflation(sim) 
# Tests of residuals check out

# To report effect size and CI in main text and table
hfr_effect.raw <- tidy(hemlock_glmm.raw.zi, effects = "fixed", conf.int = TRUE, conf.level = 0.95) %>% 
  filter(term == "trtgirdled")

write_csv(hfr_effect.raw, "Datasets/Effect sizes/Raw/hfr_effect.raw.csv")

## Visualization
# Get predicted values and CIs for each treatment
preds <- ggpredict(hemlock_glmm.raw.zi, terms = "trt")

ggplot() +
  geom_jitter(data = gird_log.saps, aes(x = reorder(trt, dens_ha_hemlock), y = dens_ha_hemlock),
              color = "#6B6C58",
              width = 0.15, 
              alpha = 0.6) +
  geom_line(data = preds, 
            aes(x = x, y = predicted, group = group)) +
  geom_errorbar(
    data = preds,
    aes(x = x, ymin = conf.low, ymax = conf.high),
    width = 0) +
  geom_point(data = preds,
             aes(x = x, y = predicted),
             size = 3,
             color = "black") +
  geom_point(data = preds,
             aes(x = x, y = predicted),
             size = 3,
             color = "#6B6C58") +
  labs(x = "Dead hemlock status") +
  theme_classic(base_size = 14)

# Log-transformed for visual purposes
ggplot(gird_log.saps, aes(x = reorder(trt, dens_ha_hemlock), y = dens_ha_hemlock + 1)) +
  geom_jitter(color = "#6B6C58",
              width = 0.15, 
              alpha = 0.6) +
  geom_line(data = preds, aes(x = x, y = predicted + 1, group = group)) +
  geom_pointrange(data = preds, aes(x = x, y = predicted + 1, ymin = conf.low + 1, ymax = conf.high + 1),
                  inherit.aes = FALSE,
                  size = 0.75, 
                  color = "black", 
                  fill = "#6B6C58", 
                  shape = 21) +
  scale_x_discrete(labels = c("girdled" = "Standing", "logged" = "Removed")) +
  scale_y_continuous(trans = "log10",
                     name = expression("Hemlock sapling density (no./ha/yr, log"[10]*")")) +
  labs(x = "Dead hemlock status") +
  theme_classic(base_size = 14)

## Z-score standardized model
# Scale the response variable
gird_log.saps$dens.z <- scale(gird_log.saps$dens_ha_hemlock)

# Fit the zero-inflated model using the z-scored response, and gaussian distribution (Tweedie not needed with scaled response)
hemlock_glmm.z <- glmmTMB(
  dens.z ~ trt + (1 | block/plot) + (1 | year),
  dispformula = ~ trt,  # Allow residual variance to differ by treatment to account for heteroscedasticity
  family = gaussian(),
  data = gird_log.saps
)
summary(hemlock_glmm.z)

# Diagnostics
res.z <- simulateResiduals(hemlock_glmm.z)
plot(res.z) # Tests are non-significant


## Visualization
# Get predicted values and CIs for each treatment
preds.z <- ggpredict(hemlock_glmm.z, terms = "trt")

ggplot() +
  geom_jitter(data = gird_log.saps, 
              aes(x = reorder(trt, dens.z), y = dens.z),
              color = "darkgrey",
              width = 0.15, 
              alpha = 0.6) +
  geom_point(data = preds.z, 
                  aes(x = x, y = predicted),
                  size = 3) +
  geom_errorbar(data = preds.z, 
                aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high),
                width = 0) +
  geom_line(data = preds.z, 
            aes(x = x, y = predicted, group = group)) +
  scale_x_discrete(labels = c("girdled" = "Standing", "logged" = "Removed")) +
  scale_y_continuous(limits = c(-2, 3.4)) +
  labs(x = "Dead hemlock status",
       y = "Sapling density (Z-score)") +
  theme_classic(base_size = 14)


## Extract effect size and 95% CI using broom.mixed
hfr_effect.z <- tidy(hemlock_sap.glmm.z, effects = "fixed", conf.int = TRUE) %>%
  filter(term == "trtgirdled")

print(hfr_effect.z)

write_csv(hfr_effect.z, "Datasets/Effect sizes/Standardized/hfr_effect.z.csv")
