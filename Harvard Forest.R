## Harvard Forest LTER Hemlock removal experiment ##

# Data include live and dead trees of different species and life stages

## Questions for Audrey:

# 1) What are the 'Hemlock' and 'Hardwood' treatments?
# 2) Could the 'girdled' and 'logged' treatments be grouped so as to have a 'gradient' of deadwood cover?
# 3) If yes to 2), is Hemlock or total sapling density ~ deadwood a sensible relationship to explore?


library(tidyverse)
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

# GLMM of sapling density as a function of treatment
# Hemlock only

gird_log.saps$trt <- relevel(gird_log.saps$trt, ref = "logged") # Make sure the logged treatment is the reference level

hemlock_sap.glmm.tweedie <- glmmTMB(dens_ha_hemlock ~ trt + (1 | block/plot) + (1 | year),
                                    family = tweedie(link = "log"), 
                                    data = gird_log.saps)
summary(hemlock_sap.glmm.tweedie)

hemlock_sap.glmm.tweedie.zi <- glmmTMB(dens_ha_hemlock ~ trt + (1 | block/plot) + (1 | year),
                                    family = tweedie(link = "log"),
                                    ziformula = ~1,
                                    data = gird_log.saps)
summary(hemlock_sap.glmm.tweedie.zi)
tidy(hemlock_sap.glmm.tweedie.zi, effects = "fixed", conf.int = TRUE, conf.level = 0.95)


AIC(hemlock_sap.glmm.tweedie, hemlock_sap.glmm.tweedie.zi) # Zero-inflated model has lower AIC score

# Diagnostics
sim <- simulateResiduals(hemlock_sap.glmm.tweedie.zi)
plot(sim)
testZeroInflation(sim) 
# Tests of resiudals seem to check out

## Visualization
# Get predicted values and CIs for each treatment
preds <- ggpredict(hemlock_sap.glmm.tweedie.zi, terms = "trt")

# Base plot: raw data
ggplot(gird_log.saps, aes(x = reorder(trt, dens_ha_hemlock), y = dens_ha_hemlock + 1)) +
  geom_jitter(width = 0.1, 
              alpha = 0.6, 
              size = 2) +
  geom_pointrange(data = preds, aes(x = x, y = predicted + 1, ymin = conf.low + 1, ymax = conf.high + 1),
                  inherit.aes = FALSE,
                  size = 1, 
                  color = "black", 
                  fill = "black", 
                  shape = 21) +
  geom_line(data = preds, aes(x = x, y = predicted + 1, group = group)) +
  scale_y_continuous(trans = "log10",
                     name = "Sapling density (no./ha, log + 1)") +
  labs(x = "Dead hemlock status") +
  theme_classic(base_size = 14)

## Z-score standardized model
# Scale the response variable
gird_log.saps$dens_z <- scale(gird_log.saps$dens_ha_hemlock)

# Fit the zero-inflated model using the z-scored response, and gaussian distribution (Tweedie not needed with scaled response)
hemlock_sap.glmm.z <- glmmTMB(
  dens_z ~ trt + (1 | block/plot) + (1 | year),
  family = gaussian(),
  data = gird_log.saps
)
summary(hemlock_sap.glmm.z)

# Extract effect size and 95% CI using broom.mixed
hemlock_effect <- tidy(hemlock_sap.glmm.z, effects = "fixed", conf.int = TRUE) %>%
  filter(term == "trtgirdled")

print(hemlock_effect)

write_csv(hemlock_effect, "Datasets/Effect sizes/hfr.effect_size.csv")





#### Exploratory analyses ####

# Dataframe for just girdled and logged treatments and adult trees
gird_log.trees <- pred.resp %>% 
  filter(trt %in% c("girdled", "logged"),
         stratum == "Tree") %>% 
  select(-c(14,16,17,19))

## Exploratory graphs
# Time series
ggplot(gird_log.trees, aes(x = year, y = biomass_gm2_hemlock, color = trt)) +
  geom_point() +
  geom_line() +
  labs(x = "Year",
       y = "Biomass of Adult Hemlock") +
  facet_wrap(~block) +
  theme_minimal()

# Adult tree biomass ~ deadwood
ggplot(gird_log.trees, aes(x = volm3ha_TotalDead, y = biomass_gm2_hemlock, color = trt)) +
  geom_point()






















## Saplings ~ deadwood
ggplot(gird_log.saps, aes(x = volm3ha_TotalDead, y = dens_ha_hemlock, color = trt)) +
  geom_point()

ggplot(gird_log.saps, aes(x = massgm2_TotalDead, y = dens_ha_hemlock, color = trt)) +
  geom_point()

## All species of saplings
ggplot(gird_log.saps, aes(x = massgm2_TotalDead, y = dens_ha_total, color = trt)) +
  geom_point()

ggplot(gird_log.saps, aes(x = volm3ha_TotalDead, y = dens_ha_total, color = trt)) +
  geom_point() 


#### Initial input of dead wood -------
#another data frame to see how response rate is affected by INITIAL input of legacy material (let's pick 2007 since that is when the girdled trees would all have died)
predictor07<-subset(predictorsW, predictorsW$year == 2007)
predictor07<-predictor07[,c(1:3, 5:12)]
pred.respINIT<- predictor07 %>% inner_join( responseW, 
                                            by=c('plot', 'trt', 'block'))

pred.resp$trt <- ordered(pred.resp$trt, levels = c("hemlock", "girdled", "logged", "hardwood"))
