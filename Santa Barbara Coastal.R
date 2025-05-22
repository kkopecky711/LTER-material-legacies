# Code for anlaysis of Santa Barbara Coastal LTER data
# Analysis of the inlfuence of dead giant kelp holdfasts on demographic processes in live macroalgae



# Packages
library(tidyverse)
library(janitor)
library(glmmTMB)
library(ggeffects) 

# Working directory
setwd('/Users/kako4300/Library/CloudStorage/OneDrive-UCB-O365/Projects/LTER material legacy synthesis')

#### Annual kelp forest monitoring data ####
# Read in datafile and filter for only algae observations 
sbc.annual <- read_csv("Datasets/Santa Barbara Coastal/Annual_Quad_Swath_All_Years_20240823.csv") %>% 
  clean_names() %>% 
  filter(group == "ALGAE") 

# Create dataframe with only juvenile macrocystis
sbc.macro <- sbc.annual %>%
  filter(scientific_name == "Macrocystis pyrifera",
         count > -1) %>% 
  mutate(obvs_id = paste(year, site, transect, quad, side, sep = "_")) %>% 
  select(-c(2,3,12:23))

# Create dataframe with only dead macrocystis holdfasts
sbc.dead_macro <- read_csv("Datasets/Santa Barbara Coastal/SBS_Cover_All_Years_20240823.csv") %>%
  clean_names() %>% 
  filter(sp_code == "DMH") %>% 
  mutate(obvs_id = paste(year, site, transect, quad, side, sep = "_")) %>% 
  select(c(obvs_id, percent_cover)) %>% 
  rename(dead_cover = percent_cover)

# Merge juvenile macro and dead holdfast dataframes
sbc.juv_dead <- merge(sbc.macro, sbc.dead_macro, by = "obvs_id")
 
# Exploratory graph of juvenile macrocystis count ~ % cover of dead holdfasts
ggplot(sbc.juv_dead, aes(x = dead_cover, y = count)) +
  geom_jitter() +
  geom_smooth(method = "lm") +
  theme_minimal()

# Mixed effects model of juv count ~ dead cover, site as random effect
sbc.juv_dead.glmm <- glmmTMB(count ~ dead_cover + (1|site),
                             family = gaussian(),
                             data = sbc.juv_dead)
summary(sbc.juv_dead.glmm)

plot(ggpredict(sbc.juv_dead.glmm, terms = ~dead_cover))
sbc.juv_dead.predictions <- as.data.frame(ggpredict(sbc.juv_dead.glmm, terms = ~dead_cover))

# Visualization of GLMM
sbc.juv_dead %>%
  ggplot()+
  geom_ribbon(data = sbc.juv_dead.predictions, 
              aes(x = x, y = predicted, ymin =conf.low, ymax = conf.high), 
              alpha =0.5,
              fill = "#3B4F8E")+
  geom_line(data = sbc.juv_dead.predictions, 
            aes(x = x, y = predicted))+
  geom_jitter(aes(x = dead_cover, y = count),
              alpha = 0.6,
              size = 2,
              width = 0.5) +
  labs(x = "Percent cover of dead holdfasts",
       y = "Density of juvenile kelp") +
  theme_classic(base_size = 12)

#### Long Term Experiment (LTE): long term giant kelp removal data ####
# Read in datafile and filter for only algae observations 
sbc.lte <- read_csv("Datasets/Santa Barbara Coastal/LTE_Cover_All_Years_20240501.csv") %>% 
  clean_names() %>% 
  filter(group == "ALGAE",
         taxon_genus == "Macrocystis")

# # Create dataframe with only juvenile macrocystis
# sbc.macro <- sbc.annual %>%
#   filter(scientific_name == "Macrocystis pyrifera",
#          count > -1) %>% 
#   mutate(obvs_id = paste(year, site, transect, quad, side, sep = "_")) %>% 
#   select(-c(2,3,12:23))

# Create dataframe with only dead macrocystis holdfasts
lte.dead_macro <- sbc.lte %>%
  filter(sp_code == "DMH",
         percent_cover > 0) %>% 
  mutate(obvs_id = paste(year, site, transect, quad, side, sep = "_")) %>% 
  #select(c(obvs_id, percent_cover)) %>% 
  rename(dead_cover = percent_cover)

ggplot(lte.dead_macro, aes(x = dead_cover)) +
  geom_histogram(binwidth = 5,
                 color = 'black') +
  facet_wrap(~treatment) +
  theme_minimal()

# Merge juvenile macro and dead holdfast dataframes
sbc.juv_dead <- merge(sbc.macro, sbc.dead_macro, by = "obvs_id")

# Exploratory graph of juvenile macrocystis count ~ % cover of dead holdfasts
ggplot(sbc.juv_dead, aes(x = dead_cover, y = count)) +
  geom_jitter() +
  geom_smooth(method = "lm") +
  theme_minimal()

# Mixed effects model of juv count ~ dead cover, site as random effect
sbc.juv_dead.glmm <- glmmTMB(count ~ dead_cover + (1|site),
                             family = gaussian(),
                             data = sbc.juv_dead)
summary(sbc.juv_dead.glmm)

plot(ggpredict(sbc.juv_dead.glmm, terms = ~dead_cover))
sbc.juv_dead.predictions <- as.data.frame(ggpredict(sbc.juv_dead.glmm, terms = ~dead_cover))

# Visualization of GLMM
sbc.juv_dead %>%
  ggplot()+
  geom_ribbon(data = sbc.juv_dead.predictions, 
              aes(x = x, y = predicted, ymin =conf.low, ymax = conf.high), 
              alpha =0.5,
              fill = "#3B4F8E")+
  geom_line(data = sbc.juv_dead.predictions, 
            aes(x = x, y = predicted))+
  geom_jitter(aes(x = dead_cover, y = count),
              alpha = 0.6,
              size = 2,
              width = 0.5) +
  labs(x = "Percent cover of dead holdfasts",
       y = "Density of juvenile kelp") +
  theme_classic(base_size = 12)
