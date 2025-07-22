#### Niwot ridge alpine grasslang biomass data

# Packages
library(tidyverse)
library(janitor)
library(glmmTMB)
library(ggeffects)
library(DHARMa)
library(broom.mixed)

# Working directory
setwd('/Users/kako4300/Library/CloudStorage/OneDrive-UCB-O365/Projects/LTER material legacy synthesis/Datasets/Niwot Ridge/')

# Read in and clean data
nwt.raw <- read_csv("biomnitr.mf.data.csv") %>% 
  clean_names()

# Create wide dataframe that separates live and dead
nwt.live_dead <- nwt.raw %>%
  select(-c(tissue_n, n_content))

nwt.live_dead.wide <- nwt.live_dead %>% 
  mutate(year = year(date)) %>% 
  filter(above_below == "aboveground") %>% 
  pivot_wider(
    id_cols = c(year, date, sample, meadow, site_num),
    names_from = fraction,
    values_from = biomass,
    names_prefix = "biomass_"
  )

