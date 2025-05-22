#### Konza Prairie LTER ####

## Notes
# Explore relationship of previous year dead biomass on current year productivity
# Should not group forbs and grass; grass could be considered the foundation species, and forbs neither produce legacies in the same way, nor are they annuals

# Look at relationship between dead biomass and biomass of forbs and woody 

# potential relationships between dead biomass and abiotic conditions (soil temperature, soil moisture, N content)

# Look into % cover data as well

# Watershed = fire treatment: 004 means burned every four years, 020 = every 20 years, b = statistical replicate
# Soiltype: tu = lowland, thick, moist; fl = upland, thin, rocky, drier

# Consider trying a running mean approach

setwd('/Users/kako4300/Library/CloudStorage/OneDrive-UCB-O365/Projects/LTER material legacy synthesis')

knz <- read_csv("Datasets/Konza Prairie/PAB011.csv") %>% 
  clean_names() %>% 
  mutate(watershed = as.factor(watershed))

hist(knz$lvgrass)
hist(log(knz$lvgrass))

hist(knz$pryrdead)
hist(log(knz$pryrdead))

# Exploratory graph of livegrass biomass ~ previous year dead grass biomass, separated by burn frequency treatment
ggplot(knz, aes(x = pryrdead, y = lvgrass, color = watershed)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", 
              aes(color = watershed)) +
  #scale_x_continuous(expand = c(0.01,1)) +
  # scale_y_continuous(expand = c(0.01,1)) +
  labs(x = "Previous year dead biomass (g/m^2)",
       y = "Current year live biomass (g/m^2)") +
  theme_classic()

#### Four year burns only -----
# Remove watersheds with annual and 20 year burns; remove all observations from burn years within 4-year watershed; remove observations with > 100 g/m^2 of woody biomass
knz.four_year <- knz %>% 
  filter(watershed == "004b",
         # pryrdead > 0,
         # lvgrass > 0,
         woody < 100) %>% 
  select(recyear, soiltype, transect, plotnum, lvgrass, pryrdead) %>% 
  mutate(recyear = as.numeric(recyear),
         soiltype = as.factor(soiltype)) %>% 
  rename(year = recyear)

hist(knz.four_year$lvgrass)
hist(log(knz.four_year$lvgrass))

hist(knz.four_year$pryrdead)
hist(log(knz.four_year$pryrdead))


# Create dataframe for years when burning was conducted
# Load burn year data and filter for only four-year watershed
burns.four_year <- read_csv("Datasets/Konza Prairie/KFH011.csv") %>% 
  clean_names() %>% 
  filter(watershed == "4B",
         code != "004A") %>% 
  select(year) %>% 
  distinct(year, .keep_all = TRUE) %>% 
  mutate(burn_year = "yes")

# Create complete set of burned and unburned years
all_years <- tibble(year = seq(1991, max(knz.four_year$year)))

# Full join with all years and replace NAs with "no"
burn_data_filled <- all_years %>%
  left_join(burns.four_year, by = "year") %>%
  mutate(burn_year = replace_na(burn_year, "no"))

# Merge burn data with biomass data
knz.four_year <- knz.four_year %>% 
  left_join(burn_data_filled, by = "year") 

knz.four_year <- knz.four_year %>% 
  mutate(burn_year = as.factor(burn_year))

# Create numerical variable for years since burn
knz.four_year <- knz.four_year %>% 
  arrange(year) %>%
  mutate(last_burn = ifelse(burn_year == "yes", year, NA)) %>%
  fill(last_burn, .direction = "down") %>%  
  mutate(years_since_burn = year - last_burn,
         years_since_burn = as.factor(years_since_burn)) %>% 
  filter(lvgrass > 0,
         pryrdead > 0) %>% 
  drop_na()

# GLMM of live biomass ~ previous year's dead biomass
knz.four_year.glmm <- glmmTMB(lvgrass ~ pryrdead*years_since_burn + (1|soiltype) + (1|year) + (1|transect),
                              family = "gaussian",
                              data = knz.four_year)
summary(knz.four_year.glmm)

predictions.knz.four_year <- ggpredict(knz.four_year.glmm, terms = ~ pryrdead*years_since_burn)
plot(predictions.knz.four_year)
predictions.knz.four_year <- as.data.frame(predictions.knz.four_year)

# Model visualization
ggplot()+
  geom_point(data = knz.four_year, 
             aes(x = pryrdead, 
                 y = lvgrass,
                 color = years_since_burn),
             alpha = 0.6,
             size = 2) +
  geom_ribbon(data = predictions.knz.four_year, 
              aes(x = x, 
                  y = predicted, 
                  ymin =conf.low, 
                  ymax = conf.high,
                  group = group,
                  fill = group), 
              #fill = "#3B4F8E",
              alpha =0.5)+
  geom_line(data = predictions.knz.four_year, 
            aes(x = x , 
                y = predicted, 
                group = group,
                color = group)) +
  scale_x_continuous(expand = c(0.01,1)) +
  scale_y_continuous(#limits = c(0,160), 
    expand = c(0.01,1)) +
  labs(x = Previous~year~dead~biomass~(g/m^2),
       y = Current~year~live~biomass~(g/m^2)) +
  #facet_wrap(~years_since_burn) +
  theme_classic()

# Model visualization
ggplot()+
  geom_point(data = knz.four_year, 
             aes(x = pryrdead, 
                 y = lvgrass),
             alpha = 0.6,
             size = 2) +
  geom_ribbon(data = predictions.knz.four_year, 
              aes(x = x, 
                  y = predicted, 
                  ymin =conf.low, 
                  ymax = conf.high,
                  group = group), 
              #fill = "#3B4F8E",
              alpha =0.5)+
  geom_line(data = predictions.knz.four_year, 
            aes(x = x , 
                y = predicted, 
                group = group)) +
  scale_x_continuous(expand = c(0.01,1)) +
  scale_y_continuous(#limits = c(0,160), 
    expand = c(0.01,1)) +
  labs(x = Previous~year~dead~biomass~(g/m^2),
       y = Current~year~live~biomass~(g/m^2)) +
  #facet_wrap(~years_since_burn) +
  theme_classic()

#### 20 year burns only ----
# Remove watersheds with annual and 20 year burns; remove all observations from burn years within 4-year watershed; remove observations with > 100 g/m^2 of woody biomass
knz.20_year <- knz %>% 
  mutate(woody = replace(woody, is.na(woody), 0)) %>% 
  filter(watershed == "020b",
         # pryrdead > 0,
         # lvgrass > 0,
         woody <= 100) %>% 
  select(recyear, soiltype, transect, plotnum, lvgrass, pryrdead) %>% 
  mutate(recyear = as.numeric(recyear),
         soiltype = as.factor(soiltype)) %>% 
  rename(year = recyear)

hist(knz.20_year$lvgrass)
hist(log(knz.20_year$lvgrass))

hist(knz.20_year$pryrdead)
hist(log(knz.20_year$pryrdead))


# Create dataframe for years when burning was conducted
# Load burn year data and filter for only four-year watershed
burns.20_year <- read_csv("Datasets/Konza Prairie/KFH011.csv") %>% 
  clean_names() %>% 
  filter(watershed == "20B") %>% # May want to include 020C, and R20A and R20B, but only after 2000 when 20-year burn frequency was initiated
  select(year) %>% 
  distinct(year, .keep_all = TRUE) %>% 
  mutate(burn_year = "yes")

# Create complete set of burned and unburned years
all_years <- tibble(year = seq(1975, max(knz.20_year$year)))

# Full join with all years and replace NAs with "no"
burn_data_filled <- all_years %>%
  left_join(burns.20_year, by = "year") %>%
  mutate(burn_year = replace_na(burn_year, "no"))

# Merge burn data with biomass data
knz.20_year <- knz.20_year %>% 
  left_join(burn_data_filled, by = "year") 

knz.20_year <- knz.20_year %>% 
  mutate(burn_year = as.factor(burn_year))

# Create numerical variable for years since burn
knz.20_year.yrs_since <- knz.20_year %>% 
  arrange(year) %>%
  mutate(last_burn = case_when(year < 1991 ~ 1975,
                               year >= 1991 & year < 2017 ~ 1991,
                               TRUE ~ 2017)) %>%
  mutate(years_since_burn = year - last_burn) %>% 
  filter(lvgrass > 0,
         pryrdead > 0,
         years_since_burn > 19) %>% 
  drop_na()

# GLMM of live biomass ~ previous year's dead biomass
knz.20_year.glmm <- glmmTMB(lvgrass ~ pryrdead,
                              family = "gaussian",
                              data = knz.20_year.yrs_since)
summary(knz.20_year.glmm)

predictions.knz.20_year <- ggpredict(knz.20_year.glmm, terms = ~ pryrdead)
plot(predictions.knz.20_year)
predictions.knz.20_year <- as.data.frame(predictions.knz.20_year)

# Model visualization
ggplot()+
  geom_point(data = knz.20_year.yrs_since, 
             aes(x = pryrdead, 
                 y = lvgrass),
             alpha = 0.6,
             size = 2) +
  geom_ribbon(data = predictions.knz.20_year, 
              aes(x = x, 
                  y = predicted, 
                  ymin =conf.low, 
                  ymax = conf.high,
                  group = group,
                  fill = group), 
              #fill = "#3B4F8E",
              alpha =0.5)+
  geom_line(data = predictions.knz.20_year, 
            aes(x = x , 
                y = predicted, 
                group = group,
                color = group)) +
  scale_x_continuous(expand = c(0.01,1)) +
  scale_y_continuous(#limits = c(0,160), 
    expand = c(0.01,1)) +
  labs(x = Previous~year~dead~biomass~(g/m^2),
       y = Current~year~live~biomass~(g/m^2)) +
  #facet_wrap(~years_since_burn) +
  theme_classic()
