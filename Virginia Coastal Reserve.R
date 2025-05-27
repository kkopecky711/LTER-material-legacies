#### Virginia Coastal Reserve oyster reef data ####

# Use only reference sites (i.e., no restoration sites)

# Plot juv density against day of year 

# lump seasons; conservative approach could be to limit analysis to summer-fall because juvenile oysters grow quickly enough to enter the 'adult' life stage by summer

# plot live oyster against dead


# Working directory
setwd('/Users/kako4300/Library/CloudStorage/OneDrive-UCB-O365/Projects/LTER material legacy synthesis')

vcr.raw <- read_csv("Datasets/Virginia Coastal Reserve/Oyster_Count_Data_2022_03_10.csv")

# Filter original data for just live and dead oysters
vcr.juv_dead <- vcr.raw %>% 
  clean_names() %>%
  filter(species == "Box Adult Oyster" | species == "Spat Oyster",
         restoration == "Reference") %>% 
  mutate(year = substr(date, nchar(date) - 1, nchar(date)),
         year = paste0("20", year),
         date = mdy(date), # Convert the date column to Date format
         month = month(date), # Extract the month
         species = case_when(species == "Spat Oyster" ~ "Juvenile",
                             TRUE ~ "Dead"),
         season = case_when( # Assign seasons based on month
           month %in% c(12, 1, 2) ~ "Winter",
           month %in% c(3, 4, 5) ~ "Spring",
           month %in% c(6, 7, 8) ~ "Summer",
           month %in% c(9, 10, 11) ~ "Fall",
           TRUE ~ NA_character_),
         obvs_id = paste(year, season, site, sample, sep = "_")) %>%
  select(obvs_id, site, season, species, species_count) %>% 
  pivot_wider(names_from = "species", 
              values_from = "species_count", 
              values_fn = sum)
  
# Exploratory graph of Juvenile oyster density ~ dead oyster density in previous year
ggplot(vcr.juv_dead, aes(x = Dead, y = Juvenile, color = season)) +
  geom_point(alpha = 0.5,
             aes(color = season)) +
  geom_smooth(method = "lm",
              aes(group = season)) +
  #scale_x_continuous(expand = c(0.01,0)) +
  #scale_y_continuous(expand = c(0.01,0)) +
  labs(x = "Dead oyster density",
       y = "Juvenile oyster density") +
  theme_classic()

## Analyses
hist(vcr.juv_dead$Juvenile) # Highly skewed toward zero
mean(vcr.juv_dead$Juvenile)
var(vcr.juv_dead$Juvenile) # Variance >> mean ==> overdispersed, use negative binomial distribution

# GLMM with negative binomial distribution (due to overdispersion) and zero-inflation
mod_zinb <- glmmTMB(
  Juvenile ~ Dead + (1 | site) + (1 | season),
  ziformula = ~1,
  data = vcr.juv_dead,
  family = nbinom2
)

# Standard negative binomial model (no zero inflation)
mod_nb <- glmmTMB(
  Juvenile ~ Dead + (1 | site) + (1 | season),
  data = vcr.juv_dead,
  family = nbinom2
)

# Compare AIC
AIC(mod_nb, mod_zinb) # zero-inflated has better fit

# Simulate residuals
simres <- simulateResiduals(mod_zinb)

# Plot diagnostics
plot(simres)

## Visualization
preds <- ggpredict(
  mod_zinb,
  terms = "Dead")

# Plot predicted values with 95% confidence intervals
ggplot(preds, aes(x = x, y = predicted)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), 
              alpha = 0.3,
              fill = "#3B4F8E") +
  geom_jitter(data = vcr.juv_dead, aes(x = Dead, y = Juvenile),
              alpha = 0.6,
              size = 2,
              width = 0.1) +
  labs(
    x = Dead~oyster~density~(no./0.25~m^2),
    y = "Juvenile oyster density (no./m²)"
  ) +
  theme_classic(base_size = 14)

#### Normalized model for ecosystem comparison ----
# Create a copy of the data
vcr.juv_dead.z <- vcr.juv_dead

# Z-score standardize the response and predictor
vcr.juv_dead.z$juv_density_z <- scale(vcr.juv_dead.z$Juvenile)[, 1]
vcr.juv_dead.z$dead_density_z  <- scale(vcr.juv_dead.z$Dead)[, 1]

zinb_z_gaussian <- glmmTMB(
  juv_density_z ~ dead_density_z + (1 | site) + (1 | season),
  ziformula = ~1,
  family = gaussian,
  data = vcr.juv_dead.z
)

library(broom.mixed)

effect_z <- tidy(zinb_z_gaussian, effects = "fixed", conf.int = TRUE) |>
  filter(term == "dead_density_z")

print(effect_z)
