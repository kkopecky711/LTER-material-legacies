#### Cross-system analyses ####

library(tidyverse)

# Working directory
setwd('/Users/kako4300/Library/CloudStorage/OneDrive-UCB-O365/Projects/LTER material legacy synthesis/Datasets/Effect sizes')

# Read in effect size datafiles
mcr_effect <- read_csv("Standardized/mcr_effect.z.csv")

bnz_effect <- read_csv("Standardized/bnz_effect.z.csv")

vcr_effect <- read_csv("Standardized/vcr_effect.z.csv")

gce_effect <- read_csv("Standardized/gce_effect.z.csv")

luq_effect <- read_csv("Standardized/luq_effect.z.csv")

fce_effect <- read_csv("Standardized/fce_effect.z.csv")

songs_effect <- read_csv("Standardized/songs_effect.z.csv")

hfr_effect <- read_csv("Standardized/hfr_effect.z.csv")

knz_effect <- read_csv("Standardized/knz_effect.z.csv")

hja_effect <- read_csv("Standardized/hja_effect.z.csv")


# Combine into one dataframe, filter out unneeded values, add ecosystem names
effect_sizes <- rbind(mcr_effect, bnz_effect, vcr_effect, gce_effect, luq_effect, fce_effect, songs_effect, hfr_effect, knz_effect, hja_effect)

effect_sizes <- effect_sizes %>% 
  mutate(ecosystem = c("Coral reef", "Boreal forest",  "Oyster reef", "Salt marsh", "Tropical forest", "Mangrove forest", "Kelp forest", "Temperate forest", "Tallgrass prairie", "Temperate rainforest"),
         significance = ifelse(p.value < 0.05, "Significant", "Not Significant"))

# Add extra category for marine vs terrestrial systems and transparency dependent on significance
effect_sizes <- effect_sizes %>%
  mutate(system = case_when(ecosystem %in% c("Coral reef", "Kelp forest", "Oyster reef", "Salt marsh", "Mangrove forest") ~ "Marine",
      TRUE ~ "Terrestrial"),
    alpha_level = ifelse(significance == "Significant", 1, 0.4))

## Forest plot with system type and transparency based on significance 
ggplot(effect_sizes, aes(x = estimate, y = reorder(ecosystem, estimate), color = system)) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high,
                     alpha = alpha_level),
                 height = 0,
                 linewidth = 1) +  
  geom_point(size = 3,
             color = "black") +
  geom_point(size = 2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Effect size (± 95% CI)",
       y = "Ecosystem type") +
  scale_alpha_identity() +  # uses raw alpha values
  scale_color_manual(values = c("Marine" = "#20618D", "Terrestrial" = "#6B6C58"),
                     name = "") +
  scale_x_continuous(expand = c(0.09, 0)) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

# Big Sur color palette
bigsur2 = c("#20618D", "#91AAC4", "#6B6C58", "#464724", "#83932D", "#CAB89F")

## Absolute values

# Add variables for absolute values of effects sizes and CIs

effect_sizes <- effect_sizes %>% 
  mutate(estimate.abs = abs(estimate),
         conf_low.abs = if_else(estimate < 0, abs(conf.low), conf.low), # to ensure lower CIs of positive relationships aren't artificially shrunk
         conf_high.abs = abs(conf.high))

# Plot absolute values
ggplot(effect_sizes, aes(x = estimate.abs, y = reorder(ecosystem, estimate.abs), color = system)) +
  geom_errorbarh(aes(xmin = conf_low.abs, xmax = conf_high.abs,
                     alpha = alpha_level),
                 height = 0,
                 linewidth = 1) +  
  geom_point(size = 3,
             color = "black") +
  geom_point(size = 2) +
  labs(x = "Effect size (± 95% CI), abs. value",
       y = "Ecosystem type") +
  scale_alpha_identity() +  # uses raw alpha values
  scale_color_manual(values = c("Marine" = "#20618D", "Terrestrial" = "#6B6C58"),
                     name = "") +
  scale_x_continuous(expand = c(0.09, 0)) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")
