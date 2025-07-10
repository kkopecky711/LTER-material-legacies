#### Cross-system analyses ####

library(tidyverse)

# Working directory
setwd('/Users/kako4300/Library/CloudStorage/OneDrive-UCB-O365/Projects/LTER material legacy synthesis')

# Read in effect size datafiles
mcr_effect <- read_csv("Datasets/Effect sizes/mcr.effect_size.csv")

bnz_effect <- read_csv("Datasets/Effect sizes/bnz.effect_size.csv")

vcr_effect <- read_csv("Datasets/Effect sizes/vcr.effect_size.csv")

gce_effect <- read_csv("Datasets/Effect sizes/gce.effect_size.csv")

luq_effect <- read_csv("Datasets/Effect sizes/luq.effect_size.csv")

fce_effect <- read_csv("Datasets/Effect sizes/fce.effect_size.csv")

songs_effect <- read_csv("Datasets/Effect sizes/songs.effect_size.csv")

hfr_effect <- read_csv("Datasets/Effect sizes/hfr.effect_size.csv")

knz_effect <- read_csv("Datasets/Effect sizes/knz.effect_size.csv")


# Combine into one dataframe, filter out unneeded values, add ecosystem names
effect_sizes <- rbind(mcr_effect, bnz_effect, vcr_effect, gce_effect, luq_effect, fce_effect, songs_effect, hfr_effect, knz_effect)

effect_sizes <- effect_sizes %>% 
  filter(term != "(Intercept)") %>% 
  mutate(ecosystem = c("Coral reef", "Boreal forest",  "Oyster reef", "Saltmarsh", "Tropical forest", "Mangrove forest", "Kelp forest", "Temperate forest", "Tallgrass prairie"),
         significance = ifelse(p.value < 0.05, "Significant", "Not Significant"))

# Add extra category for marine vs terrestrial systems and transparency dependent on significance
effect_sizes <- effect_sizes %>%
  mutate(
    system = case_when(
      ecosystem %in% c("Coral reef", "Kelp forest", "Oyster reef", "Saltmarsh", "Mangrove forest") ~ "Marine",
      TRUE ~ "Terrestrial"
    ),
    alpha_level = ifelse(significance == "Significant", 1, 0.4)
  )


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
  scale_x_continuous(expand = c(0.13, 0)) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

# Big Sur color palette
bigsur2 = c("#20618D", "#91AAC4", "#6B6C58", "#464724", "#83932D", "#CAB89F")

