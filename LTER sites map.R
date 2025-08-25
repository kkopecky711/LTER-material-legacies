#### Map of LTER sites used in material legacy synthesis project ####

# Load libraries
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(cowplot)
library(ggrepel)

# Load world map
world <- ne_countries(scale = "large", returnclass = "sf")

# LTER site data with 3-letter codes
lter_sites <- data.frame(
  site = c("MCR", "SONGS", "BNZ", "KNZ", "HFR", "VCR", "GCE", "FCE", "LUQ", "AND"),
  lon = c(
    -149.8260, -117.3745, -147.8839, -96.5610, -72.1760,
    -75.6850, -81.2889, -80.6800, -65.8201, -122.1762
  ),
  lat = c(
    -17.4920, 33.2502, 64.7925, 39.1006, 42.5378,
    37.4151, 31.3722, 25.3800, 18.3392, 44.2332
  )
)

# Separate MCR and the rest
mcr_site <- subset(lter_sites, site == "MCR")
main_sites <- subset(lter_sites, site != "MCR")

# Big Sur color palette
bigsur2 = c("#20618D", "#91AAC4", "#6B6C58", "#464724", "#83932D", "#CAB89F")

# ---- MAIN MAP (Excludes MCR) ----
main_map <- ggplot(data = world) +
  geom_sf(fill = "gray90", color = "gray40") +
  geom_point(data = main_sites, 
             aes(x = lon, y = lat), 
             color = "black", 
             size = 2, 
             shape = 24) +
  geom_point(data = main_sites, 
             aes(x = lon, y = lat), 
             color = "#20618D", 
             size = 1.5, 
             shape = 17) +
  geom_label(data = main_sites, aes(x = lon, y = lat, label = site),
             vjust = -0.4, 
             hjust = 0.5, 
             size = 3,
             label.size = 0.2, 
             fill = "white", 
             color = "black") +
  coord_sf(xlim = c(-152.2, -57), 
           ylim = c(15, 70), 
           expand = FALSE) +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "lightblue"),
    panel.grid.major = element_line(color = "white", linewidth = 0.1)) +
  labs(x = "Longitude",
       y = "Latitude")

# ---- INSET MAP for MCR ----
inset_map <- ggplot(data = world) +
  geom_sf(fill = "gray90", color = "gray40") +
  geom_point(data = mcr_site, 
             aes(x = lon, y = lat), 
             color = "black", 
             size = 2, 
             shape = 24) +
  geom_point(data = mcr_site, 
             aes(x = lon, y = lat), 
             color = "#20618D", 
             size = 1.5, 
             shape = 17) +
  geom_label(data = mcr_site, aes(x = lon, y = lat, label = site),
             vjust = -0.4, 
             hjust = 0.5, 
             size = 3,
             label.size = 0.2, 
             fill = "white", 
             color = "black") +
  coord_sf(xlim = c(-150.2, -149), 
           ylim = c(-18, -17), 
           expand = FALSE) +
  scale_x_continuous(breaks = c(-150, -149)) +
  scale_y_continuous(breaks = c(-17, -17.5, -18)) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "lightblue"),
    panel.grid.major = element_line(color = "white", linewidth = 0.1),
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 7),
    plot.title = element_text(size = 9)
  ) +
  labs(x = "",
       y = "")

# ---- Combine using cowplot ----
final_map <- ggdraw() +
  draw_plot(main_map) +
  draw_plot(inset_map, x = 0.12, y = 0.1, width = 0.3, height = 0.3)

# Display combined map
print(final_map)

