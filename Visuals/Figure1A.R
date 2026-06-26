################################################################################
#################### Mapping the samples into Ethiopia  ########################
################################################################################

# Isabela Gerdes Gyuricza - Parr Lab
# 6/23/2026

# Cleaning up workspace
rm(list = ls())
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)

# Setting up working directory
setwd("Ethiopia_2024_Resurgence_Investigation_repo") 
# Update to your local path

# Load libraries
library(sf)
library(tidyverse)

# Loading processed metadata
metadata  <- read.csv("Data/processed_metadata.csv") 

districts <- st_read("Data/Ethiopia_shape_files/Ethiopia_Buffere.shp") 

districts_joined <- districts %>%
  right_join(metadata, by = c("ADM3_EN" = "District")) %>% 
  # There are two samples that don't have information about District, removing them..
  filter(!is.na(ADM3_EN)) %>% 
  group_by(Region, ADM3_EN) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(
    total_samples = sum(n),
    prop = n / total_samples
  )

districts_coords <- districts_joined %>%
  st_centroid() %>%
  mutate(
    lon = st_coordinates(.)[,1],
    lat = st_coordinates(.)[,2]
  )

# Getting one label position per ADM1 region from the ADM3 shapefile
region_labels <- districts %>%
  group_by(ADM1_EN) %>%
  summarise(geometry = st_union(geometry), .groups = "drop") %>%
  st_centroid() %>%
  mutate(
    lon = st_coordinates(.)[,1],
    lat = st_coordinates(.)[,2],
    ADM1_EN = case_when(
      ADM1_EN %in% c("Addis Ababa", "Harari", "Somali", "SNNP") ~ NA_character_,
      TRUE ~ ADM1_EN
    )
  ) %>%
  filter(!is.na(ADM1_EN))

# Creating region-level background from ADM3 shapefile
eth_map_background <- districts %>%
  group_by(ADM1_EN) %>%
  summarise(geometry = st_union(geometry), .groups = "drop")

# Plotting map
pdf("sample_region_number.pdf", width = 6, height = 6)

ggplot() +
  geom_sf(data = eth_map_background, fill = "lightgrey", color = "grey33", linewidth = 0.3) +
  geom_point(
    data = districts_coords,
    aes(x = lon, y = lat, size = n, color = Region)
  ) +
  geom_label(
    data = region_labels,
    aes(x = lon, y = lat, label = ADM1_EN),
    fill = "white",
    color = "black",
    size = 3,
    linewidth = 0.2
  ) +
   scale_size_continuous(name = "Number of participants") +
  scale_color_manual(
    name = "Region",
    values = c(
      "Afar" = "#1b9e77",
      "Amhara" = "#d95f02",
      "Benishangul Gumz" = "#a6761d",
      "Central Ethiopia" = "#7570b3",
      "Dire Dawa" = "blue",
      "Gambela" = "#e6ab02",
      "Oromia" = "red3",
      "Sidama" = "black",
      "South Ethiopia" = "#e7298a",
      "South West Ethiopia" = "#fdc086",
      "Tigray" = "#66a61e"
    )
  ) +
  theme_minimal() +
  coord_sf() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.key.size = unit(0.3, "cm"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11)
  )

dev.off()