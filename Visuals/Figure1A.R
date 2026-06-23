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

# Read the shapefile to have the background of Ethiopia map with the boundaries by regions.
eth_map_background <- st_read("Data/eth_adm_csa_bofedb_2021_shp/eth_admbnda_adm1_csa_bofedb_2021.shp") 

districts <- st_read("Data/eth_adm_csa_bofedb_2021_shp/eth_admbnda_adm3_csa_bofedb_2021.shp") 

districts_joined <- districts %>%
  right_join(metadata, by = c("ADM3_EN" = "District")) %>% 
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

region_labels <- eth_map_background %>%
  st_centroid() %>%
  cbind(st_coordinates(.)) %>% 
  # Only labeling the regions of interest
  mutate(ADM1_EN = case_when(ADM1_EN == "Addis Ababa" ~ NA,
                             ADM1_EN == "Harari" ~ NA,
                             ADM1_EN == "Somali" ~ NA,
                             .default = ADM1_EN))

# Plotting
pdf("sample_region_number.pdf", width = 6, height = 6)

ggplot() +
  geom_sf(data = eth_map_background, fill = "lightgrey", color = "grey33") +
  geom_point(data = districts_coords,
             aes(x = lon, y = lat, size = n, color = Region)) +
  geom_label(
    data = region_labels,
    aes(x = X, y = Y, label = ADM1_EN),
    fill = "white",
    color = "black",
    size = 3,
    label.size = 0.2
  ) +
  scale_size_continuous(name = "Number of participants") +
  scale_color_manual(name = "Region",
                     values = c('#1b9e77','#d95f02','#a6761d','#e7298a',
                                '#66a61e','#e6ab02','#7570b3','black',
                                "red3","blue")) + 
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