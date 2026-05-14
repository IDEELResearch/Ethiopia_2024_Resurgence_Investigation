################################################################################
#################### Mapping the samples into Ethiopia  ########################
################################################################################

# Isabela Gerdes Gyuricza - Parr Lab
# 7/23/2025

# Cleaning up workspace
rm(list = ls())
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)

# Setting up working directory
setwd("/Users/isabelagyuricza/OneDrive - University of North Carolina at Chapel Hill/IDEEL_PhD")

# Load libraries
library(sf)
library(tidyverse)

# Debugging metadata!!
# metadata_debug <- read.csv("Outbreak_EPHI/metadata_outbreak_results/outbreakmetadatacomplete.csv") %>% 
#   group_by(Region, District, Sentinel_Site) %>% 
#   tally()
# 
# write.csv(metadata_debug, file = "Outbreak_EPHI/Data/metadata_to_debug.csv",
#           row.names = FALSE)

# Loading debugged metadata from Ashenafi
metadata_debug <- read.csv("Outbreak_EPHI/Assefa_Ethiopia_2026_repo/Data/metadata_cleaned.csv", 
                           na.strings = "") %>% 
  select(ADM3_EN, District) %>% 
  filter(!is.na(District)) %>% 
  unique()

# Read metadata
metadata <- read.csv("Outbreak_EPHI/Assefa_Ethiopia_2026_repo/Data/outbreakmetadatacomplete.csv") %>%
  mutate(Region = case_when(Region == "DIREDAWA" ~ "Dire Dawa",
                            Region == "Benishangul-Gumuz " ~ "Benishangul Gumz",
                            Region == "TIGRAY" ~ "Tigray",
                            Region == "AMHARA" ~ "Amhara",
                            Region == "CENTRAL-ETHIOPIA" ~ "SNNP",
                            Region == "AFAR" ~ "Afar",
                            Region == "OROMIA" ~ "Oromia",
                            Region == "S-WEST-ETHIOPIA" ~ "South West Ethiopia",
                            Region == "SOUTH-ETHIOPIA" ~ "SNNP",
                            Region == "SIDAMA" ~ "Sidama",
                            Region == "GAMBELA" ~ "Gambela")) %>% 
  left_join(metadata_debug, by = "District") %>% 
  # Fixing manually some of the district names
  mutate(ADM3_EN = case_when(Sentinel_Site == "YECHILA P/HOSPITAL" ~ "Abergele (TG)",
                             Sentinel_Site == "ABOL H/CENTER" ~ "Gog",
                             ADM3_EN == "Shone town" ~ "Shone Town",
                             ADM3_EN == "South Gondar" ~ "Gondar town",
                             ADM3_EN == "South Wello" ~ "Argoba",
                             .default = ADM3_EN))

# Save processed metadata
# write.csv(metadata, 
#           file = "Outbreak_EPHI/Data/processed_metadata.csv",
#           row.names = FALSE)

# Read the shapefile to have the background of Ethiopia map with the boundaries by regions.
eth_map_background <- st_read("Outbreak_EPHI/Data/eth_adm_csa_bofedb_2021_shp/eth_admbnda_adm1_csa_bofedb_2021.shp") 

districts <- st_read("Outbreak_EPHI/Data/eth_adm_csa_bofedb_2021_shp/eth_admbnda_adm3_csa_bofedb_2021.shp") 

districts_joined <- districts %>%
  right_join(metadata, by = "ADM3_EN") %>% 
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

# Plotting
pdf("Outbreak_EPHI/metadata_outbreak_results/sample_region_proportion.pdf",
    width = 6, height = 6)

ggplot() +
  geom_sf(data = eth_map_background, fill = "lightgrey", color = "grey33") +
  geom_point(data = districts_coords,
             aes(x = lon, y = lat, size = n, color = Region)) +
  scale_size_continuous(name = "Number of samples") +
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

png("Outbreak_EPHI/metadata_outbreak_results/sample_region_proportion.png",
    res = 600, width = 6, height = 5, units = "in")

ggplot() +
  geom_sf(data = eth_map_background, fill = "lightgrey", color = "grey33") +
  geom_point(data = districts_coords,
             aes(x = lon, y = lat, size = n, color = Region)) +
  scale_size_continuous(name = "Number of samples") +
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
