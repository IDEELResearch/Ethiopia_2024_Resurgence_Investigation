################################################################################
##############E###### Processing and cleaning metadata  ########################
################################################################################

# Isabela Gerdes Gyuricza - Parr Lab
# 6/01/2026

# Cleaning up workspace
rm(list = ls())
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)

# Setting up working directory
setwd("Ethiopia_2024_Resurgence_Investigation_repo") 
# Update to your local path

# Load libraries
library(tidyverse)

# Loading cleaned metadata from Ashenafi
metadata_clean <- read.csv("Data/metadata_cleaned.csv", 
                           na.strings = "") %>% 
  select(ADM3_EN, District) %>% 
  filter(!is.na(District)) %>% 
  unique()

# Read metadata
metadata <- read.csv("Data/outbreakmetadatacomplete.csv") %>%
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
  left_join(metadata_clean, by = "District") %>% 
  # Fixing manually some of the district names
  mutate(ADM3_EN = case_when(Sentinel_Site == "YECHILA P/HOSPITAL" ~ "Abergele (TG)",
                             Sentinel_Site == "ABOL H/CENTER" ~ "Gog",
                             ADM3_EN == "Shone town" ~ "Shone Town",
                             ADM3_EN == "South Gondar" ~ "Gondar town",
                             ADM3_EN == "South Wello" ~ "Argoba",
                             .default = ADM3_EN)) %>% 
  select(-District) %>% 
  rename(District = ADM3_EN) %>% 
  select(Barcode, Region, District, Sentinel_Site, Age, Sex)

# There are 2 samples that were sequenced and don't have metadata. Adding their
# region manually because we know from the sample names where they came from.
# The rest of the metadata will be NA

missing_samples <- tibble(
  Barcode = c("AMHA485", "DILE6100"),
  Region = c("Amhara", "Dire Dawa")
)

metadata <- bind_rows(metadata, missing_samples)

rm(missing_samples)

# Save processed metadata
write.csv(metadata,
          file = "Data/processed_metadata.csv",
          row.names = FALSE)
