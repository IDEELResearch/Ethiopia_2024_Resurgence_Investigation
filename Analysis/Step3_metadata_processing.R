#################################################################################
##############E###### Processing and cleaning metadata  ########################
################################################################################

# Isabela Gerdes Gyuricza - Parr Lab
# 6/25/2026

# Cleaning up workspace
rm(list = ls())
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)

# Setting up working directory
setwd("Ethiopia_2024_Resurgence_Investigation_repo") 
# Update to your local path

# Load libraries
library(tidyverse)

metadata <- read.csv("Data/OB_Metadata_final.csv") %>% 
  rename(Region = Regions,
         Sentinel_Site = Sentinel.site.name) %>% 
  mutate(Region = case_when(Region %in% c("DIREDAWA","DIRE DAWA") ~ "Dire Dawa",
                            Region == "B/G/R/S" ~ "Benishangul Gumz",
                            Region == "TIGRAY" ~ "Tigray",
                            Region == "AMHARA" ~ "Amhara",
                            Region == "C/ETHIOPIA" ~ "Central Ethiopia",
                            Region == "AFAR" ~ "Afar",
                            Region == "OROMIA" ~ "Oromia",
                            Region == "S/WEST" ~ "South West Ethiopia",
                            Region %in% c("S/ETHIOPIA","S/ETHIOIA","S/ETHIOIPA","S/ETH") ~ "South Ethiopia",
                            Region %in% c("SIDAMA","Sidama ") ~ "Sidama",
                            Region == "GAMBELA" ~ "Gambela",
                            .default = Region))

# Loading district information from shape file
district_info <- read.csv("Data/district_info.csv") 

metadata_combined <- metadata %>% 
  left_join(district_info, by = c("Sentinel_Site", "Region")) %>% 
  mutate(District = case_when(Sentinel_Site == "Jiga h/c" ~ "Jabi Tehnan",
                              Sentinel_Site %in% c("Sheraro h/c","Sheraro Health Center") ~ "Sheraro town",
                              Sentinel_Site == "Gara Riketa h/c" ~ "Hawassa Zuria",
                              Sentinel_Site == "Dila Health Center" ~ "Dila town",
                              Sentinel_Site == "Worer Health Center" ~ "Ambira",
                              Sentinel_Site == "Woreta Health Center" ~ "Woreta",
                              .default = District)) %>%
  mutate(Sex = ifelse(Sex == 1, "Male","Female")) %>% 
  select(Barcode, Region, District, Sentinel_Site, Age, Sex)

# There are 2 samples that were sequenced and don't have metadata. Adding their
# region manually because we know from the sample names where they came from.
# The rest of the metadata will be NA

missing_samples <- tibble(
  Barcode = c("AMHA485", "DILE6100"),
  Region = c("Amhara", "Dire Dawa")
)

metadata_combined <- bind_rows(metadata_combined, missing_samples)

rm(missing_samples)

# Save processed metadata
write.csv(metadata_combined,
          file = "Data/processed_metadata.csv",
          row.names = FALSE)

################################################################################
# Looking at PCR prevalences 

pcr <- read.csv("Data/qpcr_results.csv")

pcr_prop <- pcr %>% 
  mutate(
    Genotype = case_when(
      Pf_Positive == 0 & Pv_Positive == 0 ~ "Negative",
      Pf_Positive == 0 & Pv_Positive == 1 ~ "Vivax",
      Pf_Positive == 1 & Pv_Positive == 0 ~ "Falciparum",
      Pf_Positive == 1 & Pv_Positive == 1 ~ "Mixed",
      .default = NA_character_
    )
  ) %>% 
  count(Genotype, name = "n") %>% 
  mutate(
    proportion = n / sum(n)
  )