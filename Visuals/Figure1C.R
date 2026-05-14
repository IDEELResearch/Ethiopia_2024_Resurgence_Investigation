################################################################################
############### Plotting hrp2/3 deletion prevalence Ethiopia  ##################
################################################################################

# Isabela Gerdes Gyuricza - Parr Lab
# 7/24/2025

# Cleaning up workspace
rm(list = ls())
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)

# Setting up working directory
setwd("/Users/isabelagyuricza/OneDrive - University of North Carolina at Chapel Hill/IDEEL_PhD")

# Load libraries
library(tidyverse)

# Read genotype data
deletion <- read.csv("Outbreak_EPHI/Assefa_Ethiopia_2026_repo/Data/outbreak_pfhrp2_3_combined_3.csv") %>% 
  select(Sample.Name, Pfhrp2_call, Pfhrp3_call) %>% 
  mutate(Pfhrp2_call = gsub("hr2","hrp2",Pfhrp2_call),
         Pfhrp3_call = gsub("hr3","hrp3",Pfhrp3_call)) %>% 
  mutate(Pfhrp2_call = ifelse(grepl("Possible|N/A",Pfhrp2_call), NA, Pfhrp2_call),
         Pfhrp3_call = ifelse(grepl("Possible|N/A",Pfhrp3_call), NA, Pfhrp3_call)) %>% 
  mutate(Pfhrp2_call = gsub("Pf","", Pfhrp2_call),
         Pfhrp3_call = gsub("Pf","", Pfhrp3_call)) %>% 
  mutate(Genotype = paste(Pfhrp2_call, Pfhrp3_call, sep = "/")) %>% 
  mutate(Genotype = ifelse(grepl("NA",Genotype),NA,Genotype))

deletion_summary <- deletion %>%
  filter(!is.na(Genotype)) %>%
  group_by(Genotype) %>%
  summarise(
    x = n()
  ) %>%
  ungroup() %>%
  mutate(
    n = sum(x),
    proportion = x / n,
    se = sqrt((proportion * (1 - proportion)) / n),
    ci_low = proportion - 1.96 * se,
    ci_high = proportion + 1.96 * se
  )

# Plotting
pdf("Outbreak_EPHI/metadata_outbreak_results/deletion_prevalence.pdf",
    width = 4, height = 3)

deletion_summary %>% 
  ggplot(aes(x = Genotype, y = proportion)) +
  geom_bar(stat = "identity", fill = "black", alpha = 0.5) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), 
                width = 0.2, 
                size = 0.4, 
                color = "gray20") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Genotype", y = "Frequency (±95% CI)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()
