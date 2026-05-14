################################################################################
############### Plotting drug-resistance prevalence Ethiopia  ##################
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

# Read metadata
# metadata <- read.csv("Outbreak_EPHI/metadata_outbreak_results/outbreakmetadatacomplete.csv") %>%
#   mutate(Region = case_when(Region == "DIREDAWA" ~ "Dire Dawa",
#                             Region == "Benishangul-Gumuz " ~ "Benishangul Gumz",
#                             Region == "TIGRAY" ~ "Tigray",
#                             Region == "AMHARA" ~ "Amhara",
#                             Region == "CENTRAL-ETHIOPIA" ~ "SNNP",
#                             Region == "AFAR" ~ "Afar",
#                             Region == "OROMIA" ~ "Oromia",
#                             Region == "S-WEST-ETHIOPIA" ~ "South West Ethiopia",
#                             Region == "SOUTH-ETHIOPIA" ~ "SNNP",
#                             Region == "SIDAMA" ~ "Sidama",
#                             Region == "GAMBELA" ~ "Gambela")) 

# Read genotype data
geno <- read.csv("Outbreak_EPHI/Assefa_Ethiopia_2026_repo/Data/outbreak_metadatacoigenotype-data.csv") %>% 
  select(study_id, Region, starts_with(c("crt","dhps","dhfr","mdr1","k13","Mutation_622I"))) %>% 
  rename(k13_622I = Mutation_622I) %>% 
  mutate(k13_622I = case_when(k13_622I == "Mutant" ~ 1,
                              k13_622I == "Wildtype" ~ 0,
                              .default = NA))

geno_summary <- geno %>% 
  # mutate(Region = case_when(Region == "DIRE DAWA" ~ "Dire Dawa",
  #                           Region == "TIGRAY" ~ "Tigray",
  #                           Region == "AMHARA" ~ "Amhara",
  #                           Region == "AFAR" ~ "Afar",
  #                           Region == "S_ETHIOPIA" ~ "SNNP",
  #                           Region == "SIDAMA" ~ "Sidama",
  #                           Region == "GAMBELA" ~ "Gambela")) %>%
  pivot_longer(cols = -c(study_id, Region),
               names_to = "Mutation",
               values_to = "Prevalence") %>%
  filter(!is.na(Prevalence)) %>%
  # We don't care about mixed infections here, so changing to MUTANT (1) and 
  # wild-type (0)
  mutate(Prevalence = ifelse(Prevalence > 0, 1, 0)) %>% 
  group_by(Mutation) %>%
  summarise(
    n = n(),
    x = sum(Prevalence),  # sum of proportions approximates total number of
    #mutation-carrying samples
    mean_prevalence = x/n,
    se = sqrt((mean_prevalence * (1 - mean_prevalence)) / n),
    ci_low = mean_prevalence - 1.96 * se,
    ci_high = mean_prevalence + 1.96 * se
  ) %>%
  ungroup() %>% 
  mutate(Drug = case_when(
    grepl("dhps",Mutation) ~ "sulfadoxine",
    grepl("mdr1",Mutation) ~ "lumefantrine",
    grepl("crt",Mutation) ~ "chloroquine",
    grepl("dhfr",Mutation) ~ "pyrimethamine",
    grepl("k13",Mutation) ~ "artemisinin"))  %>% 
  mutate(Drug = factor(Drug, levels = c("lumefantrine",
                                        "pyrimethamine",
                                        "sulfadoxine",
                                        "chloroquine",
                                        "artemisinin")))

# Define the desired prefix order
prefix_order <- c("mdr1", "dhfr", "dhps", "crt", "k13")

# Create an auxiliary column with the prefix for sorting
geno_summary <- geno_summary %>%
  mutate(Mutation_prefix = case_when(
    grepl("^mdr1", Mutation) ~ "mdr1",
    grepl("^dhfr", Mutation) ~ "dhfr",
    grepl("^dhps", Mutation) ~ "dhps",
    grepl("^crt", Mutation) ~ "crt",
    grepl("^k13", Mutation) ~ "k13")) %>%
  mutate(Mutation_prefix = factor(Mutation_prefix, levels = prefix_order)) %>%
  arrange(Mutation_prefix, Mutation) %>%
  mutate(Mutation = factor(Mutation, levels = unique(Mutation))) %>%
  select(-Mutation_prefix)

# Plotting
pdf("Outbreak_EPHI/metadata_outbreak_results/mutation_prevalence.pdf",
    width = 8, height = 4)

geno_summary %>% 
  ggplot(aes(x = Mutation, y = mean_prevalence)) +
  geom_bar(aes(fill = Drug), 
           stat = "identity", 
           alpha = 0.6) +
  scale_fill_manual(values = c("#D8B70A", "#02401B", "#377eb8", "grey", "#F2300F")) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), 
                width = 0.2, 
                size = 0.4, 
                color = "gray20") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Mutation", y = "Frequency (±95% CI)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()

####### Only looking at k13 mutations now
Mut_summary <- geno %>%
  mutate(Region = case_when(Region == "DIRE DAWA" ~ "Dire Dawa",
                            Region == "TIGRAY" ~ "Tigray",
                            Region == "AMHARA" ~ "Amhara",
                            Region == "AFAR" ~ "Afar",
                            Region == "S_ETHIOPIA" ~ "SNNP",
                            Region == "SIDAMA" ~ "Sidama",
                            Region == "GAMBELA" ~ "Gambela")) %>%
  pivot_longer(cols = -c(study_id, Region),
               names_to = "Mutation",
               values_to = "Prevalence") %>%
  filter(!is.na(Prevalence),
         grepl("k13",Mutation)) %>%
  # We don't care about mixed infections here, so changing to MUTANT (1) and 
  # wild-type (0)
  mutate(Prevalence = ifelse(Prevalence > 0, 1, 0)) %>% 
  group_by(Mutation) %>%
  summarise(
    n = n(),
    x = sum(Prevalence),
    mean_prevalence = x/n,
    se = sqrt((mean_prevalence * (1 - mean_prevalence)) / n),
    ci_low = mean_prevalence - 1.96 * se,
    ci_high = mean_prevalence + 1.96 * se
  ) %>%
  ungroup()

# Plotting only k13 mutations
pdf("Outbreak_EPHI/metadata_outbreak_results/k13_prevalence.pdf",
    width = 8, height = 5)

Mut_summary %>%
  ggplot(aes(x = Mutation, y = mean_prevalence)) +
  geom_bar(fill = "#F2300F",
           stat = "identity",
           alpha = 0.6) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                width = 0.2,
                size = 0.4,
                color = "gray20") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Mutation", y = "Frequency (±95% CI)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()

####### Only looking at 622I now
Mut_summary <- geno %>% 
  mutate(Region = case_when(Region == "DIRE DAWA" ~ "Dire Dawa",
                            Region == "TIGRAY" ~ "Tigray",
                            Region == "AMHARA" ~ "Amhara",
                            Region == "AFAR" ~ "Afar",
                            Region == "S_ETHIOPIA" ~ "SNNP",
                            Region == "SIDAMA" ~ "Sidama",
                            Region == "GAMBELA" ~ "Gambela")) %>%
  pivot_longer(cols = -c(study_id, Region),
               names_to = "Mutation",
               values_to = "Prevalence") %>%
  filter(!is.na(Prevalence),
         Mutation == "k13_622I") %>%
  # We don't care about mixed infections here, so changing to MUTANT (1) and 
  # wild-type (0)
  mutate(Prevalence = ifelse(Prevalence > 0, 1, 0)) %>% 
  group_by(Region) %>%
  summarise(
    n = n(),
    x = sum(Prevalence),
    mean_prevalence = x/n,
    se = sqrt((mean_prevalence * (1 - mean_prevalence)) / n),
    ci_low = mean_prevalence - 1.96 * se,
    ci_high = mean_prevalence + 1.96 * se
  ) %>%
  ungroup()

# Plotting only 622I by region
pdf("Outbreak_EPHI/metadata_outbreak_results/622I_prevalence.pdf",
    width = 4, height = 3)

Mut_summary %>% 
  ggplot(aes(x = Region, y = mean_prevalence)) +
  geom_bar(fill = "#F2300F", 
           stat = "identity", 
           alpha = 0.6) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), 
                width = 0.2, 
                size = 0.4, 
                color = "gray20") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Region", y = "Frequency (±95% CI)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()