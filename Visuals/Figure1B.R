################################################################################
############### Plotting drug-resistance prevalence Ethiopia  ##################
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
library(survey)

# This script calculates IPW-weighted mutation prevalence and confidence 
# intervals using survey weights.

# Read genotype data
geno <- read.csv("Data/DRgenotypes.csv") %>% 
  # Matching case on metadata
  mutate(Sample_ID = toupper(Sample_ID))

# Read sample weights
wts <- read.csv("Data/metadata_all_with_ipw.csv") %>% 
  select(Barcode, Region, ipw)

# Join weights to genotype data
geno <- geno %>% 
  inner_join(wts, by = c("Sample_ID" = "Barcode"))

# Convert genotype data to long format
geno_long <- geno %>% 
  pivot_longer(
    cols = -c(Sample_ID, Region, ipw),
    names_to = "Mutation",
    values_to = "Prevalence"
  ) %>%
  filter(!is.na(Prevalence),
         !is.na(ipw)) %>%
  mutate(
    Prevalence = ifelse(Prevalence > 0, 1, 0)
  )

# Create weighted survey design
geno_svy <- svydesign(
  ids = ~1,
  weights = ~ipw,
  data = geno_long
)

# Estimate weighted mutation prevalence with confidence intervals
geno_summary <- svyby(
  ~Prevalence,
  ~Mutation,
  design = geno_svy,
  FUN = svymean,
  vartype = c("se", "ci"),
  na.rm = TRUE
) %>% 
  as.data.frame() %>% 
  as_tibble() %>% 
  rename(
    weighted_prevalence = Prevalence,
    se = se,
    ci_low = ci_l,
    ci_high = ci_u
  ) %>% 
  left_join(
    geno_long %>% 
      group_by(Mutation) %>% 
      summarise(
        n = n(),
        x = sum(Prevalence),
        .groups = "drop"
      ),
    by = "Mutation"
  ) %>% 
  mutate(
    Drug = case_when(
      grepl("pfdhps", Mutation) ~ "sulfadoxine",
      grepl("pfmdr1", Mutation) ~ "lumefantrine",
      grepl("pfcrt", Mutation) ~ "chloroquine",
      grepl("pfdhfr", Mutation) ~ "pyrimethamine",
      grepl("pfk13", Mutation) ~ "artemisinin"
    ),
    Drug = factor(
      Drug,
      levels = c("lumefantrine",
                 "pyrimethamine",
                 "sulfadoxine",
                 "chloroquine",
                 "artemisinin")
    ))

# geno_summary_2 <- geno %>%
#   # mutate(Region = case_when(Region == "DIRE DAWA" ~ "Dire Dawa",
#   #                           Region == "TIGRAY" ~ "Tigray",
#   #                           Region == "AMHARA" ~ "Amhara",
#   #                           Region == "AFAR" ~ "Afar",
#   #                           Region == "S_ETHIOPIA" ~ "SNNP",
#   #                           Region == "SIDAMA" ~ "Sidama",
#   #                           Region == "GAMBELA" ~ "Gambela")) %>%
#   pivot_longer(cols = -c(study_id, Region),
#                names_to = "Mutation",
#                values_to = "Prevalence") %>%
#   filter(!is.na(Prevalence)) %>%
#   # We don't care about mixed infections here, so changing to MUTANT (1) and
#   # wild-type (0)
#   mutate(Prevalence = ifelse(Prevalence > 0, 1, 0)) %>%
#   group_by(Mutation) %>%
#   summarise(
#     n = n(),
#     x = sum(Prevalence),  # sum of proportions approximates total number of
#     #mutation-carrying samples
#     mean_prevalence = x/n,
#     se = sqrt((mean_prevalence * (1 - mean_prevalence)) / n),
#     ci_low = mean_prevalence - 1.96 * se,
#     ci_high = mean_prevalence + 1.96 * se
#   ) %>%
#   ungroup() %>%
#   mutate(Drug = case_when(
#     grepl("dhps",Mutation) ~ "sulfadoxine",
#     grepl("mdr1",Mutation) ~ "lumefantrine",
#     grepl("crt",Mutation) ~ "chloroquine",
#     grepl("dhfr",Mutation) ~ "pyrimethamine",
#     grepl("k13",Mutation) ~ "artemisinin"))  %>%
#   mutate(Drug = factor(Drug, levels = c("lumefantrine",
#                                         "pyrimethamine",
#                                         "sulfadoxine",
#                                         "chloroquine",
#                                         "artemisinin"))) %>%
#   mutate(Mutation = paste0("pf",Mutation))

# Define the desired prefix order
prefix_order <- c("pfmdr1", "pfdhfr", "pfdhps", "pfcrt", "pfk13")

# Create an auxiliary column with the prefix for sorting
geno_summary <- geno_summary %>%
  mutate(Mutation_prefix = case_when(
    grepl("^pfmdr1", Mutation) ~ "pfmdr1",
    grepl("^pfdhfr", Mutation) ~ "pfdhfr",
    grepl("^pfdhps", Mutation) ~ "pfdhps",
    grepl("^pfcrt", Mutation) ~ "pfcrt",
    grepl("^pfk13", Mutation) ~ "pfk13")) %>%
  mutate(Mutation_prefix = factor(Mutation_prefix, levels = prefix_order)) %>%
  arrange(Mutation_prefix, Mutation) %>%
  mutate(Mutation = factor(Mutation, levels = unique(Mutation))) %>%
  select(-Mutation_prefix)

# Plotting
pdf("mutation_prevalence.pdf",width = 8, height = 4)

geno_summary %>% 
  ggplot(aes(x = Mutation, y = weighted_prevalence)) +
  geom_bar(aes(fill = Drug), 
           stat = "identity", 
           alpha = 0.8) +
  scale_fill_manual(values = c('#d9d9d9','#969696','#525252','black', "#F2300F")) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), 
                width = 0.2, 
                size = 0.4, 
                color = "gray20") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Mutation", y = "Frequency (CI)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()


# Only looking at IPW-weighted k13 622I prevalence by region.

Mut_long <- geno %>% 
  mutate(
    Region = case_when(
      Region == "DIRE DAWA" ~ "Dire Dawa",
      Region == "TIGRAY" ~ "Tigray",
      Region == "AMHARA" ~ "Amhara",
      Region == "AFAR" ~ "Afar",
      Region == "S_ETHIOPIA" ~ "SNNP",
      Region == "SIDAMA" ~ "Sidama",
      Region == "GAMBELA" ~ "Gambela",
      TRUE ~ Region
    )
  ) %>%
  pivot_longer(
    cols = -c(Sample_ID, Region, ipw),
    names_to = "Mutation",
    values_to = "Prevalence"
  ) %>%
  filter(
    !is.na(Prevalence),
    !is.na(ipw),
    Mutation == "pfk13_622I"
  ) %>%
  mutate(
    Prevalence = ifelse(Prevalence > 0, 1, 0)
  )

Mut_svy <- svydesign(
  ids = ~1,
  weights = ~ipw,
  data = Mut_long
)

Mut_summary <- svyby(
  ~Prevalence,
  ~Region,
  design = Mut_svy,
  FUN = svymean,
  vartype = c("se", "ci"),
  na.rm = TRUE
) %>% 
  as.data.frame() %>% 
  as_tibble() %>% 
  rename(
    mean_prevalence = Prevalence,
    se = se,
    ci_low = ci_l,
    ci_high = ci_u
  ) %>% 
  left_join(
    Mut_long %>% 
      group_by(Region) %>% 
      summarise(
        n = n(),
        x = sum(Prevalence),
        .groups = "drop"
      ),
    by = "Region"
  ) %>% 
  select(Region, n, x, mean_prevalence, se, ci_low, ci_high)

# Plotting only 622I by region
pdf("622I_prevalence.pdf", width = 5, height = 3)

Mut_summary %>% 
  ggplot(aes(x = Region, y = mean_prevalence)) +
  geom_bar(aes(fill = Region), 
           stat = "identity", 
           alpha = 0.6) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), 
                width = 0.2, 
                size = 0.4) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(name = "Region",
                     values = c('#1b9e77','#d95f02','#e7298a','#66a61e',
                                '#7570b3','black',"blue")) + 
  labs(x = "Region", y = "Frequency (CI)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Mut_summary %>% 
#   ggplot(aes(x = Region, y = mean_prevalence)) +
#   geom_bar(fill = "#F2300F", 
#            stat = "identity", 
#            alpha = 0.6) +
#   geom_errorbar(aes(ymin = ci_low, ymax = ci_high), 
#                 width = 0.2, 
#                 size = 0.4, 
#                 color = "gray20") +
#   scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
#   labs(x = "Region", y = "Frequency (±95% CI)") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()

# Gathering some numbers 
geno_number <- geno %>% 
  pivot_longer(-c(Sample_ID, Region, -ipw), names_to = "Mutation",
               values_to = "Proportion") %>% 
  mutate(Mutant = ifelse(Proportion > 0, TRUE, FALSE)) %>% 
  filter(grepl("pfk13",Mutation)) %>% 
  group_by(Mutation) %>% 
  summarise(n_mutant = sum(Mutant, na.rm = TRUE))

geno_number <- geno %>% 
  pivot_longer(-c(Sample_ID, Region, -ipw), names_to = "Mutation",
               values_to = "Proportion") %>% 
  mutate(Mutant = ifelse(Proportion > 0, TRUE, FALSE)) %>% 
  filter(grepl("pfk13",Mutation)) %>% 
  group_by(Mutation, Region) %>% 
  summarise(n_mutant = sum(Mutant, na.rm = TRUE))