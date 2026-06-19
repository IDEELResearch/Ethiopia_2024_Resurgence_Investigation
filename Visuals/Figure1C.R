################################################################################
############### Plotting hrp2/3 deletion prevalence Ethiopia  ##################
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

# Read deletion genotype data
deletion <- read.csv("Data/outbreak_pfhrp2_3_combined_3.csv") %>% 
  select(Sample_ID, Pfhrp2_call, Pfhrp3_call) %>% 
  mutate(Pfhrp2_call = gsub("hr2", "hrp2", Pfhrp2_call),
         Pfhrp3_call = gsub("hr3", "hrp3", Pfhrp3_call)) %>% 
  mutate(Pfhrp2_call = ifelse(grepl("Possible|N/A", Pfhrp2_call), NA, Pfhrp2_call),
         Pfhrp3_call = ifelse(grepl("Possible|N/A", Pfhrp3_call), NA, Pfhrp3_call)) %>% 
  mutate(Pfhrp2_call = gsub("Pf", "", Pfhrp2_call),
         Pfhrp3_call = gsub("Pf", "", Pfhrp3_call)) %>% 
  mutate(Genotype = paste(Pfhrp2_call, Pfhrp3_call, sep = "/")) %>% 
  mutate(Genotype = ifelse(grepl("NA", Genotype), NA, Genotype)) %>% 
  mutate(Genotype = case_when(Genotype == "hrp2-/hrp3-" ~ "pfhrp2-/3-",
                              Genotype == "hrp2-/hrp3+" ~ "pfhrp2-/3+",
                              Genotype == "hrp2+/hrp3-" ~ "pfhrp2+/3-",
                              Genotype == "hrp2+/hrp3+" ~ "pfhrp2+/3+",
                              .default = NA)) %>% 
  mutate(Sample_ID = toupper(Sample_ID))

# Read sample weights
wts <- read.csv("Data/metadata_all_with_ipw.csv") %>% 
  select(Barcode, ipw)

# Join weights to deletion data
deletion <- deletion %>% 
  left_join(wts, by = c("Sample_ID" = "Barcode"))

# Keep samples with genotype calls and weights
deletion_svy_data <- deletion %>% 
  filter(!is.na(Genotype),
         !is.na(ipw))

# This script calculates IPW-weighted pfhrp2/3 genotype proportions one genotype at a time.

genotype_levels <- unique(deletion_svy_data$Genotype)

deletion_summary <- map_dfr(genotype_levels, function(g) {
  
  deletion_temp <- deletion_svy_data %>% 
    mutate(is_genotype = ifelse(Genotype == g, 1, 0))
  
  deletion_temp_svy <- svydesign(
    ids = ~1,
    weights = ~ipw,
    data = deletion_temp
  )
  
  est <- svymean(
    ~is_genotype,
    design = deletion_temp_svy,
    na.rm = TRUE
  )
  
  ci <- confint(est)
  
  tibble(
    Genotype = g,
    n = nrow(deletion_temp),
    x = sum(deletion_temp$Genotype == g),
    proportion = as.numeric(est),
    se = as.numeric(SE(est)),
    ci_low = ci[1],
    ci_high = ci[2]
  )
})

# deletion_summary <- deletion %>%
#   filter(!is.na(Genotype)) %>%
#   group_by(Genotype) %>%
#   summarise(
#     x = n()
#   ) %>%
#   ungroup() %>%
#   mutate(
#     n = sum(x),
#     proportion = x / n,
#     se = sqrt((proportion * (1 - proportion)) / n),
#     ci_low = proportion - 1.96 * se,
#     ci_high = proportion + 1.96 * se
#   )

# Plotting
pdf("deletion_prevalence.pdf", width = 4, height = 3)

deletion_summary %>% 
  ggplot(aes(x = Genotype, y = proportion)) +
  geom_bar(stat = "identity", fill = "black", alpha = 0.5) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), 
                width = 0.2, 
                size = 0.4, 
                color = "gray20") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "Genotype", y = "Frequency (CI)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 10))

dev.off()
