################################################################################
############### Plotting hrp2/3 deletion prevalence Ethiopia  ##################
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

write_csv(deletion, file = "Data/processed_deletion_results.csv")

# Read sample weights
wts <- read.csv("Data/metadata_all_with_ipw.csv") %>% 
  select(Barcode, Region, ipw)

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
    weighted_prevalence = as.numeric(est),
    se = as.numeric(SE(est)),
    ci_low = ci[1],
    ci_high = ci[2]
  )
}) %>% 
  mutate(
    ci_low = pmax(ci_low, 0),
    Genotype = factor(Genotype, levels = c("pfhrp2+/3+","pfhrp2+/3-",
                                           "pfhrp2-/3+","pfhrp2-/3-"))
  )

write_csv(deletion_summary, file = "Results/deletion_prevalence_summary.csv")

# Plotting
pdf("deletion_prevalence.pdf", width = 4, height = 3)

deletion_summary %>% 
  ggplot(aes(x = Genotype, y = weighted_prevalence)) +
  geom_bar(aes(fill = Genotype), 
           stat = "identity", alpha = 0.8) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), 
                width = 0.2, 
                size = 0.4, 
                color = "gray20") +
  scale_fill_manual(values = c("black","black","black","#F2300F"),
                    guide = FALSE) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  labs(x = "Genotype", y = "Prevalence (CI)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 10))

dev.off()

# This script calculates IPW-weighted pfhrp2/3 genotype proportions by Region.

deletion_summary_region <- deletion_svy_data %>% 
  group_by(Region) %>% 
  group_modify(~ {
    
    genotype_levels <- unique(.x$Genotype)
    
    map_dfr(genotype_levels, function(g) {
      
      deletion_temp <- .x %>% 
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
        weighted_prevalence = as.numeric(est),
        se = as.numeric(SE(est)),
        ci_low = ci[1],
        ci_high = ci[2]
      )
    })
  }) %>% 
  ungroup() %>% 
  mutate(
    ci_low = pmax(ci_low, 0),
    Genotype = factor(Genotype, levels = c("pfhrp2+/3+","pfhrp2+/3-",
                                           "pfhrp2-/3+","pfhrp2-/3-"))
  )

write_csv(deletion_summary_region, 
          file = "Results/deletion_prevalence_summary_by_region.csv")


# Plotting deletion prevalence by region
pdf("double_del_prevalence_by_region.pdf", width = 6, height = 4)

deletion_summary_region %>%
  filter(Genotype == "pfhrp2-/3-") %>%
  ggplot(aes(x = Region, y = weighted_prevalence)) +
  geom_bar(aes(fill = Region),
           stat = "identity",
           alpha = 0.6) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                width = 0.2,
                size = 0.4) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  scale_fill_manual(
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
  labs(x = "Region", y = "Prevalence (CI)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()