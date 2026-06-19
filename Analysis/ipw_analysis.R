########################## INVERSE PROBABILITY WEIGHTING ########################
################################################################################

# Isabela Gerdes Gyuricza - Parr Lab
# 6/2/2026

# This script compares enrolled and sequenced participants and calculates 
# inverse probability of sequencing weights.

# Cleaning up workspace
rm(list = ls())
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)

# Setting up working directory
setwd("Ethiopia_2024_Resurgence_Investigation_repo") 
# Update to your local path

# Loading libraries needed
library(tidyverse)
library(gtsummary)
library(survey)
library(srvyr)

# Loading data
data_all <- read.csv("Data/processed_metadata.csv")

data_seq <- read.csv("Data/amino_acid_fracs.csv") %>% 
  mutate(Sequenced = TRUE) %>%
  select(Sample_ID, Sequenced) %>%
  rename(Barcode = Sample_ID)

df <- data_all %>% 
  left_join(data_seq) %>% 
  mutate(Sequenced = ifelse(!is.na(Sequenced), 1, 0),
    Sequenced = as.numeric(Sequenced)
  ) %>% 
  select(Barcode, Region, Age, Sex, Sequenced)

# Preparing variables for Table 1 and IPW analysis
df <- df %>% 
  mutate(
    Sequenced = as.numeric(Sequenced),
    Age = as.numeric(Age),
    Sex = as.factor(Sex),
    Region = as.factor(Region)
  )

# Creating Table 1 with all enrolled and sequenced participants

table_vars <- c("Age", "Sex", "Region")

# Table for all enrolled participants
table_all_enrolled <- df %>% 
  select(all_of(table_vars)) %>% 
  tbl_summary(
    statistic = list(
      Age ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    label = list(
      Age ~ "Age, median (IQR)",
      Sex ~ "Sex",
      Region ~ "Region"
    ),
    missing = "ifany",
    missing_text = "Missing"
  )

# Table for sequenced participants only
table_sequenced <- df %>% 
  filter(Sequenced == 1) %>% 
  select(all_of(table_vars)) %>% 
  tbl_summary(
    statistic = list(
      Age ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    label = list(
      Age ~ "Age, median (IQR)",
      Sex ~ "Sex",
      Region ~ "Region"
    ),
    missing = "ifany",
    missing_text = "Missing"
  )

# P-values comparing sequenced vs not sequenced participants
table_pvalues <- df %>% 
  select(Sequenced, all_of(table_vars)) %>% 
  tbl_summary(
    by = Sequenced,
    statistic = list(
      Age ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    missing = "no"
  ) %>% 
  add_p(
    test = list(
      Age ~ "wilcox.test",
      all_categorical() ~ "chisq.test"
    )
  )

# Merging all enrolled and sequenced columns with p-values
table_enrolled_sequenced <- tbl_merge(
  tbls = list(table_all_enrolled, table_sequenced),
  tab_spanner = c("All enrolled", "Sequenced")
) %>% 
  modify_table_body(
    ~ .x %>% 
      left_join(
        table_pvalues$table_body %>% 
          filter(row_type == "label") %>% 
          select(variable, row_type, p.value),
        by = c("variable", "row_type")
      )
  ) %>% 
  modify_header(
    label ~ "Variable",
    stat_0_1 ~ "All enrolled",
    stat_0_2 ~ "Sequenced",
    p.value ~ "p-value"
  ) %>% 
  modify_fmt_fun(
    p.value ~ gtsummary::style_pvalue
  ) %>% 
  bold_labels()

# Saving Table 1 as CSV
table_enrolled_sequenced_csv <- table_enrolled_sequenced %>% 
  as_tibble()

# Saving supplementary table
# write_csv(
#   table_enrolled_sequenced_csv,
#   "Data/enrolled_vs_sequenced_table.csv"
# )

# Creating propensity score model for probability of being sequenced
ps_model <- glm(
  Sequenced ~ Age + Sex + Region,
  family = binomial("logit"),
  data = df
)

# Removing rows that have NAs because they mess up the scores
df_ipw <- df %>% 
  filter(
    !is.na(Sequenced),
    !is.na(Age),
    !is.na(Sex),
    !is.na(Region)
  )

ps_model <- glm(
  Sequenced ~ Age + Sex + Region,
  family = binomial("logit"),
  data = df_ipw
)

df_ipw <- df_ipw %>% 
  mutate(
    ps = predict(ps_model, df_ipw, type = "response"),
    p_sequenced = mean(Sequenced == 1, na.rm = TRUE),
    ipw = case_when(
      Sequenced == 1 ~ p_sequenced / ps,
      Sequenced == 0 ~ (1 - p_sequenced) / (1 - ps),
      TRUE ~ NA_real_
    )
  )

# Saving dataframe with inverse probability weights
write_csv(
  df_ipw,
  "Data/metadata_all_with_ipw.csv"
)