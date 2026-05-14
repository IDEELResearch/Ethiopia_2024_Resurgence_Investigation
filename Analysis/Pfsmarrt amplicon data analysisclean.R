################################################################################
# Pfsmarrt: Plasmodium falciparum Amplicon Data analysis
# Script by: Abebe A. Fola
# Date: 05/13/2026
# Description: This script processes AA fraction data, merges with metadata, 
#              performs COI (Complexity of Infection) analysis and Identity By Descent (IBD) estimation using 
#              isoRelate,
################################################################################

# --- 1. Setup & Environment ---
# setwd("C:/Users/afola/Desktop/malariaresugdata") # Update to your local path
message(paste("Current Working Directory:", getwd()))

# Load Libraries
library(tidyverse)
library(dplyr)
library(tidyr)
library(scales)
library(knitr)--
library(isoRelate)
library(igraph)
library(ggraph)
library(tidygraph)
library(cowplot)
library(ggplot2)
library(RColorBrewer)
library(wesanderson)
library(rhandsontable)
library(McCOILR)

# --- 2. Data Preprocessing & Genotype Selection ---
# Load AA fraction data (from SeekDeep output)
# Ref: https://github.com/bailey-lab/seekdeep_illumina_snakemake

df_raw <- read.csv("amino_acid_fracs.csv", check.names = FALSE)

# Define target drug resistance loci mapping
# Format: "New_Name" = "Original_Name"
dr_mapping <- c(
  "sample"      = "sample",
  "crt_76T"     = "pfcrt_76T",
  "dhf_51I"     = "dhfr-51-59_51I",
  "dhfr_59R"    = "dhfr-51-59_59R",
  "dhfr_108N"   = "dhfr-108-164_108N",
  "dhps_437G"   = "dhps-436-437_437G",
  "dhps_540E"   = "dhps-540_540E",
  "dhps_581G"   = "dhps-581_581G",
  "mdr1_N86"    = "mdr1-86_86N",
  "mdr1_184F"   = "mdr1-184_184F",
  "mdr1_D1246"  = "mdr1-1246_1246D",
  "k13_433D"    = "k13-a_433D",
  "k13_441L"    = "k13-a_441L",
  "k13_574L"    = "k13-f_574L",
  "k13_622I"    = "k13-f_622I",
  "k13_658E"    = "k13-g_658E",
  "k13_675V"    = "k13-g_675V"
)

# Subset, Rename, and Round
# Note: 1 = pure mutant, 0 = wildtype, decimals = mixed, NA = missing
df_dr <- df_raw[, dr_mapping]
colnames(df_dr) <- names(dr_mapping)
df_dr[, -1] <- round(df_dr[, -1], 2) 

# Save cleaned genotype file

write.csv(df_dr, "DRgenotypes.csv", row.names = FALSE) # Use these data for down stream prevalence estimation and visualiztion


# --- 3. Metadata Integration ---
metadata <- read.csv("Outbreak_Metadata_cleaned_short.csv") 

# --- 4. COI Estimation (THE REAL McCOIL) ---
# Prepare data for McCOIL categorical tool. Only use diversity markers and dont included drug resistance
data_mccoil_raw <- read.table("COIFROMATED05.txt", sep = "", header = TRUE) # generated from AA tables 

# Apply 0.5 threshold for mixed calls
data_mccoil_formatted <- data_mccoil_raw %>%
  mutate(across(-1, ~ ifelse(. > 0 & . < 1, 0.5, .)))

write.table(data_mccoil_formatted, "coiformated1.txt", sep = "\t", row.names = FALSE)

# Run McCOIL Categorical
# Note: Ensure a 'COI' directory exists for the output path
if(!dir.exists("COI")) dir.create("COI")

set.seed(2024)
out_cat <- McCOIL_categorical(
  data_mccoil_formatted[,-1], 
  maxCOI = 25, 
  threshold_ind = 20, 
  threshold_site = 20,
  totalrun = 1000, 
  burnin = 100, 
  M0 = 15, 
  e1 = 0.05, 
  e2 = 0.05,
  err_method = 3,
  path = "COI", 
  output = "Overall_COI.tsv"
)

# --- 5. IBD Analysis  ---



#  Data Preprocessing & isoRelate Formatting ---
# Load formatted AA table to isorelate format (each locus should be chromname:genomic coordinte see example Pf_3D7_chr04:115544)
# All loci included (both drug resistance and diversiyt markers to increase chance detecting estimated ibd across Pf genome)
barcode_raw <- read_csv("ibd_formated1.csv")

# Standardize values: NA -> -1 (missing), Decimals -> 0.5 (mixed), 0/1 (pure)
barcode_clean <- barcode_raw %>%
  mutate(across(where(is.numeric), ~ case_when(
    is.na(.)      ~ -1,
    . > 0 & . < 1 ~ 0.5,
    TRUE          ~ .
  )))

write_csv(barcode_clean, "ibd1_formated1_cleaned.csv")

# Prepare MAP file (Chr and Position)
col_names <- colnames(barcode_clean)[-1]
recomb_rate <- 18000 # P. falciparum: 1 event per 18kb

map_df <- data.frame(col_names) %>%
  separate(col_names, into = c("chr", "pos"), sep = ":", remove = FALSE) %>%
  mutate(
    pos = as.numeric(pos),
    pos_cM = pos / recomb_rate
  ) %>%
  select(chr, snp_id = col_names, pos_cM, pos_bp = pos)

# Prepare PED file (PLINK format)
# Note: isoRelate requires doubling the alleles for diploid-style input
ped_matrix <- matrix(nrow = nrow(barcode_clean), ncol = 2 * (ncol(barcode_clean) - 1) + 6)

for (i in 1:nrow(barcode_clean)) {
  genotypes <- unlist(barcode_clean[i, -1])
  ped_matrix[i, seq(7, ncol(ped_matrix), 2)] <- genotypes
  ped_matrix[i, seq(8, ncol(ped_matrix), 2)] <- genotypes
}

# Encode for isoRelate/PLINK: 0 -> 2 (Alt), -1 -> 0 (Missing)
ped_matrix[ped_matrix == 0] <- 2
ped_matrix[ped_matrix == -1] <- 0

# Add Pedigree metadata
ped_df <- as.data.frame(ped_matrix)
ped_df[, 1] <- barcode_clean[[1]] # FID
ped_df[, 2] <- barcode_clean[[1]] # IID
ped_df[, 3:4] <- 0                # PID/MID
ped_df[, 5]   <- 1                # MOI/Sex
ped_df[, 6]   <- -9               # Phenotype

colnames(ped_df)[1:6] <- c("fid", "iid", "pid", "mid", "moi", "aff")

# Save combined pedmap object
barcode_pedmap <- list(ped = ped_df, map = map_df)
saveRDS(barcode_pedmap, "barcode_clean_pedmap.rds")

# IBD Estimation ---
# Convert to isoRelate genotypes
eth_genotypes <- getGenotypes(
  ped.map = barcode_pedmap,
  maf = 0.0,
  isolate.max.missing = 0.1, 
  snp.max.missing = 0.1,
  input.map.distance = "M"
)

# Run IBD estimation
eth_parameters <- getIBDparameters(ped.genotypes = eth_genotypes, number.cores = 4)
write.csv(eth_parameters, "ibd_output_parameters.csv", row.names = FALSE)

# Metadata Integration ---
# Merging IBD results with COI, Drug Resistance (622I), and HRP2/3 Deletion data
metadata <- read.csv("isorelateformtedcoihrp23_finalMOD.csv")

# Create pairwise metadata merge
mtdt_p1 <- metadata %>% rename_with(~paste0(., "_p1"))
mtdt_p2 <- metadata %>% rename_with(~paste0(., "_p2"))

ibd_merged <- eth_parameters %>%
  left_join(mtdt_p1, by = c("p1" = "study_id_p1")) %>%
  left_join(mtdt_p2, by = c("p2" = "study_id_p2"))

write.csv(ibd_merged, "ibd_outbreakSamples.csv", row.names = FALSE)
# this output used for  Visualizing IBD Distributions and Relatedness network 
