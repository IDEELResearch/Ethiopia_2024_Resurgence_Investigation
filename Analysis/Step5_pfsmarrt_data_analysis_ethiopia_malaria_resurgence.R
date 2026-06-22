################################################################################
# Pfsmarrt: Plasmodium falciparum Amplicon Data analysis
# Script by: Abebe A. Fola
# Date: 05/13/2026
# Description: This script processes AA fraction data, merges with metadata, 
#              performs COI (Complexity of Infection) analysis and Identity By Descent (IBD) estimation using 
#              isoRelate,
################################################################################

# --- 1. Setup & Environment ---
setwd("Ethiopia_2024_Resurgence_Investigation_repo") 
# Update to your local path

message(paste("Current Working Directory:", getwd()))

# Load Libraries
library(tidyverse)
#v.2.0.0
library(isoRelate)
#v0.1.0
library(McCOILR)
#v1.3.1

# --- 2. Data Preprocessing & Genotype Selection ---
# Load AA fraction data (from SeekDeep output)
# Ref: https://github.com/bailey-lab/seekdeep_illumina_snakemake

df_raw <- read.csv("Data/amino_acid_fracs.csv", check.names = FALSE)

# Define target drug resistance loci mapping
# Format: "New_Name" = "Original_Name"
dr_mapping <- c(
  "Sample_ID"      = "Sample_ID",
  "pfcrt_76T"     = "pfcrt_76T",
  "pfdhfr_51I"     = "dhfr-51-59_51I",
  "pfdhfr_59R"    = "dhfr-51-59_59R",
  "pfdhfr_108N"   = "dhfr-108-164_108N",
  "pfdhps_437G"   = "dhps-436-437_437G",
  "pfdhps_540E"   = "dhps-540_540E",
  "pfdhps_581G"   = "dhps-581_581G",
  "pfmdr1_N86"    = "mdr1-86_86N",
  "pfmdr1_184F"   = "mdr1-184_184F",
  "pfmdr1_D1246"  = "mdr1-1246_1246D",
  "pfk13_433D"    = "k13-a_433D",
  "pfk13_441L"    = "k13-a_441L",
  "pfk13_574L"    = "k13-f_574L",
  "pfk13_622I"    = "k13-f_622I",
  "pfk13_658E"    = "k13-g_658E",
  "pfk13_675V"    = "k13-g_675V"
)

# Subset, Rename, and Round
# Note: 1 = pure mutant, 0 = wildtype, decimals = mixed, NA = missing
df_dr <- df_raw[, dr_mapping]
colnames(df_dr) <- names(dr_mapping)
df_dr[, -1] <- round(df_dr[, -1], 2) 

# Save cleaned genotype file

write.csv(df_dr, "Data/DRgenotypes.csv", row.names = FALSE) 
# Use these data for down stream prevalence estimation and visualization

# --- 4. COI Estimation (THE REAL McCOIL) ---
# Prepare data for McCOIL categorical tool. Only use diversity markers and dont 
# included drug resistance
data_mccoil_raw <- read.table("Data/COIFROMATED05.txt", fileEncoding = "UTF-16",
                              header = TRUE)
# generated from AA tables 

# Apply 0.5 threshold for mixed calls
data_mccoil_formatted <- data_mccoil_raw %>%
  mutate(across(-1, ~ ifelse(. > 0 & . < 1, 0.5, .)))

# Transforming column into samples names
data_mccoil_formatted <- data_mccoil_formatted %>% 
  column_to_rownames("Sample_ID")

# Run McCOIL Categorical
# Note: Ensure a 'COI' directory exists for the output path
if(!dir.exists("COI")) dir.create("COI")

set.seed(2024)
out_cat <- McCOIL_categorical(
  data_mccoil_formatted, 
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
# load data
barcode_raw <- read_csv("Data/ibd_formated1.csv")

# Replace NA and values < 1 with 0, only for numeric columns
barcode_clean <- barcode_raw %>%
  mutate(across(where(is.numeric), ~ ifelse(is.na(.) | . < 1, 0, .)))

col_names <- colnames(barcode_clean)
chr <- NULL
pos <- NULL

# Extract chr and position from SNP column names
for (i in 2:length(col_names)) {
  col_name_i <- unlist(strsplit(col_names[i], ":"))
  chr <- c(chr, col_name_i[1])
  pos <- c(pos, as.numeric(col_name_i[2]))
}

rm(i, col_name_i)

# Recombination rate: 1 event per 18,000 bp
recomb_rate <- 18000

# Create map file
my_map <- data.frame(
  chr = chr,
  snp_id = col_names[2:length(col_names)],
  pos_cM = pos / recomb_rate,
  pos_bp = pos,
  stringsAsFactors = FALSE
)

rm(chr, col_names, pos, recomb_rate)

# Convert first two columns to character (for PLINK format compliance)
my_map$chr <- as.character(my_map$chr)
my_map$snp_id <- as.character(my_map$snp_id)

# Create PED matrix: nrow = samples, ncol = 2*number of SNPs + 6 pedigree columns
my_ped <- matrix(nrow = nrow(barcode_clean), ncol = 2 * (ncol(barcode_clean) - 1) + 6)

# Fill genotypes (duplicate each SNP value for diploid format)
for (i in 1:nrow(barcode_clean)) {
  my_ped[i, seq(7, ncol(my_ped), 2)] <- unlist(barcode_clean[i, 2:ncol(barcode_clean)])
  my_ped[i, seq(8, ncol(my_ped), 2)] <- unlist(barcode_clean[i, 2:ncol(barcode_clean)])
}
rm(i)

# Replace values for PLINK-style encoding
my_ped[my_ped == 0] <- 2  # Treat "0" as alternate allele
my_ped[my_ped == -1] <- 0  # Missing

# Add pedigree columns
my_ped[, 1] <- barcode_clean[[1]]  # Family ID
my_ped[, 2] <- barcode_clean[[1]]  # Sample ID
my_ped[, 3] <- 0  # Paternal ID
my_ped[, 4] <- 0  # Maternal ID
my_ped[, 5] <- 1  # Sex (1 = male; change if needed)
my_ped[, 6] <- -9 # Phenotype

# Turn into a ibd_final frame if desired
my_ped <- as.data.frame(my_ped, stringsAsFactors = FALSE)

# Preview outputs
head(my_map)
head(my_ped[, 1:10])


# change ped and map to ibd_final frames
my_ped <- data.frame(my_ped)
for (i in c(3:ncol(my_ped))) {
  my_ped[,i] <- as.numeric(as.character(my_ped[,i])) # numeric
}
rm(i)

for (i in 1:2) {
  my_ped[,i] <- as.character(my_ped[,i]) # characters
}
rm(i)

colnames(my_ped)[1:6] <- c("fid", "iid", "pid", "mid", "moi", "aff")

# create pedmap list
barcode_clean_pedmap <- list(my_ped, my_map)

rm(my_ped, my_map)

# lets look at the ibd_final
str(barcode_clean_pedmap)

eth_genotypes <- getGenotypes(ped.map = barcode_clean_pedmap,
                              reference.ped.map = NULL,
                              maf = 0.0000,
                              isolate.max.missing = 0.0000, # save all samples
                              snp.max.missing = 0.0000, # save all loci
                              chromosomes = NULL,
                              input.map.distance = "M",
                              reference.map.distance = "cM")

# estimate ibd1 parameters
eth_parameters <- getIBDparameters(ped.genotypes = eth_genotypes,
                                   number.cores = 4)

# Saving data for plotting
eth_parameters <- eth_parameters %>% 
  select(iid1, iid2, ibd1) %>% 
  rename(p1 = iid1, 
         p2 = iid2)

write.csv(eth_parameters, "Data/ibd_outbreakSamples_final.csv", 
          row.names = FALSE)
# this output used for  Visualizing IBD Distributions and Relatedness network 
