################################################################################
# Pfsmarrt: Plasmodium falciparum Amplicon Data analysis
# Script by: Abebe A. Fola
# Date: 06/23/2026
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

df_raw <- read.table("Data/samp_amino_acid_fracs.tsv", 
                     sep = "\t", 
                     header = TRUE)

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

rm(list = ls())
# --- 4. COI Estimation (THE REAL McCOIL) ---
# Loading aminoacid data with information about ref and alt
aa_file <- read_tsv("Data/renamed_AAChangesInfo.tsv")

aa_file_freq <- aa_file %>%
  select(sample, CoveredBy, reference_AA_pos, alternate_AA, alternate_AA_freq, coverage_AA_cnt) %>%
  filter(grepl("heome|ama1", CoveredBy)) %>% # Keeping only heome and ama1 for COI 
  unite("aa_tmp", reference_AA_pos, alternate_AA, sep = "", remove = TRUE) %>%
  unite("aa", CoveredBy, aa_tmp, sep = "_", remove = TRUE) %>%
  # Formatting column with position + aminoacid name
  mutate(
    alternate_AA_freq = case_when(
      coverage_AA_cnt == 0 ~ -1, # positions that are not covered are missing data
      alternate_AA_freq > 0 & alternate_AA_freq < 1 ~ 0.5, # freq between 0-1 are heterozygous
      TRUE ~ alternate_AA_freq # the rest remains the same
    )
  ) %>% 
  select(-coverage_AA_cnt)

# Making sure the classes of the dataframe are correct for coi
aa_file_coi <- aa_file_freq %>%
  pivot_wider(names_from = aa, values_from = alternate_AA_freq) %>%
  mutate(across(-sample, ~replace_na(., -1))) %>%
  mutate(across(-sample, as.numeric)) %>%
  as.data.frame()

# Removing sample column and transforming into rownames
aa_file_coi_final = aa_file_coi[,-1]
rownames(aa_file_coi_final) = aa_file_coi[,1]

# Run McCOIL Categorical
# Note: Ensure a 'COI' directory exists for the output path
if(!dir.exists("Results/COI")) dir.create("Results/COI")

set.seed(2024)
out_cat <- McCOIL_categorical(
  aa_file_coi_final, 
  maxCOI = 25, 
  threshold_ind = 20, 
  threshold_site = 20,
  totalrun = 1000, 
  burnin = 100, 
  M0 = 15, 
  e1 = 0.05, 
  e2 = 0.05,
  err_method = 3,
  path = "Results/COI", 
  output = "Overall_COI.tsv"
)

coi <- read.table("Results/COI/Overall_COI.tsv_summary.txt", 
                  header = TRUE)

# Rounding McCOIL posterior median estimates to integer COI values
coi <- coi %>%
  filter(CorP == "C") %>% 
  mutate(COI = round(median))

# Summarizing COI values
coi_summary <- coi %>%
  count(COI, name = "n") %>%
  arrange(COI)

# Calculating monoclonal and polyclonal proportions
n_total <- sum(coi_summary$n)
mono_pct <- sum(coi_summary$n[coi_summary$COI == 1]) / n_total * 100
poly_pct <- 100 - mono_pct

# Plotting COI distribution
# ggplot(coi_summary,
#        aes(x = factor(COI),
#            y = n,
#            fill = ifelse(COI == 1, "Monoclonal", "Polyclonal"))) +
#   geom_col(width = 0.9) +
#   scale_fill_manual(
#     values = c("Monoclonal" = "black",
#                "Polyclonal" = "darkgrey"),
#     labels = c(
#       paste0("Monoclonal (", round(mono_pct, 1), "%)"),
#       paste0("Polyclonal (", round(poly_pct, 1), "%)")
#     ),
#     name = NULL
#   ) +
#   labs(
#     x = "Complexity of Infection (COI)",
#     y = "Number of Samples"
#   ) +
#   theme_bw() +
#   theme(
#     axis.title = element_text(face = "bold", size = 16),
#     axis.text = element_text(size = 13),
#     legend.position = "top",
#     legend.text = element_text(size = 14)
#   )

rm(list = ls())
# --- 5. IBD Analysis  ---
#  Data Preprocessing & isoRelate Formatting ---

# Loading amino acid frequency matrix formatted as samples by SNPs
barcode_clean <- read_csv("Data/AA_frequency_genomic_SNP_matrix.csv", 
                          show_col_types = FALSE) %>% 
  mutate(sample = toupper(sample))

# Converting amino acid frequencies to allele calls (presence/absence)
barcode_clean <- barcode_clean %>%
  mutate(
    across(
      where(is.numeric),
      ~ case_when(
        is.na(.) ~ -1,
        . < 1 ~ 0,
        TRUE ~ 1
      )
    )
  )

# Extracting chromosome and position from SNP column names
snp_info <- tibble(
  snp_id = colnames(barcode_clean)[-1],
  chr = str_extract(snp_id, "^[^:]+"),
  pos_bp = as.numeric(str_extract(snp_id, "(?<=:)\\d+"))
) %>%
  arrange(chr, pos_bp)

# Reordering SNP columns to match map order
barcode_clean <- barcode_clean %>%
  select(1, all_of(snp_info$snp_id))

# Creating map file for isoRelate
recomb_rate <- 18000

my_map <- snp_info %>%
  mutate(
    pos_cM = pos_bp / recomb_rate,
    chr = as.character(chr),
    snp_id = as.character(snp_id),
    pos_cM = as.numeric(pos_cM),
    pos_bp = as.numeric(pos_bp)
  ) %>%
  select(chr, snp_id, pos_cM, pos_bp) %>%
  as.data.frame(stringsAsFactors = FALSE)

# Creating genotype matrix
geno_matrix <- barcode_clean %>%
  select(-1) %>%
  as.matrix()

# Create PED matrix: nrow = samples, ncol = 2*number of SNPs + 6 pedigree columns
my_ped <- matrix(nrow = nrow(barcode_clean), ncol = 2 * (ncol(barcode_clean) - 1) + 6)

# Fill genotypes (duplicate each SNP value for diploid format)
for (i in 1:nrow(barcode_clean)) {
  my_ped[i, seq(7, ncol(my_ped), 2)] <- unlist(barcode_clean[i, 2:ncol(barcode_clean)])
  my_ped[i, seq(8, ncol(my_ped), 2)] <- unlist(barcode_clean[i, 2:ncol(barcode_clean)])
}
rm(i)

# Replace values for PLINK-style encoding
my_ped[my_ped == 1] <- 3   # Temporary alternate
my_ped[my_ped == 0] <- 2   # Reference
my_ped[my_ped == -1] <- 0  # Missing
my_ped[my_ped == 3] <- 1   # Alternate

# Add pedigree columns
my_ped[, 1] <- barcode_clean[[1]]  # Family ID
my_ped[, 2] <- barcode_clean[[1]]  # Sample ID
my_ped[, 3] <- 0  # Paternal ID
my_ped[, 4] <- 0  # Maternal ID
my_ped[, 5] <- 1  # Sex (1 = male; change if needed)
my_ped[, 6] <- -9 # Phenotype

# Turn into a ibd_final frame if desired
my_ped <- as.data.frame(my_ped, stringsAsFactors = FALSE)

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

# Running isoRelate genotype formatting
eth_genotypes <- getGenotypes(
  ped.map = barcode_clean_pedmap,
  reference.ped.map = NULL,
  maf = 0.0000,
  isolate.max.missing = 0.5,
  snp.max.missing = 0.5,
  chromosomes = NULL,
  input.map.distance = "cM",
  reference.map.distance = "cM"
)

# Estimating IBD parameters
eth_parameters <- getIBDparameters(
  ped.genotypes = eth_genotypes,
  number.cores = 4
)

# Saving IBD parameter output
write_csv(eth_parameters, "Data/ibd_output_parameters.csv")