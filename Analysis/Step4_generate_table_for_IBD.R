#################### AA frequency table to genotype table for IBD ##############
################################################################################

# Isabela Gerdes Gyuricza - Parr Lab
# 6/23/2026

# Cleaning up workspace
rm(list = ls())
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)

# Setting up working directory
setwd("Ethiopia_2024_Resurgence_Investigation_repo") 
# Update to your local path

# Loading libraries needed
library(tidyverse)

# Loading amino acid frequency table
aa_freq <- read_tsv("Data/renamed_AAChangesInfo.tsv", show_col_types = FALSE) %>%
  # Everything that did not have coverage, we would like it to be NA
  mutate(
    alternate_AA_freq = ifelse(
      coverage_AA_cnt == 0,
      NA,
      alternate_AA_freq
    )
  )

# Loading GFF file
gff <- read_tsv(
  "Data/Pf3D7.gff",
  comment = "#",
  col_names = c(
    "chrom", "source", "type", "start", "end",
    "score", "strand", "phase", "attributes"
  ),
  show_col_types = FALSE
)

# Extracting values from the GFF attributes column
extract_gff_attribute <- function(attributes, key) {
  str_match(attributes, paste0(key, "=([^;]+)"))[, 2]
}

# Extracting gene ID and parent ID from the GFF attributes
gff_clean <- gff %>%
  mutate(
    ID = extract_gff_attribute(attributes, "ID"),
    Parent = extract_gff_attribute(attributes, "Parent")
  )

# Keeping protein coding genes with their coordinates
genes <- gff_clean %>%
  filter(type == "protein_coding_gene") %>%
  select(
    gene_id = ID,
    chrom,
    gene_start = start,
    gene_end = end,
    strand
  )

aa_freq %>% 
  left_join(genes, by = c("Gene_ID" = "gene_id"))

# Keeping CDS rows and assigning them to genes using Parent
cds <- gff_clean %>%
  filter(type == "CDS") %>%
  mutate(
    gene_id = str_remove(Parent, "\\.\\d+$")
  ) %>%
  select(
    gene_id,
    chrom,
    start,
    end,
    strand
  ) %>%
  filter(gene_id %in% genes$gene_id)

# Creating one row per coding base
cds_positions <- cds %>%
  rowwise() %>%
  mutate(
    genomic_pos = list(
      if (strand == "+") {
        start:end
      } else {
        end:start
      }
    )
  ) %>%
  unnest(genomic_pos) %>%
  ungroup() %>%
  group_by(gene_id) %>%
  arrange(
    if_else(strand == "+", genomic_pos, -genomic_pos),
    .by_group = TRUE
  ) %>%
  mutate(cds_pos = row_number()) %>%
  ungroup()

# Matching amino acid positions to the first base of each codon
aa_to_genome <- aa_freq %>%
  distinct(Gene_ID, reference_AA_pos) %>%
  mutate(
    codon_start_cds_pos = (reference_AA_pos - 1) * 3 + 1
  ) %>%
  left_join(
    cds_positions,
    by = c(
      "Gene_ID" = "gene_id",
      "codon_start_cds_pos" = "cds_pos"
    )
  ) %>%
  mutate(
    genome_snp = paste0(chrom, ":", genomic_pos)
  ) %>%
  select(
    Gene_ID,
    reference_AA_pos,
    genome_snp
  )

# Adding genomic SNP position to the amino acid frequency table
aa_freq_mapped <- aa_freq %>%
  left_join(
    aa_to_genome,
    by = c("Gene_ID", "reference_AA_pos")
  )

aa_freq_mapped <- aa_freq_mapped %>%
  mutate(
    genome_snp_alt = paste0(genome_snp, "_", reference_AA, ">", alternate_AA)
  )

aa_snp_matrix <- aa_freq_mapped %>%
  filter(!is.na(genome_snp_alt)) %>%
  select(sample, genome_snp_alt, alternate_AA_freq) %>%
  pivot_wider(
    names_from = genome_snp_alt,
    values_from = alternate_AA_freq,
    values_fill = NA
  )

# Saving output
write_csv(aa_snp_matrix, "Data/AA_frequency_genomic_SNP_matrix.csv")
