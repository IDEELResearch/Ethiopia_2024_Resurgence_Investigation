# Ethiopia_2024_Resurgence_Investigation

This repository hosts the data files and analysis scripts required to reproduce the findings of the manuscript *“Resurgence of Malaria in Ethiopia (2024).*”

It includes the complete workflow for estimating drug resistance prevalence, generating geographic distributions, and visualizing population genetics. Access to the repository is unrestricted.

**Repository Contents:**
- `/Data`: Processed data files as input files.
- `/Analysis`: Scripts for running SeekDeep, performing inverse probability weighting, and calculating drug resistance prevalence .
- `/Visuals`: Scripts for summarizing data and generating figures.

**Additional bioinformatics pipeline:**
- Primary processing: Pfsmarrter amplicon sequencing data were processed using [SeekDeep](https://github.com/bailey-lab/SeekDeep) to ensure high-quality variant calling.
- Variant conversion: Output files were converted using the [haplotype_variant_calling](https://github.com/bailey-lab/haplotype_variant_calling) pipeline into:
  - Amino acid tables: For antimalarial drug resistance marker prevalence.
  - VCF files: For downstream population genetic analyses.
- Genetic analysis: Identity-by-descent (IBD) estimation was performed via isoRelate R package (v0.1.0).

Raw FASTQ files are available through the SRA under accession PRJNA1465284.

