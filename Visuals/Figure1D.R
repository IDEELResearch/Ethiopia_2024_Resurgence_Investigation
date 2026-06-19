################################################################################
#################### Plotting IBD network in Ethiopia  #########################
################################################################################

# Isabela Gerdes Gyuricza - Parr Lab
# 6/01/2026

# Cleaning up workspace
rm(list = ls())
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)

# Setting up working directory
setwd("Ethiopia_2024_Resurgence_Investigation_repo") 
# Update to your local path

# Loading libraries
library(tidyverse)
library(tidygraph)
library(ggraph)

# Loading IBD results
ibd <- read.csv("Data/ibd_outbreakSamples_final.csv")

annot <- read_csv("Data/processed_metadata.csv")

# Including 622I genotype information
geno <- read_csv("Data/DRgenotypes.csv") %>% 
  select(Sample_ID, pfk13_622I) %>% 
  mutate(`pfk13_622I genotype` = ifelse(pfk13_622I > 0, "Mutant",
                             ifelse(pfk13_622I == 0, "Wildtype",
                                    NA))) %>% 
  select(Sample_ID, `pfk13_622I genotype`)

annot <- annot %>% 
  left_join(geno)

# Now, let's focus on these clones (IBD >= 0.9)
ibd_90 <- ibd %>% 
  filter(ibd1 >= 0.9)

graph <- as_tbl_graph(ibd_90, directed = FALSE)

graph_c <- graph %>%
  activate(nodes) %>%
  left_join(annot, 
             by = c("name" = "Sample_ID")) %>% 
  filter(!is.na(`pfk13_622I genotype`))

pdf("09_IBD_MLE_network_by_region.pdf", width = 8, height = 6)

set.seed(1)
ggraph(graph_c,layout='igraph',
       algorithm ='fr') +
  geom_edge_link(color = "grey") +
  geom_node_point(aes(color = Region),
                  size = 4,
                  alpha = 0.6,
                  show.legend = TRUE,
                  position = position_jitter(width = 0.1,
                                             height = 0.1,
                                             seed = 1)) +
  scale_color_manual(name = "Region",
                     values = c('#1b9e77','#d95f02','#e7298a',
                                '#7570b3','black',
                                "blue")) + 
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, size = 15))

dev.off()

pdf("09_IBD_MLE_network_by_622I.pdf", width = 8, height = 6)

set.seed(1)
ggraph(graph_c,layout='igraph',
       algorithm ='fr') +
  geom_edge_link(color = "grey") +
  geom_node_point(aes(color = `pfk13_622I genotype`),
                  size = 4,
                  alpha = 0.6,
                  show.legend = TRUE,
                  position = position_jitter(width = 0.1,
                                             height = 0.1,
                                             seed = 1)) +
  scale_color_manual(values = c("#F2300F","black")) + 
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, size = 15))

dev.off()
