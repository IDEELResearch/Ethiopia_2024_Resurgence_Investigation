################################################################################
#################### Plotting IBD network in Ethiopia  #########################
################################################################################

# Isabela Gerdes Gyuricza - Parr Lab
# 7/30/2025

# Cleaning up workspace
rm(list = ls())
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)

# Setting up working directory
setwd("/Users/isabelagyuricza/OneDrive - University of North Carolina at Chapel Hill/IDEEL_PhD")

# Cleaning up workspace
rm(list = ls())
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)

# Loading libraries
library(tidyverse)
library(tidygraph)
library(ggraph)

# Loading IBD results
ibd <- read.csv("Outbreak_EPHI/Assefa_Ethiopia_2026_repo/Data/ibd_outbreakSamples.csv") %>% 
  select(p1, p2, ibd1) %>% 
  rename(malecotf = ibd1) %>% 
  unique()

# Loading annotation
annot <- read.csv("Outbreak_EPHI/metadata_outbreak_results/outbreak_metadatacoigenotype-data.csv") %>% 
  select(study_id, Region, Mutation_622I) %>% 
  mutate(Region = case_when(Region == "DIRE DAWA" ~ "Dire Dawa",
                            Region == "Benishangul-Gumuz " ~ "Benishangul Gumz",
                            Region == "TIGRAY" ~ "Tigray",
                            Region == "AMHARA" ~ "Amhara",
                            Region == "CENTRAL-ETHIOPIA" ~ "SNNP",
                            Region == "AFAR" ~ "Afar",
                            Region == "OROMIA" ~ "Oromia",
                            Region == "S-WEST-ETHIOPIA" ~ "South West Ethiopia",
                            Region == "S_ETHIOPIA" ~ "SNNP",
                            Region == "SIDAMA" ~ "Sidama",
                            Region == "GAMBELA" ~ "Gambela")) %>% 
  unique()

# Now, let's focus on these clones (IBD >= 0.9)
ibd_90 <- ibd %>% 
  filter(malecotf >= 0.9)

graph <- as_tbl_graph(ibd_90, directed = FALSE)

graph_c <- graph %>%
  activate(nodes) %>%
  left_join(annot, 
             by = c("name" = "study_id")) %>% 
  filter(!is.na(Region),
         !is.na(Mutation_622I))

pdf("Outbreak_EPHI/metadata_outbreak_results/09_IBD_MLE_network_by_region.pdf",
    width = 8, height = 6)

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

pdf("Outbreak_EPHI/metadata_outbreak_results/09_IBD_MLE_network_by_622I.pdf",
    width = 8, height = 6)

set.seed(1)
ggraph(graph_c,layout='igraph',
       algorithm ='fr') +
  geom_edge_link(color = "grey") +
  geom_node_point(aes(color = Mutation_622I),
                  size = 4,
                  alpha = 0.6,
                  show.legend = TRUE,
                  position = position_jitter(width = 0.1,
                                             height = 0.1,
                                             seed = 1)) +
  scale_color_manual(name = "622I genotype",
                     values = c("red","black")) + 
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, size = 15))


dev.off()
