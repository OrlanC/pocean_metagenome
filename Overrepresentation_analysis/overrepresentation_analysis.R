library(phyloseq)
library(ggplot2)
#library(dplyr)
library(scales)
library(RColorBrewer)
library(ggpubr)
library(readr)
library(patchwork)
library(tidyverse)

# Read the biom file
raw_biom <- import_biom("~/dc_workshop/pocean_metagenome/data/pacific_ocean.biom")

#Give name to taxonomic categories
raw_biom@tax_table@.Data<-substring(raw_biom@tax_table@.Data,4)
colnames(raw_biom@tax_table@.Data)<- c("Kingdom","Phylum", "Class", "Order","Family", "Genus", "Species")

# Read metadata
metadatos<-read.csv("~/dc_workshop/pocean_metagenome/data/metadata_sampling.csv") 
row.names(metadatos) <- metadatos$Samples
metadatos$Zone<-factor(metadatos$Zone, levels = c("Epipelagic", "Mesopelagic", "Abyssopelagic"))

# Generate otu_table
otu<-otu_table(raw_biom@otu_table@.Data, taxa_are_rows = TRUE) 
colnames(otu) <-metadatos$Samples

# Generate taxonomy table
tax<-tax_table(raw_biom@tax_table@.Data)

# Generate metadata table
metad<-sample_data(metadatos)

# It is the complete phyloseq object with the 3 previous tables
global<-phyloseq(otu, tax, metad)

# Put the complete species names in the column for Species (instead of the epithet only)
global@tax_table@.Data[,7]<-paste(global@tax_table@.Data[,6], global@tax_table@.Data[,7], sep= "_")


# Agglomerate at Phylum level
phylum_global<- tax_glom(global, "Phylum")

# Absolute abundance
abs_phyl_global_df<-psmelt(phylum_global) # convert to data frame

# Filter OTU table
abs_phyl_global_df_E <- abs_phyl_global_df %>%
  filter(abs_phyl_global_df$Zone == "Epipelagic") %>%
  select(OTU, Sample, Abundance) %>%
  pivot_wider(names_from = OTU, values_from = Abundance)

abs_phyl_global_df_M <- abs_phyl_global_df %>%
  filter(abs_phyl_global_df$Zone == "Mesopelagic")%>%
  select(OTU, Sample, Abundance) %>%
  pivot_wider(names_from = OTU, values_from = Abundance)

abs_phyl_global_df_A <- abs_phyl_global_df %>%
  filter(abs_phyl_global_df$Zone == "Abyssopelagic")%>%
  select(OTU, Sample, Abundance) %>%
  pivot_wider(names_from = OTU, values_from = Abundance)



#Relative abundance
relative_phylum_global<- transform_sample_counts(phylum_global, function(x) x / sum(x) ) #Transform to relative abundance
rel_phyl_global_df<-psmelt(relative_phylum_global) # convert to data frame

rel_phyl_global_df$Phylum<- as.character(rel_phyl_global_df$Phylum) # convert to character type
rel_phyl_global_df$Phylum[rel_phyl_global_df$Abundance < 0.01] <- "Fila < 1% de abundancia" # Agglomerate in one name all the phyla with little abundance
rel_phyl_global_df$Phylum<- as.factor(rel_phyl_global_df$Phylum) # convert to factor type
rel_phyl_global_df$Samples<- factor(rel_phyl_global_df$Samples, levels = c("SRR5788415", "SRR5788416", "SRR5788417", "SRR5788420", "SRR5788421", "SRR5788422"))  # Put Sample Type levels in chosen object to coincide with corresponding sample ID order

#Make the abundance plot by samples
phylum_colors<- brewer.pal(length(levels(rel_phyl_global_df$Phylum)),"Dark2") 
global_rel_phyl_plot<-ggplot(rel_phyl_global_df, aes(x=Samples, y=Abundance, fill=Phylum)) +
  geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values=phylum_colors) + 
  labs(x= "Muestras", y = "Abundancia relativa") +
  facet_grid(~Zone, scales = "free_x")

ggsave("abund_rel_phyl_global.svg", plot = global_rel_phyl_plot, width = 6.5, height = 3.65, units = "in", dpi = 400)
global_rel_phyl_plot
