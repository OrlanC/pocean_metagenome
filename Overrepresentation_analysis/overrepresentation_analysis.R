## Co-ocurrence networks with Ocean data

#This script is for microbial co-ocurrence networks in 3 different enviroments. 

### Load the libraries

library(phyloseq)
library(ggplot2)
library(scales)
library(RColorBrewer)
library(ggpubr)
library(readr)
library(patchwork)
library(plyr)
library(tidyverse)
library(cluster)
library(factoextra)
library(ggcorrplot)
library(corrplot)


### Read the biom file


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



## Make the abundance plots

#Make the object

# Agglomerate at Phylum level
phylum_global<- tax_glom(global, "Phylum")

#Relative abundance
relative_phylum_global<- transform_sample_counts(phylum_global, function(x) x / sum(x) ) 
rel_phyl_global_df<-psmelt(relative_phylum_global) # convert to data frame
rel_phyl_global_df$Phylum<- as.character(rel_phyl_global_df$Phylum) # convert to character type

# Agglomerate phylum with less than 0.01 abundance
rel_phyl_global_df$Phylum[rel_phyl_global_df$Abundance < 0.01] <- "Fila < 1% de abundancia" 

rel_phyl_global_df$Phylum<- as.factor(rel_phyl_global_df$Phylum) # convert to factor type
rel_phyl_global_df$Samples<- factor(rel_phyl_global_df$Samples, levels = c("SRR5788415", "SRR5788416", "SRR5788417", "SRR5788420", "SRR5788421", "SRR5788422"))  # Put Samples names in levels


#Make the abundance plot by samples.

phylum_colors<- brewer.pal(length(levels(rel_phyl_global_df$Phylum)),"Dark2") 
global_rel_phyl_plot<-ggplot(rel_phyl_global_df, aes(x=Samples, y=Abundance, fill=Phylum)) +
  geom_bar(aes(), stat="identity", position="stack") +
  scale_fill_manual(values=phylum_colors) + 
  labs(x= "Samples", y = "Relative abundance") +
  facet_grid(~Zone, scales = "free_x")+
  ggtitle("Relative abundance by Zones")+
  theme(axis.title = element_text(size = 14), axis.text.x = element_text(size = 8), plot.title = element_text(size = 20))

ggsave("~/dc_workshop/pocean_metagenome/images/abund_rel_phyl_global.svg", plot = global_rel_phyl_plot, width = 6.5, height = 3.65, units = "in", dpi = 400)
global_rel_phyl_plot


## Filter OTUS tables with absolute abundance by Zone

#Make the objects
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


#Create the otutable directory


#mkdir otutable


#Save otu tables without row names.
#write.csv(abs_phyl_global_df_E,paste("~/dc_workshop/pocean_metagenome/Overrepresentation_analysis/otutable/","Epipelagic",".csv",sep=""), row.names= FALSE)
#write.csv(abs_phyl_global_df_M,paste("~/dc_workshop/pocean_metagenome/Overrepresentation_analysis/otutable/","Mesopelagic",".csv",sep=""), row.names= FALSE)
#write.csv(abs_phyl_global_df_A,paste("~/dc_workshop/pocean_metagenome/Overrepresentation_analysis/otutable/","Abyssopelagic",".csv",sep=""), row.names= FALSE)


## Clusters

#Calculate the `ordinate` objet with Jaccard distance.


meta.ord <- ordinate(physeq = relative_phylum_global, method = "NMDS", 
                     distance = "jaccard")


#Plot the clusters.


rel_phylum_cluster_Jaccard <- plot_ordination(physeq = relative_phylum_global, ordination = meta.ord, color="Zone")+
  geom_text(aes(label=Samples), size=3)

ggsave("cluster_jaccard.svg", plot = rel_phylum_cluster_Jaccard, width = 6.5, height = 3.65, units = "in", dpi = 400)
rel_phylum_cluster_Jaccard


dist = "bray"
ord_meths = c("DCA", "CCA", "RDA", "DPCoA", "NMDS", "MDS", "PCoA")
plist = llply(as.list(ord_meths), function(i, physeq, dist){
  ordi = ordinate(physeq, method=i, distance=dist)
  plot_ordination(physeq, ordi, "Samples", color="Zone")
}, physeq=relative_phylum_global, dist)




## Cluster analysis using `cluster` and `factoExtra` packages.

#Filter the data set by samples and Phylum.

abs_filter <- daply(
  .data = rel_phyl_global_df,
  .variables = c("Samples", "Phylum"),
  .fun = function(x) sum(x$Abundance)
)
abs_filter <- as.data.frame(abs_filter)
abs_filter[is.na(abs_filter)] = 0


#Compute the distance matrix with `pearson` distance.


#Computing correlation based distance: pearson, spearman or kendall
dist_corr <- get_dist(abs_filter, method = "pearson")


#Plot distance matrix.

fviz_dist(dist.obj = dist_corr, lab_size = 8)


#Clusters with different methods


set.seed(12345)
hc_single <- hclust(d=dist_corr, method = "single")
hc_average <- hclust(d=dist_corr, method = "average")
hc_complete <- hclust(d=dist_corr, method = "complete")
hc_ward <- hclust(d=dist_corr, method = "ward.D2")
hc_median <- hclust(d=dist_corr, method = "median")
hc_centroid <- hclust(d=dist_corr, method = "centroid")


#Plot clusters with different methods.


fviz_dend(x = hc_average, k=3, 
          cex = 0.7,
          main = "",
          xlab = "Samples",
          ylab = "Distance",
          sub = "",
          horiz = TRUE) +
  geom_hline(yintercept = 0.01, linetype = "dashed") 


fviz_dend(x = hc_ward, k=3,
          cex = 0.7,
          main = "",
          xlab = "Samples",
          ylab = "Distance",
          sub = "",
          horiz = TRUE)+
  geom_hline(yintercept = 0.02, linetype = "dashed")




fviz_dend(x = hc_single, k=3,
          cex = 0.7,
          main = "",
          xlab = "Samples",
          ylab = "Distance",
          sub = "",
          horiz = TRUE)+
  geom_hline(yintercept = 0.01, linetype = "dashed")



fviz_dend(x = hc_complete, k=3,
          cex = 0.7,
          main = "",
          xlab = "Samples",
          ylab = "Distance",
          sub = "",horiz = TRUE)+
  geom_hline(yintercept = 0.018, linetype = "dashed")





fviz_dend(x = hc_median, k=3,
          cex = 0.7,
          main = "",
          xlab = "Samples",
          ylab = "Distance",
          sub = "",horiz = TRUE)+
  geom_hline(yintercept = 0.01, linetype = "dashed")




fviz_dend(x = hc_centroid, k=3,
          cex = 0.7,
          main = "",
          xlab = "Samples",
          ylab = "Distance",
          sub = "",horiz = TRUE)+
  geom_hline(yintercept = 0.01, linetype = "dashed")



#Compare methods.

# Methods
m <- c( "average", "single", "complete", "ward.D2", "median", "centroid")
names(m) <- c( "average", "single", "complete", "ward.D2", "median", "centroid")
# Compute the correlation coefficient 
coef_cor <- function(x) {
  cor(x=dist_corr, cophenetic(hclust(d=dist_corr, method = x)))
}
# Comparative table
coef_tabla <- map_dbl(m, coef_cor) 
coef_tabla


#PCA proyection.


fviz_cluster(object = list(data=abs_filter, cluster=cutree(hc_average, k=3)),
             ellipse.type = "convex", repel = TRUE, show.clust.cent = FALSE,
             labelsize = 8)  +
  labs(title = "Clusters + PCA proyection",
       subtitle = "Pearson distance, Linkage average, K=3") +
  theme_bw() +
  theme(legend.position = "bottom")


#Groups obtained.


grupos_average <-sort(cutree(hc_average, 3))
grupos_average


## Corplot


cor_base <- cor(abs_filter, use="complete.obs", method = "spearman")

cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
# matrix of the p-value of the correlation
p.mat <- cor.mtest(cor_base)

corr_plot <- corrplot(cor_base, type = "upper", order = "hclust", tl.col = "black", tl.srt = 45,
         p.mat = p.mat, sig.level = 0.05)



## References

#[1] Ma, B., Wang, Y., Ye, S. et al. Earth microbial co-occurrence network reveals interconnection pattern across microbiomes. Microbiome 8, 82 (2020). https://doi.org/10.1186/s40168-020-00857-2

#[2] C. Zirion. Zamia-bacterial-BGCs. https://czirion.github.io/Zamia-bacterial-BGCs/

#[3] Data processing and visualization for metagenomics. https://carpentries-incubator.github.io/metagenomics/
