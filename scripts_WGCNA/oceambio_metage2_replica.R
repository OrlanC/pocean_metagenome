##Kraken2
setwd("C:/Users/Orlando Camargo/Desktop/")
library("phyloseq")
library("ggplot2")
library("readr")
library("patchwork")

merged_metagenome_op<- import_biom("pacific_ocean.biom")
class(merged_metagenome_op)
View(merged_metagenome_op)
class(merged_metagenome_op@tax_table@.Data)
merged_metagenome_op@tax_table@.Data<-substring(merged_metagenome_op@tax_table@.Data,4)
colnames(merged_metagenome_op@tax_table@.Data)<- c("Kingdom","Phylum", "Class", "Order","Family", "Genus", "Species")
unique(merged_metagenome_op@tax_table@.Data[, "Phylum"])
merged_metagenome_op<-subset_taxa()
sum(merged_metagenome_op@tax_table@.Data[,"Phylum"] == "Firmicutes")
merged_metagenome_op <- subset_taxa(merged_metagenome_op, Kingdom == "Bacteria")
merged_metagenome_op
sample_sums(merged_metagenome_op)
summary(merged_metagenome_op@otu_table@.Data)
plot_richness(physeq = merged_metagenome_op, 
              measures = c("Observed","Chao1","Shannon"))
merged_metagenome_op@tax_table@.Data

###Next
## <- indica asignación de caracteres
View(merged_metagenome_op@otu_table@.Data)
sample_sums(x = merged_metagenome_op)
summary(merged_metagenome_op@otu_table@.Data)# reads por OTUs

deep <- data.frame(Samples = sample_names(merged_metagenome_op),
                   Reads = sample_sums(merged_metagenome_op))
deep

ggplot(data = deep, mapping = aes(x = Samples,y = Reads)) +
  geom_col()

summary(merged_metagenome_op@tax_table@.Data== "")

merged_metagenome_op <- subset_taxa(merged_metagenome_op, Genus != "")

#sampleType<-c("SRR5788415", "SRR5788416", "SRR5788417", "SRR5788420", "SRR5788421", "SRR5788422")


summary(merged_metagenome_op@tax_table@.Data== "")

head(merged_metagenome_op@otu_table@.Data)

percentages = transform_sample_counts(merged_metagenome_op, function(x) x*100 / sum(x) )

percentages[sam_data]<-sampleType

head(percentages@otu_table@.Data)

meta.ord <- ordinate(physeq = percentages, method = "NMDS", distance = "bray")

#Sacar los puntos de distancia en meta.ord
meta.ord[["points"]]

#plot sin nombres ni meta-datos
plot_ordination(physeq = percentages, ordination = meta.ord, color="sample_names(percentages)")

#plot con nombres y metadatos

Zone= c("Epipelagic", "Epipelagic", "Epipelagic", "Epipelagic", "Mesopelagic", "Abyssopelagic")
Depth = c("50 m","75 m","100 m","204 m","1023 m","5100 m")
Sample = c("SRR5788415","SRR5788416","SRR5788417","SRR5788420","SRR5788421","SRR5788422")
MDS1= c(-0.06813988,-0.0613852,-0.17871905,0.23363768,0.18848639,-0.11387995)
MDS2=c(-0.057142155, -0.061742261,-0.001971943,-0.01097201,0.031550929,0.10027744)

nmds_plot<-data.frame(Zone, Depth,Sample,MDS1,MDS2)
nmds_plot

p<-ggplot(nmds_plot, aes(x = MDS1, y = MDS2)) + 
  geom_point(size = 4, aes( shape = Zone, colour = Depth))+ 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "MDS1", colour = "Depth", y = "MDS2", shape = "Zone")  + 
  scale_color_brewer(palette="Dark2")


p


#summary (plot_ordination(physeq = percentages, ordination = meta.ord, color="percentages@otu_table@.Data("SRR5788415", "SRR5788416", "SRR5788417", "SRR5788420", "SRR5788421", "SRR5788422")"))

#data(GlobalPatterns)
#GP = prune_taxa(names(sort(taxa_sums(GlobalPatterns), TRUE)[1:50]), GlobalPatterns)
#gp_bray_pcoa = ordinate(GP, "CCA", "bray")
#plot_ordination(GP, gp_bray_pcoa, "samples", color="sampleType")

#### Correlation analysis ####
total = median(sample_sums(merged_metagenome_op))
percentage_01 <- filter_taxa(merged_metagenome_op, function(x) sum(x > total*0.001) > 0, TRUE)

plot_net(percentage_005, distance = "bray", type = "taxa", maxdist = 0.7, color="Class", point_label="Genus",shape="Phylum")

View()
plot_net(merged_metagenome_op, distance = "(A+B-2*J)/(A+B)", type = "taxa", maxdist = 0.8, color="Class", point_label="Genu")


####

glom <- tax_glom(percentages, taxrank = 'Phylum')
View(glom@tax_table@.Data)
percent <- psmelt(glom)
str(percent)
head(percent)

raw <- tax_glom(physeq = merged_metagenome_op, taxrank = "Phylum")
raw.data <- psmelt(raw)
str(raw.data)

library(tidyr)
samples_new<- percent$Sample%>%gsub(pattern="SRR5788416",replacement = "75m")%>%gsub(pattern="SRR5788415",replacement = "100m")%>%gsub(pattern="SRR5788417",replacement = "50m")%>%gsub(pattern="SRR5788420",replacement = "5100m")%>%gsub(pattern="SRR5788421",replacement = "1023m")%>%gsub(pattern="SRR5788422",replacement = "204m")
samples_new

percent$Sample<-samples_new
raw.data$Sample<-samples_new


raw.plot <- ggplot(data=raw.data, aes(x=Sample, y=Abundance, fill=Phylum))+ 
  geom_bar(aes(), stat="identity", position="stack")

rel.plot <- ggplot(data=percent, aes(x=Sample, y=Abundance, fill=Phylum))+ 
  geom_bar(aes(), stat="identity", position="stack")
raw.plot | rel.plot

percent$Phylum <- as.character(percent$Phylum)
percent$Phylum[percent$Abundance < 0.5] <- "Phyla < 0.5% abund."
unique(percent$Phylum)
rel.plot <- ggplot(data=percent, aes(x=Sample, y=Abundance, fill=Phylum))+ 
  geom_bar(aes(), stat="identity", position="stack")
raw.plot | rel.plot

merged_metagenomes_bacter <- subset_taxa(merged_metagenome_op, Kingdom == "Bacteria")

cyanos <- subset_taxa(merged_metagenome_op, Phylum == "Cyanobacteria")
unique(cyanos@tax_table@.Data[,2])

cyanos  = transform_sample_counts(cyanos, function(x) x*100 / sum(x) )
glom <- tax_glom(cyanos, taxrank = "Genus")
g.cyanos <- psmelt(glom)
sample_new_cyanos<- g.cyanos$Sample%>%gsub(pattern="SRR5788416",replacement = "75m")%>%gsub(pattern="SRR5788415",replacement = "100m")%>%gsub(pattern="SRR5788417",replacement = "50m")%>%gsub(pattern="SRR5788420",replacement = "5100m")%>%gsub(pattern="SRR5788421",replacement = "1023m")%>%gsub(pattern="SRR5788422",replacement = "204m")
g.cyanos$Sample<-sample_new_cyanos
g.cyanos$Genus[g.cyanos$Abundance < 5] <- "Genera < 5.0 abund"
p.cyanos <- ggplot(data=g.cyanos, aes(x=Sample, y=Abundance, fill=Genus))+ 
  geom_bar(aes(), stat="identity", position="stack")
p.cyanos

proteo <- subset_taxa(merged_metagenome_op, Phylum == "Proteobacteria")
unique(proteo@tax_table@.Data[,2])

proteo  = transform_sample_counts(proteo, function(x) x*100 / sum(x) )
glom <- tax_glom(proteo, taxrank = "Genus")
g.proteo <- psmelt(glom)
sample_new_proteo<- g.proteo$Sample%>%gsub(pattern="SRR5788416",replacement = "75m")%>%gsub(pattern="SRR5788415",replacement = "100m")%>%gsub(pattern="SRR5788417",replacement = "50m")%>%gsub(pattern="SRR5788420",replacement = "5100m")%>%gsub(pattern="SRR5788421",replacement = "1023m")%>%gsub(pattern="SRR5788422",replacement = "204m")
g.proteo$Sample<-sample_new_proteo
g.proteo$Genus[g.proteo$Abundance < 2] <- "Genera < 2.0 abund"
p.proteo <- ggplot(data=g.proteo, aes(x=Sample, y=Abundance, fill=Genus))+ 
  geom_bar(aes(), stat="identity", position="stack")
p.proteo

p.cyanos|p.proteo

firm <- subset_taxa(merged_metagenome_op, Phylum == "Firmicutes")
unique(firm@tax_table@.Data[,2])

firm  = transform_sample_counts(firm, function(x) x*100 / sum(x) )
glom <- tax_glom(firm, taxrank = "Genus")
g.firm <- psmelt(glom)
sample_new_firm<- g.firm$Sample%>%gsub(pattern="SRR5788416",replacement = "75m")%>%gsub(pattern="SRR5788415",replacement = "100m")%>%gsub(pattern="SRR5788417",replacement = "50m")%>%gsub(pattern="SRR5788420",replacement = "5100m")%>%gsub(pattern="SRR5788421",replacement = "1023m")%>%gsub(pattern="SRR5788422",replacement = "204m")
g.firm$Sample<-sample_new_firm
g.firm$Genus[g.firm$Abundance < 2] <- "Genera < 2.0 abund"
p.firm <- ggplot(data=g.firm, aes(x=Sample, y=Abundance, fill=Genus))+ 
  geom_bar(aes(), stat="identity", position="stack")
p.firm

actino <- subset_taxa(merged_metagenome_op, Phylum == "Actinobacteria")
unique(actino@tax_table@.Data[,2])

actino  = transform_sample_counts(actino, function(x) x*100 / sum(x) )
glom <- tax_glom(actino, taxrank = "Genus")
g.actino <- psmelt(glom)
sample_new_actino<- g.actino$Sample%>%gsub(pattern="SRR5788416",replacement = "75m")%>%gsub(pattern="SRR5788415",replacement = "100m")%>%gsub(pattern="SRR5788417",replacement = "50m")%>%gsub(pattern="SRR5788420",replacement = "5100m")%>%gsub(pattern="SRR5788421",replacement = "1023m")%>%gsub(pattern="SRR5788422",replacement = "204m")
g.actino$Sample<-sample_new_actino
g.actino$Genus[g.actino$Abundance < 2] <- "Genera < 2.0 abund"
p.actino <- ggplot(data=g.actino, aes(x=Sample, y=Abundance, fill=Genus))+ 
  geom_bar(aes(), stat="identity", position="stack")
p.actino

p.firm|p.actino

########################################################
##Kraken2
setwd("~/dc_workshop/taxonomy/")
library("phyloseq")
library("ggplot2")
library("readr")
library("patchwork")
merged_metagenomes<- import_biom("cuatroc.bio")
class(merged_metagenomes)
View(merged_metagenomes)
class(merged_metagenomes@tax_table@.Data)
merged_metagenomes@tax_table@.Data<-substring(merged_metagenomes@tax_table@.Data,4)
colnames(merged_metagenomes@tax_table@.Data)<- c("Kingdom","Phylum", "Class", "Order","Family", "Genus", "Species")
unique(merged_metagenomes@tax_table@.Data[, "Phylum"])
merged_metagenomes<-subset_taxa()
sum(merged_metagenomes@tax_table@.Data[,"Phylum"] == "Firmicutes")
merged_metagenomes <- subset_taxa(merged_metagenomes, Kingdom == "Bacteria")
merged_metagenomes
sample_sums(merged_metagenomes)
summary(merged_metagenomes@otu_table@.Data)
plot_richness(physeq = merged_metagenomes, 
              measures = c("Observed","Chao1","Shannon"))

###Next
## <- indica asignación de caracteres
View(merged_metagenomes@otu_table@.Data)
sample_sums(x = merged_metagenomes)
summary(merged_metagenomes@otu_table@.Data)# reads por OTUs

deep <- data.frame(Samples = sample_names(merged_metagenomes),
                   Reads = sample_sums(merged_metagenomes))
deep

ggplot(data = deep, mapping = aes(x = Samples,y = Reads, fill=drv) +
         geom_col()
       
       ggplot(data = deep, mapping = aes(x = Samples,y = Reads)) +
         geom_point()
       
       ggplot(data = deep, mapping = aes(x = Samples,y = Reads)) +
         geom_count()
       
       summary(merged_metagenomes@tax_table@.Data== "")
       merged_metagenomes <- subset_taxa(merged_metagenomes, Genus != "")
       summary(merged_metagenomes@tax_table@.Data== "")
       percentages  = transform_sample_counts(merged_metagenomes, function(x) x*100 / sum(x) )
       head(percentages@otu_table@.Data)
       
       meta.ord <- ordinate(physeq = percentages, method = "NMDS", 
                            distance = "bray")
       head(meta.ord)
       plot_ordination(physeq = percentages, ordination = meta.ord)
       glom <- tax_glom(percentages, taxrank = 'Phylum')
       View(glom@tax_table@.Data)
       
       #pasamos a df para hacer un graf ggplot
       percentages <- psmelt(glom)
       str(percentages)
       
       raw <- tax_glom(physeq = merged_metagenomes, taxrank = "Phylum")
       raw.data <- psmelt(raw)
       str(raw.data)
       
       raw.plot <- ggplot(data=raw.data, aes(x=Sample, y=Abundance, fill=Phylum))+ 
         geom_bar(aes(), stat="identity", position="stack")
       rel.plot <- ggplot(data=percentages, aes(x=Sample, y=Abundance, fill=Phylum))+ 
         geom_bar(aes(), stat="identity", position="stack")
       raw.plot | rel.plot
       
       percentages$Phylum <- as.character(percentages$Phylum)
       percentages$Phylum[percentages$Abundance < 0.5] <- "Phyla < 0.5% abund."
       unique(percentages$Phylum)
       
       rel.plot <- ggplot(data=percentages, aes(x=Sample, y=Abundance, fill=Phylum))+ 
         geom_bar(aes(), stat="identity", position="stack")
       raw.plot | rel.plot
       
       merged_metagenomes <- subset_taxa(merged_metagenomes, Kingdom == "Bacteria")
       
       cyanos <- subset_taxa(merged_metagenomes, Phylum == "Cyanobacteria")
       unique(cyanos@tax_table@.Data[,2])
       
       cyanos  = transform_sample_counts(cyanos, function(x) x*100 / sum(x) )
       glom <- tax_glom(cyanos, taxrank = "Genus")
       g.cyanos <- psmelt(glom)
       g.cyanos$Genus[g.cyanos$Abundance < 10] <- "Genera < 10.0 abund"
       p.cyanos <- ggplot(data=g.cyanos, aes(x=Sample, y=Abundance, fill=Genus))+ 
         geom_bar(aes(), stat="identity", position="stack")
       p.cyanos
       
       firmicutes <- subset_taxa(merged_metagenomes, Phylum == "Firmicutes")
       unique(firmicutes@tax_table@.Data[,2])
       firmicutes  = transform_sample_counts(firmicutes, function(x) x*100 / sum(x) )
       glom <- tax_glom(firmicutes, taxrank = "Genus")
       g.firmicutes <- psmelt(glom)
       g.firmicutes$Genus[g.firmicutes$Abundance < 5] <- "Genera < 10.0 abund"
       p.firmicutes <- ggplot(data=g.firmicutes, aes(x=Sample, y=Abundance, fill=Genus))+ 
         geom_bar(aes(), stat="identity", position="stack")
       p.firmicutes
       