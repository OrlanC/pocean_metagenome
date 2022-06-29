#construcción de red teniendo en cuenta la paquetería WGCNA 
setwd("/Volumes/ADATA_cam/curso_morelia_27_06_22/matriz_oceanos/")
library(phyloseq)
library(WGCNA)
library(DESeq2)
library("limma")
library("edgeR")
library("statmod")
library("igraph")


datos <- import_biom("pacific_ocean.biom")


data<-datos@otu_table
datos <-as.data.frame(data)
datos


otus = datos[rowSums(datos >= 1) >=1,]

red_oceanos <- as.matrix(otus)

#normalizamos 

vst <- vst(red_oceanos)

library(WGCNA)
#Esta opción es muy importante para que WGCNA pueda trabajar:
options(stringsAsFactors = FALSE)



red <- (vst)



datExpr= t(red[,])
colnames(datExpr)= row.names(red)
rownames(datExpr)=colnames(vst)
dim(datExpr)


gsg=goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK 

if (!gsg$allOK)
{if (sum(!gsg$goodGenes)>0)
  printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse= ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse=", ")))
  datExpr= datExpr[gsg$goodSamples, gsg$goodGenes]
}


sampleTree = hclust
sampleTree = hclust(dist(datExpr), method = "average")


#Visualizamos como se comportan las bibliotecas
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

abline(h = 90, col = "red");

clust = cutreeStatic(sampleTree, cutHeight = 100, minSize = 10)
table(clust)



keepSamples = (clust==0)
datExpr = datExpr[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)


sampleTree2 = hclust(dist(datExpr), method = "average")
plot(sampleTree2, main = "Sample dendrogram", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)



powers = c(c(1:10), seq(from = 11, to=20, by=1))
# Llame a la función de análisis de topología de red
sft = pickSoftThreshold(datExpr, powerVector = powers,verbose = 5)
sft
###parametros de la imagen (plot)
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Índice de ajuste de topología sin escala en función de la potencia de umbral suave
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.70,col="red")


# Conectividad media en función de la potencia de umbral suave ---Tengo mis dudas-----
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


cor <- WGCNA::cor


?blockwiseModules
net = blockwiseModules(datExpr, 
                       #similarity matrix options 
                       corType = "pearson",
                       #opciones de matriz de adyacencia
                       power = 13,
                       networkType = "signed",
                       #opciones TOM
                       TOMType = "signed",
                       #Opciones de identificación de modulos
                       minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "ocean",
                       verbose = 3)

table(net$colors)

mergedColors = labels2colors(net$colors)
mergedColors
table(mergedColors)

colors <- net$colors

gen_hub <- chooseTopHubInEachModule(datExpr, 
                                    colors,
                                    omitColors = "grey",
                                    power = 10, 
                                    type = "signed")

gen_hub 

#?blockwiseModules
sizeGrWindow(12, 9)
mergedColors = labels2colors(net$colors)
mergedColors
table(mergedColors)
#Red: dacamos nuestro dendograma de la red
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];

#Cytoscape [2] permite al usuario ingresar un archivo de borde y un archivo de nodo, 
#lo que le permite especificar, por ejemplo, los pesos de los enlaces y los colores de los nodos. 
#Aquí demostramos la salida de dos módulos, el rojo y el marrón, a Cytoscape.

TOM = TOMsimilarityFromExpr(datExpr, power = 13, corType = "pearson", networkType = "signed");
?TOMsimilarityFromExpr

#?TOMsimilarityFromExpr
probes = colnames(datExpr)

# Exporte la red a archivos de lista de nodos y borde que Cytoscape puede leer

cyt = exportNetworkToCytoscape(TOM,
                               edgeFile = paste("CytoscapeInput-edges_corte_0.2", ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes_corte_0.2", ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.2,
                               nodeNames = probes,
                               nodeAttr = moduleColors);


