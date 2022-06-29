#script tutorial
#install.packages("devtools")
#install.packages("BiocManager")
install.packages("WGCNA")
install.packages("qgraph")
install.packages("NetCoMi")

library(WGCNA)
# Install NetCoMi
devtools::install_github("stefpeschel/NetCoMi", 
                         dependencies = c("Depends", "Imports", "LinkingTo"),
                         repos = c("https://cloud.r-project.org/",
                                   BiocManager::repositories()))

devtools::install_github("GraceYoon/SPRING", force = TRUE)
devtools::install_github("zdk123/SpiecEasi", force = TRUE)

install.packages("remotes")
remotes::install_github("stefpeschel/NetCoMi", force = TRUE)
library(NetCoMi)

#Tratamos los datos del tutorial 

data("amgut1.filt")
data("amgut2.filt.phy")


#Algunas consideraciones de los comandos a tomar en cuenta:
#Especificaciones de los comandos:
#>netConstruct().   sirve para filtrar los datos de la red
#>filtSamp.        Solo se incluyen muestras con un número total de lecturas ejemplo: de al menos 1000
#>filtTax.        El comando manada a llamar solo al #NUMERO DE TAXONES CON MAYOR FRECUENCIA EJEMPLO: MAYOR A 50
#>measure.         define la medida de asociación o disimilitud, que es "spring"en nuestro caso.
#>measurePar. nlambday rep.numse establecen en 10 para reducir el tiempo de ejecución,
# pero deberían ser más altos para datos reales.
#>La normalización, así como el manejo de cero, se realizan internamente en SPRING(). 
#Por lo tanto, establecemos normMethody zeroMethodpara "none".
#>Además, configuramos porque sparsMethod = "none"
#devuelve una red dispersa en la que no es necesario ningún paso de dispersión adicional.
#>dissFunc.     Este comando es para seleccionar si la red va ser firmada, no firmada o hibrida (signed, unsigned, hybrid)
#en el ejemplo si realiza dissFunc = "signed" ------firmado
#El argumento: verbose se establece en 3 para que netConstruct() se impriman todos los mensajes generados por funciones externas verbose = 3

net_single <- netConstruct(amgut1.filt,
                           filtTax = "highestFreq",
                           filtTaxPar = list(highestFreq = 50),
                           filtSamp = "totalReads",
                           filtSampPar = list(totalReads = 1000),
                           measure = "spring",
                           measurePar = list(nlambda=10, 
                                             rep.num=10),
                           normMethod = "none", 
                           zeroMethod = "none",
                           sparsMethod = "none", 
                           dissFunc = "signed",
                           verbose = 3,
                           seed = 123456)

#Analizando la red construida
#netAnalyze().   sirve para analizar la red construida, algunas metricas que da el resultado es:
#comandos:
#centrLCC.   --------se establece en el TRUE significado de que las centralidades se calculan solo para los nodos en el componente conectado más grande (LCC).
#cluster_fast_greedy(). ------Los clústeres se identifican mediante la optimización de la modularidad codiciosa, toma en cuenta la función igraph
#hubPar.  ----- esta función nos da los nodos hub o nodos concentradores. Los concentradores son nodos con un valor de centralidad de vector propio por encima 
#del cuantil empírico del 95 % de todas las centralidades de vector propio en la red  
#weightDeg y normDegestán estan configurados para FALSE que el grado de un nodo se defina simplemente como el número de nodos que son adyacentes al nodo.


props_single <- netAnalyze(net_single, 
                           centrLCC = TRUE,
                           clustMethod = "cluster_fast_greedy",
                           hubPar = "eigenvector",
                           weightDeg = FALSE, normDeg = FALSE)

#?summary.microNetProps
summary(props_single, numbNodes = 5L)

#>los datos de salida son los siguientes
## 
## Component sizes
## ```````````````          
## size: 48 1
##    #:  1 2
## ______________________________
## Global network properties
## `````````````````````````
## Largest connected component (LCC):
##                                  
## Relative LCC size         0.96000
## Clustering coefficient    0.33595
## Modularity                0.53407
## Positive edge percentage 88.34951
## Edge density              0.09131
## Natural connectivity      0.02856
## Vertex connectivity       1.00000
## Edge connectivity         1.00000
## Average dissimilarity*    0.97035
## Average path length**     2.36877
## 
## Whole network:
##                                  
## Number of components      3.00000
## Clustering coefficient    0.33595
## Modularity                0.53407
## Positive edge percentage 88.34951
## Edge density              0.08408
## Natural connectivity      0.02714
## 
##  *Dissimilarity = 1 - edge weight
## **Path length: Units with average dissimilarity
## 
## ______________________________
## Clusters
## - In the whole network
## - Algorithm: cluster_fast_greedy
## ```````````````````````````````` 
##                     
## name: 0  1  2  3 4 5
##    #: 2 12 17 10 5 4
## 
## ______________________________
## Hubs
## - In alphabetical/numerical order
## - Based on empirical quantiles of centralities
## ```````````````````````````````````````````````       
##  190597
##  288134
##  311477
## 
## ______________________________
## Centrality measures
## - In decreasing order
## - Centrality of disconnected components is zero
## ````````````````````````````````````````````````
## Degree (unnormalized):
##          
## 288134 10
## 190597  9
## 311477  9
## 188236  8
## 199487  8
## 
## Betweenness centrality (normalized):
##               
## 302160 0.31360
## 268332 0.24144
## 259569 0.23404
## 470973 0.21462
## 119010 0.19704
## 
## Closeness centrality (normalized):
##               
## 288134 0.68431
## 311477 0.68417
## 199487 0.68108
## 302160 0.67528
## 188236 0.66867
## 
## Eigenvector centrality (normalized):
##               
## 288134 1.00000
## 311477 0.94406
## 190597 0.90806
## 199487 0.85436
## 188236 0.72730


#Visualización de la red 

#lo que se pretende realizar en este paso visualizar nuestra red de manera visual, utilizando nodos y aristas
#Usamos los grupos determinados como colores de nodo y 
#escalamos los tamaños de nodo de acuerdo con la centralidad del vector propio del nodo.

# help page
?plot.microNetProps

sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

p <- plot(props_single, 
          nodeColor = "cluster", 
          nodeSize = "eigenvector",
          title1 = "Network on OTU level with SPRING associations", 
          showTitle = TRUE,
          cexTitle = 2.3)

legend(0.7, 1.1, cex = 2.2, title = "estimated association:",
       legend = c("+","-"), lty = 1, lwd = 3, col = c("#009900","red"), 
       bty = "n", horiz = TRUE)

p$q1$Arguments$cut

#75% 
#0.3367287 
#

#Red única con correlación de Pearson como medida de asociación
net_single2 <- netConstruct(amgut2.filt.phy,  
                            measure = "pearson",
                            normMethod = "clr", 
                            zeroMethod = "multRepl",
                            sparsMethod = "threshold", 
                            thresh = 0.3,
                            verbose = 3)
######.   salida: 

## 2 rows with zero sum removed.

## 138 taxa and 294 samples remaining.

## 
## Zero treatment:

## Execute multRepl() ...

## Warning in (function (X, label = NULL, dl = NULL, frac = 0.65, imp.missing = FALSE, : Row(s) containing more than 80% zeros/unobserved values were found (check it out using zPatterns).
##                   (You can use the z.warning argument to modify the warning threshold).

## Done.
## 
## Normalization:
## Execute clr(){SpiecEasi} ... Done.
## 
## Calculate 'pearson' associations ... Done.
## 
## Sparsify associations via 'threshold' ... Done.


#Análisis y visualización de redes
props_single2 <- netAnalyze(net_single2, clustMethod = "cluster_fast_greedy")

plot(props_single2, 
     nodeColor = "cluster", 
     nodeSize = "eigenvector",
     title1 = "Network on OTU level with Pearson correlations", 
     showTitle = TRUE,
     cexTitle = 2.3)

legend(0.7, 1.1, cex = 2.2, title = "estimated correlation:", 
       legend = c("+","-"), lty = 1, lwd = 3, col = c("#009900","red"), 
       bty = "n", horiz = TRUE)



#mejorando la visualización de la red para que se vea mas estetica y la infoamci`ón sea mas facil de procesar_
plot(props_single2, 
     nodeColor = "cluster", 
     nodeSize = "eigenvector",
     repulsion = 0.8,
     rmSingles = TRUE,
     labelScale = FALSE,
     cexLabels = 1.6,
     nodeSizeSpread = 3,
     cexNodes = 2,
     title1 = "Network on OTU level with Pearson correlations", 
     showTitle = TRUE,
     cexTitle = 2.3)

legend(0.7, 1.1, cex = 2.2, title = "estimated correlation:",
       legend = c("+","-"), lty = 1, lwd = 3, col = c("#009900","red"),
       bty = "n", horiz = TRUE)










#####---------------comandos para visualizar otro ejemplo de redes 
data("soilrep")

soil_warm_yes <- phyloseq::subset_samples(soilrep, warmed == "yes")
soil_warm_no  <- phyloseq::subset_samples(soilrep, warmed == "no")

net_seas_p <- netConstruct(soil_warm_yes, soil_warm_no,
                           filtTax = "highestVar",
                           filtTaxPar = list(highestVar = 500),
                           zeroMethod = "pseudo",
                           normMethod = "clr",
                           measure = "pearson",
                           verbose = 0)

netprops1 <- netAnalyze(net_seas_p, clustMethod = "cluster_fast_greedy")

nclust <- as.numeric(max(names(table(netprops1$clustering$clust1))))

col <- c(topo.colors(nclust), rainbow(6))

plot(netprops1, 
     sameLayout = TRUE, 
     layoutGroup = "union", 
     colorVec = col,
     borderCol = "gray40", 
     nodeSize = "degree", 
     cexNodes = 0.9, 
     nodeSizeSpread = 3, 
     edgeTranspLow = 80, 
     edgeTranspHigh = 50,
     groupNames = c("Warming", "Non-warming"), 
     showTitle = TRUE, 
     cexTitle = 2.8,
     mar = c(1,1,3,1), 
     repulsion = 0.9, 
     labels = FALSE, 
     rmSingles = "inboth",
     nodeFilter = "clustMin", 
     nodeFilterPar = 10, 
     nodeTransp = 50, 
     hubTransp = 30)



########construcción de redes 