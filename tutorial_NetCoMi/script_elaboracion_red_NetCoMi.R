#script del tutorial "NetCoMi"

#NetCoMi: es un paquete de R que integra métodos existentes para la elaboración paso a paso, 
#para construir y analizar redes de asociación microbianas individuales (partiendo de datos de secuenciación metagenómica), 
#así como para cuantificar las diferencias de la red. 

#Este es un script modificado y adaptado del tutorial publicado por  Peschel el al., 2021."NetCoMi: network construction and comparison for microbiome data in R"

##################################### Construcción de redes y comparación de datos de microbiomas en R ######################################

#El fundamento de esta paquetería es realizar redes de asociación microbianas,
#tomando en cuenta la abundancia del taxón OTU (NODO)
#El pathline de manera global es el sig.:
#1. Tomar en cuenta la abundancia del taxón OTU. Por correlación se empieza a construir 
#la red: tomamos en cuenta algun tipo de correlación a utilizar (Pearson, Bicor, Sperman).
#2. Construimos la matriz de adyacencia (usando una medida de asociación del punto 1).
#fundamento: Se calcula una matriz de asociación con entradas que expresan la relación entre pares de taxones. 
#3. Análisis de la red: este paso es importante ya que permite caracterizar la topología de la red: Camino más corto,
#Medidas de centralidad: Grado, betweenness, closeness, eigenvector centrality, hubs.
#Fundamento: Las metrícas de la red nos ayudaran a conocer como se comportan las comunidades microbianas en las condiciones 
#que queremos evaluar o el tipo de muestra que se tiene,  con ello nos da como resultado aquellos taxones clave que juegan 
#un papel importante dentro de la comunidad.

########################################### script tutorial ###############################################################
install.packages("devtools")
install.packages("BiocManager")
install.packages("WGCNA")
install.packages("qgraph")
install.packages("remotes")
remotes::install_github("stefpeschel/NetCoMi")
#Tiempo de ejecución max. 3 min


# Para instalar y descaragar la matrices del tutoril, de aplica el sig. comando:
devtools::install_github("stefpeschel/NetCoMi", 
                         dependencies = c("Depends", "Imports", "LinkingTo"),
                         repos = c("https://cloud.r-project.org/",
                                   BiocManager::repositories()))
#Tiempo de ejecución máximo 10 min, Nota: la versión de R puede influir en que no se instalen todas las librerias 
#que trae el tutorial, por lo cual, se recomienda que algunas librerias se instalen manualmente cuando marque
#algun tipo de error

#Algunas paqueterias se deben de instalar con los siguientes comandos
devtools::install_github("GraceYoon/SPRING", force = TRUE)
devtools::install_github("zdk123/SpiecEasi", force = TRUE)


#Llamamos a las librerias
library(WGCNA)
library(NetCoMi)


#Cargamos los datos del tutorial: en este caso los datos que se cargan son las tablas de las muestras con los datos 
#de abundancia de los OTUS

data("amgut1.filt")
data("amgut2.filt.phy")


#Algunas consideraciones de los comandos a tomar en cuenta :
#Especificaciones de los comandos:
#>netConstruct().   sirve para filtrar los datos de la red
#>filtSamp.        Solo se incluyen muestras con un número total de lecturas ejemplo: de al menos 1000
#>filtTax.        El comando manda a llamar solo al #NUMERO DE TAXONES CON MAYOR FRECUENCIA EJEMPLO: MAYOR A 50
#>measure.         define la medida de asociación o disimilitud, que es "spring" en nuestro caso.
#>measurePar. nlambday rep.num se establecen en 10 para reducir el tiempo de ejecución,
# pero deberan ser más altos para datos reales.
#>La normalización, así como el manejo de cero, se realizan internamente en SPRING(). 
#Por lo tanto, establecemos normMethody zeroMethodpara "none".
#>Además, configuramos porque sparsMethod = "none"
#devuelve una red dispersa en la que no es necesario ningún paso de dispersión adicional.
#>dissFunc.     Este comando es para seleccionar si la red va ser firmada, no firmada o hibrida (signed, unsigned, hybrid)
#en el ejemplo se realiza dissFunc = "signed" ------firmado
#El argumento: "verbose" se establece en 3 para que netConstruct() se impriman todos los mensajes generados por funciones externas verbose = 3

####################Construcción de la red##################################
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
#Tiempo de ejecución aprox. 1 min 
#Nota: con datos externos al del tutorial el tiempo de ejecución puede variar



########################################Analizando la red construida##########################################
#Una ves que se tiene la red, se analiza la red de acuerdo a los valores de las metricas de la red.
#La función a utilizar para determinar las metrícas es la sig,:
#netAnalyze().   sirve para analizar la red construida

#función de los comandos:
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


#Una ves que determinamos las metrícas de la red, visualizamos los resultados con el siguiente comando:
#numbNodes    es la función para que te muestre la información de los primeros 5 nodos 
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


################################################Visualización de la red############################################## 

#lo que se pretende realizar en este paso es visualizar nuestra red de manera visual, utilizando nodos y aristas
#Usamos los grupos determinados como colores de nodo y 
#escalamos los tamaños de nodo de acuerdo con la centralidad del vector propio del nodo.

png("Network on OTU level with SPRINHG associations.png")
p <- plot(props_single, 
          nodeColor = "cluster", 
          nodeSize = "eigenvector",
          title1 = "Network on OTU level with SPRING associations", 
          showTitle = TRUE,
          cexTitle = 1.3)

legend(0.7, 1, cex = 0.5, title = "estimated association:",
       legend = c("+","-"), lty = 1, lwd = 3, col = c("#009900","red"), 
       bty = "n", horiz = TRUE)
dev.off()


p$q1$Arguments$cut


#Resultado de salida
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
######Resultado de salida: 

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


###################################Análisis y visualización de redes###########################################
props_single2 <- netAnalyze(net_single2, clustMethod = "cluster_fast_greedy")

png("Network on OTU level with Pearson correlations.png")
plot(props_single2, 
     nodeColor = "cluster", 
     nodeSize = "eigenvector",
     title1 = "Network on OTU level with Pearson correlations", 
     showTitle = TRUE,
     cexTitle = 1.3)

legend(0.7, 1, cex = 0.5, title = "estimated correlation:", 
       legend = c("+","-"), lty = 1, lwd = 3, col = c("#009900","red"), 
       bty = "n", horiz = TRUE)
dev.off()


#mejorando la visualización de la red para que se vea mas estetica y la infoamción sea mas facil de procesar
png("Network on OTU level with Pearson correlations_mejorado.png")
plot(props_single2, 
     nodeColor = "cluster", 
     nodeSize = "eigenvector",
     repulsion = 0.8,
     rmSingles = TRUE,
     labelScale = FALSE,
     cexLabels = 1,
     nodeSizeSpread = 3,
     cexNodes = 2,
     title1 = "Network on OTU level with Pearson correlations", 
     showTitle = TRUE,
     cexTitle = 1.3)

legend(0.7, 1, cex = 0.6, title = "estimated correlation:",
       legend = c("+","-"), lty = 1, lwd = 3, col = c("#009900","red"),
       bty = "n", horiz = TRUE)
dev.off()