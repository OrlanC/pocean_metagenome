
#date: 27/06/2022
#name : Mario Iván Alemán Duarte
#objective: TO construct a co-ocrrence network using a .biom extension file
############################################################################


#set working directory 

#workingDir = "WRITE YOUR OWN PATH HERE";
workingDir = "/home/mario_aleman/Escritorio/OCEAN";

setwd(workingDir); 

#############################################################
#Package installation                                       #
#############################################################

#instal "BiocManager" if is not installed before
install.packages("BiocManager", dependencies = TRUE)

#install WGNA, the package to construct the network 
BiocManager::install("WGCNA", dependencies = TRUE) 

############################################################
#Load packages
###########################################################

# Load the WGCNA package
library(WGCNA);


############################################################
#seetings 
##############################################################

# The following setting is important, do not omit.

options(stringsAsFactors = FALSE);
