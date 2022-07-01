#####
setwd("C:/Users/Orlando Camargo/Desktop/")
#BiocManager::install("rhdf5")
#BiocManager::install("phyloseq")
library("phyloseq")
library("ggplot2")
library("readr")
library("patchwork")

merged_metagenome_op<- import_biom("pacific_ocean.biom")
mydata<-merged_metagenome_op@otu_table
mydata<-as.data.frame(mydata)
head(mydata)
data<-mydata
rownames(data)<-1:nrow(data)
head(data)

########### code of presencen and abscen of species#############################

data$first<-rep(1, each=nrow(data))
data<- data[ , c("first",names(data)[names(data) != "first"])]

data$last<-rep(1, each=nrow(data))

y=as.single(c())
df_n<-data.frame(row.names = rownames(data))
for (i in colnames(data)){
  print(paste("running:",i,"..."))
  df_n[,i]<-y
  for(z in rownames(df_n)){
    if(data[z,i]>=5){ #este número es arbitrario podemos usar porcentajes de abudancia
      y[z]=1
    }else{
      y[z]=0
    }
  }
}

df_n<-df_n[,-c(1:2)]
colnames(df_n)<-colnames(mydata)

head(df_n)

##############

library(combinat)
combn()
