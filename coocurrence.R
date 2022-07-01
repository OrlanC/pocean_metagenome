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

#####################################
######## base code occurrence #######
#####################################
#### Checkerboard
##concept script1
#ordered data
x1<-c(1,1,1,0,0,0,0,0,1,1)
x2<-c(1,1,1,0,1,1,1,1,0,0)
data<-data.frame(Isla1=x1, Isla2=x2)
data<-t(data)
#select no-ocurrent to Cij
co<-data[1,]!=data[2,] #compare
data[, co, drop=FALSE] #drop False and select True
drop<-data[,co, drop=FALSE]
##concept script2
#non-ordered data
x3<-c(0,0,0,1,1,1,0,1,0,1)
x4<-c(1,1,1,0,1,0,0,0,1,0)
x3.1<-c(5,2,2,10,25,80,75,50,30,49)
x4.1<-c(5,2,0,100,2,35,75,0,30,50)
data2<-data.frame(Isla1=x3, Isla2=x4)
data2<-t(data2)
co1<-data2[1,]!=data2[2,]
co1
drop1<-data2[,co1, drop=FALSE]
####
#estimate Cij
cij<-(length(drop[1,])-sum(drop[1,]))*(length(drop[2,])-sum(drop[2,]))
cij

cij1<-(length(drop1[1,])-sum(drop1[1,]))*(length(drop1[2,])-sum(drop1[2,]))
cij1

C<-sum(cij,cij1)/(4*(4-1)/2)
C

g<-cor.test(x3.1,x4.1, method="spearman")

dacars<-mtcars
cortest<-cor.test(dacars$wt, dacars$mpg, method = "pearson")
cortest

###############################################################################
for (i in colnames(mydata)) {
  print(paste("adding column:", i))
  abund<-mydata[,i]
  #print(abund)
  y=c()
  #dframe[rownames(mydata),i]<-y
  for (z in abund){
    if (z >= 5) {
      y[z]=1
    }else{
      y[z]=0
    }
  }
}
########### code of presencen and abscen of species#############################

x1<-c(4,3,8,11,4,3,1,20)
x2<-c(1,10,9,11,0,1,7,34)
x3<-c(1,1,1,1,7,1,10,6)
dataframe<-data.frame(a=x1, b=x2, c=x3)
df<-dataframe

df$first<-rep(1, each=nrow(df))
df<- df[ , c("first",names(df)[names(df) != "first"])]

df$last<-rep(1, each=nrow(df))

y=as.single(c())
df_n<-data.frame(row.names = rownames(df))
for (i in colnames(df)){
  df_n[,i]<-y
  for(z in rownames(df_n)){
    if(df[z,i]>=5){
      y[z]=1
    }else{
      y[z]=0
    }
  }
}

df_n<-df_n[,-c(1:2)]
colnames(df_n)<-colnames(dataframe)

### Permutation ###
#install.packages("combinat")
library("combinat")
comn<-combn(1:6, 2)

comn<-as.data.frame(comn)
head(comn)

per<-comn[,1]

for (i in colnames(comn)){
  print(comn[,i])
  perm<-comn[,i]
  df_n[c(min(rownames(comn)),max(rownames(comn))),i]
  
  
}
######################################################
x1<-c(4,3,8,11,4,3,1,20)
x2<-c(1,10,9,11,0,1,7,34)
x3<-c(1,1,1,1,7,1,10,6)
dataframe<-data.frame(a=x1, b=x2, c=x3)
df<-dataframe

y=as.single(c())
for (z in rownames(df)){
  print(z)
  if (df[z,1]>=5){
    y[z]=1
  }else{
    y[z]=0
  }
}
y=as.single(c())

#df_n[,i]<-y
y=as.single(c())
#list_1[[i]]<-y

#list_1<-list()
y=as.single(c())

df$first<-rep(1, each=nrow(df))
df<- df[ , c("first",    # Reorder data frame
             names(df)[names(df) != "first"])]
df$last<-rep(1, each=nrow(df))

y=as.single(c())
df_n<-data.frame(row.names = rownames(df))
for (i in colnames(df)){
  df_n[,i]<-y
  for(z in rownames(df_n)){
    if(df[z,i]>=5){
      y[z]=1
    }else{
      y[z]=0
    }
  }
}

df_n<-df_n[,-c(1:2)]
colnames(df_n)<-colnames(dataframe)

########################################################


#y=c()
#p[y]=y
if (z >= 5) {
  y[z]=1
}else{
  y[z]=0
}
}

y=c()
for (z in (x4)){
  print(z)
  #y=c()
  #p[y]=y
  if (z >= 5) {
    y[z]=1
  }else{
    y[z]=0
  }
}




for (x in abund){
  if x >=5{
    y[x]=1
  }else{
    y[x]=0
    
  }
  
}

df=list_data[[i]]
dframe[rownames(df), i]<-df[,1]
}


