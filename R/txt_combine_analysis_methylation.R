library(ChAMP)
library(sva)
library(stringr)
library(FDb.InfiniumMethylation.hg19)
library(RnBeads)
library(RnBeads.hg19)

#load .idat files
myLoad_98203<-champ.load(directory='/Volumes/storage_ext/methylation/_MDD_450K/GSE98203_idat/GSE98203_RAW',arraytype="450K")
myNorm_98203<-champ.norm(myLoad_98203$beta,method="BMIQ",arraytype='EPIC',cores=6)

#load series matrix files
myLoad_41826 =  rnb.read.geo(
  accession = '/Volumes/storage_ext/methylation/_MDD_450K/GSE41826/GSE41826_RAW/GSE41826_series_matrix.txt.gz',
  verbose = logger.isinitialized(),
  destdir = '~/Volumes/storage_ext/methylation/_MDD_450K/GSE41826')
pd_41826<-read.csv('/Volumes/storage_ext/methylation/_MDD_450K/GSE41826/GSE41826_RAW/GSE41826_meta.csv')
myNorm_41826<-champ.norm(beta=myLoad_41826@meth.sites[],method="BMIQ",arraytype='450K',cores=6)

myLoad_88890 =  rnb.read.geo(
  accession = '/Volumes/storage_ext/methylation/_MDD_450K/GSE88890/GSE88890_RAW/GSE88890_series_matrix.txt.gz',
  verbose = logger.isinitialized(),
  destdir = '~/Volumes/storage_ext/methylation/_MDD_450K/GSE88890')
pd_88890<-read.csv('/Volumes/storage_ext/methylation/_MDD_450K/GSE88890/GSE88890_RAW/GSE88890_meta.csv')
myNorm_88890<-champ.norm(beta=myLoad_88890@meth.sites[],method="BMIQ",arraytype='450K',cores=6)

#adjust batch
myCombat_98203<-champ.runCombat(beta=myNorm_98203,pd=myLoad_98203$pd,variablename="Sample_Group",batchname=c("Slide","Array"))
myCombat_41826<-champ.runCombat(beta=myNorm_41826,pd=pd_41826,variablename="Sample_Group",batchname=c("Sentrix_ID","Sentrix_Position"))
myCombat_88890<-champ.runCombat(beta=myNorm_88890,pd=pd_88890,variablename="Sample_Group",batchname=c("Sentrix_ID","Sentrix_Position"))

#adjust covariates with regression

#Find unknown confounders 
library(sva) 
mod = model.matrix(~Sample_Group + Sex + Age, data=myLoad_98203$pd)
mod0 = mod[,-c(2)] 
n.sv<-num.sv(as.matrix(myCombat_98203), mod, method = "be", vfilter = NULL, B = 20,seed = NULL)
n=n.sv
sva = svaseq(myCombat_98203,mod,mod0, n.sv =n)$sv 
colnames(sva)<-paste("sv", 1:n, sep = "")
pd_sva<-cbind(myLoad_98203$pd,sva)

#Correct known and unknown confounders
design = model.matrix(~Sample_Group +Sex + Age+sv1+sv2+sv3+sv4+sv5+sv6+sv7+sv8+sv9+sv10+sv11+sv12+sv13, data=pd_sva)
Y<-myCombat_98203
X<-as.matrix(design)
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
myRegress_98203= Y - t(X[,c(3:ncol(X))] %*% beta[c(3:nrow(beta)),])

#repeat for myCombat_41826 and myCombat_88890

#join data and adjust samples names
a<-as.data.frame(myRegress_98203)
b<-as.data.frame(myRegress_41826)
c<-as.data.frame(myRegress_88890)

df=transform(merge(a,b,by=0), row.names=Row.names, Row.names=NULL)
MDD_450K=transform(merge(df,c,by=0), row.names=Row.names, Row.names=NULL)

colnames(MDD_450K) = gsub(".x", "", colnames(MDD_450K)) 
colnames(BD_EPIC) = gsub(".y", "", colnames(MDD_450K)   
                         
new_column_names= paste('V', as.character(1:n_samples))
new_column_names= str_replace(new_column_names, ' ', '')
colnames(MDD_450K) =  new_column_names
                         
#adjust batch from different studies

pd_MDD<-read.csv('/Volumes/storage_ext/methylation/_MDD_450/MDD_450_combines_meta.csv')
myCombat<-champ.runCombat(beta=as.matrix(MDD_450),pd=pd_MDD,variablename="Sample_Group",batchname=c("Series_ID"))

