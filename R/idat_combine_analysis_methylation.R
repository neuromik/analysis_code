library(ChAMP)
library(sva)
library(stringr)

myLoad_112179<-champ.load(directory='/Volumes/storage_ext/methylation/_BD_EPIC_idat/GSE112179_idat/GSE112179_RAW',arraytype="EPIC")
myLoad_129428<-champ.load(directory='/Volumes/storage_ext/methylation/_BD_EPIC_idat/GSE129428_idat/GSE129428_RAW',arraytype="EPIC")
myLoad_191200<-champ.load(directory='/Volumes/storage_ext/methylation/_BD_EPIC_idat/GSE191200_idat/GSE191200_RAW',arraytype="EPIC")

myNorm_112179<-champ.norm(myLoad_112179$beta,method="BMIQ",arraytype='EPIC',cores=6)
myNorm_129428<-champ.norm(myLoad_129428$beta,method="BMIQ",arraytype='EPIC',cores=6)
myNorm_191200<-champ.norm(myLoad_191200$beta,method="BMIQ",arraytype='EPIC',cores=6)

myCombat_112179<-champ.runCombat(beta=myNorm_112179,pd=myLoad_112179$pd,variablename="Sample_Group",batchname=c("Slide","Array"))
myCombat_129428<-champ.runCombat(beta=myNorm_129428,pd=myLoad_129428$pd,variablename="Sample_Group",batchname=c("Slide","Array"))
myCombat_191200<-champ.runCombat(beta=myNorm_191200,pd=myLoad_191200$pd,variablename="Sample_Group",batchname=c("Slide","Array")) #does not work

#Find unknown confounders 
library(sva) 
mod = model.matrix(~Sample_Group + Sex + Age, data=myLoad_112179$pd)
mod0 = mod[,-c(2)] 
n.sv<-num.sv(as.matrix(myCombat_112179), mod, method = "be", vfilter = NULL, B = 20,seed = NULL)
n=n.sv
sva = svaseq(myCombat_112179,mod,mod0, n.sv =n)$sv 
colnames(sva)<-paste("sv", 1:n, sep = "")
pd_sva<-cbind(myLoad_112179$pd,sva)

#Correct known and unknown confounders
design = model.matrix(~Sample_Group +Sex + Age+sv1+sv2+sv3+sv4+sv5+sv6+sv7+sv8+sv9+sv10+sv11+sv12+sv13, data=pd_sva)
Y<-myCombat_112179
X<-as.matrix(design)
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
myRegress_112179= Y - t(X[,c(3:ncol(X))] %*% beta[c(3:nrow(beta)),])

#repeat for myCombat_129428 and myCombat_191200

#join data and adjust names
a<-as.data.frame(myRegress_112179)
b<-as.data.frame(myRegress_129428)
c<-as.data.frame(myRegress_191200)

df=transform(merge(a,b,by=0), row.names=Row.names, Row.names=NULL)
BD_EPIC=transform(merge(df,c,by=0), row.names=Row.names, Row.names=NULL)
colnames(BD_EPIC) = gsub(".x", "", colnames(BD_EPIC)) 
colnames(BD_EPIC) = gsub(".y", "", colnames(BD_EPIC)   
new_column_names= paste('V', as.character(1:176))
new_column_names= str_replace(new_column_names, ' ', '')
colnames(BD_EPIC) =  new_column_names

#correct batch effect from multiple studies                    
pd_BD<-read.csv('/Volumes/storage_ext/methylation/_BD_EPIC_idat/BD_EPIC_meta.csv')
myCombat<-champ.runCombat(beta=as.matrix(BD_EPIC),pd=pd_BD,variablename="Sample_Group",batchname=c("Series_ID"))

