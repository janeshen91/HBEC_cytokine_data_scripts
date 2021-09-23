c5<-read.csv(file="c5.go.bp.v7.4.symbols.gmt",sep="\t",header=F)
c5t<-t(c5)
c5tdf<-as.data.frame(c5t)
namesC5<-c5tdf[1,]
c5tdf<-c5tdf[-c(1,2),]
colnames(c5tdf)<-namesC5

combined<-c5tdf[,1]
for(i in c(2:dim(c5tdf)[1])){
	combined<-c(combined,c5tdf[,i])

}
combinedUnique<-unique(combined)
##how to figure out what is the 



combinedUnique1<- combinedUnique %in% c5tdf[,1] 
for(i in c(6822:dim(c5tdf)[2])){
combinedUnique2<- combinedUnique %in% c5tdf[,i] 
combinedUnique1<-cbind(combinedUnique1,combinedUnique2)
}


sumsOfColumns<-colSums(combinedUnique1)

rownames(combinedUnique1)<-combinedUnique
colnames(combinedUnique1)<-namesC5
combinedUnique1<-as.data.frame(combinedUnique1)

save(combinedUnique1,file="BP_C5.RData")

grepGO<-grep("GO",namesC5,value=FALSE,invert=TRUE)
combinedUnique2<-combinedUnique1[,-grepGO]
for(i in grepGO){
print(i)
combinedUnique2[,i]<-combinedUnique1[,i-1]+combinedUnique1[,i]
}
names2<-colnames(combinedUnique2)
grepGO<-grep("GO",names2,value=FALSE,invert=TRUE)
save(combinedUnique2,file="BP_C5ii.RData")
###on server

library(iDEA)

load("BP_C5.RData")

#idea <- CreateiDEAObject(summary_data, annotation_data, max_var_beta = 100, min_precent_annot = 0.0025, num_core=10)

files<-list.files("./input")
for(i in c(1:9)){
filename<-paste("./input/",files[i],sep="")
file=read.csv(file=filename,sep="\t")
idea <- CreateiDEAObject(file, combinedUnique1, max_var_beta = 100, min_precent_annot = 0.0025, num_core=30)
idea <- iDEA.fit(idea,
                 fit_noGS=FALSE,
                 init_beta=NULL, 
                 init_tau=c(-2,0.5),
                 min_degene=5,
                 em_iter=15,
                 mcmc_iter=1000, 
                 fit.tol=1e-5,
                 modelVariant = F,
                 verbose=TRUE)
idea <- iDEA.louis(idea) ## 
filenameOutput<-paste("./outputBP/",files[i],sep="")
write.table(idea@gsea,file=filenameOutput)
}