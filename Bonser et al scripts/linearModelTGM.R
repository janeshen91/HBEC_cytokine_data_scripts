##trying to figure out what these are: gene expression of something?

library(ggplot2)
library(cowplot)
library(ggpubr)

##read in  MAST data
MASTMeta<-read.csv(file="./MASTAsthmaAndHealth_MetaData.csv")
MASTlog2<-read.csv(file="./MASTAsthmaAndHealth_TMMVoomNormalized_Log2Expression.csv")

MASTraw<-read.csv(file=".//MAST_allRawCounts.csv")
sum(colnames(MASTraw) %in% MASTMeta$MAST_ID)
sum(colnames(MASTlog2) %in% MASTMeta$MAST_ID)

##read in spira data
SpiraMeta<-read.csv(file="./SpiraMetadataWithGeneSignatures.csv")
SpiraNormalized<-read.csv(file="./Spira_RMANormalizedData.csv")
SpiraSmoker<-read.csv(file="./SpiraFormersmokersonly_RMANormalizedData.csv")
sum(colnames(SpiraNormalized) %in% SpiraMeta$X)
sum(colnames(SpiraSmoker) %in% SpiraMeta$X)
SpiraMetaSmoker<-subset(SpiraMeta,SpiraMeta$characteristics..SmokingStatus=="EX")
sum(colnames(SpiraSmoker) %in% SpiraMetaSmoker$X)



##reading in genes that are upregulated by IL13
UpGenesIL13vsUT<-read.csv(file="../LukeWalterData/Up_gene_list.condition-IL13vsUT.txt",header=F)

rows<-MASTlog2$X %in% UpGenesIL13vsUT$V1
onnlyUPGenes<-MASTlog2[rows,]
onnlyUPGenesT<-t(onnlyUPGenes)
colnames(onnlyUPGenesT)<-onnlyUPGenesT[1,]

onnlyUPGenesT<-onnlyUPGenesT[-1,]
onnlyUPGenesT<-as.data.frame(onnlyUPGenesT)
onnlyUPGenesT$sampleNames<-rownames((onnlyUPGenesT))

#merge the upregulated genes with teh MAST data
UpGenesMetaMast<-merge(onnlyUPGenesT,MASTMeta[,c("MAST_ID","Disease_Status","threegenemean")],by.x="sampleNames",by.y="MAST_ID")


##subset only the asthma individuals
subsetDataAsthma<-subset(UpGenesMetaMast,UpGenesMetaMast$Disease_Status=="asthma")


save(UpGenesMetaMast,file="UpGenesMetaMast.RData")

##getting ready for doing linear regression, by setting up a coefficient dataframe
coefficientsOfUpMast<-UpGenesMetaMast
coefficientsOfUpMast<-t(UpGenesMetaMast)[,c(1:6)]
coefficientsOfUpMast<-as.data.frame(coefficientsOfUpMast)

##doing correlation and linear model of TGM and gene expression and age and sex
for(i in c(2:(numGenes+1))){
  print(i)
  corResult<-cor.test(as.numeric(UpGenesMetaMast[,i]), as.numeric(UpGenesMetaMast$threegenemean), method=c("pearson"))
  coefficientsOfUpMast[i,1]<-corResult$estimate
  coefficientsOfUpMast[i,2]<-corResult$p.value
  lmiAll<-lm(as.numeric(UpGenesMetaMast[,i]) ~  threegenemean+age+sex, data=UpGenesMetaMast)
  sumCoefsAll<-summary(lmiAll)
  coefficientsOfUpMast[i,3]<-sumCoefsAll$coefficients[2,1]
  coefficientsOfUpMast[i,4]<-sumCoefsAll$coefficients[2,4]
  subsetDataAsthma<-subset(UpGenesMetaMast,UpGenesMetaMast$Disease_Status=="asthma")
  lmi<-lm(as.numeric(subsetDataAsthma[,i]) ~  subsetDataAsthma$threegenemean+age+sex, data=subsetDataAsthma)
  sumCoefs<-summary(lmi)
  coefficientsOfUpMast[i,5]<-sumCoefs$coefficients[2,1]
  coefficientsOfUpMast[i,6]<-sumCoefs$coefficients[2,4]
}
colnames(coefficientsOfUpMast)<-c("Pearson's R", "P value","Coefficient of Linear Model All Samples","P value","Coefficient of Linear Model Asthma Samples","Pvalue")
fileNameOfTable<-paste("correlationPvalues",DEfiles[j],".csv",sep="")
numRows<-dim(coefficientsOfUpMast)[1]
numRows3<-(numRows-3)
coefficientsOfUpMast$geneNames<-rownames(coefficientsOfUpMast)
write.table(coefficientsOfUpMast[-c(1,numRows3:numRows),],file=fileNameOfTable,sep=",",quote=F,row.names=FALSE)

}

##Reading in Spira data and doing a similar linear model
head(SpiraMeta)
head(SpiraNormalized)
UpGenesIL13vsUT<-read.csv(file="/Users/janeshen/Documents/StephData//LukeWalterData/Up_gene_list.condition-IL13vsUT.txt",header=F)
#see if the upgenes are in the mast??
rows<-SpiraNormalized$X %in% UpGenesIL13vsUT$V1
onnlyUPGenes<-SpiraNormalized[rows,]
onnlyUPGenesT<-t(onnlyUPGenes)
colnames(onnlyUPGenesT)<-onnlyUPGenesT[1,]

onnlyUPGenesT<-onnlyUPGenesT[-1,]
onnlyUPGenesT<-as.data.frame(onnlyUPGenesT)
onnlyUPGenesT$sampleNames<-rownames((onnlyUPGenesT))


UpGenesMetaSpira<-merge(onnlyUPGenesT,SpiraMeta[,c("X","characteristics..HistoryOfAsthma","threegenemean","characteristics..Age","gender")],by.x="sampleNames",by.y="X")
coefficientsOfUpSpira<-UpGenesMetaSpira
coefficientsOfUpSpira<-t(UpGenesMetaSpira)[,c(1:6)]
coefficientsOfUpSpira<-as.data.frame(coefficientsOfUpSpira)
for(i in c(2:27)){
corResult<-cor.test(as.numeric(UpGenesMetaSpira[,i]), as.numeric(UpGenesMetaSpira$threegenemean), method=c("pearson"))
coefficientsOfUpSpira[i,1]<-corResult$estimate
coefficientsOfUpSpira[i,2]<-corResult$p.value
lmiAll<-lm(as.numeric(UpGenesMetaSpira[,i]) ~  threegenemean+gender+characteristics..Age, data=UpGenesMetaSpira)
sumCoefsAll<-summary(lmiAll)
coefficientsOfUpSpira[i,3]<-sumCoefsAll$coefficients[2,1]
coefficientsOfUpSpira[i,4]<-sumCoefsAll$coefficients[2,4]
subsetDataAsthma<-subset(UpGenesMetaSpira,UpGenesMetaSpira$characteristics..HistoryOfAsthma=="yes")
lmi<-lm(as.numeric(subsetDataAsthma[,i]) ~  threegenemean+gender+characteristics..Age, data=subsetDataAsthma)
sumCoefs<-summary(lmi)
coefficientsOfUpSpira[i,5]<-sumCoefs$coefficients[2,1]
coefficientsOfUpSpira[i,6]<-sumCoefs$coefficients[2,4]

}
colnames(coefficientsOfUpSpira)<-c("Pearson's R", "P value","Coefficient of Linear Model All Samples","P value","Coefficient of Linear Model Asthma Samples","Pvalue")

coefficientsOfUpSpira$gene<-rownames(coefficientsOfUpSpira)
write.table(coefficientsOfUpSpira[-1,],file="../SpiraCorrelations/coefficientsOfUpSpiraIL13.txt",row.names=F,quote=F)

