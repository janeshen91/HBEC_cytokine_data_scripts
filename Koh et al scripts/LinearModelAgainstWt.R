#analyzing KD1
#again with just coutns over wildtype?
#wildtype is ID770
setwd("~/Documents/MPRAflow/OYB2 Feb 2021/KD1")
KD1<-read.csv(file="KD1 individual CRS.txt",sep="\t")
head(KD1)
for(i in c(1:25)){
  name<-paste("pos",i*8-7,"-",i*8,sep="")

  for(j in c(1:24)){
    k<-(i-1)*24+j
    KD1[k,4]<-name
  }
}
#WHATEVER, i give up on naming them
#IL13
  IL13avg<-read.csv(file="./KD1May25/IL13/allreps.tsv",sep="\t")
head(IL13avg)
IL13avgii<-merge(IL13avg,KD1,by.x="name",by.y="")
IL13avgiii<-substr(IL13avg$name,start=3,stop=7)
IL13avg$Pos<-IL13avgii
IL13avgiii<-merge(IL13avg,KD1,by.x="Pos",by.y="Number")


ggplot(IL13avg,aes(x=Pos,y=log2))+geom_point()
#i guess we could just see what pops up when we test for significnace
ggplot(IL13avgiii,aes(y=log2,x=Notes))+geom_boxplot()+theme(axis.text.x = element_text(angle = 90))

##get pvalues by getting barcode distributions
#replicate1
rep1DNA<-read.csv(file="./KD1May25/IL13/1/IL13_1_DNA_counts.tsv",sep="\t",header=F)
head(rep1DNA)
rep1RNA<-read.csv(file="./KD1May25/IL13/1/IL13_1_RNA_counts.tsv",sep="\t",header=F)

setwd("~/Documents/MPRAflow/OYB2 Feb 2021/KD1/newAssocVariableLengtMNM0/")
filteredKD1<-read.csv(file="filteredCoordsToBarcodeFrac05.csv",header=F)
beforeFiltering<-read.csv(file="beforefilteringCoordsToBarcodes.csv",header=F)
head(beforeFiltering)
wt<-subset(beforeFiltering,beforeFiltering$V1=="ID770")
#merge with DNA and RNA
#rep1DNAii<-merge(rep1DNA,filteredKD1,by.x = "V2", by.y="V1")
rep1DNAID<-apply(rep1DNA,1,function(x){substr(x,start=9,stop=23)})
rep1DNAcount<-apply(rep1DNA,1,function(x){str_trim(substr(x,start=1,stop=8))})
rep1DNAcountii<-cbind(rep1DNAID,rep1DNAcount)

rep1RNAID<-apply(rep1RNA,1,function(x){substr(x,start=9,stop=23)})
rep1RNAcount<-apply(rep1RNA,1,function(x){str_trim(substr(x,start=1,stop=8))})
rep1RNAcountii<-cbind(rep1RNAID,rep1RNAcount)

rep1RNARatio<-merge(rep1RNAcountii,rep1DNAcountii,by.x="rep1RNAID",by.y="rep1DNAID")
rep1RNARatiodf<-as.data.frame(rep1RNARatio)
rep1RNARatioii<-merge(rep1RNARatiodf,filteredKD1,by.x="rep1RNAID",by.y="V2")
rep1RNARatioii$ratio<-as.numeric(rep1RNARatioii$rep1RNAcount)/as.numeric(rep1RNARatioii$rep1DNAcount)
allIDsKD1<-unique(rep1RNARatioii$V1)
df<-rep1RNARatioii[c(1:length(allIDsKD1)),c(1,2)]
y<-subset(rep1RNARatioii,rep1RNARatioii$V1=="ID770")
for(id in c(1:length(allIDsKD1))){
  x<-subset(rep1RNARatioii,rep1RNARatioii$V1==rep1RNARatioii$V1[id])
  x<-x$ratio
  testResult<-wilcox.test(x, y, alternative = "two.sided")
  df[id,1]<-allIDsKD1[id]
  df[id,2]<-testResult$p.value
}
#df
colnames(df)<-c("ID","pvalue2sided")
KD1$ID<-paste("ID",KD1$Number,sep="")
dfKD1rep1<-merge(df,KD1,by="ID")
dfKD1rep1$log10<--log10(as.numeric(dfKD1rep1$pvalue2sided))
ggplot(dfKD1rep1,aes(x=Number,y=log10))+geom_point()
rep1RNAIDiv<-merge(df,rep1RNARatioii,by.x ="barcode",by.y="rep1RNAID")
