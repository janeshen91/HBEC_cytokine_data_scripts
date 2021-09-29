library(Seurat)
library(MAST)
library(ggplot2)
library(dplyr)
library(data.table)
library(ggridges)
library(stringr)
colnames(soTreatmentBD2@meta.data)
table(soTreatmentBD2@meta.data$cellType)
table(soTreatmentBD2@meta.data$cellType_Treatment)
#subset to only the 3 celltypes and 3 treatments
soThreeCellTypes<-subset(soTreatmentBD2, idents = c("Basal", "Secretory","Ciliated"))
Idents(soThreeCellTypes)<-soThreeCellTypes@meta.data$Treatment
soThreeCellTypes2<-subset(soThreeCellTypes, idents = c("UT", "IL17","IL13", "IFN"))
table(soThreeCellTypes2@meta.data$Treatment)
Idents(soThreeCellTypes2)
rm(soThreeCellTypes)
###
#downsamples
table(soThreeCellTypes2@meta.data$cellType,soThreeCellTypes2@meta.data$Treatment)
#BASAL: 1197, SECRETORY: 753, CILIATED: 405
basalNumbers<-1144 #IL13
secrectoryNumbers<-753 #IL17
CiliatedNumbers<-405 #IL17
metadata<-soThreeCellTypes2@meta.data[,c("cellType","Treatment","cellType_Treatment")]
#how to name them so that I can sample
cellsToSample<-list()

#can we just randomly* get like 50 samples or something?

metadata$RowNum<-seq(from=1,to=11580,by=1)
smallSubset400<-sample(metadata$RowNum,size=400)
new_df <- metadata %>% group_by(cellType_Treatment) %>% sample_n(50)
new_df$subsampled<-"TRUE"
new_df2<-merge(metadata,new_df[,c(4,5)],by="RowNum",all.x=T)
rownames(new_df2)<-rownames(metadata)
AddMetaData(soThreeCellTypes2[,"subsampled"],new_df2,col.name="subsampled")
soThreeCellTypes2@meta.data$subsampled<-new_df2$subsampled
head(soThreeCellTypes2@meta.data)

###
geneNames = rownames(as.matrix(soTreatmentBD2@assays$RNA@data))
names(geneNames) = geneNames
geneNames = as.data.frame(geneNames)




dataMat = as.matrix(soThreeCellTypes2@assays$RNA@counts)

dataMat = (dataMat/colSums((dataMat)))*1e6
dataMat = log2(dataMat + 1)

# cellToTestList = levels(soThreeCellTypes2@meta.data$cellType)[c(2,4,6)]
# cellList = rownames(soThreeCellTypes2@meta.data[(soThreeCellTypes2@meta.data$cellType %in% cellToTestList),])
cellList<-rownames(subset(new_df2,new_df2$subsampled==TRUE))

conds = soThreeCellTypes2@meta.data[cellList,]
conds$cellType = factor(conds$cellType,levels=cellToTestList)

dataMat = dataMat[,cellList]
geneNames<-as.data.frame(geneNames)

scaRaw <- FromMatrix(exprsArray =  dataMat, 
                     cData = conds, 
                     fData = geneNames)

sca = scaRaw

scaSample <- sca[sample(which(freq(sca)>.1), 20),]
flat <- as(scaSample, 'data.table')
ggplot(flat, aes(x=value))+geom_density() +facet_wrap(~geneNames, scale='free_y')
range(assay(sca))

##problem with range?

thres <- thresholdSCRNACountMatrix(assay(sca), cutbins=c(.01,.1,.2,.5,1,2,3,4,5,6,7,20), min_per_bin = 30)
thres <- thresholdSCRNACountMatrix(assay(sca), min_per_bin = 30)
par(mfrow=c(5,4))
plot(thres)
##

freq_expressed <- 0.01
FCTHRESHOLD <- log2(1.25)

freqSca = freq(sca)
head(freqSca)
expressed_genes <- freqSca > freq_expressed
sca <- sca[expressed_genes,]

cond<-factor(colData(sca)$Treatment)
cond<-relevel(cond,"UT")
colData(sca)$Treatment<-cond

cond<-factor(colData(sca)$cellType)
cond<-relevel(cond,"Basal")
colData(sca)$cellType<-cond

cond<-factor(colData(sca)$IL13)
cond<-relevel(cond,"No")
colData(sca)$IL13<-cond

cond<-factor(colData(sca)$IL17)
cond<-relevel(cond,"No")
colData(sca)$IL17<-cond

cond<-factor(colData(sca)$IFN)
cond<-relevel(cond,"No")
colData(sca)$IFN <-cond


cdr2 <-colSums(assay(sca)>0)
qplot(x=cdr2, y=colData(sca)$nGeneOn) + xlab('New CDR') + ylab('Old CDR')
##can't find "nGeneOn"
colData(sca)$cngeneson <- scale(cdr2)
gc()
options(mc.cores=5)

zlmMod <- zlm(~ cellType + IL13 + IL17 + IFN + IL13*IFN + 
                IL13*cellType + IL17*cellType + IFN*cellType  +
                IL13*IFN*cellType + batch +
                Ind + cngeneson, sca, parallel = T)

sca@colData$cellType_Treatment = factor(paste(sca@colData$cellType, sca@colData$Treatment,sep="_"))
sca@colData$cellType_Treatment = relevel(sca@colData$cellType_Treatment,"Basal_UT")

zlmCond <- zlm(~ cellType_Treatment + batch + Ind + cngeneson, 
               sca, parallel = T)

mastTables = c()
mastTablesDE = c()
dir.create("DE/DEfiles/mastAllTerms",recursive = T)
dir.create("DE/DEfiles/mastAllTerms/All", recursive = T)
dir.create("DE/DEfiles/mastAllTerms/DEonly")

sigcounts <- as.data.frame(setNames(replicate(2,numeric(0), simplify = F),letters[0:2]))

condsNames = colnames(zlmMod@coefD)
#condsNames = condsNames[c(2:12,25:56)]
#condsNames = condsNames[c(26:37)]
condsNames = condsNames[c(19:24)]
for (i in condsNames) {
  summaryTable <- summary(zlmMod, doLRT=i)
  summaryTableDt <- summaryTable$datatable
  fcHurdle <- merge(summaryTableDt[contrast==i & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                    summaryTableDt[contrast==i & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid')
  fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  setorder(fcHurdle, fdr)
  
  fcHurdleSig <- fcHurdle[fdr<.1 & abs(coef)>FCTHRESHOLD]
  setorder(fcHurdleSig, fdr)
  
  mastTables[[i]] <- fcHurdle
  mastTablesDE[[i]] = fcHurdleSig
  
  fdrPass = nrow(fcHurdleSig)
  rawPass = nrow(fcHurdle[(fcHurdle$`Pr(>Chisq)` < 0.1),])
  sigcounts <- rbind(sigcounts,c(fdrPass,rawPass))
  
  print(i)
  print(c(fdrPass,rawPass))
  print(head(fcHurdleSig))
  
  write.table(fcHurdle,file = paste("./DE/DEfiles/mastAllTerms/All/",i,".all.full.txt",sep=""),row.names = F,col.names = T,sep = "\t",quote=F)
  write.table(fcHurdleSig,file = paste("./DE/DEfiles/mastAllTerms/DEonly/",i,".DE.full.txt",sep=""),row.names = F,col.names = T,sep = "\t",quote=F)
}


###um, so maybe the de etc. model is sort of worked out, let's see how we can find coverage for these datasets



rownames(soTreatmentBD2@assays$RNA)

countsPerGene<-rowSums(soTreatmentBD2@assays$RNA)


countsPerGene<-rowSums(soTreatmentBD2[["RNA"]]@counts)
#but we want counts distribution over one celltype/treatment 
treatmentConditions<-levels(as.factor(soThreeCellTypes2@meta.data$cellType_Treatment))
##save a list of counts of genes for each condition
GeneCounts<-list()
for (i in treatmentConditions){
  oneCondition<-subset(soThreeCellTypes2,subset=cellType_Treatment==i)

  countsPerGene<-colSums(oneCondition[["RNA"]]@counts)
  GeneCounts[[i]]<-countsPerGene
}
##plot all the basals together??

par(mfrow=c(2,2))
for(i in c(1:4)){
hist(log10(GeneCounts[[i]]),main=names(GeneCounts)[i],xlab="log10 counts of number of cells",ylim = c(0,700),breaks = seq(from=3,to=5,by=0.25))
}

par(mfrow=c(2,2))
for(i in c(5:8)){
  hist(log10(GeneCounts[[i]]),main=names(GeneCounts)[i],xlab="log10 counts of number of cells",ylim = c(0,400),breaks = seq(from=3,to=5,by=0.25))
}

par(mfrow=c(2,2))
for(i in c(9:12)){
  hist(log10(GeneCounts[[i]]),main=names(GeneCounts)[i],xlab="log10 counts of number of cells",ylim = c(0,500),breaks = seq(from=3,to=5,by=0.25))
}


###so we'll first downsample to thenumber of genes
###and then downsample th erad count by ranking and then downsampleing
###ranodmly?
celltypeTreatmentLevels<-levels(as.factor(metadata$cellType_Treatment))
metadata$cellNames<-rownames(metadata)
for (i in c(1:12)){
  print(i)
  Basal_IFN<-subset(metadata,metadata$cellType_Treatment==celltypeTreatmentLevels[i])
  
  if(i<=4){
    print("basal")
    sampledBasal<-sample_n(Basal_IFN,size=basalNumbers)
  }
  else if(i<=8){
    print("ciliated")
    sampledBasal<-sample_n(Basal_IFN,size=CiliatedNumbers)
  }
  else {
    print("Secretory")
    sampledBasal<-sample_n(Basal_IFN,size=secrectoryNumbers)
  }
  cellsToSample[[celltypeTreatmentLevels[i]]]<-sampledBasal
}

##now we have which rownumbers of cells to sample
##we then want like figure out what the counts in each cell are

###
###this loop goes throught the 12 celltype/treatment conditions and has a list of cell names and counts
orderedGeneCounts<-list()
for(i in c(1:12)){
basal1<-as.data.frame(GeneCounts[[i]])
basal1$cellNames<-rownames(basal1)
colnames(basal1)[1]<-"countsPerCell"
basal1<-basal1[order(as.numeric(basal1$countsPerCell)),]

sampledBasal<-cellsToSample[[i]]
basal1<-basal1 %>% arrange(countsPerCell) %>% inner_join (sampledBasal,by="cellNames") 
basal1$rank<-rownames(basal1)
orderedGeneCounts[[i]]<-basal1
}
#ok so we want something that's like, per celltype/condition
#CellNames,counts per cell, cellRowNames, cellChosenInFirstRoundTorF, RankInSecondRound, downsizeGoal,
#for every cell type, we now need to go through and 


dataMatAll<-as.matrix(soThreeCellTypes2[["RNA"]]@counts)
subListOfPools<-list()
#for(i in seq(from=1,to=12,by=4)){#i is the celltype/treatment condition
  #skip in 4
i=1
  basal1i<-merge(orderedGeneCounts[[i]][,c("cellNames","countsPerCell","rank")],orderedGeneCounts[[i+1]][,c("cellNames","countsPerCell","rank")],by="rank")
  basal2i<-merge(orderedGeneCounts[[i+2]][,c("cellNames","countsPerCell","rank")],orderedGeneCounts[[i+3]][,c("cellNames","countsPerCell","rank")],by="rank")
  basal3i<-merge(basal1i,basal2i,by="rank")
  minPerRow<-apply(cbind(basal3i$countsPerCell.x.x,basal3i$countsPerCell.x.y,basal3i$countsPerCell.y.x,basal3i$countsPerCell.y.y),1,min)
  basal3i$min<-minPerRow
  basal3iCellNames<-c(basal3i$cellNames.x.x,basal3i$cellNames.x.y,basal3i$cellNames.y.x,basal3i$cellNames.y.y)
  columns<-colnames((soThreeCellTypes2[["RNA"]]@counts)) %in% basal3iCellNames
  subsetedCells<-dataMatAll[,columns]
  #for each cell, figure out how many counts we're supposed to have and then downsample from there
  basal3iCellNamesAndCounts<-basal3i[,c("cellNames.x.x","cellNames.x.y","cellNames.y.x","cellNames.y.y","min")]
  for(j in c(1:dim(subsetedCells)[2])){
    #how do we get the number of reads we need
    #if Reads=number of reads we need
    #print(j)
    cellName<-colnames(subsetedCells)[j]
    matchesOf<-apply(basal3iCellNamesAndCounts[,c(1:4)],2,function(x){match(cellName,x)})
    rowNum<-matchesOf[!is.na(matchesOf)]
    numCounts<-basal3iCellNamesAndCounts[rowNum,5]
    #we want to iterate over each column and then each row
    #but only use the appropriate row names
    #and subsample over each column
    listOfPools<-list()
    poolOfStrings<-""
    for(k in c(1:dim(subsetedCells)[1])){
      if(subsetedCells[k,j]>0){
      #print(k)
      repeated<-rep(rownames(subsetedCells)[k],subsetedCells[k,j])
      #print(repeated)
      poolOfStrings<-c(poolOfStrings,repeated)
      }
      #print(head(poolOfStrings))
   }
    #listOfPools[[j]]<-poolOfStrings
    #print(j)
    #print(cellName)
    
    subListOfPools[[cellName]]<-sample(poolOfStrings,size=numCounts)
} #end for j in going through columns of cells
    
    
#}#end going over one to 12 conditions

##ok, should we go over one condition (ie i from 1 to4) and then assess whether it changed the way we wanted it to?
##we also need to go from subListOfPools back to numbers in a matrix
#ie frist get unique
#make a dataframe of columns=number of cells in that conditions and rows=number of genes
#

  dataFrameOfWordCounts<-as.data.frame(subsetedCells)
  geneNames<-rownames(dataFrameOfWordCounts)
  
  for(i in dim(dataFrameOfWordCounts)[1]){
    dataFrameOfWordCounts[i,]<-rep(0,times=dim(dataFrameOfWordCounts)[2])
  }

#ok, so dataframeofwordcounts right now do have only the cells that are in this condition (ie i = 1 to 4)
  
for (j in names(subListOfPools)){
  print(j)
  
  uniqueStrings<-unique(subListOfPools[[j]])
  for(s in uniqueStrings){
  print(s)
  if (s !=""){
  rowNum<-match(s,geneNames)
  colNum<-match(j,colnames(dataFrameOfWordCounts))
  dataFrameOfWordCounts[rowNum,colNum]<-sum(str_count(subListOfPools[[j]],pattern=s))
  }
  }
}
  
##ok, so how do we look at whether the dataFrameOfWordcounts is the new correctly reduced/downsampled counts of genes?
  #we need to match the cells to the conditions and  then plot them

  head(metadata)
  downSampledBasal<-CreateSeuratObject(counts=dataFrameOfWordCounts,meta.data = metadata)
  downSampledBasal
  save(dataFrameOfWordCounts,file="dataFrameOfWordCountsBasal.RData")
  
  treatmentConditions<-levels(as.factor(downSampledBasal@meta.data$cellType_Treatment))
  ##save a list of counts of genes for each condition
  GeneCounts<-list()
  for (i in treatmentConditions){
    oneCondition<-subset(downSampledBasal,subset=cellType_Treatment==i)
    
    countsPerGene<-colSums(oneCondition[["RNA"]]@counts)
    GeneCounts[[i]]<-countsPerGene
  }
  ##plot all the basals together??
  
  par(mfrow=c(2,2))
  for(i in c(1:4)){
    hist(log10(GeneCounts[[i]]),main=names(GeneCounts)[i],xlab="log10 counts of number of cells",ylim = c(0,700),breaks = seq(from=3,to=5,by=0.25))
  }

basalNumbers<-1144 #IL13
secrectoryNumbers<-753 #IL17
CiliatedNumbers<-405 #IL17
basalNumbers+secrectoryNumbers+CiliatedNumbers
##2302 cell sin total after sampling
#how to um sample from within that?
#we need to go and get that cell and then get the 


##ok, let's check it is the same number of cells for each rank
orderedGeneCounts<-list()
for(i in c(1:4)){
  basal1<-as.data.frame(GeneCounts[[i]])
  basal1$cellNames<-rownames(basal1)
  colnames(basal1)[1]<-"countsPerCell"
  basal1<-basal1[order(as.numeric(basal1$countsPerCell)),]
  
  sampledBasal<-cellsToSample[[i]]
  basal1<-basal1 %>% arrange(countsPerCell) %>% inner_join (sampledBasal,by="cellNames") 
  basal1$rank<-rownames(basal1)
  orderedGeneCounts[[i]]<-basal1
}
i=1
basal1i<-merge(orderedGeneCounts[[i]][,c("cellNames","countsPerCell","rank")],orderedGeneCounts[[i+1]][,c("cellNames","countsPerCell","rank")],by="rank")
basal2i<-merge(orderedGeneCounts[[i+2]][,c("cellNames","countsPerCell","rank")],orderedGeneCounts[[i+3]][,c("cellNames","countsPerCell","rank")],by="rank")
basal3i<-merge(basal1i,basal2i,by="rank")


#hwo to get a smaller subset of cells and genes and see if the script is doing what we think it is doing
