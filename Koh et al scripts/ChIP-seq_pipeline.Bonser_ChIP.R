library("cowplot")
library("stringr")
library(biomaRt)
library(dplyr)
library(DESeq2)
library(genefilter)
library("devtools")
library("BiocParallel")
library(Rtsne)
#register(MulticoreParam(4))
library(gage)
library(gageData)
data(kegg.sets.hs)
data(sigmet.idx.hs)
kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]
library(R.utils)
sourceDirectory("~/data/BatchArray-source/src")
source("~/data/BatchArray-source/load.R")
library(gplots)
library(ggplot2)
library(limma)
project <- c("Bonser_ChIP_cleanedPeaks")
dir <- paste(paste("~/data",project, sep = "/"),"Analysis_with_old_peakFix", sep = "/")
dir.create(dir)
setwd(dir)
library(gplots)
library(ggplot2)
library(limma)
library("calibrate")
library("rJava")
library("RDAVIDWebService")
library("pathview")
library(ComplexHeatmap)
library(dendextend)
library(circlize)


dir.create("DE")
dir.create("QC")
dir.create("QC/ChIPnorms/")
dir.create("DE/beds")
dir.create("DE/broadPeaks")
dir.create("DE/MAplots/")
dir.create("DE/VolcanoPlots/")
dir.create("DE/MMplots/")
dir.create("DE/pathways/")
dir.create("DE/Heatmaps/")
dir.create("DE/DEfiles/")
dir.create("DE/DAVID/")

file = "../counts.ChIP.HUMAN.txt"

chip_gene_pairs_TimeSeries <- read.table("../../Bonser_ChIP_02_2018_pilot/fixed_pilot_counts/Partitioned_peaks.closestTimeSeries.bed",header=F, sep="\t")

data <- read.table(file, header=T, row.names=1)

common <- read.table("../FinalMergedPeaks/Partitioned_peaks.cleaned.common.bed",sep = "\t")
all <- read.table("../FinalMergedPeaks/Partitioned_peaks.cleaned.bed",sep="\t")
print("got all files")


if (file  == "../counts.ChIP.HUMAN.txt"){
annot <- read.table("~/data/human_ens_GRCh38_83_10x_annot.extended.v2.txt", sep="\t", quote="", header=T, row.names=1, stringsAsFactors=F, fill=T)
	data(go.sets.hs)
	data(go.subs.hs)
	gobpsets = go.sets.hs[go.subs.hs$BP]
	gomfsets = go.sets.hs[go.subs.hs$MF]
	data(kegg.sets.hs)
	data(sigmet.idx.hs)
	kegg.sets = kegg.sets.hs[sigmet.idx.hs]
	spec = "hsa"
}
if (file  == "../counts.STARMOUSE.txt"){
	annot <- read.table("~/data/mouse_ens_GRCm38_annot.extended.txt", sep="\t", quote="", header=T, row.names=1, stringsAsFactors=F, fill=T)
	data(go.sets.mm)
	data(go.subs.mm)
	gobpsets = go.sets.mm[go.subs.mm$BP]
	gomfsets = go.sets.mm[go.subs.mm$MF]
	data(kegg.sets.mm)
	data(sigmet.idx.mm)
	kegg.sets = kegg.sets.mm[sigmet.idx.mm]
	spec = "mmu"
}
if (file  == "../counts.STARMOUSEMM9.txt" | file == "../counts.scRNA.STARMOUSEMM9.txt"){
	annot <- read.table("~/data/mouse_ens_GRCm37.annot.extended.txt", sep="\t", quote="", header=T, row.names=1, stringsAsFactors=F, fill=T)
	data(go.sets.mm)
	data(go.subs.mm)
	gobpsets = go.sets.mm[go.subs.mm$BP]
	gomfsets = go.sets.mm[go.subs.mm$MF]
	data(kegg.sets.mm)
	data(sigmet.idx.mm)
	kegg.sets = kegg.sets.mm[sigmet.idx.mm]
	spec = "mmu"
}
if (file  == "../counts.STARZEBRAFISH.txt"){
	annot <- read.table("~/data/zfish_ens_GRCz10_annot.extended.txt", sep="\t", quote="", header=T, row.names=1, stringsAsFactors=F, fill=T)
}

#dataOld <- read.table("../../Bonser_ChIP_02_2018_pilot/fixed_pilot_counts/counts.ChIP.HUMAN.txt", header=T, row.names=1)
#colnames(dataOld) <- c("D28_13_K27ac_repHBES","D28_13_K27ac_repHBET","D28_UT_K27ac_repHBES","D28_UT_K27ac_repHBET")

#dataAll <- merge(data, dataOld, by="row.names")
#rownames(dataAll) <-dataAll[,1]
#dataAll <-dataAll[,-c(1)]

dataAll = data[,-c(1,10,22,35)]

colnames(dataAll)[c(1,32)] = c("D28_13_K27ac_repHBET","D28_UT_K27ac_repHBET")

#rownames(annot) <- annot$Ensembl_ID

#data <- dataAll[,-c(8,20,33,35)]

data = dataAll

targets <- as.data.frame(as.matrix(colnames(data)))

conds <- cbind(targets, do.call(rbind, strsplit(as.character(targets$V1),'_')))
rownames(conds) <- conds[,c(1)]
conds <- conds[,-c(1)]

colnames(conds) <- c("Time","Treatment","ChIP","Ind")
conds$Treatment <- gsub("13","IL13",conds$Treatment)
conds$Treatment <- gsub("17","IL17",conds$Treatment)

conds$Time <- gsub("D0","D00",conds$Time)
conds$Time <- gsub("D3","D03",conds$Time)
conds$Time <- gsub("D9","D09",conds$Time)
conds$Time <- gsub("D28","D24",conds$Time)



conds$condition <- paste(conds$Time, conds$Treatment, sep="_")

conds$Treatment <- as.factor(conds$Treatment)
conds[,"Treatment"] <- factor(conds[,"Treatment"], levels=rev(levels(conds[,"Treatment"])))
#conds$class <- as.factor(conds$class)
#conds[,"class"] <- factor(conds[,"class"], levels=rev(levels(conds[,"class"])))
#
#conds$condition <- as.factor(paste(conds$genotype, conds$class, sep="_"))
#
#
conds$Time <- as.factor(conds$Time)
conds$condition <- as.factor(conds$condition)
conds[,"condition"] <- factor(conds[,"condition"], levels=levels(conds[,"condition"])[c(1,2,3,8,7,6,5,4)])

conds = conds[c(2:31,1,32),]
data = data[,c(2:31,1,32)]

conds$IL13 <- c(rep("N",4), rep("Y",3), rep("N",4), rep("Y",4),rep("N",15),"Y","N")
conds$IL17 <- c(rep("N",7), rep("Y",4),rep("N",21))
conds$IFN <- c(rep("N",11), rep("Y",7),rep("N",14))

colnames(data) = gsub("D28","D24",colnames(data))
rownames(conds) = gsub("D28","D24",rownames(conds))

dataFiltered <- data[ rowSums( data >= 1 ) >= 4, ]
dataFiltered <- dataFiltered[ rowSums(dataFiltered) > 50, ]

all <- all[which(all$V4 %in% rownames(dataFiltered)),]

write.table(conds, file=paste("conditions.", project, ".txt", sep=""), sep="\t", quote=F, row.names=T, col.names=NA)



ddsModel <- DESeqDataSetFromMatrix(countData = dataFiltered,
							  colData = conds,
							  design = ~  Ind + Time + IL13 + IL17 + IFN + IL13*IFN)


ddsConds <- DESeqDataSetFromMatrix(countData = dataFiltered,
							  colData = conds,
							  design = ~  Ind + condition)


ddsModel <- DESeq(ddsModel)


#RldOld <- rlogTransformation(ddsModel, blind=FALSE)
normData <- t(t(dataFiltered) / sizeFactors(ddsModel))
N = ncol(normData)
normData <- addmargins(normData, FUN = list(Total = sum), quiet = TRUE)
n<-dim(normData)[1]
normData<-normData[1:(n-1),]
TAll = (as.data.frame(normData)$Total)/N
print("got norms")

rownames(all) <- all$V4
commonName <- common$V4
allName <- all$V4
normDataCommon <- as.data.frame(normData)[commonName,]

Tcommon = (as.data.frame(normDataCommon)$Total)/N

J = N + 1
K = N - 1

l = 1
normMatrix <- as.matrix(normData)[,-1:-J]
for(i in names(dataFiltered)) {
  S = as.data.frame(normDataCommon)[i]
  M = log((S+1)/(Tcommon+1), 2)
  A = log(Tcommon+1,2)
  mat = cbind(M,A)
  #mat = mat[!rowSums(!is.finite(mat)),]
  Mi = as.data.frame(mat)[1]
  Ai = as.data.frame(mat)[2]
  Mname = colnames(mat)[1]
  Aname = colnames(mat)[2]
  b <- rlm(as.matrix(Mi)~as.matrix(Ai))$coefficients
  c = lm(as.matrix(Mi)~as.matrix(Ai))$coefficients
  pdf(paste("QC/ChIPnorms/", i, ".model.MAplot_before_rescaling.pdf", sep=""))
  p = ggplot(mat, aes(x = eval(parse(text=Aname)),y = eval(parse(text=Mname))), alpha = I(1/100), xlim=c(0,20), ylim=c(-4,4)) + geom_point()
  print(p + geom_abline(intercept = b[1], slope = b[2], col = "blue") + geom_abline(intercept = c[1], slope = c[2], col="green"))
  dev.off()
  SAll <- as.data.frame(normData)[i]
  SAllLog <- log(SAll+1,2)
  S1Log <- (2-b[2])*(SAllLog)/(2+b[2]) - 2*b[1]/(2+b[2])
  S1raw = (2^S1Log)
  ratio = S1raw/SAll
  normFactor = 1/ratio * as.data.frame.list(sizeFactors(ddsModel))[i][1,]
  normMatrix <- cbind(normMatrix, normFactor)
  l= l + 1
}



normMatrix[normMatrix == 0] <- 1
ddsModel <- DESeqDataSetFromMatrix(countData = dataFiltered, colData = conds, design = ~ Ind + Time + IL13 + IL17 + IFN + IL13*IFN)
normalizationFactors(ddsModel) <- as.matrix(normMatrix)
ddsModel <- DESeq(ddsModel)

#Rld <- rlogTransformation(ddsModel, blind=FALSE)


##### round 2 MA plot generation #####

normCounts <- DESeq2::counts(ddsModel, normalized=T)
N = ncol(normCounts)
J = N + 1
K = N - 1
normCounts <- addmargins(normCounts, FUN = list(Total = sum), quiet = TRUE)
n<-dim(normCounts)[1]
normCounts<-normCounts[1:(n-1),]
TAll = (as.data.frame(normCounts)$Total)/N

#common <- read.table(args[4],sep = "\t")
#all <- read.table(args[5],sep="\t")
#rownames(all) <- all$V4
#commonName <- common$V4

normCountsAll <- as.data.frame(normCounts)[allName,]

Taverage = (as.data.frame(normCountsAll)$Total)/N
normCountsAll <- normCountsAll[,-J]


for(i in names(normCountsAll)) {
  S = as.data.frame(normCountsAll)[i]
  M = log((S+1)/(Taverage+1), 2)
  A = log(Taverage+1,2)
  mat = cbind(M,A)
  #mat = mat[!rowSums(!is.finite(mat)),]
  Mi = as.data.frame(mat)[1]
  Ai = as.data.frame(mat)[2]
  Mname = colnames(mat)[1]
  Aname = colnames(mat)[2]
  b <- rlm(as.matrix(Mi)~as.matrix(Ai))$coefficients
  c = lm(as.matrix(Mi)~as.matrix(Ai))$coefficients
  pdf(paste("QC/ChIPnorms/", i, ".model.MAplot_after_rescaling.pdf", sep=""))
  p = ggplot(mat, aes(x = eval(parse(text=Aname)),y = eval(parse(text=Mname))), alpha = I(1/100), xlim=c(0,20), ylim=c(-4,4)) + geom_point()
  print(p + geom_abline(intercept = b[1], slope = b[2], col = "blue") + geom_abline(intercept = c[1], slope = c[2], col="green"))
  dev.off()
  #SAll <- as.data.frame(normData)[i]
  #SAllLog <- log(SAll+1,2)
  #S1Log <- (2-b[2])*(SAllLog)/(2+b[2]) - 2*b[1]/(2+b[2])
  #S1raw = (2^S1Log)
  #ratio = S1raw/SAll
  #normFactor = 1/ratio * as.data.frame.list(sizeFactors(ddsModel))[i][1,]
  #normMatrix <- cbind(normMatrix, normFactor)
}



ddsConds <- DESeq(ddsConds)


#RldOld <- rlogTransformation(ddsConds, blind=FALSE)

normData <- t(t(dataFiltered) / sizeFactors(ddsConds))
N = ncol(normData)
normData <- addmargins(normData, FUN = list(Total = sum), quiet = TRUE)
n<-dim(normData)[1]
normData<-normData[1:(n-1),]
TAll = (as.data.frame(normData)$Total)/N
print("got norms")
rownames(all) <- all$V4
commonName <- common$V4
allName <- all$V4
normDataCommon <- as.data.frame(normData)[commonName,]

Tcommon = (as.data.frame(normDataCommon)$Total)/N

J = N + 1
K = N - 1

l = 1
normMatrix <- as.matrix(normData)[,-1:-J]
for(i in names(dataFiltered)) {
  S = as.data.frame(normDataCommon)[i]
  M = log((S+1)/(Tcommon+1), 2)
  A = log(Tcommon+1,2)
  mat = cbind(M,A)
  #mat = mat[!rowSums(!is.fini  te(mat)),]
  Mi = as.data.frame(mat)[1]
  Ai = as.data.frame(mat)[2]
  Mname = colnames(mat)[1]
  Aname = colnames(mat)[2]
  b <- rlm(as.matrix(Mi)~as.matrix(Ai))$coefficients
  c = lm(as.matrix(Mi)~as.matrix(Ai))$coefficients
  pdf(paste("QC/ChIPnorms/", i, ".conds.MAplot_before_rescaling.pdf", sep=""))
  p = ggplot(mat, aes(x = eval(parse(text=Aname)),y = eval(parse(text=Mname))), alpha = I(1/100), xlim=c(0,20), ylim=c(-4,4)) + geom_point()
  print(p + geom_abline(intercept = b[1], slope = b[2], col = "blue") + geom_abline(intercept = c[1], slope = c[2], col="green"))
  dev.off()
  SAll <- as.data.frame(normData)[i]
  SAllLog <- log(SAll+1,2)
  S1Log <- (2-b[2])*(SAllLog)/(2+b[2]) - 2*b[1]/(2+b[2])
  S1raw = (2^S1Log)
  ratio = S1raw/SAll
  normFactor = 1/ratio * as.data.frame.list(sizeFactors(ddsConds))[i][1,]
  normMatrix <- cbind(normMatrix, normFactor)
  l= l + 1
}



normMatrix[normMatrix == 0] <- 1
ddsConds <- DESeqDataSetFromMatrix(countData = dataFiltered, colData = conds, design = ~  Ind + condition)
normalizationFactors(ddsConds) <- as.matrix(normMatrix)
ddsConds <- DESeq(ddsConds)


#### YOU CUT HERE ####


#Rld <- rlogTransformation(ddsConds, blind=FALSE)


##### round 2 MA plot generation #####

normCounts <- DESeq2::counts(ddsConds, normalized=T)
N = ncol(normCounts)
J = N + 1
K = N - 1
normCounts <- addmargins(normCounts, FUN = list(Total = sum), quiet = TRUE)
n<-dim(normCounts)[1]
normCounts<-normCounts[1:(n-1),]
TAll = (as.data.frame(normCounts)$Total)/N

#common <- read.table(args[4],sep = "\t")
#all <- read.table(args[5],sep="\t")
#rownames(all) <- all$V4
#commonName <- common$V4

normCountsAll <- as.data.frame(normCounts)[allName,]

Taverage = (as.data.frame(normCountsAll)$Total)/N
normCountsAll <- normCountsAll[,-J]


for(i in names(normCountsAll)) {
  S = as.data.frame(normCountsAll)[i]
  M = log((S+1)/(Taverage+1), 2)
  A = log(Taverage+1,2)
  mat = cbind(M,A)
  #mat = mat[!rowSums(!is.finite(mat)),]
  Mi = as.data.frame(mat)[1]
  Ai = as.data.frame(mat)[2]
  Mname = colnames(mat)[1]
  Aname = colnames(mat)[2]
  b <- rlm(as.matrix(Mi)~as.matrix(Ai))$coefficients
  c = lm(as.matrix(Mi)~as.matrix(Ai))$coefficients
  pdf(paste("QC/ChIPnorms/", i, ".conds.MAplot_after_rescaling.pdf", sep=""))
  p = ggplot(mat, aes(x = eval(parse(text=Aname)),y = eval(parse(text=Mname))), alpha = I(1/100), xlim=c(0,20), ylim=c(-4,4)) + geom_point()
  print(p + geom_abline(intercept = b[1], slope = b[2], col = "blue") + geom_abline(intercept = c[1], slope = c[2], col="green"))
  dev.off()
  #SAll <- as.data.frame(normData)[i]
  #SAllLog <- log(SAll+1,2)
  #S1Log <- (2-b[2])*(SAllLog)/(2+b[2]) - 2*b[1]/(2+b[2])
  #S1raw = (2^S1Log)
  #ratio = S1raw/SAll
  #normFactor = 1/ratio * as.data.frame.list(sizeFactors(ddsConds))[i][1,]
  #normMatrix <- cbind(normMatrix, normFactor)
}


norm_counts<- round(DESeq2::counts(ddsModel, normalized=TRUE),1)
allgenes <- norm_counts
stabilizedBlind <- assay(varianceStabilizingTransformation(ddsModel, blind=TRUE))
stabilized <- assay(varianceStabilizingTransformation(ddsModel, blind=FALSE))
stabilizedBR <- removeBatchEffect(as.data.frame(stabilized), batch=conds$Ind)
allgenes <- merge(allgenes, stabilized, by="row.names", suffixes=c(".norm",".stabilized"))
row.names(allgenes) <- allgenes$Row.names
allgenes <- allgenes[,-c(1)]


rv <- rowVars(stabilizedBlind)
select <- order(rv, decreasing=T)[seq_len(5000)]
pca <- prcomp(t(stabilizedBlind[select,]))

list = c(1,2,3,4,5)
for (i in 2:length(list)) {
  i
}
for (i in 2:length(list)) {                 
  pcaobject <- as.data.frame(pca$x[,c(1,i)])
  pcaobject$condition <- conds$condition
  pcaobject$Individual <- conds$Ind
  
  comps <- colnames(pcaobject)[1:2]
  labels <- as.vector(pca$sdev^2/sum(pca$sdev^2))
  xlabel <- round(labels[1]*100,1)
  ylabel <- round(labels[i]*100,1)
  #print(ylabel)
  #print(paste("PC",i,sep=""))
  ggplot2::ggplot(pcaobject, aes_string(x = comps[1], y = comps[2])) +
    geom_point(aes(color=condition, shape=Individual), size=6) +
    xlab(paste0("PC1: ", xlabel, "% variance")) +
    ylab(paste0("PC",i,": ", ylabel, "% variance")) +
    theme_bw()
  ggsave(paste("QC/PCAplot.PC1vPC",i,".pdf",sep=""), height=7, width=7, units=c("in"))}
for (i in 3:length(list)) {                 
  pcaobject <- as.data.frame(pca$x[,c(2,i)])
  pcaobject$condition <- conds$condition
  pcaobject$Individual <- conds$Ind
  
  comps <- colnames(pcaobject)[1:2]
  labels <- as.vector(pca$sdev^2/sum(pca$sdev^2))
  xlabel <- round(labels[2]*100,1)
  ylabel <- round(labels[i]*100,1)
  #print(ylabel)
  #print(paste("PC",i,sep=""))
  ggplot2::ggplot(pcaobject, aes_string(x = comps[1], y = comps[2])) +
    geom_point(aes(color=condition, shape=Individual), size=6) +
    xlab(paste0("PC2: ", xlabel, "% variance")) +
    ylab(paste0("PC",i,": ", ylabel, "% variance")) +
    theme_bw()
  ggsave(paste("QC/PCAplot.PC2vPC",i,".pdf",sep=""), height=7, width=7, units=c("in"))}
for (i in 4:length(list)) {                 
  pcaobject <- as.data.frame(pca$x[,c(3,i)])
  #print(i)
  pcaobject$condition <- conds$condition
  pcaobject$Individual <- conds$Ind
  
  comps <- colnames(pcaobject)[1:2]
  labels <- as.vector(pca$sdev^2/sum(pca$sdev^2))
  xlabel <- round(labels[3]*100,1)
  ylabel <- round(labels[i]*100,1)
  #print(ylabel)
  #print(paste("PC",i,sep=""))
  ggplot2::ggplot(pcaobject, aes_string(x = comps[1], y = comps[2])) +
    geom_point(aes(color=condition, shape=Individual), size=6) +
    xlab(paste0("PC3: ", xlabel, "% variance")) +
    ylab(paste0("PC",i,": ", ylabel, "% variance")) +
    theme_bw()
  ggsave(paste("QC/PCAplot.PC3vPC",i,".pdf",sep=""), height=7, width=7, units=c("in"))}


#ggplot(pcaobject, aes(pcaobject[,c(1)], pcaobject[,c(2)], color=GenotypeL)) +
#	geom_point(size=3) +
#	xlab(paste0("PC1: ", xlabel, "% variance")) +
#	ylab(paste0("PC",i,": ", ylabel, "% variance"))
#ggsave(paste("QC/PCAplot.PC3vPC",i,".pdf",sep=""), height=7, width=7, units=c("in"))



M <- stabilizedBlind[select,]

list = c(3,4,5,6,8,10,12,15,20,30,40,50,75,100)
for (i in list){
  tsne_out <- Rtsne::Rtsne(t(M), perplexity = i, max.iter=100000)
  df_to_plot <- as.data.frame(tsne_out$Y)
  comps <- colnames(df_to_plot)[1:2]
  df_to_plot$condition <- conds$condition
  df_to_plot$Individual <- conds$Ind
  
  ggplot2::ggplot(df_to_plot, aes_string(x = comps[1], y = comps[2])) +
    geom_point(aes(color=condition, shape=Individual), size=6) +
    xlab("Dimension 1") +
    ylab("Dimension 2") +
    scale_shape_manual(values=c(15,16,17,18,8)) + 
    stat_ellipse(aes(color=condition, fill=condition), geom = "polygon", alpha = 0.3, linetype=0,level = .8) +
    theme_bw()
  ggsave(paste("QC/Tsne.plot.", i, ".elipse.pdf", sep=""), height=7, width=7, units=c("in"))
}

list = c(3,4,5,6,8,10,12,15,20,30,40,50,75,100)
for (i in list){
  tsne_out <- Rtsne::Rtsne(t(M), perplexity = i, max.iter=100000)
  df_to_plot <- as.data.frame(tsne_out$Y)
  comps <- colnames(df_to_plot)[1:2]
  df_to_plot$condition <- conds$condition
  df_to_plot$Individual <- conds$Ind
  
  ggplot2::ggplot(df_to_plot, aes_string(x = comps[1], y = comps[2])) +
    geom_point(aes(color=condition, shape=Individual), size=6) +
    xlab("Dimension 1") +
    ylab("Dimension 2") +
    scale_shape_manual(values=c(15,16,17,18,8)) + 
    #stat_ellipse(aes(color=condition, fill=condition), geom = "polygon", alpha = 0.3, linetype=0,level = .8) +
    theme_bw()
  ggsave(paste("QC/Tsne.plot.", i, ".pdf", sep=""), height=7, width=7, units=c("in"))
}

###


rvBR <- rowVars(stabilizedBR)
select <- order(rvBR, decreasing=T)[seq_len(5000)]
pca <- prcomp(t(stabilizedBR[select,]))

list = c(1,2,3,4,5)
for (i in 2:length(list)) {
  i
}
for (i in 2:length(list)) {                 
  pcaobject <- as.data.frame(pca$x[,c(1,i)])
  pcaobject$condition <- conds$condition
  pcaobject$Individual <- conds$Ind
  
  comps <- colnames(pcaobject)[1:2]
  labels <- as.vector(pca$sdev^2/sum(pca$sdev^2))
  xlabel <- round(labels[1]*100,1)
  ylabel <- round(labels[i]*100,1)
  #print(ylabel)
  #print(paste("PC",i,sep=""))
  ggplot2::ggplot(pcaobject, aes_string(x = comps[1], y = comps[2])) +
    geom_point(aes(color=condition, shape=Individual), size=6) +
    xlab(paste0("PC1: ", xlabel, "% variance")) +
    ylab(paste0("PC",i,": ", ylabel, "% variance")) +
    theme_bw()
  ggsave(paste("QC/PCAplot.PC1vPC",i,".IndBR.pdf",sep=""), height=7, width=7, units=c("in"))}
for (i in 3:length(list)) {                 
  pcaobject <- as.data.frame(pca$x[,c(2,i)])
  pcaobject$condition <- conds$condition
  pcaobject$Individual <- conds$Ind
  
  comps <- colnames(pcaobject)[1:2]
  labels <- as.vector(pca$sdev^2/sum(pca$sdev^2))
  xlabel <- round(labels[2]*100,1)
  ylabel <- round(labels[i]*100,1)
  #print(ylabel)
  #print(paste("PC",i,sep=""))
  ggplot2::ggplot(pcaobject, aes_string(x = comps[1], y = comps[2])) +
    geom_point(aes(color=condition, shape=Individual), size=6) +
    xlab(paste0("PC2: ", xlabel, "% variance")) +
    ylab(paste0("PC",i,": ", ylabel, "% variance")) +
    theme_bw()
  ggsave(paste("QC/PCAplot.PC2vPC",i,".IndBR.pdf",sep=""), height=7, width=7, units=c("in"))}
for (i in 4:length(list)) {                 
  pcaobject <- as.data.frame(pca$x[,c(3,i)])
  #print(i)
  pcaobject$condition <- conds$condition
  pcaobject$Individual <- conds$Ind
  
  comps <- colnames(pcaobject)[1:2]
  labels <- as.vector(pca$sdev^2/sum(pca$sdev^2))
  xlabel <- round(labels[3]*100,1)
  ylabel <- round(labels[i]*100,1)
  #print(ylabel)
  #print(paste("PC",i,sep=""))
  ggplot2::ggplot(pcaobject, aes_string(x = comps[1], y = comps[2])) +
    geom_point(aes(color=condition, shape=Individual), size=6) +
    xlab(paste0("PC3: ", xlabel, "% variance")) +
    ylab(paste0("PC",i,": ", ylabel, "% variance")) +
    theme_bw()
  ggsave(paste("QC/PCAplot.PC3vPC",i,".IndBR.pdf",sep=""), height=7, width=7, units=c("in"))}


#ggplot(pcaobject, aes(pcaobject[,c(1)], pcaobject[,c(2)], color=GenotypeL)) +
#	geom_point(size=3) +
#	xlab(paste0("PC1: ", xlabel, "% variance")) +
#	ylab(paste0("PC",i,": ", ylabel, "% variance"))
#ggsave(paste("QC/PCAplot.PC3vPC",i,".pdf",sep=""), height=7, width=7, units=c("in"))



M <- stabilizedBR[select,]

list = c(3,4,5,6,8,10,12,15,20,30,40,50,75,100)
for (i in list){
  tsne_out <- Rtsne::Rtsne(t(M), perplexity = i, max.iter=100000)
  df_to_plot <- as.data.frame(tsne_out$Y)
  comps <- colnames(df_to_plot)[1:2]
  df_to_plot$condition <- conds$condition
  df_to_plot$Individual <- conds$Ind
  
  ggplot2::ggplot(df_to_plot, aes_string(x = comps[1], y = comps[2])) +
    geom_point(aes(color=condition, shape=Individual), size=6) +
    xlab("Dimension 1") +
    ylab("Dimension 2") +
    stat_ellipse(aes(color=condition, fill=condition), geom = "polygon", alpha = 0.3, linetype=0,level = .8) +
    scale_shape_manual(values=c(15,16,17,18,8)) + 
    theme_bw()
  ggsave(paste("QC/Tsne.plot.", i, ".IndBR.elipse.pdf", sep=""), height=7, width=7, units=c("in"))
}


M <- stabilizedBR[select,]

list = c(3,4,5,6,8,10,12,15,20,30,40,50,75,100)
for (i in list){
  tsne_out <- Rtsne::Rtsne(t(M), perplexity = i, max.iter=100000)
  df_to_plot <- as.data.frame(tsne_out$Y)
  comps <- colnames(df_to_plot)[1:2]
  df_to_plot$condition <- conds$condition
  df_to_plot$Individual <- conds$Ind
  
  ggplot2::ggplot(df_to_plot, aes_string(x = comps[1], y = comps[2])) +
    geom_point(aes(color=condition, shape=Individual), size=6) +
    xlab("Dimension 1") +
    ylab("Dimension 2") +
    #stat_ellipse(aes(color=condition, fill=condition), geom = "polygon", alpha = 0.3, linetype=0,level = .8) +
    scale_shape_manual(values=c(15,16,17,18,8)) + 
    theme_bw()
  ggsave(paste("QC/Tsne.plot.", i, ".IndBR.pdf", sep=""), height=7, width=7, units=c("in"))
}



####



select <- order(rv, decreasing=T)[seq_len(5000)]
M <- stabilizedBlind[select,]
M <- M - rowMeans(M, na.rm=T)
#pairs.agilent(stabilizedBlind, filename="QC/pairs_log2normRC.png", main="Pairwise comparison of log2 normalized read counts", cex=2)


pdf("QC/boxplot.pdf")
par(mar=c(7,5,5,5))
boxplot(stabilizedBlind, las=2, col="blue", main="Log2 Normalized Read Count by Sample", pch=19, cex=.5)
dev.off()


DE.heatmap(M, dendrogram="column", height=1024, fact=conds[,c("Time","Treatment","Ind"), drop=F], lhei=c(1,.5,.5,.5,6,1), main="Heatmap of log2 FC rel. global average for top 1000 genes by variance", file.prefix="QC/QC_heatmap", breaks=.2)

select <- order(rvBR, decreasing=T)[seq_len(5000)]
M <- stabilizedBR[select,]
M <- M - rowMeans(M, na.rm=T)
DE.heatmap(M, dendrogram="column", height=1024, fact=conds[,c("Time","Treatment","Ind"), drop=F], lhei=c(1,.5,.5,.5,6,1), main="Heatmap of log2 FC rel. global average for top 3000 genes by variance", file.prefix="QC/QC_heatmap.IndBR", breaks=.2)


###############################
#### Pairwise comparisons #####


#condsOrder <- conds[order(conds[,"Time"],decreasing = T),]
#condsOrder <- condsOrder[order(condsOrder[,"Treatment"],decreasing = T),]

comps <- c("condition-D24_UTvsD00_UT","condition-D09_UTvsD00_UT","condition-D03_UTvsD00_UT",
           "condition-D24_UTvsD03_UT","condition-D09_UTvsD03_UT","condition-D24_UTvsD09_UT",
           "condition-D24_COvsD24_UT","condition-D24_IL17vsD24_UT",
           "condition-D24_IL13vsD24_UT","condition-D24_IFNvsD24_UT")

#comps_1 <- resultsNames(dds)[(grep("Int",resultsNames(dds),invert=T))]

sigcounts <- as.data.frame(setNames(replicate(2,numeric(0), simplify = F),letters[0:2]))
syms <- c()

fdrcutoff <- 0.1
rawpcutoff <- 0.1
lfc_cutoff <- 1


print("Running pairwise comparisons")
for (i in rev(comps)) {
	print(i)
	dir.create(paste("DE/DAVID/",i,sep=""))
	c <- unlist(strsplit(i,"-"))[1]
	t <- unlist(strsplit(i,"-"))[2]
	f <- unlist(strsplit(t,"vs"))[1]
	s <- unlist(strsplit(t,"vs"))[2]
	res <- results(ddsConds, contrast=c(c,f,s), alpha=fdrcutoff, cooksCutoff =.9, test="Wald")
	#res <- results(dds, contrast=list(comps_1[c(5,i)]), alpha=fdrcutoff, cooksCutoff =.9)
	fdr <- res[which(res$padj < fdrcutoff & abs(res$log2FoldChange) > 0),]
	rawp <- res[which(res$pvalue < rawpcutoff),]
	sigcounts <- rbind(sigcounts,c(nrow(fdr), nrow(rawp)))
	syms <- c(syms,rownames(fdr))
	
	
	
	res2 <- as.data.frame(res[,c(2,5,6)])
	res2$fc <- 2**res2[,1]
	res2 <- res2[,c(1,4,2,3)]
	
	colnames(res2) <- paste(i, c("log2FC", "FC","RawP", "FDR"), sep=".")
	allgenes <- cbind(as.data.frame(res2), allgenes)
	}




##########################################
#### DE tests for individual effects #####



comps_1 <- resultsNames(ddsModel)[(grep("Int",resultsNames(ddsModel),invert=T))]

#i <- comps_1[4]

print("Individual effects")
for (i in rev(comps_1)) {
	print(i)
	dir.create(paste("DE/DAVID/",i,sep=""))
	#c <- unlist(strsplit(i,"-"))[1]
	#t <- unlist(strsplit(i,"-"))[2]
	#f <- unlist(strsplit(t,"vs"))[1]
	#s <- unlist(strsplit(t,"vs"))[2]
	#res <- results(ddsI, contrast=c(c,f,s), alpha=fdrcutoff, cooksCutoff =.9, test="Wald")
	res <- results(ddsModel, name=i, alpha=fdrcutoff, cooksCutoff =.9)
	fdr <- res[which(res$padj < fdrcutoff & abs(res$log2FoldChange) > 0),]
	rawp <- res[which(res$pvalue < rawpcutoff),]
	sigcounts <- rbind(sigcounts,c(nrow(fdr), nrow(rawp)))
	syms <- c(syms,rownames(fdr))
	
	
	
	res2 <- as.data.frame(res[,c(2,5,6)])
	res2$fc <- 2**res2[,1]
	res2 <- res2[,c(1,4,2,3)]
	colnames(res2) <- paste(i, c("log2FC", "FC","RawP", "FDR"), sep=".")
	allgenes <- cbind(as.data.frame(res2), allgenes)
}

baseMean <- res$baseMean
baseMean <- as.matrix(baseMean)
row.names(baseMean) <- row.names(res)
colnames(baseMean) <- c("Mean_Normalized_Counts")
allgenes <- merge(allgenes, baseMean, by="row.names", all.y=T)


syms <- unique(syms)
#genes <- annot[syms,]$Gene
M <- stabilized[syms,]
#M <- M[,c(1:12,19:24,31:36)]
#M <-subset(M,,names)
M <- M - rowMeans(M, na.rm=T)


#sigcounts <- sigcounts[nrow(sigcounts):1,]
#syms <- unique(syms)
#genes <- annot[syms,]$Gene

#stabilizedBR <- removeBatchEffect(as.data.frame(stabilized), batch=c("A", "B", "B", "A", "B", "B", "A", "B", "B", "A", "B", "B", "A", "B", "B", "A", "B", "B", "A", "B", "B", "B", "B", "B"))

#M <- stabilizedBR[syms,]
#M <- M[,c(1:12,19:24,31:36)]
#M <-subset(M,,names)
#M <- M - rowMeans(M, na.rm=T)



rownames(allgenes) <- allgenes[,"Row.names"]
allgenes <- allgenes[,-c(1)]
#allgenes1 <- allgenes[,c(ncol(allgenes), ncol(allgenes)-1, 2:(ncol(allgenes)-2))]
#allgenes2 <- merge(annot, allgenes1, by="Ensembl_ID", all.y=T)
#allgenes3 <- merge(allgenes2, fpkmAve, by.x = "Ensembl_ID", by.y = "row.names", all.x=T)
norm_counts_ord <- norm_counts[order(rownames(norm_counts), decreasing=F),]

rownames(sigcounts) <- c(rev(comps),rev(comps_1))
colnames(sigcounts) <- c(paste("FDR < ",fdrcutoff), paste("RawP <", rawpcutoff))
#sigcounts <- sigcounts[order(rev(rownames(sigcounts))),]

dir.create(paste("DE/DEfiles/RawP",as.character(rawpcutoff), sep=""))
dir.create(paste("DE/DEfiles/FDR",as.character(fdrcutoff), sep=""))


colnames(norm_counts_ord) = gsub("D28","D24",colnames(norm_counts_ord))
colnames(norm_counts_ord) = gsub("_13_","_IL13_",colnames(norm_counts_ord))
colnames(norm_counts_ord) = gsub("_17_","_IL17_",colnames(norm_counts_ord))
colnames(norm_counts_ord) = gsub("D9","D09",colnames(norm_counts_ord))
colnames(norm_counts_ord) = gsub("D3","D03",colnames(norm_counts_ord))
colnames(norm_counts_ord) = gsub("D0","D00",colnames(norm_counts_ord))

colnames(all) = c("chr","start","stop","name")
allgenes = merge(all,  allgenes, by.x = "name",by.y='row.names')

for (i in c(comps,comps_1)) {
  a <- allgenes[,grep(i, colnames(allgenes))]
  b <- allgenes[,c(1:3)]
  j <- strsplit(i,"-")[[1]][2]
  n <- paste(strsplit(j,"vs")[[1]][1],"_", sep="")
  m <- paste(strsplit(j,"vs")[[1]][2],"_" , sep="")
  c <- norm_counts_ord[,grep(n,colnames(norm_counts_ord))]
  d <- norm_counts_ord[,grep(m,colnames(norm_counts_ord))]
  defile <- cbind(b,a,c,d)
  defile[,4] <- as.numeric(defile[,4])
  defile[,7] <- as.numeric(defile[,7])
  defile <- defile[order(defile[,4], decreasing=T),]
  pos <- grep("FDR", colnames(defile))
  defile <- defile[which(defile[,pos] < fdrcutoff & abs(defile[,4]) > 0),]
  write.table(defile, file=paste(paste("DE/DEfiles/FDR",as.character(fdrcutoff), sep=""), "/",i,".xls", sep=""), quote=F, row.names=F, sep="\t")}

for (i in c(comps,comps_1)) {
  a <- allgenes[,grep(i, colnames(allgenes))]
  b <- allgenes[,c(1:3)]
  j <- strsplit(i,"-")[[1]][2]
  n <- paste(strsplit(j,"vs")[[1]][1],"_", sep="")
  m <- paste(strsplit(j,"vs")[[1]][2],"_" , sep="")
  c <- norm_counts_ord[,grep(n,colnames(norm_counts_ord))]
  d <- norm_counts_ord[,grep(m,colnames(norm_counts_ord))]
  defile <- cbind(b,a,c,d)
  defile[,4] <- as.numeric(defile[,4])
  defile[,7] <- as.numeric(defile[,7])
  defile <- defile[order(defile[,4], decreasing=T),]
  pos <- grep("RawP", colnames(defile))
  defile <- defile[which(defile[,pos] < rawpcutoff),]
  write.table(defile, file=paste(paste("DE/DEfiles/RawP",as.character(rawpcutoff), sep=""), "/",i,".xls", sep=""), quote=F, row.names=F, sep="\t")} 

#rownames(sigcounts) <- names
#colnames(sigcounts) <- c(paste("FDR < ",fdrcutoff), paste("RawP <", rawpcutoff))
rownames(allgenes) = allgenes[,1]
allgenes2 <- merge(all, allgenes, by = "row.names")
rownames(allgenes2) <- allgenes2[,1]
colnames(allgenes2)[1:5] <- c("name","chr","start","stop","name")
allgenes2 <- allgenes2[,-c(1)] 

colnames(allgenes2)[4] <- "peak_ID"
allgenes2 <- allgenes2[,-c(5:20)]
last <- dim(allgenes2)[2]
allgenes2 <- allgenes2[,c(1:4,last,5:(last-1))]

select <- order(rv, decreasing=T)[seq_len(3000)]
v <- allgenes[select,][,c(1:3)]
write.table(v, file=paste(paste("DE/Top_3000_var_genes_from_PCA",project, sep = "."),"txt", sep="."), sep="\t", quote=F, row.names=F)
write.table(allgenes2, file=paste(paste("DE/All_genes",project, sep = "."),"txt", sep="."), sep="\t", quote=F, row.names=T, col.names = NA)
write.table(sigcounts, file=paste(paste("DE/DEcounts",project, sep = "."),"txt", sep="."), sep="\t", quote=F, row.names=T, col.names=NA)

save.image(file=paste(project,".RData",sep=""))

DE.heatmap(M, dendrogram="both", fact=conds[,c("Treatment","Time","Ind"), drop=F], lhei=c(1,.25,.25,.25,5,.75), height = 1200, width=1200, main=paste("DE genes - Union of all Pairwise Comparisons - FDR < ", fdrcutoff), file.prefix="DE/Heatmaps/DE_heatmap_nosymbols", breaks=.2, Symbols=NULL, min.col = "blue", mid.col = "white", max.col = "red")
DE.heatmap(M, dendrogram="both", fact=conds[,c("Treatment","Time"), drop=F], lhei=c(1,.25,.25,5,.5),  width=1200,  main=paste("DE genes - Union of all Comparisons - FDR < ", fdrcutoff), file.prefix="DE/Heatmaps/DE_heatmap_symbols", breaks=.2, Symbols=syms, min.col = "blue", mid.col = "white", max.col = "red")

M <- stabilizedBR[syms,]
#M <- M[,c(1:12,19:24,31:36)]
#M <-subset(M,,names)
M <- M - rowMeans(M, na.rm=T)

DE.heatmap(M, dendrogram="both", fact=conds[,c("Treatment","Time","Ind"), drop=F], lhei=c(1,.25,.25,.25,5,.75), height = 1200, width=1200, main=paste("DE genes - Union of all Pairwise Comparisons - FDR < ", fdrcutoff), file.prefix="DE/Heatmaps/DE_heatmap_nosymbols", breaks=.2, Symbols=NULL, min.col = "blue", mid.col = "white", max.col = "red")
DE.heatmap(M, dendrogram="both", fact=conds[,c("Treatment","Time"), drop=F], lhei=c(1,.25,.25,5,.5),  width=1200,  main=paste("DE genes - Union of all Comparisons - FDR < ", fdrcutoff), file.prefix="DE/Heatmaps/DE_heatmap_symbols", breaks=.2, Symbols=syms, min.col = "blue", mid.col = "white", max.col = "red")


###############################################
#### Generate heatmaps for all comparisons ####
###############################################

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}



sigcounts <- sigcounts[nrow(sigcounts):1,]
c = 1
for (i in grep(".FDR", colnames(allgenes2))) {
	#j <- colnames(allgene2s)[i]
	name <- gsub(".FDR","", colnames(allgenes2)[i])
	list <- rownames(allgenes2[which(allgenes2[,i] < fdrcutoff),])
	if (length(list) > 2){
	  x = (((length(list) - 30)/30) * 700) + 1000
      x <- max(1000, x)
      y = round((((length(list) - 30)/30) * 4) + 6)
      y <- max(6, y)
	  dir.create(paste("DE/Heatmaps/",name,sep=""))
	  if (as.list(strsplit(name, "-|_|vs"))[[1]][c(1)] == "condition") {
	  	group <- as.list(strsplit(name, "-|vs"))[[1]][c(2,3)]
	  	print(group)
	  	condsSelect <- conds[which(conds$condition %in% group),]
	  	M <- stabilizedBR[,colnames(stabilizedBR) %in% rownames(condsSelect)]
	  	M <- M[list,]
	  	M <- M - rowMeans(M, na.rm=T)
	  	#genes <- annot[list,]$Gene
	  	if (nrow(M) > 0) {
	  		DE.heatmap(M, dendrogram="both", fact=condsSelect[,c("Treatment","Time","Ind"), drop=F], lhei=c(1,.5,.5,.5,6,1), height = 1200, width=1200, main=paste("DE genes - from ", name, " - FDR < ", fdrcutoff), file.prefix=paste("DE/Heatmaps/",name,"/DE_heatmap_nosymbols.",name,".only",sep=""), breaks=.1, Symbols=NULL, min.col = "blue", mid.col = "white", max.col = "red")
	  		#DE.heatmap(M, dendrogram="both", fact=conds[,c("Treatment","Time","Ind"), drop=F], lhei=c(1,.5,.5,.5,y,1),  width=1200,  height=x, main=paste("DE genes - from ", name, " - FDR < ", fdrcutoff), file.prefix=paste("DE/Heatmaps/",name,"/DE_heatmap_symbols",name,".only.withGenes",sep=""), breaks=.2, Symbols=genes, min.col = "blue", mid.col = "white", max.col = "red")
	  	}
	  }
	  if (as.list(strsplit(name, "-|_|vs"))[[1]][c(1)] != "condition" && grepl("_vs_",name) == TRUE){
	  	cond <- as.list(strsplit(name, "_"))[[1]][c(1)]
	  	group <- as.list(strsplit(name, "_"))[[1]][c(2,4)]
	  	condsSelect <- conds[which(conds[,cond] %in% group),]
	  	M <- stabilizedBR[,colnames(stabilizedBR) %in% rownames(condsSelect)]
	  	M <- M[list,]
	  	M <- M - rowMeans(M, na.rm=T)
	  	#genes <- annot[list,]$Gene
	  	if (nrow(M) > 0) {
	  		DE.heatmap(M, dendrogram="both", fact=condsSelect[,c("Treatment","Time","Ind"), drop=F], lhei=c(1,.5,.5,.5,6,1), height = 1200, width=1200, main=paste("DE genes - from ", name, " - FDR < ", fdrcutoff), file.prefix=paste("DE/Heatmaps/",name,"/DE_heatmap_nosymbols.",name,".only",sep=""), breaks=.1, Symbols=NULL, min.col = "blue", mid.col = "white", max.col = "red")
	  		#DE.heatmap(M, dendrogram="both", fact=conds[,c("Treatment","Time","Ind"), drop=F], lhei=c(1,.5,.5,.5,y,1),  width=1200,  height=x,  main=paste("DE genes - from ", name, " - FDR < ", fdrcutoff), file.prefix=paste("DE/Heatmaps/",name,"/DE_heatmap_symbols",name,".only.withGenes",sep=""), breaks=.2, Symbols=genes, min.col = "blue", mid.col = "white", max.col = "red")
	  	}
	  }
	  M <- stabilizedBR[list,]
	  M <- M - rowMeans(M, na.rm=T)
	  #genes <- annot[list,]$Gene
	  if (nrow(M) > 0) {
	  	DE.heatmap(M, dendrogram="both", fact=conds[,c("Treatment","Time","Ind"), drop=F], lhei=c(1,.5,.5,.5,6,1), height = 1200, width=1200, main=paste("DE genes - from ", name, " - FDR < ", fdrcutoff), file.prefix=paste("DE/Heatmaps/",name,"/DE_heatmap_nosymbols.",name,".allSamples",sep=""), breaks=.1, Symbols=NULL, min.col = "blue", mid.col = "white", max.col = "red")
	  	#DE.heatmap(M, dendrogram="both", fact=conds[,c("Treatment","Time","Ind"), drop=F], lhei=c(1,.5,.5,.5,y,1),  width=1200,  height=x,  main=paste("DE genes - from ", name, " - FDR < ", fdrcutoff), file.prefix=paste("DE/Heatmaps/",name,"/DE_heatmap_symbols",name,".allSamples.withGenes",sep=""), breaks=.2, Symbols=genes, min.col = "blue", mid.col = "white", max.col = "red")
	  }
	}
	c = c + 1
}

colnames(allgenes2) <- gsub("D28","D24",colnames(allgenes2))

bigListTime <- c()
bigListTreat <- c()
bigListInd <- c()

comps <- gsub("-",".",comps)
indList <- c(comps_1[1:4])
timeList <- c(comps_1[5:7], comps[1:6])
treatmentList <- c(comps_1[8:11], comps[7:10])

timeList <- paste(timeList,".FDR",sep="")
treatmentList <- paste(treatmentList,".FDR",sep="")
indList <- paste(indList,".FDR",sep="")

for (i in timeList) {
  #j <- colnames(allgenes)[i]
  #name <- gsub(".FDR","", colnames(allgenes)[i])
  #j = i - 3
  j = gsub("FDR","log2FC", i)
  k = gsub("FDR","RawP", i)
  shortList <- allgenes2[,c(j,k,i,"Mean_Normalized_Counts")]
  shortList <- shortList[which(shortList[,4] > 20),]
  shortList <- shortList[,c(1:3)]
  shortList <- shortList[order(shortList[,2]),]
  shortHead <- head(shortList, n=50)
  shortHead$abs <- abs(shortHead[,1])
  shortHead <- shortHead[order(shortHead[,4], decreasing = T),]
  list <- rownames(head(shortHead, n=20))
  bigListTime <- c(bigListTime, list)
}

for (i in treatmentList) {
  #j <- colnames(allgenes)[i]
  #name <- gsub(".FDR","", colnames(allgenes)[i])
  j = gsub("FDR","log2FC", i)
  k = gsub("FDR","RawP", i)
  shortList <- allgenes2[,c(j,k,i,"Mean_Normalized_Counts")]
  shortList <- shortList[which(shortList[,4] > 20),]
  shortList <- shortList[,c(1:3)]
  #shortList <- allgenes[,c(j,k,i)]
  shortList <- shortList[order(shortList[,2]),]
  shortHead <- head(shortList, n=50)
  shortHead$abs <- abs(shortHead[,1])
  shortHead <- shortHead[order(shortHead[,4], decreasing = T),]
  list <- rownames(head(shortHead, n=20))
  bigListTreat <- c(bigListTreat, list)
}

#for (i in indList) {
#  #j <- colnames(allgenes)[i]
#  #name <- gsub(".FDR","", colnames(allgenes)[i])
#  j = gsub("FDR","log2FC", i)
#  k = gsub("FDR","RawP", i)
#  shortList <- allgenes2[,c(j,k,i,"Mean_Normalized_Counts")]
#  shortList <- shortList[which(shortList[,4] > 20),]
#  shortList <- shortList[,c(1:3)]
#  #shortList <- allgenes[,c(j,k,i)]
#  shortList <- shortList[order(shortList[,2]),]
#  shortHead <- head(shortList, n=50)
#  shortHead$abs <- abs(shortHead[,1])
#  shortHead <- shortHead[order(shortHead[,4], decreasing = T),]
#  list <- rownames(head(shortHead, n=20))
#  bigListInd <- c(bigListInd, list)
#}

bigListTime <- unique(bigListTime)
M <- stabilizedBR[bigListTime,]
condsTime <- conds[which(conds$Treatment == "UT"),]
TimeList <- rownames(condsTime)
M <- M[,which(colnames(M) %in% TimeList)]

#M <- M[,c(1:12,19:24,31:36)]
#M <-subset(M,,names)
M <- M - rowMeans(M, na.rm=T)

DE.heatmap(M, dendrogram="both", fact=condsTime[,c("Time","Ind"), drop=F], lhei=c(1,.5,.5,5,.75), height = 1200, width=1200, main=paste("DE genes - Union of all Pairwise Comparisons - FDR < ", fdrcutoff), file.prefix="DE/Heatmaps/DE_heatmap_nosymbols.Time", breaks=.02, Symbols=NULL, min.col = "blue", mid.col = "white", max.col = "red")

bigListTreat <- unique(bigListTreat)
M <- stabilizedBR[bigListTreat,]
condsTreat <- conds[which(conds$Time == "D24"),]
TreatmentList <- rownames(condsTreat)
M <- M[,which(colnames(M) %in% TreatmentList)]
#M <- M[,c(1:12,19:24,31:36)]
#M <-subset(M,,names)
M <- M - rowMeans(M, na.rm=T)

DE.heatmap(M, dendrogram="both", fact=condsTreat[,c("Treatment","Ind"), drop=F], lhei=c(1,.5,.5,5,.75), height = 1200, width=1200, main=paste("DE genes - Union of all Pairwise Comparisons - FDR < ", fdrcutoff), file.prefix="DE/Heatmaps/DE_heatmap_nosymbols.Treatment", breaks=.02, Symbols=NULL, min.col = "blue", mid.col = "white", max.col = "red")


bigListInd <- unique(bigListInd)
M <- stabilized[bigListInd,]
condsInd <- conds[which(conds$Time == "D28" & conds$Treatment == "UT"),]
IndList <- rownames(condsInd)
M <- M[,which(colnames(M) %in% IndList)]
#M <- M[,c(1:12,19:24,31:36)]
#M <-subset(M,,names)
M <- M - rowMeans(M, na.rm=T)

DE.heatmap(M, dendrogram="both", fact=condsInd[,c("Ind"), drop=F], lhei=c(1,.5,5,.75), height = 1200, width=1200, main=paste("DE genes - Union of all Pairwise Comparisons - FDR < ", fdrcutoff), file.prefix="DE/Heatmaps/DE_heatmap_nosymbols.Ind", breaks=.02, Symbols=NULL, min.col = "blue", mid.col = "white", max.col = "red")

bigListInd <- unique(bigListInd)
M <- stabilized[bigListInd,]

#M <- M[,c(1:12,19:24,31:36)]
#M <-subset(M,,names)
M <- M - rowMeans(M, na.rm=T)
DE.heatmap(M, dendrogram="both", fact=conds[,c("Treatment","Time","Ind"), drop=F], lhei=c(1,.5,.5,.5,5,.75), height = 1200, width=1200, main=paste("DE genes - Union of all Pairwise Comparisons - FDR < ", fdrcutoff), file.prefix="DE/Heatmaps/DE_heatmap_nosymbols.Ind_all", breaks=.02, Symbols=NULL, min.col = "blue", mid.col = "white", max.col = "red")



sigcountsTime <- as.data.frame(setNames(replicate(2,numeric(0), simplify = F),letters[0:2]))
symsUp <- c()
symsDown <- c()


timeList <- c(comps_1[5:7], comps[1:6])

timeList <- paste(timeList,".FDR",sep="")

for (i in timeList) {
  #j <- colnames(allgenes)[i]
  #name <- gsub(".FDR","", colnames(allgenes)[i])
  #j = i - 3
  list <- c()
  j = gsub("FDR","log2FC", i)
  k = gsub("FDR","RawP", i)
  shortList <- allgenes2[,c(j,k,i,"Mean_Normalized_Counts")]
  shortList <- shortList[which(shortList[,3] <= 0.01),]
  shortListUp <- shortList[which(shortList[,1] > 0),]
  shortListDown <- shortList[which(shortList[,1] < 0),]
  upList <- rownames(shortListUp)
  downList <- rownames(shortListDown)
  upCount <- length(upList)
  downCount <- length(downList)
  tmp <- c(upCount,downCount)
  symsUp <- c(symsUp, upList)
  symsDown <- c(symsDown, downList)
  sigcountsTime <- rbind(sigcountsTime, tmp)
  #  shortList$abs <- abs(shortList[,1])
}

symsUp <- unique(symsUp)
symsDown <- unique(symsDown)

allCount <- c(length(symsUp),length(symsDown))

sigcountsTime <- rbind(sigcountsTime, allCount)
rownames(sigcountsTime) <- c(rownames(sigcountsTime)[-c(10)],"Total_Unique")

write.table(sigcountsTime, file=paste(paste("DE/DEcounts.Time",project, sep = "."),"txt", sep="."), sep="\t", quote=F, row.names=T, col.names=NA)


bigListTime <- c()
bigListTimeAll <- c()


timeList <- c(comps_1[5:7], comps[1:6])

timeList <- paste(timeList,".FDR",sep="")

for (i in timeList) {
  #j <- colnames(allgenes)[i]
  #name <- gsub(".FDR","", colnames(allgenes)[i])
  #j = i - 3
  list <- c()
  j = gsub("FDR","log2FC", i)
  k = gsub("FDR","RawP", i)
  shortList <- allgenes2[,c(j,k,i,"Mean_Normalized_Counts")]
  shortList <- shortList[which(shortList[,4] > 10),]
  shortList <- shortList[which(shortList[,2] <= 0.1),]
  shortList0.01 <- shortList[which(shortList[,2] <= 0.01),]
  p0.01List <- rownames(shortList0.01)
  bigListTimeAll <- c(bigListTimeAll, p0.01List)
#  shortList$abs <- abs(shortList[,1])
  shortList <- shortList[order(shortList[,2]),]
  list <- rownames(head(shortList, n=50))
  bigListTime <- c(bigListTime, list)
  
}

bigListTime <- unique(bigListTime)
length(bigListTime)

bigListTimeAll <- unique(bigListTimeAll)
length(bigListTimeAll)



condsTimeAll <- conds[which(conds$Treatment == "UT"),]

rownames(condsTimeAll) <- gsub("_K27ac","",rownames(condsTimeAll))
rownames(condsTimeAll) <- gsub("rep","",rownames(condsTimeAll))
colnames(condsTimeAll) <- gsub("Ind","Subject", colnames(condsTimeAll))
condsTimeAll$Subject <- gsub("rep","",condsTimeAll$Subject)

colnames(Mtime) <- gsub("D28","D24",colnames(Mtime))
rownames(condsTimeAll) <- gsub("D28","D24",rownames(condsTimeAll))
colnames(stabilized) <- gsub("_K27ac","",colnames(stabilized))
colnames(stabilized) <- gsub("rep","",colnames(stabilized))
colnames(stabilized) <- gsub("D28","D24",colnames(stabilized))
condsTimeAll <- condsTimeAll[-c(17),]

stabilized_noT<- stabilized[,which(colnames(stabilized) %in% rownames(condsTimeAll))]
stabilizedBR_noT <- removeBatchEffect(stabilized_noT, batch = condsTimeAll$Subject)



M <- stabilizedBR_noT[bigListTimeAll,]
M <- M - rowMeans(M)

allgenes2_deTime <- allgenes2[(allgenes2$peak_ID %in% bigListTimeAll),]
allgenes2_deTime <- allgenes2_deTime[,c(1:73)]
write.table(allgenes2_deTime, file="DE/All_genes.Bonser_ChIP_02_2018.deTime.bed", sep="\t", quote=F, row.names = F,col.names = T)

Mtime <- M[,which(colnames(M) %in% rownames(condsTimeAll))]

TimeCols <- brewer.pal(4, "Reds")
IndCols <- brewer.pal(4,"Accent")

condsTimeAllSmall  <- as.data.frame(condsTimeAll[,c("Time","Subject")])
colnames(condsTimeAllSmall ) <- c("Time","Subject")
rownames(condsTimeAllSmall) <- rownames(condsTimeAll)


haTime = HeatmapAnnotation(df = condsTimeAllSmall,
                            col = list(Time = c("D00" = TimeCols[1],"D03" =TimeCols[2], 
                                                "D09" = TimeCols[3], "D24" = TimeCols[4]),
                                       Subject = c("HBE1243" = IndCols[1],"HBE1326" = IndCols[2],
                                                   "HBE1328" = IndCols[3], "HBE1333" = IndCols[4])),
                           annotation_legend_param = list(Time = list(title = "\nTime", title_gp = gpar(fontsize = 8), 
                                                                      labels_gp = gpar(fontsize = 6)),
                                                          Subject = list(title = "\nSubject",title_gp = gpar(fontsize = 8), 
                                                                      labels_gp = gpar(fontsize = 6))),
                                                          #peak_cluster = list(title = "Peak Cluster",title_gp = gpar(fontsize = 8), 
                                                          #               labels_gp = gpar(fontsize = 6))),
                           annotation_height = c(.5,.5),
                           show_legend = TRUE)


day_dist = dist(t(Mtime), method='euclidean')
hc_day = hclust(day_dist, method='complete')
dend_day <- as.dendrogram(hc_day)
#dend_day <- rotate(x = dend_cell,order = c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))

peak_dist = dist(Mtime, method='euclidean')
hc_peak = hclust(peak_dist, method='complete')
dend <- as.dendrogram(hc_peak)

hc_peak_clusters <- cutree(hc_peak, h=max(hc_peak$height/2.45))
hc_peak_clusters <- as.data.frame(hc_peak_clusters)

clusterColors <- gg_color_hue(length(unique(hc_peak_clusters$hc_peak_clusters)))
names(clusterColors) <-  unique(hc_peak_clusters$hc_peak_clusters)
colnames(hc_peak_clusters) <- "peak_cluster"

#pdf("dude.pdf",height = 18)
#print(plot(dend))
#dev.off()

ha_row = rowAnnotation(hc_peak_clusters,
                       col = list(peak_cluster = clusterColors), width = unit(.5, "cm"),
                       annotation_legend_param = list(peak_cluster = list(title = "Peak Cluster", title_gp = gpar(fontsize = 8), 
                                                                  labels_gp = gpar(fontsize = 6)))
                       )


htTime = Heatmap(Mtime, name = "Log2 Normalized \nEnrichment", 
                  top_annotation = haTime, 
                  col = colorRamp2(c(-.25, 0, .25), c("purple", "black", "orange")), 
                  cluster_columns = dend_day, cluster_rows = dend, show_heatmap_legend = TRUE,
                  show_column_names = F, show_row_names = F, column_dend_reorder = F,
                  row_dend_width = unit(2, "cm"),
                  heatmap_legend_param = list(title = "Log2 Normalized \nEnrichment", 
                                              title_gp = gpar(fontsize = 8), 
                                             labels_gp = gpar(fontsize = 6)))

png("tmp.png",width = 5, height = 5,units = "in",res = 300)
print(htTime + ha_row)
dev.off()


pdf("TimeSeries_DE.all.v4.rowAnnot.2.45.pdf",width = 5, height = 5)
print(htTime + ha_row)
dev.off()

png("TimeSeries_DE.all.v4.rowAnnot.2.45.png",width = 5, height = 5,units = "in",res = 300)
print(htTime + ha_row)
dev.off()

timeMatrix <- c()
hc_peak_clusters$dup <- hc_peak_clusters$peak_cluster
clusterMeans <- c()
for (i in unique(hc_peak_clusters$peak_cluster)) {
  peaksInCluster <- hc_peak_clusters[(hc_peak_clusters$peak_cluster == i),]
  stabilizedInCluster <- stabilizedBR_noT[which(rownames(stabilizedBR_noT) %in% rownames(peaksInCluster)),]
  if (nrow(peaksInCluster) > 1) {
    timeMeans <- stabilizedInCluster[,FALSE]
    for (j in c("D00","D03","D09","D24")){
      condsSmall <- conds[which(conds$Time == j & conds$Treatment == "UT"),]
      rownames(condsSmall) <- gsub("D28","D24",rownames(condsSmall))
      rownames(condsSmall) <- gsub("_K27ac","",rownames(condsSmall))
      rownames(condsSmall) <- gsub("rep","",rownames(condsSmall))
      colnames(stabilizedInCluster) <- gsub("_K27ac","",colnames(stabilizedInCluster))
      colnames(stabilizedInCluster) <- gsub("rep","",colnames(stabilizedInCluster))
      colnames(stabilizedInCluster) <- gsub("D28","D24",colnames(stabilizedInCluster))
      stabilizedInClusterSmall <- stabilizedInCluster[,which(colnames(stabilizedInCluster) %in% rownames(condsSmall))]
      timeMean <- rowMeans(stabilizedInClusterSmall)
      timeMeans <- cbind(timeMeans, as.data.frame(timeMean))
    }
    colnames(timeMeans) <- c("D00","D03","D09","D24")
    timeMeans <- timeMeans - rowMeans(timeMeans)
    clusterMean <- colMeans(timeMeans)
    print(clusterMean)
  } else {
    timeMeans <- c()
    for (j in c("D00","D03","D09","D24")){
      condsSmall <- conds[which(conds$Time == j & conds$Treatment == "UT"),]
      rownames(condsSmall) <- gsub("D28","D24",rownames(condsSmall))
      rownames(condsSmall) <- gsub("_K27ac","",rownames(condsSmall))
      rownames(condsSmall) <- gsub("rep","",rownames(condsSmall))
      names(stabilizedInCluster) <- gsub("_K27ac","",names(stabilizedInCluster))
      names(stabilizedInCluster) <- gsub("rep","",names(stabilizedInCluster))
      names(stabilizedInCluster) <- gsub("D28","D24",names(stabilizedInCluster))
      stabilizedInClusterSmall <- stabilizedInCluster[which(names(stabilizedInCluster) %in% rownames(condsSmall))]
      timeMean <- mean(stabilizedInClusterSmall)
      timeMeans <- c(timeMeans, timeMean)
    }
    timeMeans <- t(as.data.frame(timeMeans))
    colnames(timeMeans) <- c("D00","D03","D09","D24")
    timeMeans <- timeMeans - rowMeans(timeMeans)
    clusterMean <- timeMeans
    print(clusterMean)
  }
  clusterMeans <- rbind(clusterMeans, clusterMean)
}
rownames(clusterMeans) <- paste("cluster_",unique(hc_peak_clusters$peak_cluster),sep="")

htTime = Heatmap(clusterMeans, name = "Log2 Normalized \nEnrichment", 
                 col = colorRamp2(c(-.2, 0, .2), c("purple", "black", "orange")), 
                 cluster_columns = F, show_heatmap_legend = TRUE,
                 cluster_rows = T,
                 show_column_names = T, show_row_names = T, column_dend_reorder = F,
                 heatmap_legend_param = list(title = "Log2 Normalized \nEnrichment", 
                                             title_gp = gpar(fontsize = 8), 
                                             labels_gp = gpar(fontsize = 6)))
ChIP_clusterMeans <- clusterMeans

write.table(ChIP_clusterMeans, file = "DE/ChIP_clusterMeans_byDay.txt", row.names = T, col.names = NA, sep="\t", quote=F)

ChIP_clusterMeans <- read.table("DE/ChIP_clusterMeans_byDay.txt",header=T,sep="\t",row.names = 1)
scRNA_clusterMeans <- read.table("DE/scRNA_clusterMeans.txt",header=T,sep="\t",row.names = 1)
positiveMarkersPerCluster <- readRDS("../Analysis_with_old/DE/positiveMarkersPerCluster.rds")
deClusterObjectPos <- readRDS("../Analysis_with_old/DE/deCluster.pos.rds")
corM <- read.table("DE/scRNA_ChIP.correlation_Matrix.txt",header=T,sep = "\t",row.names = 1)



scRNA_dist = dist(scRNA_clusterMeans, method='euclidean')
hc_scRNA = hclust(scRNA_dist, method='complete')
dend_scRNA <- as.dendrogram(hc_scRNA)

ChIP_dist = dist(ChIP_clusterMeans, method='euclidean')
hc_ChIP = hclust(ChIP_dist, method='complete')
dend_ChIP <- as.dendrogram(hc_ChIP)

corM <- cor(t(ChIP_clusterMeans),t(scRNA_clusterMeans))

cor_ChIP_dist = dist(corM, method='euclidean')
hc_cor_ChIP = hclust(cor_ChIP_dist, method='complete')
dend_cor_ChIP <- as.dendrogram(hc_cor_ChIP)

cor_scRNA_dist = dist(t(corM), method='euclidean')
hc_cor_scRNA = hclust(cor_scRNA_dist, method='complete')
dend_cor_scRNA <- as.dendrogram(hc_cor_scRNA)


ht_scRNA = Heatmap(t(scRNA_clusterMeans), name = "Log2 Normalized \nEnrichment", 
                   col = colorRamp2(c(-.2, 0, .2), c("purple", "black", "orange")), 
                   cluster_columns = dend_cor_scRNA, show_heatmap_legend = TRUE,
                   cluster_rows = F, 
                   show_column_names = F, show_row_names = T,column_dend_reorder = F,
                   heatmap_legend_param = list(title = "Z-score", 
                                               title_gp = gpar(fontsize = 10), 
                                               labels_gp = gpar(fontsize = 8)))

ht_ChIP = Heatmap(ChIP_clusterMeans, name = "Log2 Normalized \nEnrichment", 
                  col = colorRamp2(c(-.2, 0, .2), c("purple", "black", "orange")), 
                  cluster_columns = F, show_heatmap_legend = F,
                  cluster_rows = dend_cor_ChIP,
                  show_column_names = T, show_row_names = F, column_dend_reorder = F,
                  heatmap_legend_param = list(title = "Z-score", 
                                              title_gp = gpar(fontsize = 10), 
                                              labels_gp = gpar(fontsize = 8)))

ht_cor = Heatmap(corM, name = "Log2 Normalized \nEnrichment", 
                 col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")), 
                 cluster_columns = dend_cor_scRNA, show_heatmap_legend = TRUE,
                 cluster_rows = dend_cor_ChIP, show_row_dend = F, show_column_dend = F,
                 show_column_names = T, show_row_names = T,column_dend_reorder = F,
                 heatmap_legend_param = list(title = "Correlation \nCoefficient", 
                                             title_gp = gpar(fontsize = 10), 
                                             labels_gp = gpar(fontsize = 8)))
png("scRNA_ChIP_cluster_correlation.png",height = 6, width = 6, res=300,units = "in")
print(ht_cor)
dev.off()

png("scRNA_cluster_meanByDay.png",height = 6, width = 6, res=300,units = "in")
print(ht_scRNA)
dev.off()

png("ChIP_cluster_meanByDay.png",height = 6, width = 6, res=300,units = "in")
print(ht_ChIP)
dev.off()
#positiveMarkersPerCluster
ChIP_cluster_genes <- c()
for (i in rownames(corM)){
  v <- corM[i,]
  v <- v[order(v,decreasing = T)]
  bestV <- v[1]
  cutV <- bestV*.9
  tmpList <- c()
  for (j in names(v)){
    if (v[j] > cutV){
      name <- paste(j,".markerList",sep="")
      ChIP_cluster_genes[[i]] <- c()
      tmpList <- c(tmpList,positiveMarkersPerCluster[[name]])
    }
  }
  tmpList <- unique(tmpList)
  ChIP_cluster_genes[[i]] <- c(tmpList)
}
#geneID_to_geneName <- read.table("../../Erle_Bonser_scRNA_12_12_2017_round2/geneID_to_geneName.txt",sep="\t", header = F)
#txdb <- makeTxDbFromGFF("../../Erle_Bonser_scRNA_12_12_2017_round2/genes.gtf",format="gtf")

NearestGene <- function(peakFile, annotFile,maxDistance){
  NearestGene_To_Peak <- c()
  for (i in 1:nrow(peakFile)){
    newLine  <- c(0)
    chr <- as.vector(peakFile[i,1])
    start <- as.vector(peakFile[i,2])
    stop <- as.vector(peakFile[i,3])
    name <- as.vector(peakFile[i,4])
    chrAnnot <- annotFile[(annotFile$Chromosome_Name == chr),]
    if (nrow(chrAnnot) >= 1 ){
      chrAnnot$TSS <- "NA"
      print(head(chrAnnot))
      chrAnnot$TSS[(chrAnnot$Strand == -1)] <- annotFile$End_Position
      chrAnnot$TSS[(chrAnnot$Strand ==  1)] <- chrAnnot$Start_Position
      chrAnnot$TSS <- as.numeric(chrAnnot$TSS)
      chrAnnot$Distance <- min(c(abs(chrAnnot$TSS - start),abs(chrAnnot$TSS - stop)))
      chrAnnot <- chrAnnot[order(chrAnnot$Distance),]
      if (chrAnnot$Distance <= maxDistance){
        best <- chrAnnot[1,c(4,5,6,1,2,7,15)]
        best <- as.vector(as.matrix(best))
      } else { best <- c(rep("NA",7))}
    }
    newLine <- c(chr,start,stop,name,best)
                 #best$Chromosome_Name,best$Start_Position,best$End_Position,
                 #best$Ensembl_ID,best$Gene,best$Strand,best$Distance)
    NearestGene_To_Peak <- rbind(NearestGene_To_Peak,newLine)
  }
  rownames(NearestGene_To_Peak) <- NearestGene_To_Peak[,4]
  NearestGene_To_Peak <- as.data.frame(NearestGene_To_Peak)
  colnames(NearestGene_To_Peak) <- c("peak_chr","peak_start","peak_stop","peak_ID","gene_chr","gene_start",
                                     "gene_stop","Ensembl_ID","Gene","Strand","Distance")
  return(NearestGene_To_Peak)
}

tadFile <- read.table("../../Erle_Bonser_scRNA_12_12_2017_round2/IMR90_Rao_2014-raw_TADs.txt",header=F)
annot <- read.delim("~/data/human_ens_GRCh38_83_10x_annot.extended.txt",sep="\t",header = T,row.names = 1)
tadFile$V1 <- gsub("chr","",tadFile$V1)
colnames(tadFile) <- c("chr","start","stop")
#corM
hc_peak_clusters$peak_cluster <- factor(hc_peak_clusters$peak_cluster)
peak_clusters <- hc_peak_clusters$peak_cluster
names(peak_clusters) <- rownames(hc_peak_clusters)

TSSes <- c()
print("getting TSSes")
annotTmp <- annot
for (l in rownames(annotTmp)){
  line <- annotTmp[l,]
  if (line$Strand == -1){
    TSSes <- c(TSSes, line$End_Position)
  }
  if (line$Strand == 1){
    TSSes <- c(TSSes, line$Start_Position)
  }
}
print("Got TSSes")
annotTmp$TSS <- TSSes

annot <- annotTmp

annotTmp <- NULL

#for testing
peakClusters <- peak_clusters
allPeaks <- all
annotFile <- annot
tads <- tadFile
tads <- rbind(tads,c("Y","1","2"))
clusterCorMatrix <- corM
maxDistance <- 1000000
corCutOff <- 0.5

NearestGeneFit <- function(peakClusters,allPeaks, annotFile,tads,clusterCorMatrix,maxDistance,corCutOff){
  BestGenePeakPairs <- c() # this will be a list of matrices, one for each chip cluster, with the peak-gene pairs
  annotTmp <- annotFile
  for (h in levels(peakClusters)){
    individualPeakCluster <- peakClusters[(peakClusters == h)]
    individualPeakNames <- names(individualPeakCluster)
    NearestGene_To_Peak <- c()
    ChIP_clust_name <- paste("cluster_",h,sep="")
    print(paste("starting ",ChIP_clust_name,sep=""))
    BestGenePeakPairs[[ChIP_clust_name]] <- c()
    allPairs <- c()
    for (l in rownames(annotTmp)){
      geneName <-  as.vector(annotTmp[l,"Gene"])
      geneStart <- annotTmp[l,"Start_Position"]
      geneStop <- annotTmp[l,"End_Position"]
      geneStrand <- annotTmp[l,"Strand"]
      geneCors <- c()
      bestCluster <- c()
      bestPair <- c()
      for (m in names(positiveMarkersPerCluster)){
        geneList <- positiveMarkersPerCluster[[m]]
        if (geneName %in% geneList) {
          geneCors <- c(geneCors,m)
        } else {geneCors <- c(geneCors)}
      }
      if (length(geneCors) > 1){
        geneCors <- gsub(".markerList","",geneCors)
        geneClusterCors <- clusterCorMatrix[ChIP_clust_name,geneCors]
        geneClusterCors <- geneClusterCors[order(geneClusterCors,decreasing = T)]
        bestCor <- as.numeric(geneClusterCors[1])
        bestCluster <- names(geneClusterCors[1])
        if (bestCor >= corCutOff){
          bestPair <- c(geneName,bestCluster,bestCor)
        } else {
          bestPair <- c(geneName,bestCluster,0)
        }
      }
      if (length(geneCors) == 1){
        geneCors <- gsub(".markerList","",geneCors)
        geneClusterCors <- clusterCorMatrix[ChIP_clust_name,geneCors]
        bestCor <- as.numeric(geneClusterCors[1])
        bestCluster <- geneCors
        if (bestCor >= corCutOff){
          bestPair <- c(geneName,bestCluster,bestCor)
        } else {
          bestPair <- c(geneName,bestCluster,0)
        }
      }
      if (length(geneCors) == 0) { 
        bestPair <- c(geneName,NA,0)
      }
      allPairs <- rbind(allPairs,bestPair)
    }
    print("Got gene correlations")
    colnames(allPairs) <- c("Gene","best_cluster","R")
    rownames(allPairs) <- allPairs[,1]
    allPairs <- as.data.frame(allPairs)
    allPairs$R  <- as.numeric(levels(allPairs$R))[allPairs$R]
    allPairs <- allPairs[order(allPairs$R, decreasing = T),]
    allPairs = allPairs[!duplicated(allPairs$Gene),]
    annotTmp2 <- merge(annotTmp, allPairs, by="Gene")
    rownames(annotTmp2) <- annotTmp2$Ensembl_ID
    print("starting peak matching")
    peakBests <- c()
    # for each peak finding the best gene inside by distance, correlation and tad structure #
    for (i in individualPeakNames){
      peak <- allPeaks[(allPeaks[,4] == i),]
      chr <- as.vector(as.matrix(peak[1]))
      start <- as.vector(as.matrix(peak[2]))
      stop <- as.vector(as.matrix(peak[3]))
      name <- as.vector(as.matrix(peak[4]))
      chrAnnot <- annotTmp2[(annotTmp2$Chromosome_Name == chr),]
      chrTad <- tads[(tads$chr == chr),]
      geneTadIdx <- c()
      tadStatus.mat <- c()
      p.dists.mat <- c()
      hasTad <- FALSE
      ### erroring here because this peak is hitting two tads, need to make matrix if this happens and take max
      if (nrow(chrAnnot) >= 1 ){
        chrAnnot$Distance <- rowMin(cbind(abs(chrAnnot$TSS - as.numeric(start)),abs(chrAnnot$TSS - as.numeric(stop))))
        chrAnnot <- chrAnnot[order(chrAnnot$Distance),]
        if (nrow(chrTad) >= 1) {
          for (j in 1:nrow(chrTad)){
            p.dists <- c()
            tadStatus <- c()
            if ( ((chrTad[j,2] - 2e4) <= start && (chrTad[j,3] + 2e4) >= start) || ((chrTad[j,2] - 2e4) <= stop && (chrTad[j,3] + 2e4) >= stop)){
              hasTad <- TRUE
              for (k in 1:nrow(chrAnnot)){
                if ( ((chrTad[j,2] - 2e4) <= chrAnnot[k,14] && (chrTad[j,3] + 2e4) >= chrAnnot[k,14])){
                  geneTadIdx <- c(geneTadIdx, chrAnnot[k,1])
                  tadStatus <- c(tadStatus,1)
                  p.dist <- exp(-1*chrAnnot[k,"Distance"]*(1e-6))
                  p.dists <- c(p.dists,p.dist)
                } else { 
                  tadStatus <- c(tadStatus,0)
                  p.dist <- exp(-1*chrAnnot[k,"Distance"]*(5e-6))
                  p.dists <- c(p.dists,p.dist)
                }
              }
              p.dists <- as.matrix(p.dists)
              p.dists.mat <- cbind(p.dists.mat,p.dists)
              tadStatus <- as.matrix(tadStatus)
              tadStatus.mat <- cbind(tadStatus.mat,tadStatus)
              p.dists.maxs <- rowMax(p.dists.mat)
            }
          }
        }
        if (hasTad == FALSE){
          p.dists <- c()
          tadStatus <- c()
          p.dist <- c()
          for (p in rownames(chrAnnot)){
            tadStatus <- c(tadStatus,0)
            p.dist <- exp(-1*chrAnnot[p,"Distance"]*(2.5e-6))
            p.dists <- c(p.dists,p.dist)
          }
          p.dists.mat <- as.matrix(p.dists)
          chrAnnot$p.dist <- p.dists.mat
          chrAnnot$inTad <- tadStatus
        }
        if (hasTad == TRUE){
          chrAnnot$p.dist <- p.dists.maxs
          chrAnnot$inTad <- tadStatus
        }
        chrAnnot$SortFactor <- as.numeric(as.character(chrAnnot$R)) * as.numeric(as.character(chrAnnot$p.dist))
        chrAnnot <- chrAnnot[order(chrAnnot$SortFactor, decreasing = T),]
        maxCor <- max(chrAnnot$SortFactor)
        chrAnnot <- chrAnnot[(chrAnnot$SortFactor >= (maxCor*.8)),]
        peakBestLines <- c()
        peakInfo <- as.vector(as.matrix(peak))
        if (nrow(chrAnnot) >= 1 ){
          for (o in rownames(chrAnnot)){
            geneInfo <- as.vector(as.matrix(chrAnnot[o,c("Start_Position","End_Position","Ensembl_ID","Gene","best_cluster","R","Distance","p.dist","SortFactor")]))
            peakBestLines <- rbind(peakBestLines,c(peakInfo,geneInfo))
          }
        } else {peakBestLines <- rbind(peakBestLines,c(peakInfo,rep(NA,9)))}
      }
      peakBests <- rbind(peakBests,peakBestLines)
    }
    peakBests <- as.data.frame(peakBests)
    colnames(peakBests) <- c("chr","peak_start","peak_stop","peak_name","gene_start","gene_stop","ensembl_id","gene","gene_cluster",
                             "cluster_cor","distance","p.dist","Score")
    BestGenePeakPairs[[ChIP_clust_name]] <- peakBests
  }
  return(BestGenePeakPairs)
}


#peakClusters <- NULL
#allPeaks <- NULL
#annotFile <- NULL
#tads <- NULL
#clusterCorMatrix <- NULL
#maxDistance <- NULL
#corCutOff <- NULL

#peakClusters <- peak_clusters
#allPeaks <- all
#annotFile <- annot
#tads <- tadFile
#tads <- rbind(tads,c("Y",20001,40001))
#tads$start <- as.numeric(tads$start)
#tads$stop <- as.numeric(tads$stop)
#clusterCorMatrix <- corM
#maxDistance <- 1000000
#corCutOff <- 0.5

x <- NearestGeneFit(peakClusters,allPeaks, annotFile,tads,clusterCorMatrix,maxDistance,corCutOff)

tmpMat <- merge(as.data.frame(tmpList), geneID_to_geneName, by.x="V1", by.y="V2")
colnames(tmpMat) <- c("Gene","Ensembl_ID")
tmpAnnot <- merge(tmpMat,annot, by="Ensembl_ID")
rownames(tmpAnnot) <- tmpAnnot$Ensembl_ID
tmpAnnot <- tmpAnnot[,-c(3)]
colnames(tmpAnnot)[2] <- "Gene"
k <- as.numeric(gsub("cluster_","",i))
tmpPeaks <- hc_peak_clusters[(hc_peak_clusters$peak_cluster == k),]
tmpBed <- all[(all$V4 %in% rownames(tmpPeaks)),]
annotBed <- tmpAnnot[,c(4,5,6,1,7)]
#hc_peak_clusters
#all

a <- NearestGene(tmpBed, tmpAnnot, 1e6)

M <- stabilizedBR[bigListTimeAll,]
M <- M - rowMeans(M)

colnames(M) <- gsub("_K27ac","",colnames(M))
colnames(M) <- gsub("rep","",colnames(M))

condsTimeAll <- conds[which(conds$Treatment == "UT"),]

rownames(condsTimeAll) <- gsub("_K27ac","",rownames(condsTimeAll))
rownames(condsTimeAll) <- gsub("rep","",rownames(condsTimeAll))
colnames(condsTimeAll) <- gsub("Ind","Subject", colnames(condsTimeAll))
condsTimeAll$Subject <- gsub("rep","",condsTimeAll$Subject)


Mtime <- M[,which(colnames(M) %in% rownames(condsTimeAll))]

TimeCols <- brewer.pal(4, "Reds")
IndCols <- brewer.pal(4,"Accent")

condsTimeAllSmall  <- as.data.frame(condsTimeAll[,c("Time")])
colnames(condsTimeAllSmall ) <- c("Time")
rownames(condsTimeAllSmall) <- rownames(condsTimeAll)

haTime = HeatmapAnnotation(df = condsTimeAllSmall,
                           col = list(Time = c("D00" = TimeCols[1],"D03" =TimeCols[2], 
                                               "D09" = TimeCols[3], "D24" = TimeCols[4])),
                           annotation_legend_param = list(Time = list(title = "\nTime", title_gp = gpar(fontsize = 8), 
                                                                      labels_gp = gpar(fontsize = 6))),
                           show_legend = TRUE)

sample_dist = dist(t(Mtime), method='euclidean')
hc_samples = hclust(sample_dist, method='complete')
sample_dend <- as.dendrogram(hc_samples)

gene_dist = dist(Mtime, method='euclidean')
hc_gene = hclust(gene_dist, method='complete')
dend <- as.dendrogram(hc_gene)
mycl <- cutree(hc_gene, h=max(hc_gene$height/1.5))
#dend <- rotate(x = dend,order = c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))

peakClusters <- as.data.frame(mycl)
colnames(peakClusters) <- c("peak_clusters")

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

peakClusterCols <- gg_color_hue(length(unique(peakClusters$peak_clusters)))

names(peakClusterCols) <- unique(peakClusters$peak_clusters)

haTimeRow = HeatmapAnnotation(df = peakClusters, 
                              col = list(peak_clusters = peakClusterCols),
                       which = "row",
                       annotation_legend_param = list(peak_clusters = list(title = "Peak Clusters", 
                                                                          title_gp = gpar(fontsize = 8), 
                                                                  labels_gp = gpar(fontsize = 6))),
                       show_legend = TRUE,gap = unit(.1,"mm"))

pdf("dude.pdf",height = 18)
print(plot(dend))
dev.off()

htTime = Heatmap(Mtime, name = "Log2 Normalized \nEnrichment", 
                 top_annotation = haTime, 
                 col = colorRamp2(c(-.3, 0, .3), c("purple", "black", "orange")), 
                 cluster_columns = sample_dend, show_heatmap_legend = TRUE,
                 cluster_rows = dend,
                 show_column_names = F, show_row_names = F, column_dend_reorder = F,
                 heatmap_legend_param = list(title = "Log2 Normalized \nEnrichment", 
                                             title_gp = gpar(fontsize = 8), 
                                             labels_gp = gpar(fontsize = 6)))

pdf("TimeSeries_DE.v3.all.noBR.pdf",width = 5, height = 5)
print(htTime)
dev.off()

haTimeRow + htTime
htTime + haTimeRow 


bigListTreat <- c()


treatmentList <- c(comps_1[8:11], comps[7:10])

treatmentList <- paste(treatmentList,".FDR",sep="")

for (i in treatmentList) {
  #j <- colnames(allgenes)[i]
  #name <- gsub(".FDR","", colnames(allgenes)[i])
  #j = i - 3
  list <- c()
  j = gsub("FDR","log2FC", i)
  k = gsub("FDR","RawP", i)
  shortList <- allgenes2[,c(j,k,i,"Mean_Normalized_Counts")]
  shortList <- shortList[which(shortList[,4] > 10),]
  shortList <- shortList[which(shortList[,2] <= 0.01),]
  #  shortList$abs <- abs(shortList[,1])
  shortList <- shortList[order(shortList[,2]),]
  list <- rownames(head(shortList, n=553))
  bigListTreat <- c(bigListTreat, list)
}

bigListTreat <- unique(bigListTreat)
length(bigListTreat)

M <- stabilizedBR[bigListTreat,]
M <- M - rowMeans(M)

colnames(M) <- gsub("_K27ac","",colnames(M))
colnames(M) <- gsub("rep","",colnames(M))

condsTreatAll <- conds[which(conds$Time == "D24"),]

rownames(condsTreatAll) <- gsub("_K27ac","",rownames(condsTreatAll))
rownames(condsTreatAll) <- gsub("rep","",rownames(condsTreatAll))
colnames(condsTreatAll) <- gsub("Ind","Subject", colnames(condsTreatAll))
condsTreatAll$Subject <- gsub("rep","",condsTreatAll$Subject)


Mtreat <- M[,which(colnames(M) %in% rownames(condsTreatAll))]

TreatCols <- brewer.pal(5, "Set1")
IndCols <- brewer.pal(4,"Accent")

condsTreatAllSmall  <- as.data.frame(condsTreatAll[,c("Treatment")])
colnames(condsTreatAllSmall ) <- c("Treatment")
rownames(condsTreatAllSmall) <- rownames(condsTreatAll)

haTreat = HeatmapAnnotation(df = condsTreatAllSmall,
                           col = list(Treatment = c("IFN" = TreatCols[1],"UT" =TreatCols[2], 
                                               "IL13" = TreatCols[3], "IL17" = TreatCols[5], "CO" = TreatCols[4])),
                           annotation_legend_param = list(Treatment = list(title = "\nTreatment", title_gp = gpar(fontsize = 8), 
                                                                      labels_gp = gpar(fontsize = 6))),
                           show_legend = TRUE)

gene_dist = dist(t(Mtreat), method='euclidean')
hc_genes = hclust(gene_dist, method='complete')
dend <- as.dendrogram(hc_genes)
#dend <- rotate(x = dend,order = c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))

pdf("dude.pdf",height = 18)
print(plot(dend))
dev.off()

htTreat = Heatmap(Mtreat, name = "Log2 Normalized \nEnrichment", 
                 top_annotation = haTreat, 
                 col = colorRamp2(c(-.2, 0, .2), c("purple", "black", "orange")), 
                 cluster_columns = dend, show_heatmap_legend = TRUE,
                 show_column_names = F, show_row_names = F, column_dend_reorder = F,
                 heatmap_legend_param = list(title = "Log2 Normalized \nEnrichment", 
                                             title_gp = gpar(fontsize = 8), 
                                             labels_gp = gpar(fontsize = 6)))

png("Treatment_DE.v4.png",width = 5, height = 5, res = 300, units = "in")
print(htTreat)
dev.off()






dir.create("DE")
dir.create("QC")


colnames(all) <- c("#chr","start","stop","ID")
for (i in grep(".FDR", colnames(allgenes2))) {
  name <- gsub(".FDR","", colnames(allgenes2)[i])
  j = i - 3
  res4 <- allgenes2[,c(j:i)]
  res4 <- as.data.frame(res4)
  d <- merge(all,res4,by="row.names",all.x=TRUE)[,-1]
  d <- d[order(d[,c(8)]),]
  write.table(d, file = paste("DE/beds/",name, ".DESeq2_results.bed", sep=""), row.names=FALSE, sep="\t", quote=FALSE)
  e <- d[(d[,8] <= 0.1),]
  e = e[!(is.na(e$ID)),]
  write.table(e, file = paste("DE/beds/",name, ".DESeq2_results.DEonly.bed", sep=""), row.names=FALSE, sep="\t", quote=FALSE)
  #colnames(d) <- NULL
  dNoNA <- d[complete.cases(d),]
  broad <- dNoNA[c(1,2,3,4)]
  broad$V5 <- -10*(log10(dNoNA[,c(8)]))
  broad$V6 <- "."
  broad <- broad[c(1,2,3,4,5,6,5,5,5)]
  rownames(broad) <- broad$V4
  ResDE <- d[d[,c(8)] <= 0.1,]
  ResUp <- ResDE[ResDE[,c(5)] > 0,]
  ResDown <- ResDE[ResDE[,c(5)] < 0,]
  ResUpList <- ResUp$ID
  ResDownList <- ResDown$ID
  broadUp <- broad[which(broad$ID %in% ResUpList),]
  broadDown <- broad[which(broad$ID %in% ResDownList),]
  write.table(broad, file = paste("DE/broadPeaks/", name, ".DESeq2_results.broadPeak", sep=""), row.names=FALSE, sep="\t", quote=FALSE)
  write.table(broadUp, file = paste("DE/broadPeaks/", name, ".DESeq2_results.Up.DE.broadPeak", sep=""), row.names=FALSE, sep="\t", quote=FALSE)
  write.table(broadDown, file = paste("DE/broadPeaks/", name, ".DESeq2_results.Down.DE.broadPeak", sep=""), row.names=FALSE, sep="\t", quote=FALSE)
}





stabilizedT <- t(stabilized)
stabilizedTAnnotated <- merge(stabilizedT, conds, by="row.names")
rownames(stabilizedTAnnotated) <- stabilizedTAnnotated[,1]
stabilizedTAnnotated <- stabilizedTAnnotated[,-c(1)]


MeltStablized <- melt(stabilizedTAnnotated)
meanStabilized <- aggregate(MeltStablized[,"value"], list(condition = MeltStablized[,"condition"], peak = MeltStablized[,"variable"]), mean)

stabilizedMeans = dcast( meanStabilized , peak~condition )
rownames(stabilizedMeans) <- stabilizedMeans$peak
stabilizedMeans <- stabilizedMeans[,-c(1)]

condsMean <- conds[which(conds$Ind == "repHBE1328"),]
rownames(condsMean) <- condsMean$condition

bigListTime <- unique(bigListTime)
M <- stabilizedMeans[bigListTime,]
condsTime <- condsMean[which(condsMean$Treatment == "UT"),]
condsTime <- condsTime[order(condsTime$Time),]
TimeList <- rownames(condsTime)
M <- M[,which(colnames(M) %in% TimeList)]
M <- M - rowMeans(M, na.rm=T)

DE.heatmap(as.matrix(M), dendrogram="both", fact=condsTime[,c("Time"), drop=F], lhei=c(1,.5,5,.75), height = 1200, width=1200, main=paste("DE genes - Union of all Pairwise Comparisons - FDR < ", fdrcutoff), file.prefix="DE/Heatmaps/DE_heatmap_nosymbols.Time.Mean", breaks=.02, Symbols=NULL, min.col = "blue", mid.col = "white", max.col = "red")

bigListTreat <- unique(bigListTreat)
M <- stabilizedMeans[bigListTreat,]
condsTreat <- condsMean[which(condsMean$Time == "D28"),]
condsTreat <- condsTreat[order(condsTreat$Treatment),]
TreatmentList <- rownames(condsTreat)
M <- M[,which(colnames(M) %in% TreatmentList)]
#M <- M[,c(1:12,19:24,31:36)]
#M <-subset(M,,names)
M <- M - rowMeans(M, na.rm=T)

DE.heatmap(as.matrix(M), dendrogram="both", fact=condsTreat[,c("Treatment"), drop=F], lhei=c(1,.5,5,.75), height = 1200, width=1200, main=paste("DE genes - Union of all Pairwise Comparisons - FDR < ", fdrcutoff), file.prefix="DE/Heatmaps/DE_heatmap_nosymbols.Treatment.Mean", breaks=.02, Symbols=NULL, min.col = "blue", mid.col = "white", max.col = "red")



pdf("NFKB1.box.pdf")
ggplot(data=stabilizedForBoxAnnotated) + geom_boxplot(aes(x=as.character(genotype),y=ENSG00000109320,fill=sex)) + labs(title="Box Plot of NFKB1", x="Genotype of rs28450894", y="Variance Stabilized Counts of NFKB1")
dev.off()


pdf("NFKB1.box.bdr.pdf")
ggplot(data=stabilizedForBoxAnnotated) + geom_boxplot(aes(x=catBDR,y=ENSG00000109320,fill=sex)) + labs(title="Box Plot of NFKB1", x="BDR catagory", y="Variance Stabilized Counts of NFKB1")
dev.off()



pdf("SLC39A8.box.pdf")
ggplot(data=stabilizedForBoxAnnotated) + geom_boxplot(aes(x=as.character(genotype),y=ENSG00000138821,fill=sex)) + labs(title="Box Plot of SLC39A81", x="Genotype of rs28450894", y="Variance Stabilized Counts of SLC39A81")
dev.off()


pdf("SLC39A8.box.bdr.pdf")
ggplot(data=stabilizedForBoxAnnotated) + geom_boxplot(aes(x=catBDR,y=ENSG00000138821,fill=sex)) + labs(title="Box Plot of SLC39A8", x="BDR catagory", y="Variance Stabilized Counts of SLC39A81")
dev.off()


bed <- allgenes2[,c(4:6,1,58:97)]

write.table(bed, file=paste(paste("Stabilized.BDR.bed",project, sep = "."),"txt", sep="."), sep="\t", quote=F, row.names=F)
write.table(conds, file="conds.bdr.qtl.txt",sep="\t", quote=F, row.names=T, col.names=NA)

pdf("NFKB1.box.pdf")
x1 <- stabilizedForBoxAnnotated$ENSG00000109320[stabilizedForBoxAnnotated$genotype==0]
x2 <- stabilizedForBoxAnnotated$ENSG00000109320[stabilizedForBoxAnnotated$genotype==1]
x3 <- stabilizedForBoxAnnotated$ENSG00000109320[stabilizedForBoxAnnotated$genotype==2]
vioplot(x1, x2, x3, names=c("0", "1", "2"), col="gold")
title("Box plot of NFKB1", xlab="Genotype of rs28450894", ylab="Variance Stabilized Counts")
dev.off()

pdf("NFKB1.pdf")
x1 <- stabilizedForBoxAnnotated$ENSG00000109320[stabilizedForBoxAnnotated$genotype==0]
x2 <- stabilizedForBoxAnnotated$ENSG00000109320[stabilizedForBoxAnnotated$genotype==1]
x3 <- stabilizedForBoxAnnotated$ENSG00000109320[stabilizedForBoxAnnotated$genotype==2]
vioplot(x1, x2, x3, names=c("0", "1", "2"), col="gold")
title(main="Violin Plot of NFKB1", xlab="Genotype of rs28450894", ylab="Variance Stabilized Counts")
dev.off()

pdf("SLC39A8.pdf")
x1 <- stabilizedForBoxAnnotated$ENSG00000138821[stabilizedForBoxAnnotated$genotype==0]
x2 <- stabilizedForBoxAnnotated$ENSG00000138821[stabilizedForBoxAnnotated$genotype==1]
x3 <- stabilizedForBoxAnnotated$ENSG00000138821[stabilizedForBoxAnnotated$genotype==2]
vioplot(x1, x2, x3, names=c("0", "1", "2"), col="gold")
title(main="Violin Plot of SLC39A8", xlab="Genotype of rs28450894", ylab="Variance Stabilized Counts")
dev.off()

pdf("SLC39A8.paired.pdf")
x1 <- stabilizedForBoxAnnotated$ENSG00000138821[stabilizedForBoxAnnotated$genotype==0 && stabilizedForBoxAnnotated$sex=="M"]
x2 <- stabilizedForBoxAnnotated$ENSG00000138821[stabilizedForBoxAnnotated$genotype==1 && stabilizedForBoxAnnotated$sex=="M"]
x3 <- stabilizedForBoxAnnotated$ENSG00000138821[stabilizedForBoxAnnotated$genotype==2 && stabilizedForBoxAnnotated$sex=="M"]
x4 <- stabilizedForBoxAnnotated$ENSG00000138821[stabilizedForBoxAnnotated$genotype==0 && stabilizedForBoxAnnotated$sex=="F"]
x5 <- stabilizedForBoxAnnotated$ENSG00000138821[stabilizedForBoxAnnotated$genotype==1 && stabilizedForBoxAnnotated$sex=="F"]
x6 <- stabilizedForBoxAnnotated$ENSG00000138821[stabilizedForBoxAnnotated$genotype==2 && stabilizedForBoxAnnotated$sex=="F"]
vioplot(x1, x4, x2, x5, x3, x6, names=c("0M", "0F", "1M", "1F", "2M", "2F"), col=c("gold","darkgreen"))
title(main="Violin Plot of SLC39A8", xlab="Genotype of rs28450894", ylab="Variance Stabilized Counts")
dev.off()

pdf("SLC39A8.box.pdf")
x1 <- stabilizedForBoxAnnotated$ENSG00000138821[stabilizedForBoxAnnotated$genotype==0]
x2 <- stabilizedForBoxAnnotated$ENSG00000138821[stabilizedForBoxAnnotated$genotype==1]
x3 <- stabilizedForBoxAnnotated$ENSG00000138821[stabilizedForBoxAnnotated$genotype==2]
boxplot(x1, x2, x3, names=c("0", "1", "2"), col="gold")
title(main="Box Plot of SLC39A8", xlab="Genotype of rs28450894", ylab="Variance Stabilized Counts")
dev.off()


list <- c("Criz_YvsCriz","Criz_YvsY","Criz_YvsDMSO","CrizvsY","CrizvsDMSO","YvsDMSO")
v <- combn(list, 2, simplify = F)
for (i in v) {
a <- allgenes2[,paste(i[1], ".log2FC", sep="")]
b <- allgenes2[,paste(i[2], ".log2FC", sep="")]
a[a > 10] <- 10
a[a < -10] <- -10
b[b > 10] <- 10
b[b < -10] <- -10
de1 <- allgenes2[,paste(i[1], ".FDR", sep="")]
de2 <- allgenes2[,paste(i[2], ".FDR", sep="")]

pdf(paste("DE/MMplots/",i[1], 'VS', i[2],".pdf", sep=""), 7,7)
plot(a,b, cex=.4, pch=19, xlab=paste("log2FC:",i[1]), ylab=paste("log2FC:",i[2]), main=paste('MMplot:',i[1], 'VS', i[2], sep=" "), cex.lab=1.2, xlim=c(-6,6), ylim=c(-6,6))
points(a[which(de1 < fdrcutoff)], b[which(de1 < fdrcutoff)], col="green", cex=.4, pch=19)
points(a[which(de2 < fdrcutoff)], b[which(de2 < fdrcutoff)], col="blue", cex=.4, pch=19)
points(a[which(de1 < fdrcutoff & de2 < fdrcutoff)], b[which(de1 < fdrcutoff & de2 < fdrcutoff)], col="red", cex=.4, pch=19)
abline(h=0)
abline(v=0)
abline(a=0,b=1)
legend(-6,6, c(paste("DE in", i[1], sep=" "), paste("DE in", i[2], sep=" "), "DE in both"), pch=c(19,19,19), cex=c(rep(.8,3)), col=c("green", "blue", "red"),bty = "n")
dev.off()
pdf(paste("DE/MMplots/",i[1], 'VS', i[2],".right.pdf", sep=""), 7,7)
plot(a,b, cex=.4, pch=19, xlab=paste("log2FC:",i[1]), ylab=paste("log2FC:",i[2]), main=paste('MMplot:',i[1], 'VS', i[2], sep=" "), cex.lab=1.2, xlim=c(-6,6), ylim=c(-6,6))
points(a[which(de1 < fdrcutoff)], b[which(de1 < fdrcutoff)], col="green", cex=.4, pch=19)
points(a[which(de2 < fdrcutoff)], b[which(de2 < fdrcutoff)], col="blue", cex=.4, pch=19)
points(a[which(de1 < fdrcutoff & de2 < fdrcutoff)], b[which(de1 < fdrcutoff & de2 < fdrcutoff)], col="red", cex=.4, pch=19)
abline(h=0)
abline(v=0)
abline(a=0,b=-1)
legend(0,6, c(paste("DE in", i[1], sep=" "), paste("DE in", i[2], sep=" "), "DE in both"), pch=c(19,19,19), cex=c(rep(.8,3)), col=c("green", "blue", "red"),bty = "n")
dev.off()
}


## Heatmaps for DE genes pairwise comparisons  
DE.heatmap(M, dendrogram="both", fact=conds[,c("sex","binAgeS","genotypeL","catBDR"), drop=F], height = 1200, width=1200, lhei=c(1,.25,.25,.25,.25,6,.5), main=paste("DE genes - Union of all Pairwise and Main Effect Comparisons - FDR < ", fdrcutoff), file.prefix="DE/DE_heatmap_nosymbols.all_DE", breaks=.2, Symbols=NULL)
DE.heatmap(M, dendrogram="both", fact=conds[,c("sex","binAgeS","genotypeL","catBDR"), drop=F], width=1200, lhei=c(.1,.25,.25,.25,.25,4,.5), main=paste("DE genes - Union of all Pairwise and Main Effect Comparisons - FDR < ", fdrcutoff), file.prefix="DE/DE_heatmap_symbols.all_DE", breaks=.2, Symbols=genes)
#.interaction_plus_Criz_YvsCriz

cr = cor(M, method='spearman')
centered_data = t(scale(t(M), scale=F))
#centered_data = scale(data, scale=T) # center rows, mean substracted
gene_dist = dist(centered_data, method='euclidean')
hc_genes = hclust(gene_dist, method='complete')
hc_samples = hclust(as.dist(1-cr), method="complete") # cluster conditions
myheatcol = redgreen(255)
gene_partition_assignments <- cutree(as.hclust(hc_genes), k=6);
partition_colors = rainbow(length(unique(gene_partition_assignments)), s=.5, v=.75, start=0.1, end=0.95)
gene_colors = partition_colors[gene_partition_assignments]


ha = HeatmapAnnotation(conds, col = list(cell = c("A" = "red", "B" = "blue"), genotype = c("C" = "green", "D" = "purple")))

list(cell = c("A" = "red", "B" = "blue"), genotype = c("C" = "green", "D" = "purple"))

Heatmap(mat, name = "foo", split = paste0("cluster", gene_partition_assignments))




col_assignments <- targets$condition[c(1:12,19:24,31:36)]
col_colors = rainbow(length(unique(col_assignments)), s=.5, v=.75, start=0.1, end=0.95)
sample_colors = col_colors[col_assignments]
quantile.range <- quantile(centered_data, probs = seq(0, 1, 0.001))
palette.breaks <- seq(quantile.range["1.9%"], quantile.range["97.9%"], 0.01)
color.palette  <- colorRampPalette(c("blue", "white", "red"))(length(palette.breaks) - 1)
png("DE/DE_heatmap_nosymbols.blue_to_red.hcClust.interaction.png", height=900, width=900)
heatmap.2(centered_data, ColSideColors=sample_colors, symkey=T, symm=T, symbreaks=T, dendrogram='both', Rowv=reorder(as.dendrogram(hc_genes), 10:1), Colv=as.dendrogram(hc_samples), col = color.palette, breaks= palette.breaks, scale="none", density.info="none", trace="none", key=TRUE, keysize=1.2, cexCol=1, margins=c(12,5), labRow=FALSE)
dev.off()



## Heatmap for DE genes (rawp) from Interaction 

syms <- rownames(Interaction[which(Interaction$pval < .01),])
genes <- annot[syms,]$Gene
M <- stabilized[syms,]
M <- M - rowMeans(M, na.rm=T)

DE.heatmap(M, dendrogram="both", fact=targets[,c("condition"), drop=F], width=1200, lhei=c(1,.5,6,.5), main=paste("DE genes - Interaction - RawP < ", rawpcutoff), file.prefix="DE/DE_Interaction_symbols", breaks=.2, Symbols=genes)

