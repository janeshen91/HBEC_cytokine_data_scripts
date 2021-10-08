library(ggrepel)
library("cowplot")
library("stringr")
library(ComplexHeatmap)
library(circlize)
library(biomaRt)
library(dplyr)
library(DESeq2)
library(genefilter)
library("devtools")
library("BiocParallel")
library(Rtsne)
#register(MulticoreParam(4))
library(R.utils)
sourceDirectory("~/data/BatchArray-source/src")
source("~/data/BatchArray-source/load.R")
source("~/tools/R.RNA-seq/DAVID.R")
source("~/tools/R.RNA-seq/Pathway.analysis.R")
library(gplots)
library(ggplot2)
library(limma)
project <- c("Erle_Bonser_0529")
dir <- paste(paste("~/data",project, sep = "/"),"Analysis", sep = "/")
dir.create(dir)
setwd(dir)
library(gplots)
library(ggplot2)
library(limma)
library("calibrate")
library("rJava")
library("RDAVIDWebService")
library("pathview")
library(gage)
library(gageData)


dir.create("DE")
dir.create("QC")
dir.create("DE/MAplots/")
dir.create("DE/VolcanoPlots/")
dir.create("DE/MMplots/")
dir.create("DE/Heatmaps/")
dir.create("DE/DEfiles/")

file = "../counts.STAR.HUMAN.txt"

data <- read.table(file, header=T, row.names=1)

RunDavid <- TRUE
if (RunDavid == TRUE){
  dir.create("DE/DAVID/")
}
RunPathway <- TRUE
if (RunPathway == TRUE) {
  dir.create("DE/pathways/")
}

if (file  == "../counts.STAR.HUMAN.txt"){
  annot <- read.table("~/data/human_ens_GRCh38_annot.extended.v2.txt", sep="\t", quote="", header=T, row.names=1, stringsAsFactors=F, fill=T)
  data(go.sets.hs)
  data(go.subs.hs)
  gobpsets = go.sets.hs[go.subs.hs$BP]
  gomfsets = go.sets.hs[go.subs.hs$MF]
  data(kegg.sets.hs)
  data(sigmet.idx.hs)
  kegg.sets = kegg.sets.hs[sigmet.idx.hs]
  spec = "hsa"
  species = "hg38"
}
if (file  == "../counts.STAR.HUMANHG37.txt"){
  annot <- read.table("~/data/human_ens_GRCh37_annot.extended.txt", sep="\t", quote="", header=T, row.names=1, stringsAsFactors=F, fill=T)
  spec = "hsa"
  species = "hg19"
}
if (file  == "../counts.STAR.MOUSE.txt"){
  annot <- read.table("~/data/mouse_ens_GRCm38_annot.extended.txt", sep="\t", quote="", header=T, row.names=1, stringsAsFactors=F, fill=T)
  data(go.sets.mm)
  data(go.subs.mm)
  gobpsets = go.sets.mm[go.subs.mm$BP]
  gomfsets = go.sets.mm[go.subs.mm$MF]
  data(kegg.sets.mm)
  data(sigmet.idx.mm)
  kegg.sets = kegg.sets.mm[sigmet.idx.mm]
  spec = "mmu"
  species = "mm10"
}
if (file  == "../counts.STAR.MOUSEMM9.txt" | file == "../counts.scRNA.STARMOUSEMM9.txt"){
	annot <- read.table("~/data/mouse_ens_GRCm37.annot.extended.txt", sep="\t", quote="", header=T, row.names=1, stringsAsFactors=F, fill=T)
	spec = "mmu"
	species = "mm9"
}
if (file  == "../counts.STAR.ZEBRAFISH.txt"){
	annot <- read.table("~/data/zfish_ens_GRCz10_annot.extended.txt", sep="\t", quote="", header=T, row.names=1, stringsAsFactors=F, fill=T)
	RunPathway <- FALSE
}


#rownames(annot) <- annot$Ensembl_ID


targets <- as.data.frame(as.matrix(colnames(data)))

conds <- cbind(targets, do.call(rbind, strsplit(as.character(targets$V1),'_')))
rownames(conds) <- conds[,c(1)]
conds <- conds[,-c(1)]


condition <- c(rep("CoT",6),rep("IFNa",6),rep("IFNg",6),rep("IL13",6),rep("IL17",6),rep("UT",6))
targets <- as.data.frame(condition)
targets$IL13 <- c(rep("Y",6), rep("N",12), rep("Y",6), rep("N",12))
targets$IL17 <- c(rep("N",24), rep("Y",6), rep("N",6))
targets$IFNa <- c(rep("Y",12), rep("N",24))
targets$IFNg <- c(rep("N",12), rep("Y",6), rep("N",18))
targets$ind <- c(rep(c("A","B","C","D","E","F"),6))
rownames(targets) <- colnames(data)

conds = targets

#conds[,"genotype"] <- factor(conds[,"genotype"], levels=rev(levels(conds[,"genotype"])))
#conds$class <- as.factor(conds$class)
#conds[,"class"] <- factor(conds[,"class"], levels=rev(levels(conds[,"class"])))
#
#conds$condition <- as.factor(paste(conds$genotype, conds$class, sep="_"))
#
#
#conds[,"condition"] <- factor(conds[,"condition"], levels=levels(conds[,"condition"])[c(3,4,1,2)])
#conds$condition <- paste(conds$Location, conds$Tissue, sep="_")

write.table(conds, file=paste("conditions.", project, ".txt", sep=""), sep="\t", quote=F, row.names=T, col.names=NA)


ddsModel <- DESeqDataSetFromMatrix(countData = data,
							  colData = conds,
							  design = ~ ind + IL13 + IL17 + IFNa + IFNg + IL13*IFNa)


ddsConds <- DESeqDataSetFromMatrix(countData = data,
							  colData = conds,
							  design = ~  ind + condition)


ddsModel <- DESeq(ddsModel)
ddsConds <- DESeq(ddsConds)

norm_counts<- round(DESeq2::counts(ddsModel, normalized=TRUE),1)
allgenes <- norm_counts
stabilizedBlind <- assay(varianceStabilizingTransformation(ddsModel, blind=TRUE))
stabilized <- assay(varianceStabilizingTransformation(ddsModel, blind=FALSE))
stabilizedBR <- removeBatchEffect(stabilized, batch = conds$ind)

allgenes <- merge(allgenes, stabilized, by="row.names", suffixes=c(".norm",".stabilized"))
row.names(allgenes) <- allgenes$Row.names
allgenes <- allgenes[,-c(1)]


rv <- rowVars(stabilizedBlind)
select <- order(rv, decreasing=T)[seq_len(3000)]
pca <- prcomp(t(stabilizedBlind[select,]))

list = c(1,2,3,4,5)
for (i in 2:length(list)) {
	i
}
for (i in 2:length(list)) {                 
	pcaobject <- as.data.frame(pca$x[,c(1,i)])
	pcaobject$condition <- conds$condition
	pcaobject$subject <- conds$ind

	comps <- colnames(pcaobject)[1:2]
	labels <- as.vector(pca$sdev^2/sum(pca$sdev^2))
	xlabel <- round(labels[1]*100,1)
	ylabel <- round(labels[i]*100,1)
	#print(ylabel)
	#print(paste("PC",i,sep=""))
	ggplot2::ggplot(pcaobject, aes_string(x = comps[1], y = comps[2])) +
    	geom_point(aes(color=condition, shape=subject), size=6) +
    	xlab(paste0("PC1: ", xlabel, "% variance")) +
		ylab(paste0("PC",i,": ", ylabel, "% variance")) +
    	theme_bw()
	ggsave(paste("QC/PCAplot.PC1vPC",i,".pdf",sep=""), height=7, width=7, units=c("in"))
	ggplot2::ggplot(pcaobject, aes_string(x = comps[1], y = comps[2])) +
	  geom_point(aes(color=condition, shape=subject), size=6) +
	  xlab(paste0("PC1: ", xlabel, "% variance")) +
	  ylab(paste0("PC",i,": ", ylabel, "% variance")) +
	  theme_classic()
	ggsave(paste("QC/PCAplot.PC1vPC",i,".classic.pdf",sep=""), height=7, width=7, units=c("in"))}
for (i in 3:length(list)) {                 
	pcaobject <- as.data.frame(pca$x[,c(2,i)])
	pcaobject$condition <- conds$condition
	pcaobject$subject <- conds$ind

	comps <- colnames(pcaobject)[1:2]
	labels <- as.vector(pca$sdev^2/sum(pca$sdev^2))
	xlabel <- round(labels[2]*100,1)
	ylabel <- round(labels[i]*100,1)
	#print(ylabel)
	#print(paste("PC",i,sep=""))
	ggplot2::ggplot(pcaobject, aes_string(x = comps[1], y = comps[2])) +
	    geom_point(aes(color=condition, shape=subject), size=6) +
    	xlab(paste0("PC2: ", xlabel, "% variance")) +
		ylab(paste0("PC",i,": ", ylabel, "% variance")) +
    	theme_bw()
	ggsave(paste("QC/PCAplot.PC2vPC",i,".pdf",sep=""), height=7, width=7, units=c("in"))}
for (i in 4:length(list)) {                 
	pcaobject <- as.data.frame(pca$x[,c(3,i)])
	#print(i)
	pcaobject$condition <- conds$condition
	pcaobject$subject <- conds$ind

	comps <- colnames(pcaobject)[1:2]
	labels <- as.vector(pca$sdev^2/sum(pca$sdev^2))
	xlabel <- round(labels[3]*100,1)
	ylabel <- round(labels[i]*100,1)
	#print(ylabel)
	#print(paste("PC",i,sep=""))
	ggplot2::ggplot(pcaobject, aes_string(x = comps[1], y = comps[2])) +
	    geom_point(aes(color=condition, shape=subject), size=6) +
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
	df_to_plot$subject <- conds$ind
	#df_to_plot$Ethnicity <- conds$Ethnicity

	ggplot2::ggplot(df_to_plot, aes_string(x = comps[1], y = comps[2])) +
    	geom_point(aes(color=condition, shape=subject), size=6) +
    	xlab("Dimension 1") +
    	ylab("Dimension 2") +
    	theme_bw()
	ggsave(paste("QC/Tsne.plot.", i, ".pdf", sep=""), height=7, width=7, units=c("in"))
}

conds$subject = conds$ind

select <- order(rv, decreasing=T)[seq_len(1000)]
M <- stabilizedBlind[select,]
M <- M - rowMeans(M, na.rm=T)
#pairs.agilent(stabilizedBlind, filename="QC/pairs_log2normRC.png", main="Pairwise comparison of log2 normalized read counts", cex=2)


pdf("QC/boxplot.pdf")
par(mar=c(7,5,5,5))
boxplot(stabilizedBlind, las=2, col="blue", main="Log2 Normalized Read Count by Sample", pch=19, cex=.5)
dev.off()


DE.heatmap(M, dendrogram="column", height=1024, fact=conds[,c("condition","subject"), drop=F], lhei=c(1,.5,.5,6,1), main="Heatmap of log2 FC rel. global average for top 1000 genes by variance", file.prefix="QC/QC_heatmap", breaks=.5, margins=c(12,12,6,6))

select <- order(rv, decreasing=T)[seq_len(1000)]
M <- stabilizedBR[select,]
M <- M - rowMeans(M, na.rm=T)

DE.heatmap(M, dendrogram="column", height=1024, fact=conds[,c("condition","subject"), drop=F], lhei=c(1,.5,.5,6,1), main="Heatmap of log2 FC rel. global average for top 1000 genes by variance", file.prefix="QC/QC_heatmap.br", breaks=.5, margins=c(12,12,6,6))


rv <- rowVars(stabilizedBR)
select <- order(rv, decreasing=T)[seq_len(3000)]
pca <- prcomp(t(stabilizedBR[select,]))

list = c(1,2,3,4,5)
for (i in 2:length(list)) {
  i
}
for (i in 2:length(list)) {                 
  pcaobject <- as.data.frame(pca$x[,c(1,i)])
  pcaobject$condition <- conds$condition
  pcaobject$subject <- conds$ind
  
  comps <- colnames(pcaobject)[1:2]
  labels <- as.vector(pca$sdev^2/sum(pca$sdev^2))
  xlabel <- round(labels[1]*100,1)
  ylabel <- round(labels[i]*100,1)
  #print(ylabel)
  #print(paste("PC",i,sep=""))
  ggplot2::ggplot(pcaobject, aes_string(x = comps[1], y = comps[2])) +
    geom_point(aes(color=condition, shape=subject), size=6) +
    xlab(paste0("PC1: ", xlabel, "% variance")) +
    ylab(paste0("PC",i,": ", ylabel, "% variance")) +
    theme_bw()
  ggsave(paste("QC/PCAplot.br.PC1vPC",i,".pdf",sep=""), height=7, width=7, units=c("in"))}
for (i in 3:length(list)) {                 
  pcaobject <- as.data.frame(pca$x[,c(2,i)])
  pcaobject$condition <- conds$condition
  pcaobject$subject <- conds$ind
  
  comps <- colnames(pcaobject)[1:2]
  labels <- as.vector(pca$sdev^2/sum(pca$sdev^2))
  xlabel <- round(labels[2]*100,1)
  ylabel <- round(labels[i]*100,1)
  #print(ylabel)
  #print(paste("PC",i,sep=""))
  ggplot2::ggplot(pcaobject, aes_string(x = comps[1], y = comps[2])) +
    geom_point(aes(color=condition, shape=subject), size=6) +
    xlab(paste0("PC2: ", xlabel, "% variance")) +
    ylab(paste0("PC",i,": ", ylabel, "% variance")) +
    theme_bw()
  ggsave(paste("QC/PCAplot.br.PC2vPC",i,".pdf",sep=""), height=7, width=7, units=c("in"))}
for (i in 4:length(list)) {                 
  pcaobject <- as.data.frame(pca$x[,c(3,i)])
  #print(i)
  pcaobject$condition <- conds$condition
  pcaobject$subject <- conds$ind
  
  comps <- colnames(pcaobject)[1:2]
  labels <- as.vector(pca$sdev^2/sum(pca$sdev^2))
  xlabel <- round(labels[3]*100,1)
  ylabel <- round(labels[i]*100,1)
  #print(ylabel)
  #print(paste("PC",i,sep=""))
  ggplot2::ggplot(pcaobject, aes_string(x = comps[1], y = comps[2])) +
    geom_point(aes(color=condition, shape=subject), size=6) +
    xlab(paste0("PC3: ", xlabel, "% variance")) +
    ylab(paste0("PC",i,": ", ylabel, "% variance")) +
    theme_bw()
  ggsave(paste("QC/PCAplot.br.PC3vPC",i,".pdf",sep=""), height=7, width=7, units=c("in"))}


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
  df_to_plot$subject <- conds$ind
  #df_to_plot$Ethnicity <- conds$Ethnicity
  
  ggplot2::ggplot(df_to_plot, aes_string(x = comps[1], y = comps[2])) +
    geom_point(aes(color=condition, shape=subject), size=6) +
    xlab("Dimension 1") +
    ylab("Dimension 2") +
    theme_bw()
  ggsave(paste("QC/Tsne.plot.br.", i, ".pdf", sep=""), height=7, width=7, units=c("in"))
}


#############################
#### Multifactorial design ####
#############################


combined <- t(as.data.frame(combn(unique(conds$condition), 2, simplify = F)))
comps <- paste("condition-",paste(combined[,c(1)], combined[,c(2)], sep="vs"), sep="")
comps <- comps[order(comps)]
#comps_1 <- resultsNames(dds)[(grep("Int",resultsNames(dds),invert=T))]
comps = c("condition-CoTvsUT","condition-IL13vsUT","condition-IL17vsUT","condition-IFNavsUT","condition-IFNgvsUT")
sigcounts <- as.data.frame(setNames(replicate(2,numeric(0), simplify = F),letters[0:2]))
syms <- c()

fdrcutoff <- 0.1
rawpcutoff <- 0.01
lfc_cutoff <- 1


print("Running pairwise comparisons")
for (i in rev(comps)) {
	print(i)
	c <- unlist(strsplit(i,"-"))[1]
	t <- unlist(strsplit(i,"-"))[2]
	f <- unlist(strsplit(t,"vs"))[1]
	s <- unlist(strsplit(t,"vs"))[2]
	res <- results(ddsConds, contrast=c(c,f,s), alpha=fdrcutoff, cooksCutoff =1, test="Wald")
	#res <- results(dds, contrast=list(comps_1[c(5,i)]), alpha=fdrcutoff, cooksCutoff =.9)
	fdr <- res[which(res$padj < fdrcutoff & abs(res$log2FoldChange) > 0),]
	rawp <- res[which(res$pvalue < rawpcutoff),]
	sigcounts <- rbind(sigcounts,c(nrow(fdr), nrow(rawp)))
	syms <- c(syms,rownames(fdr))
	
	
	
	res2 <- as.data.frame(res[,c(2,5,6)])
	res2$fc <- 2**res2[,1]
	res2 <- res2[,c(1,4,2,3)]
	
	d <- merge(annot, res2, by='row.names')

	
	if (RunDavid == TRUE) {
	  DAVID.analysis(res2, fdrFilter = fdrcutoff, lfcFilter = lfc_cutoff, name = i)
	}
	
  if (RunPathway == TRUE) {
	  dirPathway <- paste("DE/pathways/",i,sep="")
	  dir.create(dirPathway)
	  dirCd <- paste(dir, dirPathway, sep="/")
	  setwd(dirCd)
	  
	  Pathway.analysis(d, species, i)
	  
	  setwd(dir)
  }
#
#
	

	
	d$padjf <- d$padj + 1e-300
	sub <- subset(d, log2FoldChange>1)
	up <- sub[order(sub$padj),][5,18]
	if ( up <= 1e-300 & !is.na(up)){
		up = 0 }
	sub <- subset(d, log2FoldChange< -1)
	down <- sub[order(sub$padj),][5,18]
	if ( down <= 1e-300 & !is.na(down)){
		down = 0 }
	if ( is.na(up)){
		up = 0.1
	}
	if ( is.na(down)){
		down = 0.1
	}

	
	upFCcut <- 2
	downFCcut <- -2
	
	
	
	
	d$color[d$log2FoldChange < upFCcut & d$log2FoldChange > downFCcut & d$padjf > fdrcutoff ]<- "aNotDE_log2FCb1"
	d$color[d$log2FoldChange < upFCcut & d$log2FoldChange > downFCcut & d$padjf <= fdrcutoff ]<- "cDE_log2FCb1"
	d$color[d$log2FoldChange >= upFCcut & d$padjf > fdrcutoff | d$log2FoldChange <= downFCcut & d$padjf > fdrcutoff ]<- "bNotDE_log2FCa1"
	d$color[d$log2FoldChange >= upFCcut & d$padjf <= fdrcutoff | d$log2FoldChange <= downFCcut & d$padjf <= fdrcutoff ] <- "dDE_log2FCa1"
	d$color[is.na(d$color)] <- "aNotDE_log2FCb1"
	
	r <- round(max(abs(na.omit(d$log2FoldChange)))) + 1
	
	
	d$lab <- d$Gene
	d$lab[d$log2FoldChange >= 0 & d$padjf >= up | is.na(d$padjf)] <- NA
	d$lab[d$log2FoldChange <= 0 & d$padjf >= down | is.na(d$padjf) ] <- NA
	
	
	d$Gene <- as.character(d$Gene)
	d$padjf <- as.numeric(d$padjf)
	d$logPadjf <- -log(d$padjf, base = 10)

	### make volcano plot
	
	p <- ggplot2::ggplot(data=d, aes(x = log2FoldChange, y = logPadjf)) +
	  geom_point(aes(color=d$color), size=3, alpha = 0.5) +
	  geom_text_repel(aes(label = d$lab), point.padding = 0.25, size = 3) + 
	  xlab(paste("log2 Fold change ( ", i, " )", sep="")) + scale_x_continuous(limits = c(-r, r)) +
	  ylab("-log10 FDR adjusted p-value") + 
	  theme(legend.key.size = unit(2, 'lines')) + 
	  labs(color=d$color) +
	  theme_bw() + theme(text = element_text(size=10), 
	     legend.text=element_text(size=8), legend.key.height = unit(2, 'lines'),
	     legend.position = "right")
	
	blackN <- nrow(d[which(d$color == "aNotDE_log2FCb1"),])
	redN <- nrow(d[which(d$color == "cDE_log2FCb1"),])
	orangeN <- nrow(d[which(d$color == "bNotDE_log2FCa1"),])
	darkredN <- nrow(d[which(d$color == "dDE_log2FCa1"),])
	
	
	vList <- c()
	cList <- c()
	
	if (blackN > 0 ) {
	  cList <- append(cList,"black")
	  vList <- append(vList, paste("Not DE & \nabsolute log2FC < ", upFCcut))
	}
	if (orangeN > 0 ) {
	  cList <- append(cList,"orange")
	  vList <- append(vList, paste("Not DE & \nabsolute log2FC > ", upFCcut))
	}
	if (redN > 0 ) {
	  cList <- append(cList,"red")
	  vList <- append(vList, paste("FDR < ", fdrcutoff, " & \nabsolute log2FC < ", upFCcut))
	}
	if (darkredN > 0 ) {
	  cList <- append(cList,"darkred")
	  vList <- append(vList, paste("FDR < ", fdrcutoff, " & \nabsolute log2FC > ", upFCcut))
	}
	
	c <- scale_color_manual(values= cList, name = "Status",
	                        labels = vList)
	
	
	p <- p + c
	
	pdf(paste("DE/VolcanoPlots/",paste(i,"pdf", sep="."), sep=""))
	print(p)
	dev.off()
	

	colnames(res2) <- paste(i, c("log2FC", "FC","RawP", "FDR"), sep=".")
	allgenes <- cbind(as.data.frame(res2), allgenes)
	png(filename=paste("DE/MAplots/",paste(i,"png", sep="."), sep=""), height=700, width=700)
	plot(log2(res$baseMean),res$log2FoldChange, pch=19, col="black", cex=.9, xlab="log2 mean read count", ylab="log2 Fold Change", main=i, cex.lab=1.4, cex.main=1.4, cex.axis=1.4, ylim=c(-r,r))
	points(log2(fdr$baseMean),fdr$log2FoldChange, pch=19, col="red", cex=.9)
	legend("topright", paste("DE FDR < ",fdrcutoff), pch=19, col="red", cex=1.4)
	abline(h=0, lwd=3)
	abline(h=1, col="blue", lty=3, lwd=3)
	abline(h=-1, col="blue", lty=3, lwd=3)
	dev.off()
	}





#### DE tests for individual effects ####





comps_1 <- resultsNames(ddsModel)[(grep("Int",resultsNames(ddsModel),invert=T))][-c(1:5)]

#i <- comps_1[4]

print("Individual effects")
for (i in rev(comps_1)) {
	print(i)
	res <- results(ddsModel, name=i, alpha=fdrcutoff, cooksCutoff =.9)
	fdr <- res[which(res$padj < fdrcutoff & abs(res$log2FoldChange) > 0),]
	rawp <- res[which(res$pvalue < rawpcutoff),]
	#sigcounts <- rbind(sigcounts,c(nrow(fdr), nrow(rawp)))
	#syms <- c(syms,rownames(fdr))
	
	
	
	res2 <- as.data.frame(res[,c(2,5,6)])
	res2$fc <- 2**res2[,1]
	res2 <- res2[,c(1,4,2,3)]
	
	d <- merge(annot, res2, by='row.names')
	
	#if (RunDavid == TRUE) {
	#  DAVID.analysis(res2, fdrFilter = fdrcutoff, lfcFilter = lfc_cutoff, name = i)
	#}
	#
	#if (RunPathway == TRUE) {
	#  dirPathway <- paste("DE/pathways/",i,sep="")
	#  dir.create(dirPathway)
	#  dirCd <- paste(dir, dirPathway, sep="/")
	#  setwd(dirCd)
	#  
	#  Pathway.analysis(d, species, i)
	#  
	#  setwd(dir)
	#}
	
	
	
	d$padjf <- d$padj + 1e-300
	sub <- subset(d, log2FoldChange>1)
	up <- sub[order(sub$padj),][5,18]
	if ( up <= 1e-300 & !is.na(up)){
	  up = 0 }
	sub <- subset(d, log2FoldChange< -1)
	down <- sub[order(sub$padj),][5,18]
	if ( down <= 1e-300 & !is.na(down)){
	  down = 0 }
	if ( is.na(up)){
	  up = 0.1
	}
	if ( is.na(down)){
	  down = 0.1
	}
	
	
	upFCcut <- 2
	downFCcut <- -2
	
	
	
	
	d$color[d$log2FoldChange < upFCcut & d$log2FoldChange > downFCcut & d$padjf > fdrcutoff ]<- "aNotDE_log2FCb1"
	d$color[d$log2FoldChange < upFCcut & d$log2FoldChange > downFCcut & d$padjf <= fdrcutoff ]<- "cDE_log2FCb1"
	d$color[d$log2FoldChange >= upFCcut & d$padjf > fdrcutoff | d$log2FoldChange <= downFCcut & d$padjf > fdrcutoff ]<- "bNotDE_log2FCa1"
	d$color[d$log2FoldChange >= upFCcut & d$padjf <= fdrcutoff | d$log2FoldChange <= downFCcut & d$padjf <= fdrcutoff ] <- "dDE_log2FCa1"
	d$color[is.na(d$color)] <- "aNotDE_log2FCb1"
	
	r <- round(max(abs(na.omit(d$log2FoldChange)))) + 1
	
	
	d$lab <- d$Gene
	d$lab[d$log2FoldChange >= 0 & d$padjf >= up | is.na(d$padjf)] <- NA
	d$lab[d$log2FoldChange <= 0 & d$padjf >= down | is.na(d$padjf) ] <- NA
	
	
	d$Gene <- as.character(d$Gene)
	d$padjf <- as.numeric(d$padjf)
	d$logPadjf <- -log(d$padjf, base = 10)
	
	### make volcano plot
	
	p <- ggplot2::ggplot(data=d, aes(x = log2FoldChange, y = logPadjf)) +
	  geom_point(aes(color=d$color), size=3, alpha = 0.5) +
	  geom_text_repel(aes(label = d$lab), point.padding = 0.25, size = 3) + 
	  xlab(paste("log2 Fold change ( ", i, " )", sep="")) + scale_x_continuous(limits = c(-r, r)) +
	  ylab("-log10 FDR adjusted p-value") + 
	  theme(legend.key.size = unit(2, 'lines')) + 
	  labs(color=d$color) +
	  theme_bw() + theme(text = element_text(size=10), 
	                     legend.text=element_text(size=8), legend.key.height = unit(2, 'lines'),
	                     legend.position = "right")
	
	blackN <- nrow(d[which(d$color == "aNotDE_log2FCb1"),])
	redN <- nrow(d[which(d$color == "cDE_log2FCb1"),])
	orangeN <- nrow(d[which(d$color == "bNotDE_log2FCa1"),])
	darkredN <- nrow(d[which(d$color == "dDE_log2FCa1"),])
	
	
	vList <- c()
	cList <- c()
	
	if (blackN > 0 ) {
	  cList <- append(cList,"black")
	  vList <- append(vList, paste("Not DE & \nabsolute log2FC < ", upFCcut))
	}
	if (orangeN > 0 ) {
	  cList <- append(cList,"orange")
	  vList <- append(vList, paste("Not DE & \nabsolute log2FC > ", upFCcut))
	}
	if (redN > 0 ) {
	  cList <- append(cList,"red")
	  vList <- append(vList, paste("FDR < ", fdrcutoff, " & \nabsolute log2FC < ", upFCcut))
	}
	if (darkredN > 0 ) {
	  cList <- append(cList,"darkred")
	  vList <- append(vList, paste("FDR < ", fdrcutoff, " & \nabsolute log2FC > ", upFCcut))
	}
	
	c <- scale_color_manual(values= cList, name = "Status",
	                        labels = vList)
	
	
	p <- p + c
	
	pdf(paste("DE/VolcanoPlots/",paste(i,"pdf", sep="."), sep=""))
	print(p)
	dev.off()
	
	colnames(res2) <- paste(i, c("log2FC", "FC","RawP", "FDR"), sep=".")
	#allgenes <- cbind(as.data.frame(res2), allgenes)
	png(filename=paste("DE/MAplots/",paste(i,"png", sep="."), sep=""), height=700, width=700)
	plot(log2(res$baseMean),res$log2FoldChange, pch=19, col="black", cex=.9, xlab="log2 mean read count", ylab="log2 Fold Change", main=i, cex.lab=1.4, cex.main=1.4, cex.axis=1.4, ylim=c(-r,r))
	points(log2(fdr$baseMean),fdr$log2FoldChange, pch=19, col="red", cex=.9)
	legend("topright", paste("DE FDR < ",fdrcutoff), pch=19, col="red", cex=1.4)
	abline(h=0, lwd=3)
	abline(h=1, col="blue", lty=3, lwd=3)
	abline(h=-1, col="blue", lty=3, lwd=3)
	dev.off()
}

baseMean <- res$baseMean
baseMean <- as.matrix(baseMean)
row.names(baseMean) <- row.names(res)
colnames(baseMean) <- c("Mean_Normalized_Counts")
allgenes <- merge(allgenes, baseMean, by="row.names", all.y=T)


syms <- unique(syms)
genes <- annot[syms,]$Gene
M <- stabilizedBR[syms,]
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



allgenes$Ensembl_ID <- allgenes[,"Row.names"]
allgenes1 <- allgenes[,c(ncol(allgenes), ncol(allgenes)-1, 2:(ncol(allgenes)-2))]
allgenes2 <- merge(annot, allgenes1, by="Ensembl_ID", all.y=T)
#d <- merge(allgenes2, fpkmAve, by.x = "Ensembl_ID", by.y = "row.names", all.x=T)
norm_counts_ord <- norm_counts[order(rownames(norm_counts), decreasing=F),]

rownames(sigcounts) <- c(rev(comps),rev(comps_1))
colnames(sigcounts) <- c(paste("FDR < ",fdrcutoff), paste("RawP <", rawpcutoff))
#sigcounts <- sigcounts[order(rev(rownames(sigcounts))),]

dir.create(paste("DE/DEfiles/RawP",as.character(rawpcutoff), sep=""))
dir.create(paste("DE/DEfiles/FDR",as.character(fdrcutoff), sep=""))


for (i in comps) {
  a <- allgenes2[,grep(i, colnames(allgenes2))]
  b <- allgenes2[,c(1:3)]
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
  
for (i in comps) {
  a <- allgenes2[,grep(i, colnames(allgenes2))]
  b <- allgenes2[,c(1:3)]
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


select <- order(rv, decreasing=T)[seq_len(3000)]
v <- allgenes2[select,][,c(1:3)]
write.table(v, file=paste(paste("DE/Top_3000_var_genes_from_PCA",project, sep = "."),"txt", sep="."), sep="\t", quote=F, row.names=F)
write.table(allgenes2, file=paste(paste("DE/All_genes",project, sep = "."),"txt", sep="."), sep="\t", quote=F, row.names=F)
write.table(sigcounts, file=paste(paste("DE/DEcounts",project, sep = "."),"txt", sep="."), sep="\t", quote=F, row.names=T, col.names=NA)

save.image(file=paste(project,".RData",sep=""))

DE.heatmap(M, dendrogram="both", fact=conds[,c("condition"), drop=F], zlim =c(-2,2), lhei=c(1,.25,5,1), height = 1200, width=1200, main=paste("DE genes - Union of all Pairwise Comparisons - FDR < ", fdrcutoff), file.prefix="DE/Heatmaps/DE_heatmap_nosymbols", breaks=.2, Symbols=NULL, min.col = "blue", mid.col = "white", max.col = "red", margins = c(12,12,12,12))
DE.heatmap(M, dendrogram="both", fact=conds[,c("condition"), drop=F], zlim =c(-2,2), lhei=c(.5,.125,10,.5),  width=1200,  main=paste("DE genes - Union of all Comparisons - FDR < ", fdrcutoff), file.prefix="DE/Heatmaps/DE_heatmap_symbols", breaks=.2, Symbols=genes, min.col = "blue", mid.col = "white", max.col = "red", margins = c(12,12,12,12))

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}



condsSmall <- as.data.frame(conds[,c("condition")])
colnames(condsSmall) <- "condition"
rownames(condsSmall) <- rownames(conds)

cols <- gg_color_hue(6)
names(cols) = levels(condsSmall$condition)

ha = HeatmapAnnotation(df = condsSmall, 
                            col = list(condition = cols),
                            show_legend = T)

ht = Heatmap(M, name = "Normalized \nExpression", 
                top_annotation = ha, 
                col = colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red")), 
                cluster_rows = T, cluster_columns = T, show_heatmap_legend = T,
                show_column_names = T, show_row_names = F)

png("DE/Heatmaps/DE_heatmap_nosymbols_complex.png",height = 12, width = 11,res = 150,units = "in")
print(ht)
dev.off()


#### COVID-19 ####

gg_color_hue <- function(n) {
  hues = seq(25, 365, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

condsSmall <- as.data.frame(conds[,c("condition")])
colnames(condsSmall) <- "condition"
rownames(condsSmall) <- rownames(conds)

cols <- gg_color_hue(6)
names(cols) = levels(condsSmall$condition)

SARS_COV_2_table = read.table("../SARS_COV_2_interactions_paper.txt",sep="\t",header=T)

SARS_COV_2_table$Source[SARS_COV_2_table$Gene == "FURIN"] =  "Manual"
SARS_COV_2_table$Source[SARS_COV_2_table$Gene == "PCSK4"] =  "Manual"
SARS_COV_2_table$Source[SARS_COV_2_table$Gene == "PCSK7"] =  "Manual"

dim(SARS_COV_2_table[SARS_COV_2_table$Source == "Krogan",])

#ATP5MG <>  ATP5L
#CEP43 <>   FGFR1OP
#MTARC1 <>  MARC1
#POGLUT2 <> KDELC1
#POGLUT3 <> KDELC2
#TLE5 <>    AES

SARS_COV_2_table$Gene = factor(SARS_COV_2_table$Gene,levels=c(levels(SARS_COV_2_table$Gene),"ATP5L","FGFR1OP","MARC1","KDELC1","KDELC2","AES"))

SARS_COV_2_table$Gene[SARS_COV_2_table$Gene == "ATP5MG"] =  "ATP5L"
SARS_COV_2_table$Gene[SARS_COV_2_table$Gene == "CEP43"] =   "FGFR1OP"
SARS_COV_2_table$Gene[SARS_COV_2_table$Gene == "MTARC1"] =  "MARC1"
SARS_COV_2_table$Gene[SARS_COV_2_table$Gene == "POGLUT2"] = "KDELC1"
SARS_COV_2_table$Gene[SARS_COV_2_table$Gene == "POGLUT3"] = "KDELC2"
SARS_COV_2_table$Gene[SARS_COV_2_table$Gene == "TLE5"] =    "AES"

SARS_COV_2_table$Gene = factor(SARS_COV_2_table$Gene)

SARS_COV_2_table[,c("Virus","Viral_Protein")]  = str_split_fixed(SARS_COV_2_table$Bait, " ", 2)
SARS_COV_allgenes = merge(SARS_COV_2_table, allgenes2, by="Gene")
SARS_COV_allgenes_reduc = SARS_COV_allgenes[,c(1,3,4,5,6,7,8,18,seq(39,58))]
colnames(SARS_COV_allgenes_reduc)[5] = "Bait"
write.table(SARS_COV_allgenes_reduc,
            file = "~/Box/Erle group presentations/LRB/SARS-CoV-2 paper/Bulk_RNA_Krogen_plus_manual_results.txt",
            sep="\t",quote=F,row.names = F)
write.table(SARS_COV_allgenes,
            file = "~/Box/Erle group presentations/LRB/SARS-CoV-2 paper/Bulk_RNA_Krogen_plus_manual_results.full.txt",
            sep="\t",quote=F,row.names = F)

SARS_COV_2_list = c(levels(SARS_COV_2_table$Gene))

SARS_COV_2_table_Krogan = SARS_COV_2_table[SARS_COV_2_table$Source == "Krogan",]

dataPerMillion = sweep(data, 2, (colSums(data)/1e6), FUN = '/')

dataPerMillion_geneMeans = rowMeans(dataPerMillion)
#dataPerMillion_geneMeans_gt1 = dataPerMillion_geneMeans[dataPerMillion_geneMeans >= 1]

dataPerMillion_binary = dataPerMillion
dataPerMillion_binary[dataPerMillion_binary >=1 ] = 1
dataPerMillion_binary[dataPerMillion_binary <1 ] = 0
dataPerMillion_binary_50 = dataPerMillion_binary[rowSums(dataPerMillion_binary) >= 16,]

#dataPerMillion_binary_50_aveGT1 = dataPerMillion_binary_50[rownames(dataPerMillion_binary_50) %in% names(dataPerMillion_geneMeans),]



allgenes_Krogan = allgenes2[allgenes2$Gene %in% SARS_COV_2_table_Krogan$Gene,]
allgenes_Krogan = allgenes_Krogan[order(allgenes_Krogan$Mean_Normalized_Counts,decreasing = T),]
allgenes_Krogan = allgenes_Krogan[!(duplicated(allgenes_Krogan$Gene)),]

dim(allgenes_Krogan)

allgenes_Krogan_g1 = allgenes_Krogan[allgenes_Krogan$Mean_Normalized_Counts >= 1,]

norm_counts_binary = norm_counts
norm_counts_binary[norm_counts_binary >= 1] = 1

dim(norm_counts_binary)

norm_counts_binary_50 = norm_counts_binary[rowSums(norm_counts_binary) >= 16,]
norm_counts_binary_50 = norm_counts_binary_50[rownames(norm_counts_binary_50) %in% allgenes_Krogan_g1$Ensembl_ID,]

dim(allgenes_Krogan[allgenes_Krogan$Ensembl_ID %in% rownames(dataPerMillion_binary_50),])

d = SARS_COV_allgenes_reduc[(SARS_COV_allgenes_reduc$Source != "Krogan"),]
dim(d[d$Ensembl_ID %in% rownames(dataPerMillion_binary_50),])

dim(SARS_COV_allgenes_reduc[SARS_COV_allgenes_reduc$Ensembl_ID %in% rownames(dataPerMillion_binary_50),])

#d = read.table("~/Downloads/media-5_v2.txt",sep="\t",header=T)
#
#d[!(d$PreyGene %in% SARS_COV_2_table_Krogan$Gene),]
#
#SARS_COV_2_table_Krogan[!(SARS_COV_2_table_Krogan$Gene %in% d$PreyGene),]

### Add furon ### 

### Ordering for heatmap UT, IL17, IFNa, IFNg, IL13, CO

## better control for expression level

## scRNA-seq celltype specific markers from the COV2 interactions

## Add some classic cell type markers

SARS_COV_2_list_DE = c()

FDR_cut = 0.05
fc_cut = log2(1.5)

for (i in comps){
  print(i)
  FDR_col = paste(i,"FDR",sep=".")
  logFC_col = paste(i,"log2FC",sep=".")
  allgenes_tmp = allgenes2[(allgenes2[,FDR_col] <= FDR_cut),]
  allgenes_tmp = allgenes_tmp[(abs(allgenes_tmp[,logFC_col]) >= fc_cut),]
  allgenes_tmp = allgenes_tmp[!(is.na(allgenes_tmp$Ensembl_ID)),]
  allgenes_tmp = allgenes_tmp[(allgenes_tmp$Gene %in% SARS_COV_2_list),]
  
  #SARS_COV_allgenes_tmp = SARS_COV_allgenes[SARS_COV_allgenes$Gene %in% SARS_COV_2_list,]
  allgenes_tmp_forSteph = SARS_COV_allgenes[(SARS_COV_allgenes$Ensembl_ID %in% allgenes_tmp$Ensembl_ID),]
  allgenes_tmp_forSteph = allgenes_tmp_forSteph[,c("Gene","Bait","Type","Source","Virus","Viral_Protein","Ensembl_ID",
                                                   paste(i,c("log2FC","FC","RawP","FDR"),sep="."))]
  print(dim(allgenes_tmp_forSteph))
  write.table(allgenes_tmp_forSteph,
              file = paste("~/Box/Erle group presentations/LRB/SARS-CoV-2 paper/Bulk_RNA_Krogen_plus_manual_results.",i,".txt",sep=""),
              sep="\t",quote=F,row.names = F)
  
  allgenes_tmp_up = allgenes_tmp[(allgenes_tmp[,logFC_col] > 0),]
  allgenes_tmp_down = allgenes_tmp[(allgenes_tmp[,logFC_col] < 0),]
  
  allgenes_tmp_up_list = allgenes_tmp_up$Gene
  allgenes_tmp_down_list = allgenes_tmp_down$Gene
  
  write.table(allgenes_tmp_up_list,
              file = paste("~/Box/Erle group presentations/LRB/SARS-CoV-2 paper/Up_gene_list.",i,".txt",sep=""),
              sep="\t",quote=F,row.names = F,col.names = F)
  
  write.table(allgenes_tmp_down_list,
              file = paste("~/Box/Erle group presentations/LRB/SARS-CoV-2 paper/Down_gene_list.",i,".txt",sep=""),
              sep="\t",quote=F,row.names = F,col.names = F)
  
  
  SARS_COV_2_list_DE = c(SARS_COV_2_list_DE,allgenes_tmp[,"Gene"])
  SARS_COV_2_list_DE = unique(SARS_COV_2_list_DE)
  print(length(unique(allgenes_tmp[,"Gene"])))
}

allgenes_SARS2 = allgenes2[allgenes2$Gene %in% SARS_COV_2_list,]
write.table(allgenes_SARS2, file="~/Box/COV_plot_and_files/allgenes_SARS2_interactions.txt",sep='\t',quote = F,row.names = F)

allgenes_SARS2_DE = allgenes2[allgenes2$Gene %in% SARS_COV_2_list_DE,]

SARS_COV_2_list_DE_Ens = allgenes_SARS2_DE$Ensembl_ID

allgenes_SARS2_tmp = merge(allgenes_SARS2, SARS_COV_2_table,by.x="Gene",by.y="Gene")

write.table(SARS_COV_allgenes,
            file = "~/Box/Erle group presentations/LRB/SARS-CoV-2 paper/Bulk_RNA_Krogen_plus_manual_results.txt",
            sep="\t",quote=F,row.names = F)


write.table(allgenes_SARS2_DE,
            file = "~/Box/Erle group presentations/LRB/SARS-CoV-2 paper/Bulk_RNA_Krogen_plus_manual_results.DE_only.txt",
            sep="\t",quote=F,row.names = F)

#### Krogen set ####

SARS_COV_2_Krogen_list_DE = c()

SARS_COV_2_table_Krogan_annot = merge(SARS_COV_2_table_Krogan,annot,by="Gene")

FDR_cut = 0.05
fc_cut = log2(1.5)

for (i in comps){
  FDR_col = paste(i,"FDR",sep=".")
  logFC_col = paste(i,"log2FC",sep=".")
  allgenes_tmp = allgenes2[(allgenes2[,FDR_col] <= FDR_cut),]
  allgenes_tmp = allgenes_tmp[(abs(allgenes_tmp[,logFC_col]) >= fc_cut),]
  allgenes_tmp = allgenes_tmp[!(is.na(allgenes_tmp$Ensembl_ID)),]
  allgenes_tmp = allgenes_tmp[(allgenes_tmp$Ensembl_ID %in% SARS_COV_2_table_Krogan_annot$Ensembl_ID),]
  print(i)
  print(nrow(allgenes_tmp))
  print((nrow(allgenes_tmp)/318)*100)
  SARS_COV_2_Krogen_list_DE = c(SARS_COV_2_Krogen_list_DE,allgenes_tmp[,"Gene"])
  SARS_COV_2_Krogen_list_DE = unique(SARS_COV_2_Krogen_list_DE)
}

allgenes_SARS2 = allgenes2[allgenes2$Gene %in% SARS_COV_2_Krogen_list_DE,]
write.table(allgenes_SARS2, file="~/Box/COV_plot_and_files/allgenes_SARS2_interactions.Krogan_only.txt",sep='\t',quote = F,row.names = F)


allgenes_SARS2_DE = allgenes2[allgenes2$Gene %in% SARS_COV_2_Krogen_list_DE,]

SARS_COV_2_list_DE_Ens = allgenes_SARS2_DE$Ensembl_ID

allgenes_SARS2_tmp = merge(allgenes_SARS2, SARS_COV_2_table,by.x="Gene",by.y="Gene")

allgenes_SARS2_annot = merge(SARS_COV_2_table,allgenes_SARS2,by.x="Gene",by.y="Gene")
write.table(allgenes_SARS2_annot, file="~/Box/COV_plot_and_files/allgenes_SARS2_interactions.Krogan_only.txt",sep='\t',quote = F,row.names = F)

allgenes_SARS2_annot_DE = allgenes_SARS2_annot[allgenes_SARS2_annot$Ensembl_ID %in% SARS_COV_2_list_DE_Ens,]
write.table(allgenes_SARS2_annot_DE, file="~/Box/COV_plot_and_files/allgenes_SARS2_interactions.Krogan_only.DE_ONLY.txt",sep='\t',quote = F,row.names = F)

a = allgenes2[,c("Ensembl_ID","Mean_Normalized_Counts")]
b = allgenes_SARS2[,c("Ensembl_ID","Mean_Normalized_Counts")]

a$source = "all_genes"
b$source = "COV2_associated"

d = rowMeans(allgenes_SARS2_tmp[,seq(84,89)])
e = rowMeans(allgenes2[,seq(84,89)])

d = as.data.frame(d)
e = as.data.frame(e)

d$gene = allgenes_SARS2_tmp$Gene
e$gene = allgenes2$Gene

d$source = "COV2_associated"
e$source = "all_genes"

colnames(d) = c("Mean_Normalized_Counts","gene","source")
colnames(e) = c("Mean_Normalized_Counts","gene","source")

c = rbind(a,b)

f = rbind(e,d)

ACE2_exp = f[f$gene == "ACE2",]
ACE2_exp = ACE2_exp[!(is.na(ACE2_exp$gene)),]
ACE2_exp = ACE2_exp[1,"Mean_Normalized_Counts"]

TMPRSS2_exp = f[f$gene == "TMPRSS2",]
TMPRSS2_exp = TMPRSS2_exp[!(is.na(TMPRSS2_exp$gene)),]
TMPRSS2_exp = TMPRSS2_exp[1,"Mean_Normalized_Counts"]

ggplot(f, aes(x=Mean_Normalized_Counts, fill=source)) + 
  geom_density(alpha=.5) + 
  #scale_y_log10() + 
  scale_x_log10() +
  geom_vline(xintercept=ACE2_exp, linetype="dashed", color = "darkgreen") +
  geom_vline(xintercept=TMPRSS2_exp, linetype="dashed", color = "navyblue") +
  annotate("text", x = (ACE2_exp/2.2), y = .7, label = "ACE2",color="darkgreen",fontface="italic") +
  annotate("text", x = (TMPRSS2_exp*2.7), y = .7, label = "TMPRSS2",color="navyblue",fontface="italic") +
  theme_classic()
ggsave(filename="~/Box/Erle group presentations/LRB/SARS-CoV-2 paper/Walter_COV_2_figures/gene_class_density.UT.pdf",width=6,height=4.5,units = "in",dpi = 300)

ggplot(f, aes(x=Mean_Normalized_Counts, fill=source)) + 
  geom_density(alpha=.5) + 
  #scale_y_log10() + 
  scale_x_log10() +
  #geom_vline(xintercept=ACE2_exp, linetype="dashed", color = "darkgreen") +
  #geom_vline(xintercept=TMPRSS2_exp, linetype="dashed", color = "navyblue") +
  #annotate("text", x = (ACE2_exp/2), y = .7, label = "ACE2",color="darkgreen",fontface="italic") +
  #annotate("text", x = (TMPRSS2_exp*2.5), y = .7, label = "TMPRSS2",color="navyblue",fontface="italic") +
  theme_classic()
ggsave(filename="~/Box/Erle group presentations/LRB/SARS-CoV-2 paper/Walter_COV_2_figures/gene_class_density.UT.unannoted.pdf",width=6,height=4.5,units = "in",dpi = 300)


#### Percentile plots ####

head(norm_counts)

ACE2_ens = "ENSG00000130234"
TMPRSS2_ens = "ENSG00000184012"

norm_counts_filtered = norm_counts[which(rownames(norm_counts) %in% rownames(dataPerMillion_binary_50)),]

geneRanks <- as.data.frame(setNames(replicate(2,numeric(0), simplify = F),letters[0:2]))
colnames(geneRanks) = c("Gene","Individual")
rankVect = c()

norm_counts_UT = norm_counts_filtered[,c(31:36)]

for (i in colnames(norm_counts_UT)){
  norm_counts_UT_tmp = norm_counts_UT[,i]
  ACE2_exp = norm_counts_UT_tmp[ACE2_ens]
  TMPRSS2_exp = norm_counts_UT_tmp[TMPRSS2_ens]
  
  ACE_rank = norm_counts_UT_tmp[norm_counts_UT_tmp >= ACE2_exp]
  ACE_rank = length(ACE_rank[!(is.na(ACE_rank))])
  ACE_perc = 1 - ACE_rank/nrow(norm_counts_filtered)
  
  TMPRSS2_rank = norm_counts_UT_tmp[norm_counts_UT_tmp >= TMPRSS2_exp]
  TMPRSS2_rank = length(ACE_rank[!(is.na(TMPRSS2_rank))])
  TMPRSS2_perc = 1 - TMPRSS2_rank/nrow(norm_counts_filtered)
  
  ACE_row = t(as.data.frame(c("ACE2",i)))
  rownames(ACE_row) = paste("ACE2",i,sep="_")
  geneRanks = rbind(geneRanks,ACE_row)
  
  TMPRSS2_row = t(as.data.frame(c("TMPRSS2",i)))
  rownames(TMPRSS2_row) = paste("TMPRSS2",i,sep="_")
  geneRanks = rbind(geneRanks,TMPRSS2_row)
  
  rankVect = c(rankVect,ACE_perc,TMPRSS2_perc)
  
}
geneRanks_2 = geneRanks
geneRanks_2$Percentile = rankVect

dx = geneRanks_2[geneRanks_2$V1 == "ACE2",]
ACE_mean = mean(dx$Percentile)

dx = geneRanks_2[geneRanks_2$V1 == "TMPRSS2",]
TMPRSS2_mean = mean(dx$Percentile)

ACE_row = t(as.data.frame(c("ACE2","meanPerc")))
rownames(ACE_row) = paste("ACE2","meanPerc",sep="_")
geneRanks = rbind(geneRanks,ACE_row)

TMPRSS2_row = t(as.data.frame(c("TMPRSS2","meanPerc")))
rownames(TMPRSS2_row) = paste("TMPRSS2","meanPerc",sep="_")
geneRanks = rbind(geneRanks,TMPRSS2_row)

rankVect = c(rankVect,ACE_mean,TMPRSS2_mean)

allAverages = rowMeans(norm_counts_UT)

ACE_rank = allAverages[allAverages  >= ACE2_exp]
ACE_rank = length(ACE_rank[!(is.na(ACE_rank))])
ACE_perc = 1 - ACE_rank/nrow(norm_counts_filtered)

TMPRSS2_rank = allAverages[allAverages  >= TMPRSS2_exp]
TMPRSS2_rank = length(ACE_rank[!(is.na(TMPRSS2_rank))])
TMPRSS2_perc = 1 - TMPRSS2_rank/nrow(norm_counts_filtered)

ACE_row = t(as.data.frame(c("ACE2","allAve")))
rownames(ACE_row) = paste("ACE2","allAve",sep="_")
geneRanks = rbind(geneRanks,ACE_row)

TMPRSS2_row = t(as.data.frame(c("TMPRSS2","allAve")))
rownames(TMPRSS2_row) = paste("TMPRSS2","allAve",sep="_")
geneRanks = rbind(geneRanks,TMPRSS2_row)

rankVect = c(rankVect,ACE_perc,TMPRSS2_perc)

colnames(geneRanks) = c("Gene","Individual")
geneRanks$Percentile = rankVect

saveRDS(geneRanks,"geneRanks.rds")

#### COVID Dot plots ####

SARS_dotMatCellTypeTreatmentScaled_withDE = read.table("~/Box/Erle group presentations/LRB/SARS-CoV-2 paper/COV2_cellType_and_cellTypeXtreatment_results.txt",
                                                       sep="\t",header = T)



SARS_dotMatCellTypeTreatmentScaled_withDE$Treatment = factor(SARS_dotMatCellTypeTreatmentScaled_withDE$Treatment, 
                                                             levels=rev(c("UT","IL17","IFN","IL13","CO")))



Gene = "ACE2"

SARS_dotMatCellTypeTreatmentScaled_withDE_Gene = SARS_dotMatCellTypeTreatmentScaled_withDE[SARS_dotMatCellTypeTreatmentScaled_withDE$Gene == Gene,]

tmpLevels = c(levels(SARS_dotMatCellTypeTreatmentScaled_withDE_Gene$Bait),"ACE2","Proteases")
tmpLevels = tmpLevels[c(28,29,27,seq(1,26))]

SARS_dotMatCellTypeTreatmentScaled_withDE$Bait = factor(SARS_dotMatCellTypeTreatmentScaled_withDE$Bait,levels = tmpLevels)

SARS_dotMatCellTypeTreatmentScaled_withDE[SARS_dotMatCellTypeTreatmentScaled_withDE$Type == "Proteases", "Bait"] <- "Proteases"
SARS_dotMatCellTypeTreatmentScaled_withDE[SARS_dotMatCellTypeTreatmentScaled_withDE$Gene == "ACE2", "Bait"] <- "ACE2"

SARS_dotMatCellTypeTreatmentScaled_withDE$cellType = factor(SARS_dotMatCellTypeTreatmentScaled_withDE$cellType,
                                                            levels=c("Proliferating_Basal","Basal","Intermediate",
                                                                     "Secretory","Preciliated","Ciliated","Ionocyte"))

SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut = SARS_dotMatCellTypeTreatmentScaled_withDE[SARS_dotMatCellTypeTreatmentScaled_withDE$cellType.p_val_adj <= 0.1,]
SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut = SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut[!(is.na(SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut$Gene)),]
length(unique(SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut$Gene))

markerMatrix <- as.data.frame(setNames(replicate(14,numeric(0), simplify = F),letters[0:14]))
geneList = c()
for (i in levels(SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut$cellType)){
  tmpData = SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut[(SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut$cellType == i),]
  tmpData = tmpData[tmpData$cellType.avg_logFC >0,]
  tmpData = tmpData[order(tmpData$cellType.p_val_adj),]
  tmpData = head(tmpData,n=20)
  tmpData = tmpData[order(tmpData$cellType.avg_logFC,decreasing = T),]
  tmpData = head(tmpData,n=5)
  markerMatrix = rbind(markerMatrix,tmpData[,seq(1:14)])
  geneList = c(geneList,as.vector(tmpData$Gene))
}

geneList = unique(geneList)

SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut_cellMarkers = SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut[SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut$Gene %in% geneList,]
SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut_cellMarkers = SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut_cellMarkers[!(is.na(SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut$Gene)),]
SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut_cellMarkers = SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut_cellMarkers[match(geneList,
                                                                                       SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut_cellMarkers$Gene),]

SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut_cellMarkers$Gene = factor(SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut_cellMarkers$Gene , levels=geneList)

markerMatrix$Gene = factor(markerMatrix$Gene,levels=unique(as.vector(markerMatrix$Gene)))

markerMatrix$cellType = factor(markerMatrix$cellType, levels=rev(levels(markerMatrix$cellType)))

ggplot(markerMatrix, aes(y = cellType, x = Gene)) +
  geom_point(aes(size = pctExp, colour=aveExp)) +
  scale_size_area() +
  scale_colour_gradientn(colours = c("purple","grey","orange"), 
                         limits=c(-1*abs(max(markerMatrix$sclAveExp)),
                                  abs(max(markerMatrix$sclAveExp)))) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1,color="black",size=12),
        axis.title=element_text(size=16), 
        axis.text.y = element_text(color="black",size=12),
        legend.text=element_text(size=10),legend.title=element_text(size=14),
        panel.grid = element_blank(),plot.title = element_text(hjust = 0.5)) +
  labs(size="Percent\nExpressed", colour="Average\nExpression") + ylab("Cell Type") +
  guides(color = guide_colorbar(order = 1),size = guide_legend(order=2))
ggsave(paste("~/Box/COV_plot_and_files//Curtain.",FileLabel,".twoColor.pdf",sep=""),height = 4,width = 7,units = "in",dpi = 300)


ggplot(markerMatrix, aes(y = cellType, x = Gene)) +
  geom_point(aes(size = pctExp, colour=aveExp)) +
  scale_size_area() +
  scale_colour_gradientn(colours = c("grey","red"), 
                         limits=c(0,abs(max(markerMatrix$aveExp)))) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1,color="black",size=12),
        axis.title=element_text(size=16), 
        axis.text.y = element_text(color="black",size=12),
        legend.text=element_text(size=10),legend.title=element_text(size=14),
        panel.grid = element_blank(),plot.title = element_text(hjust = 0.5)) +
  labs(size="Percent\nExpressed", colour="Average\nExpression") + ylab("Cell Type") +
  guides(color = guide_colorbar(order = 1),size = guide_legend(order=2))
ggsave(paste("~/Box/COV_plot_and_files//Curtain.oneColor.cellTypeMarkerGenes.pdf",sep=""),height = 4,width = 7,units = "in",dpi = 300)

#### Heatmap of marker genes ####

SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut_fc_cut = SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut[
            abs(SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut$cellType.avg_logFC) > log2(1.25),]

SARS_dotMatCellTypeTreatmentScaled_withDE_cellMarkers = SARS_dotMatCellTypeTreatmentScaled_withDE[SARS_dotMatCellTypeTreatmentScaled_withDE$Gene %in%
                                                                                                unique(SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut_fc_cut$Gene),]

SARS_dotMatCellTypeTreatmentScaled_withDE_cellMarkers_UT = SARS_dotMatCellTypeTreatmentScaled_withDE_cellMarkers[SARS_dotMatCellTypeTreatmentScaled_withDE_cellMarkers$Treatment == "UT",]

#SARS_dotMatCellTypeTreatmentScaled_withDE_cellMarkers_UT = SARS_dotMatCellTypeTreatmentScaled_withDE_cellMarkers_UT[abs(SARS_dotMatCellTypeTreatmentScaled_withDE_cellMarkers_UT$cellType.avg_logFC) > .25,]
#SARS_dotMatCellTypeTreatmentScaled_withDE_cellMarkers_UT = SARS_dotMatCellTypeTreatmentScaled_withDE_cellMarkers_UT[!(is.na(SARS_dotMatCellTypeTreatmentScaled_withDE_cellMarkers_UT$Gene)),]

mat = acast(SARS_dotMatCellTypeTreatmentScaled_withDE_cellMarkers_UT, Gene ~ cellType,value.var =  "aveExp")

mat = mat[,-c(7)]

M = mat
M = M - rowMeans(M)


condMat = SARS_dotMatCellTypeTreatmentScaled_withDE_cellMarkers_UT[!(duplicated(SARS_dotMatCellTypeTreatmentScaled_withDE_cellMarkers_UT$Gene)),]
condMat = condMat[!(is.na(condMat$Gene)),]
condMat = condMat[,c(1,2,3,4,5)]

condMat$Type = factor(condMat$Type)



cols_mat <- gg_color_hue(length(unique(condMat$Type)))
names(cols_mat) = levels(condMat$Type)

rownames(condMat) = condMat$Gene
condMat = condMat[match(rownames(M),rownames(condMat)),]

#cols_mat = cols_mat[rownames(cols_mat) %in% rownames(M),]

#cols_mat = cols_mat[match(rownames(M),)]

condMat_Short = condMat$Type
names(condMat_Short) = rownames(condMat)
condMat_Short = as.data.frame(condMat_Short)
colnames(condMat_Short) = "Type"

ha_r = rowAnnotation(df = condMat_Short, 
                     col = list(Type = cols_mat),
                     show_legend = T)


ht = Heatmap(M, name = "Average \nExpression", 
             right_annotation = ha_r,
             col = colorRamp2(c(-.75,0, .75), c("blue","white", "red")), 
             cluster_rows = T, cluster_columns = T, show_heatmap_legend = T,
             show_column_names = T, show_row_names = T,row_dend_width = unit(20, "mm"))


pdf("~/Box/Erle group presentations/LRB/SARS-CoV-2 paper/Walter_COV_2_figures/CellType_marker_Heatmap.colCluster.noIon.pdf",height = 16, width = 11)
print(ht)
dev.off()

#### Dot plot of cell markers ####


SARS_dotMatCellTypeTreatmentScaled_withDE$cellType = factor(SARS_dotMatCellTypeTreatmentScaled_withDE$cellType ,
                                                           levels=c("Proliferating_Basal","Basal","Intermediate",
                                                                    "Secretory","Preciliated","Ciliated","Ionocyte"))

SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut = SARS_dotMatCellTypeTreatmentScaled_withDE[SARS_dotMatCellTypeTreatmentScaled_withDE$cellType.p_val_adj <= 0.05,]
SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut = SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut[!(is.na(SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut$Gene)),]

SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut$cellType = factor(SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut$cellType ,
                                                                   levels=c("Proliferating_Basal","Basal","Intermediate",
                                                                            "Secretory","Preciliated","Ciliated"))



SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut = SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut[
                                                              abs(SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut$cellType.avg_logFC) >= log2(1.25),]
SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut_noION = SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut[
                                                              (SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut$cellType != "Ionocyte"),]

plotMat = SARS_dotMatCellTypeTreatmentScaled_withDE[SARS_dotMatCellTypeTreatmentScaled_withDE$Gene %in% 
                                                      SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut_noION$Gene,]
plotMat = plotMat[!(is.na(plotMat$Gene)),]

plotMat = plotMat[plotMat$Treatment == "UT",]
plotMat = plotMat[!(is.na(plotMat$Gene)),]


mat = acast(plotMat, Gene ~ cellType,value.var =  "aveExp")

mat = mat[,-c(7)] # removing Ionocytes

dend = as.dendrogram(hclust(d = dist(x = mat)))
geneOrder = dend %>% labels

#SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut_noION = SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut_noION[
#                                                                order(SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut_noION$cellType.avg_logFC,decreasing = T),]
#
#SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut_noION = SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut_noION[
#                                                                !duplicated(SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut_noION$Gene),]
#
#SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut_noION = SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut_noION[
#                                                                order(SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut_noION$cellType),]
#orderByMarker = as.vector(SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut_noION$Gene)
#orderByMarker = orderByMarker[!duplicated(orderByMarker)]

plotMatForOrder = plotMat[order(plotMat$cellType.avg_logFC,decreasing = T),]
plotMatForOrder = plotMatForOrder[(plotMatForOrder$cellType.p_val_adj <= 0.05),]
plotMatForOrder = plotMatForOrder[!(is.na(plotMatForOrder$Gene)),]
plotMatForOrder = plotMatForOrder[!duplicated(plotMatForOrder$Gene),]
plotMatForOrder = plotMatForOrder[!(is.na(plotMatForOrder$Gene)),]
plotMatForOrder = plotMatForOrder[order(plotMatForOrder$cellType),]

orderByMarker = as.vector(plotMatForOrder$Gene)
orderByMarker = orderByMarker[!duplicated(orderByMarker)]

#plotMat$Gene = factor(plotMat$Gene, levels = geneOrder)
plotMat$Gene = factor(plotMat$Gene, levels = rev(orderByMarker))

plotMat = plotMat[(plotMat$cellType != "Ionocyte"),]

plotMat$plot_sclAveExp = plotMat$sclAveExp
plotMat$plot_sclAveExp[plotMat$plot_sclAveExp >1] = 1
plotMat$plot_sclAveExp[plotMat$plot_sclAveExp < -1] = -1

p = ggplot(plotMat, aes(y = Gene, x = cellType)) +
  geom_point(aes(size = pctExp, colour=sclAveExp)) +
  scale_size_area() +
  scale_colour_gradientn(colours = c("blue","grey","red"), 
                         limits=c(-2,2)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1,color="black",size=12),
        axis.title=element_text(size=16), 
        axis.text.y = element_text(color="black",size=12),
        legend.text=element_text(size=10),legend.title=element_text(size=14),
        panel.grid = element_blank(),plot.title = element_text(hjust = 0.5)) +
  labs(size="Percent\nExpressed", colour="Average\nExpression") + ylab("Cell Type") +
  guides(color = guide_colorbar(order = 1),size = guide_legend(order=2))

p
ggsave(paste("~/Box/COV_plot_and_files/Cell_Type_SARS2_dotPlot.noBoxes.pdf",sep=""),height = 16,width = 16,units = "in",dpi = 300)

a = p +
  geom_point(data=plotMat[(plotMat$cellType.p_val_adj <= 0.1),],
             aes(size=pctExp*1.85,alpha=abs(cellType.avg_logFC)),pch=21, fill=NA, colour="black", stroke=1.25)+
  scale_y_discrete(limits=rev(orderByMarker)) +
  scale_x_discrete(limits=c("Proliferating_Basal","Basal","Intermediate",
                                "Secretory","Preciliated","Ciliated"))+
  scale_alpha_continuous("Absolute\nLog2 Fold Change\nFDR <= 0.05",breaks = c(.25, .5, 1),range=c(.25,1),limits=c(.25,1))
a

ggsave(paste("~/Box/COV_plot_and_files/Cell_Type_SARS2_dotPlot.CirclesByFC.v2.pdf",sep=""),height = 12,width = 12,units = "in",dpi = 300)

a = p +
  geom_point(data=plotMat[(plotMat$cellType.p_val_adj <= 0.1),],
             aes(size=pctExp*1.85,),pch=21, fill=NA, colour="black", stroke=1.25)+
  scale_y_discrete(limits=rev(orderByMarker)) +
  scale_x_discrete(limits=c("Proliferating_Basal","Basal","Intermediate",
                            "Secretory","Preciliated","Ciliated"))+
  scale_alpha_continuous("Absolute\nLog2 Fold Change\nFDR <= 0.05",breaks = c(.25, .5, 1),range=c(.25,1),limits=c(.25,1))
a

ggsave(paste("~/Box/COV_plot_and_files/Cell_Type_SARS2_dotPlot.Circles.v2.pdf",sep=""),height = 12,width = 12,units = "in",dpi = 300)

a = p +
  geom_point(data=plotMat[(plotMat$cellType.p_val_adj <= 0.1 & abs(plotMat$cellType.avg_logFC) >= log2(1.25)),],
             aes(size=pctExp*1.85,),pch=21, fill=NA, colour="black", stroke=1.25)+
  scale_y_discrete(limits=rev(orderByMarker)) +
  scale_x_discrete(limits=c("Proliferating_Basal","Basal","Intermediate",
                            "Secretory","Preciliated","Ciliated"))+
  scale_alpha_continuous("Absolute\nLog2 Fold Change\nFDR <= 0.1",breaks = c(.25, .5, 1),range=c(.25,1),limits=c(.25,1))
a

ggsave(paste("~/Box/COV_plot_and_files/Cell_Type_SARS2_dotPlot.CirclesWithFC_cut.pdf",sep=""),height = 12,width = 12,units = "in",dpi = 300)

a = p+
  annotate(xmin = .5, xmax = 1.5, 
           ymin = 27.5, ymax = 31.5, 
           geom = "rect", alpha = 0.2) +
  annotate(xmin = 1.5, xmax = 2.5, 
           ymin = 23.5, ymax = 27.5, 
           geom = "rect", alpha = 0.2) + 
  annotate(xmin = 2.5, xmax = 3.5, 
         ymin = 22.5, ymax = 23.5, 
         geom = "rect", alpha = 0.2) +
  annotate(xmin = 3.5, xmax = 4.5, 
           ymin = 18.5, ymax = 22.5, 
           geom = "rect", alpha = 0.2) + 
  annotate(xmin = 4.5, xmax = 5.5, 
           ymin = 13.5, ymax = 18.5, 
           geom = "rect", alpha = 0.2) +
  annotate(xmin = 4.5, xmax = 5.5, 
           ymin = 30.5, ymax = 31.5,  
           geom = "rect", alpha = 0.2) +
  annotate(xmin = 5.5, xmax = 6.5, 
           ymin = 0.5, ymax = 13.5, 
           geom = "rect", alpha = 0.2) +
  annotate(xmin = 5.5, xmax = 6.5, 
           ymin = 16.5, ymax = 17.5, 
           geom = "rect", alpha = 0.2) +
  annotate(xmin = 5.5, xmax = 6.5, 
           ymin = 19.5, ymax = 20.5, 
           geom = "rect", alpha = 0.2)
a
  
ggsave(paste("~/Box/COV_plot_and_files/Cell_Type_SARS2_dotPlot.boxes.pdf",sep=""),height = 8,width = 8,units = "in",dpi = 300)

#### Cytokine specific effects of genes,dot pltos####

SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut_Cyt = SARS_dotMatCellTypeTreatmentScaled_withDE[SARS_dotMatCellTypeTreatmentScaled_withDE$cellTypeXtreatment.p_val_adj <= 0.05,]
SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut_Cyt = SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut_Cyt[!(is.na(SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut_Cyt$Gene)),]

SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut_Cyt$cellType = factor(SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut_Cyt$cellType ,
                                                                   levels=c("Proliferating_Basal","Basal","Intermediate",
                                                                            "Secretory","Preciliated","Ciliated"))



SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut_Cyt = SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut_Cyt[
  abs(SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut_Cyt$cellTypeXtreatment.logFC) >= log2(1.5),]
SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut_Cyt_noION = SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut_Cyt[
  !(SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut_Cyt$cellType %in% c("Preciliated","Ionocyte")),]

deTableCellxTreatment <- as.data.frame(setNames(replicate(5,numeric(0), simplify = F),letters[0:5]))


for (i in levels(SARS_dotMatCellTypeTreatmentScaled_withDE$cellType)){
  cellRow = c()
  for (j in rev(levels(SARS_dotMatCellTypeTreatmentScaled_withDE$Treatment))[-1]){
    k = SARS_dotMatCellTypeTreatmentScaled_withDE[SARS_dotMatCellTypeTreatmentScaled_withDE$cellType == i,]
    k = k[k$Treatment == j,]
    k = k[k$cellTypeXtreatment.p_val_adj <= 0.05,]
    k = k[!(is.na(k$Gene)),]
    k = k[abs(k$cellTypeXtreatment.logFC) >= log2(1.25),]
    k = k[!(is.na(k$Gene)),]
    print(paste(i,j,sep=" "))
    print(length(unique(k$Gene)))
    cellRow = c(cellRow,length(unique(k$Gene)))
  }
  deTableCellxTreatment = rbind(deTableCellxTreatment ,cellRow)

}

colnames(deTableCellxTreatment) = rev(levels(SARS_dotMatCellTypeTreatmentScaled_withDE$Treatment))[-1]
rownames(deTableCellxTreatment) = levels(SARS_dotMatCellTypeTreatmentScaled_withDE$cellType)

plotMat = SARS_dotMatCellTypeTreatmentScaled_withDE[SARS_dotMatCellTypeTreatmentScaled_withDE$Gene %in% 
                                                      SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut_Cyt_noION$Gene,]

plotMat = plotMat[!(is.na(plotMat$Gene)),]

#plotMat = plotMat[plotMat$Treatment == "UT",]
plotMat = plotMat[!(is.na(plotMat$Gene)),]
plotMat = plotMat[!(plotMat$cellType %in% c("Ionocyte","Preciliated")),]
plotMat = plotMat[!(is.na(plotMat$Gene)),]

mat = acast(plotMat, Gene ~ cellType_Treatment,value.var =  "aveExp")

dend = as.dendrogram(hclust(d = dist(x = mat)))
geneOrder = dend %>% labels

#SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut_noION = SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut_noION[
#                                                                order(SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut_noION$cellType.avg_logFC,decreasing = T),]
#
#SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut_noION = SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut_noION[
#                                                                !duplicated(SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut_noION$Gene),]
#
#SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut_noION = SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut_noION[
#                                                                order(SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut_noION$cellType),]
#orderByMarker = as.vector(SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut_noION$Gene)
#orderByMarker = orderByMarker[!duplicated(orderByMarker)]

plotMatForOrder = plotMat[order(plotMat$cellType.avg_logFC,decreasing = T),]
plotMatForOrder = plotMatForOrder[(plotMatForOrder$cellType.p_val_adj <= 0.05),]
plotMatForOrder = plotMatForOrder[!(is.na(plotMatForOrder$Gene)),]
plotMatForOrder = plotMatForOrder[!duplicated(plotMatForOrder$Gene),]
plotMatForOrder = plotMatForOrder[!(is.na(plotMatForOrder$Gene)),]
plotMatForOrder = plotMatForOrder[order(plotMatForOrder$cellType),]

orderByMarker = as.vector(plotMatForOrder$Gene)
orderByMarker = orderByMarker[!duplicated(orderByMarker)]

plotMat$Gene = factor(plotMat$Gene, levels = geneOrder)
#plotMat$Gene = factor(plotMat$Gene, levels = rev(orderByMarker))

#plotMat = plotMat[(plotMat$cellType != "Ionocyte"),]

plotMat$plot_aveExp = plotMat$aveExp
plotMat$plot_aveExp[plotMat$plot_aveExp >2] = 2
#plotMat$plot_aveExp[plotMat$plot_aveExp < -1] = -1

#  geom_point(aes(size = pctExp, colour=sclAveExp),position = position_dodge(width = 0.75)) +

p = ggplot(plotMat, aes(y = Gene, x = cellType,group=Treatment)) +
  geom_point(aes(size = pctExp, alpha=plot_aveExp,color=Treatment),position = position_dodge(width = 0.75)) +
  scale_size_area() +
  #scale_colour_gradientn(colours = c("blue","grey","red"), 
  #                       limits=c(-1,1)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1,color="black",size=12),
        axis.title=element_text(size=16), 
        axis.text.y = element_text(color="black",size=12),
        legend.text=element_text(size=10),legend.title=element_text(size=14),
        panel.grid = element_blank(),plot.title = element_text(hjust = 0.5)) +
  labs(size="Percent\nExpressed", colour="Average\nExpression") + ylab("Cell Type") +
  guides(color = guide_colorbar(order = 1),size = guide_legend(order=2))

p
ggsave(paste("~/Box/COV_plot_and_files/Cell_Type_SARS2_dotPlot.noBoxes.pdf",sep=""),height = 16,width = 16,units = "in",dpi = 300)

a = p +
  geom_point(data=plotMat[(plotMat$cellTypeXtreatment.p_val_adj <= 0.05 | plotMat$Treatment == "UT"),],
             aes(size=pctExp*1.85,alpha=abs(cellTypeXtreatment.logFC),group=Treatment),position = position_dodge(width = 0.75),
             pch=21, fill=NA, colour="black", stroke=1.25)+
  #scale_y_discrete(limits=rev(orderByMarker)) +
  #scale_x_discrete(limits=c("Proliferating_Basal","Basal","Intermediate",
 #                           "Secretory","Preciliated","Ciliated"))+
  scale_alpha_continuous("Absolute\nLog2 Fold Change\nFDR <= 0.05",breaks = c(.25, .5, 1),range=c(.25,1),limits=c(.25,1))
a

ggsave(paste("~/Box/COV_plot_and_files/Cell_Type_SARS2_dotPlot.CirclesByFC.v2.pdf",sep=""),height = 12,width = 12,units = "in",dpi = 300)

a = p +
  geom_point(data=plotMat[(plotMat$cellType.p_val_adj <= 0.1),],
             aes(size=pctExp*1.85,),pch=21, fill=NA, colour="black", stroke=1.25)+
  scale_y_discrete(limits=rev(orderByMarker)) +
  scale_x_discrete(limits=c("Proliferating_Basal","Basal","Intermediate",
                            "Secretory","Preciliated","Ciliated"))+
  scale_alpha_continuous("Absolute\nLog2 Fold Change\nFDR <= 0.05",breaks = c(.25, .5, 1),range=c(.25,1),limits=c(.25,1))
a

ggsave(paste("~/Box/COV_plot_and_files/Cell_Type_SARS2_dotPlot.Circles.v2.pdf",sep=""),height = 12,width = 12,units = "in",dpi = 300)



#### Cell type specific DOT PLOTS! ####

dotMatCellTypeTreatmentScaled = read.table("~/Box/COV_plot_and_files/dot_plot_data_file.txt",
                                           row.names = 1,header=T)



BasalMarkers_Final=                 c("COL17A1","BCAM","S100A2","IGFBP6")
SecretoryMarkers_Final =            c("BPIFA1","SCGB3A1","WFDC2","TCN1")
CiliatedMarkers_Final =           c("DNAAF1","DNAH7","TUBB4B","ERICH3")
PreciliatedMarkers_Final = c("FOXN4","MCIDAS","E2F7","PLK4")
IntermediateMarkers_Final = c("KRT13","SERPINB13","SPRR1B","CLCA4")
Proliferating_BasalMarkers_Final = c("CENPF","TOP2A","MKI67","ASPM")

markerList = c("CENPF","TOP2A","MKI67","ASPM",
               "COL17A1","BCAM","S100A2","IGFBP6",
               "KRT13","SERPINB13","SPRR1B","CLCA4",
               "BPIFA1","SCGB3A1","KRT7","CEACAM6",
               #"FOXN4","MCIDAS","E2F7","PLK4",
               "DNAAF1","DNAH7","TUBB4B","ERICH3")


dotMatCellTypeTreatmentScaled_cellType = dotMatCellTypeTreatmentScaled[(dotMatCellTypeTreatmentScaled$gene %in% markerList),]

plotMat2 = dotMatCellTypeTreatmentScaled_cellType


plotMat2 = plotMat2[!(is.na(plotMat2$gene)),]

plotMat2 = plotMat2[plotMat2$Treatment == "UT",]
plotMat2 = plotMat2[!(is.na(plotMat2$gene)),]
head(plotMat2[is.na(plotMat2$cellType),])

mat = acast(plotMat2, Gene ~ cellType,value.var =  "aveExp")

mat = mat[,-c(7)] # removing Ionocytes

dend = as.dendrogram(hclust(d = dist(x = mat)))
geneOrder = dend %>% labels

#SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut_noION = SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut_noION[
#                                                                order(SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut_noION$cellType.avg_logFC,decreasing = T),]
#
#SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut_noION = SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut_noION[
#                                                                !duplicated(SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut_noION$Gene),]
#
#SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut_noION = SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut_noION[
#                                                                order(SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut_noION$cellType),]
#orderByMarker = as.vector(SARS_dotMatCellTypeTreatmentScaled_withDE_DE_cut_noION$Gene)
#orderByMarker = orderByMarker[!duplicated(orderByMarker)]



plotMat2$gene = factor(plotMat2$gene, levels = rev(markerList))
head(plotMat2[is.na(plotMat2$cellType),])
plotMat2 = plotMat2[(plotMat2$cellType != "Ionocyte"),]
plotMat2 = plotMat2[(plotMat2$cellType != "Preciliated"),]

plotMat2$cellType = factor(plotMat2$cellType, levels = c("Proliferating_Basal","Basal","Intermediate",
                                                           "Secretory",
                                                         #"Preciliated",
                                                         "Ciliated"))
plotMat2 = plotMat2[!(is.na(plotMat2$gene)),]

plotMat2$plot_aveExp = plotMat2$aveExp
plotMat2$plot_aveExp[plotMat2$plot_aveExp >3] = 3
#plotMat2$plot_aveExp[plotMat2$plot_aveExp < -1] = -1

ggplot(plotMat2, aes(y = gene, x = cellType)) +
  geom_point(aes(size = pctExp, colour=plot_aveExp)) +
  scale_size_area() +
  scale_colour_gradientn(colours = c("grey","red")) +#, 
                         #range=c(0,4)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1,color="black",size=12),
        axis.title=element_text(size=16), 
        axis.text.y = element_text(color="black",size=12),
        legend.text=element_text(size=10),legend.title=element_text(size=14),
        panel.grid = element_blank(),plot.title = element_text(hjust = 0.5)) +
  labs(size="Percent\nExpressed", colour="Average\nExpression") + ylab("Cell Type") +
  guides(color = guide_colorbar(order = 1),size = guide_legend(order=2)) +
  geom_hline(aes(yintercept=4.5),color="gray") +
  geom_hline(aes(yintercept=8.5),color="gray") +
  geom_hline(aes(yintercept=12.5),color="gray") +
  geom_hline(aes(yintercept=16.5),color="gray")

ggsave(paste("~/Box/Erle group presentations/LRB/SARS-CoV-2 paper/Walter_COV_2_figures/Cell_Type_Marker.dotPlot.pdf",sep=""),height = 8,width = 8,units = "in",dpi = 300)

ggplot(plotMat2, aes(y = gene, x = cellType)) +
  geom_point(aes(size = pctExp, colour=plot_aveExp)) +
  scale_size_area() +
  scale_colour_gradientn(colours = c("grey","red")) +#, 
  #range=c(0,4)) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1,color="black",size=12),
        axis.title=element_text(size=16), 
        axis.text.y = element_text(color="black",size=12),
        legend.text=element_text(size=10),legend.title=element_text(size=14),
        panel.grid = element_blank(),plot.title = element_text(hjust = 0.5)) +
  labs(size="Percent\nExpressed", colour="Average\nExpression") + ylab("Cell Type") +
  guides(color = guide_colorbar(order = 1),size = guide_legend(order=2)) 
ggsave(paste("~/Box/Erle group presentations/LRB/SARS-CoV-2 paper/Walter_COV_2_figures/Cell_Type_Marker.dotPlot.noHorzLines.pdf",sep=""),height = 8,width = 8,units = "in",dpi = 300)

#### IL13 response Dot plots ####




#### Old indvidual genes ####

dir.create("~/Box/COV_plot_and_files/")

Gene = "ACE2"


SARS_dotMatCellTypeTreatmentScaled_withDE_Gene = SARS_dotMatCellTypeTreatmentScaled_withDE[SARS_dotMatCellTypeTreatmentScaled_withDE$Gene == Gene,]

SARS_dotMatCellTypeTreatmentScaled_withDE_Gene = SARS_dotMatCellTypeTreatmentScaled_withDE_Gene[!(SARS_dotMatCellTypeTreatmentScaled_withDE_Gene$cellType =="Ionocyte"),]

SARS_dotMatCellTypeTreatmentScaled_withDE_Gene = SARS_dotMatCellTypeTreatmentScaled_withDE_Gene[!(is.na(SARS_dotMatCellTypeTreatmentScaled_withDE_Gene$Gene)),]

SARS_dotMatCellTypeTreatmentScaled_withDE_Gene$cellType = factor(SARS_dotMatCellTypeTreatmentScaled_withDE_Gene$cellType,levels=c("Proliferating_Basal","Basal","Intermediate",
                                                                                                                                  "Secretory","Preciliated","Ciliated"))

ggplot(SARS_dotMatCellTypeTreatmentScaled_withDE_Gene, aes(y = cellType, x = Treatment)) +
  geom_point(aes(size = pctExp, colour=sclAveExp)) +
  geom_point(data=SARS_dotMatCellTypeTreatmentScaled_withDE_Gene[(SARS_dotMatCellTypeTreatmentScaled_withDE_Gene$cellType.p_val_adj <= 0.1 & 
                                                                    SARS_dotMatCellTypeTreatmentScaled_withDE_Gene$cellType.avg_logFC >0) ,],
             aes(size=pctExp*1.85),pch=21, fill=NA, colour="black", stroke=1.5) +
  geom_point(data=SARS_dotMatCellTypeTreatmentScaled_withDE_Gene[(SARS_dotMatCellTypeTreatmentScaled_withDE_Gene$cellTypeXtreatment.p_val_adj <= 0.1 &
                                                                    SARS_dotMatCellTypeTreatmentScaled_withDE_Gene$cellTypeXtreatment.logFC >0),],
             aes(size=pctExp*1.85),pch=21, fill=NA, colour="red", stroke=1.5) +
  geom_point(data=SARS_dotMatCellTypeTreatmentScaled_withDE_Gene[(SARS_dotMatCellTypeTreatmentScaled_withDE_Gene$cellTypeXtreatment.p_val_adj <= 0.1 &
                                                                    SARS_dotMatCellTypeTreatmentScaled_withDE_Gene$cellTypeXtreatment.logFC <0),],
             aes(size=pctExp*1.85),pch=21, fill=NA, colour="royalblue2", stroke=1.5) +
  scale_size_area() +
  scale_colour_gradientn(colours = c("purple","grey","orange"), 
                         limits=c(-1*abs(max(SARS_dotMatCellTypeTreatmentScaled_withDE_Gene$sclAveExp)),
                                  abs(max(SARS_dotMatCellTypeTreatmentScaled_withDE_Gene$sclAveExp)))) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1,color="black",size=12),
        axis.title=element_text(size=16), 
        axis.text.y = element_text(color="black",size=12),
        legend.text=element_text(size=10),legend.title=element_text(size=14),
        panel.grid = element_blank(),plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(limits=c("UT","IL17","IFN","IL13","CO")) +
  scale_y_discrete(limits=rev(c("Proliferating_Basal","Basal","Intermediate",
                                "Secretory","Preciliated","Ciliated")))+
  labs(size="Percent\nExpressed", colour="Average\nExpression") + ylab("Cell Type") +
  guides(color = guide_colorbar(order = 1),size = guide_legend(order=2)) +
  ggtitle(Gene)
ggsave(paste("~/Box/Erle group presentations/LRB/SARS-CoV-2 paper/Curtain.",Gene,".twoColor.DE_noted.pdf",sep=""),height = 5,width = 7,units = "in",dpi = 300)


ggplot(SARS_dotMatCellTypeTreatmentScaled_withDE_Gene, aes(y = cellType, x = Treatment)) +
  geom_point(aes(size = pctExp, colour=aveExp)) +
  geom_point(data=SARS_dotMatCellTypeTreatmentScaled_withDE_Gene[(SARS_dotMatCellTypeTreatmentScaled_withDE_Gene$cellType.p_val_adj <= 0.1 & 
                                                                    SARS_dotMatCellTypeTreatmentScaled_withDE_Gene$cellType.avg_logFC >0) ,],
             aes(size=pctExp*1.85),pch=21, fill=NA, colour="black", stroke=1.5) +
  geom_point(data=SARS_dotMatCellTypeTreatmentScaled_withDE_Gene[(SARS_dotMatCellTypeTreatmentScaled_withDE_Gene$cellTypeXtreatment.p_val_adj <= 0.1),],
             aes(size=pctExp*1.85),pch=21, fill=NA, colour="black", stroke=1.5) +
  scale_size_area() +
  scale_colour_gradientn(colours = c("grey","red"), 
                         limits=c(0,abs(max(SARS_dotMatCellTypeTreatmentScaled_withDE_Gene$aveExp)))) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1,color="black",size=12),
        axis.title=element_text(size=16), 
        axis.text.y = element_text(color="black",size=12),
        legend.text=element_text(size=10),legend.title=element_text(size=14),
        panel.grid = element_blank(),plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(limits=c("UT","IL17","IFN","IL13","CO")) +
  scale_y_discrete(limits=rev(c("Proliferating_Basal","Basal","Intermediate",
                                "Secretory","Preciliated","Ciliated")))+
  labs(size="Percent\nExpressed", colour="Average\nExpression") + ylab("Cell Type") +
  guides(color = guide_colorbar(order = 1),size = guide_legend(order=2)) +
  ggtitle(Gene)
ggsave(paste("~/Box/Erle group presentations/LRB/SARS-CoV-2 paper/Curtain.",Gene,".oneColor.DE_noted.pdf",sep=""),height = 5,width = 7,units = "in",dpi = 300)

#### SARS Grouped genes ####

SARS_dotMatCellTypeTreatmentScaled_withDE_UT = SARS_dotMatCellTypeTreatmentScaled_withDE[SARS_dotMatCellTypeTreatmentScaled_withDE$Treatment == "UT",]


for (i in unique(SARS_dotMatCellTypeTreatmentScaled_withDE_UT$cellType)){
  tmp_data = SARS_dotMatCellTypeTreatmentScaled_withDE_UT[SARS_dotMatCellTypeTreatmentScaled_withDE_UT$cellType == i,]
  print(ggplot(tmp_data, aes(x=Bait,y=cellType.avg_logFC)) + geom_point(fill=""))
}
ggplot(SARS_dotMatCellTypeTreatmentScaled_withDE_UT, aes(x=Bait,y=cellType.avg_logFC,
                            fill=cellType,color=cellType,group = interaction(Bait,cellType))) + 
  geom_point(aes(color=cellType))

SARS_dotMatCellTypeTreatmentScaled_withDE_UT$cellType_DE = "Not_DE"
SARS_dotMatCellTypeTreatmentScaled_withDE_UT$cellType_DE[(SARS_dotMatCellTypeTreatmentScaled_withDE_UT$cellType.p_val_adj <= 0.1)] <- "DE"

SARS_dotMatCellTypeTreatmentScaled_withDE[SARS_dotMatCellTypeTreatmentScaled_withDE$Type == "Proteases", "Bait"] <- "Proteases"

SARS_dotMatCellTypeTreatmentScaled_withDE_UT = SARS_dotMatCellTypeTreatmentScaled_withDE_UT[order(SARS_dotMatCellTypeTreatmentScaled_withDE_UT$cellType.p_val_adj,
                                                                                                  decreasing = T),]

SARS_dotMatCellTypeTreatmentScaled_withDE_UT = SARS_dotMatCellTypeTreatmentScaled_withDE_UT[!(is.na(SARS_dotMatCellTypeTreatmentScaled_withDE_UT$cellType.avg_logFC)),]
#newFactors = levels(SARS_dotMatCellTypeTreatmentScaled_withDE_UT$Bait)[-c(1)]
#SARS_dotMatCellTypeTreatmentScaled_withDE_UT$Bait = factor(SARS_dotMatCellTypeTreatmentScaled_withDE_UT$Bait,levels = newFactors)

p = ggplot(SARS_dotMatCellTypeTreatmentScaled_withDE_UT, aes(x=Bait,y=cellType.avg_logFC,fill=cellType_DE)) + 
  geom_dotplot(binaxis = "y", stackdir = "center",binwidth = .15) +
  scale_fill_manual(values = c("red", "black")) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1,color="black",size=12),
        axis.title=element_text(size=16), 
        axis.text.y = element_text(color="black",size=12),
        legend.text=element_text(size=10),legend.title=element_text(size=14)) +
  ylab("log2 Fold Change")
p + facet_wrap(~cellType,ncol = 1)
ggsave(paste("~/Box/Erle group presentations/LRB/SARS-CoV-2 paper/CellType_Markers.FC_by_cellType.plot.pdf",sep=""),height = 9,width = 10,units = "in",dpi = 300)

#### Individual Genes ####


Gene = "CTSL"  #### For individual gene plots, change this.


dotMatCellTypeTreatmentScaled = read.table("~/Box/COV_plot_and_files/dot_plot_data_file.txt",
                                           row.names = 1,header=T)


dotMatCellTypeTreatmentScaled_Gene = dotMatCellTypeTreatmentScaled[dotMatCellTypeTreatmentScaled$gene == Gene,]

dir.create("~/Box/COV_plot_and_files/")

ggplot(dotMatCellTypeTreatmentScaled_Gene, aes(y = cellType, x = Treatment)) +
  geom_point(aes(size = pctExp, colour=sclAveExp)) +
  scale_size_area() +
  scale_colour_gradientn(colours = c("purple","grey","orange"), 
                         limits=c(-1*abs(max(dotMatCellTypeTreatmentScaled_Gene$sclAveExp)),
                                  abs(max(dotMatCellTypeTreatmentScaled_Gene$sclAveExp)))) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1,color="black",size=12),
        axis.title=element_text(size=16), 
        axis.text.y = element_text(color="black",size=12),
        legend.text=element_text(size=10),legend.title=element_text(size=14),
        panel.grid = element_blank(),plot.title = element_text(hjust = 0.5)) +
  labs(size="Percent\nExpressed", colour="Average\nExpression") + ylab("Cell Type") +
  guides(color = guide_colorbar(order = 1),size = guide_legend(order=2)) +
  ggtitle(Gene)
ggsave(paste("~/Box/COV_plot_and_files//Curtain.",Gene,".twoColor.pdf",sep=""),height = 5,width = 7,units = "in",dpi = 300)


ggplot(dotMatCellTypeTreatmentScaled_Gene, aes(y = cellType, x = Treatment)) +
  geom_point(aes(size = pctExp, colour=aveExp)) +
  scale_size_area() +
  scale_colour_gradientn(colours = c("grey","red"), 
                         limits=c(0,abs(max(dotMatCellTypeTreatmentScaled_Gene$aveExp)))) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1,color="black",size=12),
        axis.title=element_text(size=16), 
        axis.text.y = element_text(color="black",size=12),
        legend.text=element_text(size=10),legend.title=element_text(size=14),
        panel.grid = element_blank(),plot.title = element_text(hjust = 0.5)) +
  labs(size="Percent\nExpressed", colour="Average\nExpression") + ylab("Cell Type") +
  guides(color = guide_colorbar(order = 1),size = guide_legend(order=2)) +
  ggtitle(Gene)
ggsave(paste("~/Box/COV_plot_and_files//Curtain.",Gene,".oneColor.pdf",sep=""),height = 5,width = 7,units = "in",dpi = 300)


#### Dotplot multi individual ####


Gene = "CTSL"  #### For individual gene plots, change this.


dotMatCellTypeTreatmentIndScaled = read.table("~/Box/COV_plot_and_files/dot_plot_data_file_split_inds.txt",
                                           row.names = 1,header=T)


dotMatCellTypeTreatmentIndScaled_Gene = dotMatCellTypeTreatmentIndScaled[dotMatCellTypeTreatmentIndScaled$gene == Gene,]


ggplot(dotMatCellTypeTreatmentIndScaled_Gene, aes(y = cellType, x = Treatment, group=Subject)) +
  geom_point(aes(size = pctExp, colour=sclAveExp),position = position_dodge(width = 0.75)) +
  scale_size_area() +
  scale_colour_gradientn(colours = c("purple","grey","orange"), 
                         limits=c(-1*abs(max(dotMatCellTypeTreatmentIndScaled_Gene$sclAveExp)),
                                  abs(max(dotMatCellTypeTreatmentIndScaled_Gene$sclAveExp)))) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1,color="black",size=12),
        axis.title=element_text(size=16), 
        axis.text.y = element_text(color="black",size=12),
        legend.text=element_text(size=10),legend.title=element_text(size=14),
        panel.grid = element_blank(),plot.title = element_text(hjust = 0.5)) +
  labs(size="Percent\nExpressed", colour="Average\nExpression") + ylab("Cell Type") +
  guides(color = guide_colorbar(order = 1),size = guide_legend(order=2)) +
  ggtitle(Gene)
ggsave(paste("~/Box/COV_plot_and_files//Curtain.",Gene,".twoColor.withReps.pdf",sep=""),height = 4,width = 7,units = "in",dpi = 300)


ggplot(dotMatCellTypeTreatmentIndScaled_Gene, aes(y = cellType, x = Treatment, group=Subject)) +
  geom_point(aes(size = pctExp, colour=aveExp),position = position_dodge(width = 0.75)) +
  scale_size_area() +
  scale_colour_gradientn(colours = c("grey","red"), 
                         limits=c(0,abs(max(dotMatCellTypeTreatmentIndScaled_Gene$aveExp)))) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1,color="black",size=12),
        axis.title=element_text(size=16), 
        axis.text.y = element_text(color="black",size=12),
        legend.text=element_text(size=10),legend.title=element_text(size=14),
        panel.grid = element_blank(),plot.title = element_text(hjust = 0.5)) +
  labs(size="Percent\nExpressed", colour="Average\nExpression") + ylab("Cell Type") +
  guides(color = guide_colorbar(order = 1),size = guide_legend(order=2)) +
  ggtitle(Gene)
ggsave(paste("~/Box/COV_plot_and_files//Curtain.",Gene,".oneColor.withReps.pdf",sep=""),height = 4,width = 7,units = "in",dpi = 300)




#### Multiple Genes in UT (or other combination of celltypes) ####

#### first without replicates

#write.table(SARS_COV_2_table,"~/Box/COV_plot_and_files/SARS_COV2_Interacti",sep="\t",quote=F,row.names = F)

SARS_COV_2_table  = read.table("~/data/Erle_Bonser_0529/SARS_COV_2_interactions.txt",sep="\t",header = T)

SARS_COV_2_table_Spike = SARS_COV_2_table[SARS_COV_2_table$Viral_Protein == "Spike",]  #### Reducing this list for plotting

Genes = SARS_COV_2_table_Spike$Gene  # Can be changed to any gene list.

Treatments = c("UT") # make any combination, another example with IL13 will be below.

Genes = c("ACE2","CTSL","TMPRSS2")
dotMatCellTypeTreatmentScaled_Genes = dotMatCellTypeTreatmentScaled[dotMatCellTypeTreatmentIndScaled$gene %in% Genes,]

dotMatCellTypeTreatmentScaled_Genes_Treat = dotMatCellTypeTreatmentScaled_Genes[(dotMatCellTypeTreatmentScaled_Genes$Treatment %in% Treatments),]

FileLabel = "Spike_Genes_UT_celltypes" # make a unique file label to describe what you're doing

ggplot(dotMatCellTypeTreatmentScaled_Genes_Treat, aes(y = cellType, x = gene)) +
  geom_point(aes(size = pctExp, colour=sclAveExp)) +
  scale_size_area() +
  scale_colour_gradientn(colours = c("purple","grey","orange"), 
                         limits=c(-1*abs(max(dotMatCellTypeTreatmentScaled_Genes_Treat$sclAveExp)),
                                  abs(max(dotMatCellTypeTreatmentScaled_Genes_Treat$sclAveExp)))) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1,color="black",size=12),
        axis.title=element_text(size=16), 
        axis.text.y = element_text(color="black",size=12),
        legend.text=element_text(size=10),legend.title=element_text(size=14),
        panel.grid = element_blank(),plot.title = element_text(hjust = 0.5)) +
  labs(size="Percent\nExpressed", colour="Average\nExpression") + ylab("Cell Type") +
  guides(color = guide_colorbar(order = 1),size = guide_legend(order=2))
ggsave(paste("~/Box/COV_plot_and_files//Curtain.",FileLabel,".twoColor.pdf",sep=""),height = 4,width = 7,units = "in",dpi = 300)


ggplot(dotMatCellTypeTreatmentScaled_Genes_Treat, aes(y = cellType, x = gene)) +
  geom_point(aes(size = pctExp, colour=aveExp)) +
  scale_size_area() +
  scale_colour_gradientn(colours = c("grey","red"), 
                         limits=c(0,abs(max(dotMatCellTypeTreatmentScaled_Genes_Treat$aveExp)))) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1,color="black",size=12),
        axis.title=element_text(size=16), 
        axis.text.y = element_text(color="black",size=12),
        legend.text=element_text(size=10),legend.title=element_text(size=14),
        panel.grid = element_blank(),plot.title = element_text(hjust = 0.5)) +
  labs(size="Percent\nExpressed", colour="Average\nExpression") + ylab("Cell Type") +
  guides(color = guide_colorbar(order = 1),size = guide_legend(order=2))
ggsave(paste("~/Box/COV_plot_and_files//Curtain.",FileLabel,".oneColor.pdf",sep=""),height = 4,width = 7,units = "in",dpi = 300)


#### with replicates

dotMatCellTypeTreatmentIndScaled_Genes = dotMatCellTypeTreatmentIndScaled[dotMatCellTypeTreatmentIndScaled$gene %in% Genes,]

dotMatCellTypeTreatmentIndScaled_Genes_Treat = dotMatCellTypeTreatmentIndScaled_Genes[(dotMatCellTypeTreatmentIndScaled_Genes$Treatment %in% Treatments),]

FileLabel = "Spike_Genes_UT_celltypes" # make a unique file label to describe what you're doing

ggplot(dotMatCellTypeTreatmentIndScaled_Genes_Treat, aes(y = cellType, x = gene, group=Subject)) +
  geom_point(aes(size = pctExp, colour=sclAveExp),position = position_dodge(width = 0.75)) +
  scale_size_area() +
  scale_colour_gradientn(colours = c("purple","grey","orange"), 
                         limits=c(-1*abs(max(dotMatCellTypeTreatmentIndScaled_Genes_Treat$sclAveExp)),
                                  abs(max(dotMatCellTypeTreatmentIndScaled_Genes_Treat$sclAveExp)))) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1,color="black",size=12),
        axis.title=element_text(size=16), 
        axis.text.y = element_text(color="black",size=12),
        legend.text=element_text(size=10),legend.title=element_text(size=14),
        panel.grid = element_blank(),plot.title = element_text(hjust = 0.5)) +
  labs(size="Percent\nExpressed", colour="Average\nExpression") + ylab("Cell Type") +
  guides(color = guide_colorbar(order = 1),size = guide_legend(order=2)) 
ggsave(paste("~/Box/COV_plot_and_files//Curtain.",FileLabel,".twoColor.withReps.pdf",sep=""),height = 4,width = 7,units = "in",dpi = 300)


ggplot(dotMatCellTypeTreatmentIndScaled_Genes_Treat, aes(y = cellType, x = gene, group=Subject)) +
  geom_point(aes(size = pctExp, colour=aveExp),position = position_dodge(width = 0.75)) +
  scale_size_area() +
  scale_colour_gradientn(colours = c("grey","red"), 
                         limits=c(0,abs(max(dotMatCellTypeTreatmentIndScaled_Genes_Treat$aveExp)))) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1,color="black",size=12),
        axis.title=element_text(size=16), 
        axis.text.y = element_text(color="black",size=12),
        legend.text=element_text(size=10),legend.title=element_text(size=14),
        panel.grid = element_blank(),plot.title = element_text(hjust = 0.5)) +
  labs(size="Percent\nExpressed", colour="Average\nExpression") + ylab("Cell Type") +
  guides(color = guide_colorbar(order = 1),size = guide_legend(order=2))
ggsave(paste("~/Box/COV_plot_and_files//Curtain.",FileLabel,".oneColor.withReps.pdf",sep=""),height = 4,width = 7,units = "in",dpi = 300)

#### Now with multiple treatments ####


SARS_COV_2_table_Spike = SARS_COV_2_table[SARS_COV_2_table$Viral_Protein == "Spike",]  #### Reducing this list for plotting

Genes = SARS_COV_2_table_Spike$Gene  # Can be changed to any gene list.

Treatments = c("UT","IL13","IFN") # make any combination


dotMatCellTypeTreatmentScaled_Genes = dotMatCellTypeTreatmentScaled[dotMatCellTypeTreatmentIndScaled$gene %in% Genes,]

dotMatCellTypeTreatmentScaled_Genes_Treat = dotMatCellTypeTreatmentScaled_Genes[(dotMatCellTypeTreatmentScaled_Genes$Treatment %in% Treatments),]

FileLabel = "Spike_Genes_UT_celltypes" # make a unique file label to describe what you're doing

dotMatCellTypeTreatmentScaled_Genes_Treat$cellType_Treatment = paste(dotMatCellTypeTreatmentScaled_Genes_Treat$cellType, 
                                                                     dotMatCellTypeTreatmentScaled_Genes_Treat$Treatment,sep="_")

dotMatCellTypeTreatmentScaled_Genes_Treat$cellType_Treatment = factor(dotMatCellTypeTreatmentScaled_Genes_Treat$cellType_Treatment, 
                                                                     levels=rev(unique(dotMatCellTypeTreatmentScaled_Genes_Treat$cellType_Treatment)))

ggplot(dotMatCellTypeTreatmentScaled_Genes_Treat, aes(y = cellType_Treatment, x = gene,colour=sclAveExp)) +
  geom_point(aes(size = pctExp)) +
  scale_size_area() +
  scale_colour_gradientn(colours = c("purple","grey","orange"), 
                         limits=c(-1*abs(max(dotMatCellTypeTreatmentScaled_Genes_Treat$sclAveExp)),
                                  abs(max(dotMatCellTypeTreatmentScaled_Genes_Treat$sclAveExp)))) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1,color="black",size=12),
        axis.title=element_text(size=16), 
        axis.text.y = element_text(color="black",size=12),
        legend.text=element_text(size=10),legend.title=element_text(size=14),
        panel.grid = element_blank(),plot.title = element_text(hjust = 0.5)) +
  labs(size="Percent\nExpressed", colour="Average\nExpression") + ylab("Cell Type") +
  guides(color = guide_colorbar(order = 1),size = guide_legend(order=2))
ggsave(paste("~/Box/COV_plot_and_files/Curtain.",FileLabel,".scaledTwoColor.pdf",sep=""),
       height = 4,width = 7,units = "in",dpi = 300)


ggplot(dotMatCellTypeTreatmentScaled_Genes_Treat, aes(y = cellType_Treatment, x = gene,colour=aveExp)) +
  geom_point(aes(size = pctExp)) +
  scale_size_area() +
  scale_colour_gradientn(colours = c("grey","red"), 
                         limits=c(0,
                                  abs(max(dotMatCellTypeTreatmentScaled_Genes_Treat$aveExp)))) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1,color="black",size=12),
        axis.title=element_text(size=16), 
        axis.text.y = element_text(color="black",size=12),
        legend.text=element_text(size=10),legend.title=element_text(size=14),
        panel.grid = element_blank(),plot.title = element_text(hjust = 0.5)) +
  labs(size="Percent\nExpressed", colour="Average\nExpression") + ylab("Cell Type") +
  guides(color = guide_colorbar(order = 1),size = guide_legend(order=2))
ggsave(paste("~/Box/COV_plot_and_files//Curtain.",FileLabel,".multiGene.oneColor.pdf",sep=""),height = 4,width = 7,units = "in",dpi = 300)


#### with replicates

dotMatCellTypeTreatmentIndScaled_Genes = dotMatCellTypeTreatmentIndScaled[dotMatCellTypeTreatmentIndScaled$gene %in% Genes,]

dotMatCellTypeTreatmentIndScaled_Genes_Treat = dotMatCellTypeTreatmentIndScaled_Genes[(dotMatCellTypeTreatmentIndScaled_Genes$Treatment %in% Treatments),]

FileLabel = "Spike_Genes_UT_celltypes" # make a unique file label to describe what you're doing

ggplot(dotMatCellTypeTreatmentIndScaled_Genes_Treat, aes(y = cellType, x = gene, group=Subject)) +
  geom_point(aes(size = pctExp, colour=sclAveExp),position = position_dodge(width = 0.75)) +
  scale_size_area() +
  scale_colour_gradientn(colours = c("purple","grey","orange"), 
                         limits=c(-1*abs(max(dotMatCellTypeTreatmentIndScaled_Genes_Treat$sclAveExp)),
                                  abs(max(dotMatCellTypeTreatmentIndScaled_Genes_Treat$sclAveExp)))) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1,color="black",size=12),
        axis.title=element_text(size=16), 
        axis.text.y = element_text(color="black",size=12),
        legend.text=element_text(size=10),legend.title=element_text(size=14),
        panel.grid = element_blank(),plot.title = element_text(hjust = 0.5)) +
  labs(size="Percent\nExpressed", colour="Average\nExpression") + ylab("Cell Type") +
  guides(color = guide_colorbar(order = 1),size = guide_legend(order=2)) 
ggsave(paste("~/Box/COV_plot_and_files//Curtain.",FileLabel,".twoColor.withReps.pdf",sep=""),height = 4,width = 7,units = "in",dpi = 300)


ggplot(dotMatCellTypeTreatmentIndScaled_Genes_Treat, aes(y = cellType, x = gene)) +
  geom_point(aes(size = pctExp, colour=aveExp, group=Subject),position = position_dodge(width = 0.75)) +
  scale_size_area() +
  scale_colour_gradientn(colours = c("grey","red"), 
                         limits=c(0,abs(max(dotMatCellTypeTreatmentIndScaled_Genes_Treat$aveExp)))) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1,color="black",size=12),
        axis.title=element_text(size=16), 
        axis.text.y = element_text(color="black",size=12),
        legend.text=element_text(size=10),legend.title=element_text(size=14),
        panel.grid = element_blank(),plot.title = element_text(hjust = 0.5)) +
  labs(size="Percent\nExpressed", colour="Average\nExpression") + ylab("Cell Type") +
  guides(color = guide_colorbar(order = 1),size = guide_legend(order=2))
ggsave(paste("~/Box/COV_plot_and_files//Curtain.",FileLabel,".oneColor.withReps.pdf",sep=""),height = 4,width = 7,units = "in",dpi = 300)

#### SARS enrichments ####
#### OLD ####


#### COVID Enrichment Test ####

enrichment <- as.data.frame(setNames(replicate(3,numeric(0), simplify = F),letters[0:3]))
colnames(enrichment) = c("condition","pvalue","Obs","Exp")


SARS_COV_2_list_Ensembl = allgenes2[(allgenes2$Gene %in% SARS_COV_2_list),][,"Ensembl_ID"]
allgenes2_SARS = allgenes2[(allgenes2$Ensembl_ID %in% SARS_COV_2_list_Ensembl),]

genePerBin = c()

nTests = 1000

bins = seq(-1,5,.2)
for (i in seq(2,length(bins))){
  high = bins[i]
  low = bins[(i-1)]
  allgenesBin = allgenes2_SARS[(log10(allgenes2_SARS$Mean_Normalized_Counts) < high),]
  allgenesBin = allgenesBin[(log10(allgenesBin$Mean_Normalized_Counts) >= low),]
  print(nrow(allgenesBin))
  genePerBin = c(genePerBin,nrow(allgenesBin))
}

names(genePerBin) = bins[c(2:length(bins))]

for (i in comps){
  FDR_col = paste(i,"FDR",sep=".")
  Pval_col = paste(i,"RawP",sep=".")
  logFC_col = paste(i,"log2FC",sep=".")
  
  counter = 1
  counts = c()
  
  allgenes_tmp_DE = allgenes2[(allgenes2[,FDR_col] <= FDR_cut),]
  allgenes_tmp_DE = allgenes_tmp_DE[(allgenes_tmp_DE[,logFC_col] <= fc_cut),]
  allgenes_tmp_DE = allgenes_tmp_DE[!(is.na(allgenes_tmp_DE$Ensembl_ID)),]
  DE_count_total = nrow(allgenes_tmp_DE)
  allgenes_tmp_DE_in_list = allgenes_tmp_DE[allgenes_tmp_DE$Ensembl_ID %in% SARS_COV_2_list_Ensembl,]
  DE_count_in_list = nrow(allgenes_tmp_DE_in_list)
  
  counter = 1
  counts = c()
  
  
  for (j in seq(1,nTests)){
    genes = c()
    for (k in names(genePerBin)){
      high = as.numeric(k)
      low = high - 0.2
      nGeneInBin = as.vector(genePerBin[k])
      allgenesBin = allgenes2[(log10(allgenes2$Mean_Normalized_Counts) < high),]
      allgenesBin = allgenesBin[(log10(allgenesBin$Mean_Normalized_Counts) >= low),]
      allgenesBinShuff = allgenesBin[sample(nrow(allgenesBin)),]
      geneKs = allgenesBinShuff[c(0:nGeneInBin),"Ensembl_ID"]
      genes = c(genes,geneKs)
    }
    allgenes_tmp_genes = allgenes2[allgenes2$Ensembl_ID %in% genes,]
    allgenes_tmp_genes_DE = allgenes_tmp_genes[allgenes_tmp_genes[,FDR_col] <= FDR_cut,]
    allgenes_tmp_genes_DE = allgenes_tmp_genes_DE[(allgenes_tmp_genes_DE[,logFC_col] <= fc_cut),]
    counts = c(counts,nrow(allgenes_tmp_genes_DE))
    if (nrow(allgenes_tmp_genes_DE) >= DE_count_in_list){
      counter = counter + 1
    }
  }
  print(i)
  print(counter/(nTests + 1))
  print(DE_count_in_list)
  print(mean(counts))
  
  enrich_i = c(counter/(nTests + 1),DE_count_in_list,mean(counts))
  enrich_i = as.data.frame(enrich_i)
  enrich_i = t(as.data.frame(enrich_i))
  colnames(enrich_i) = c("pvalue","Obs","Exp")
  rownames(enrich_i) = i
  enrichment = rbind(enrichment,enrich_i)
  print(enrichment)
  
  cDF = as.data.frame(counts)
  binWidth = round((range(counts)[2] - range(counts)[1])/20,digits = 0)
  
  p = ggplot(cDF, aes(x=counts)) + 
    geom_histogram(binwidth = binWidth ,color="black", fill="white") + 
    theme_classic() + 
    geom_vline(xintercept = DE_count_in_list, linetype="dashed", color = "red", size=1.5)
  
  pdf(paste(i,"enrichment_plot.pdf",sep=""),height = 6,width=6)
  print(p)
  dev.off()
  
}

for (i in comps){
  FDR_col = paste(i,"FDR",sep=".")
  Pval_col = paste(i,"RawP",sep=".")
  logFC_col = paste(i,"log2FC",sep=".")
  
  allgenes_tmp = allgenes2[order(allgenes2$Mean_Normalized_Counts),]
  allgenes_tmp = allgenes_tmp[!(is.na(allgenes_tmp[,Pval_col])),]
  allgenes_tmp = allgenes_tmp[!(is.na(allgenes_tmp[,FDR_col])),]
  readCut = allgenes_tmp[1,"Mean_Normalized_Counts"]
  allgenes_tmp_cut = allgenes2[allgenes_tmp$Mean_Normalized_Counts >= readCut,]
  allgenes_tmp_cut = allgenes_tmp_cut[!(is.na(allgenes_tmp_cut$Gene)),]
  
  allgenes_tmp_DE = allgenes2[(allgenes2[,FDR_col] <= FDR_cut),]
  allgenes_tmp_DE = allgenes_tmp_DE[!(is.na(allgenes_tmp_DE$Ensembl_ID)),]
  DE_count_total = nrow(allgenes_tmp_DE)
  allgenes_tmp_DE_in_list = allgenes_tmp_DE[allgenes_tmp_DE$Ensembl_ID %in% SARS_COV_2_list_DE_Ens,]
  DE_count_in_list = nrow(allgenes_tmp_DE_in_list)
  
  counter = 1
  counts = c()
  
  for (j in seq(1,1000)){
    allgenes_tmp_cut_shuff = allgenes_tmp_cut[sample(nrow(allgenes_tmp_cut)),]
    allgenes_tmp_cut_shuff_de = head(allgenes_tmp_cut_shuff, n=DE_count_total)
    allgenes_tmp_cut_shuff_de_in_list = allgenes_tmp_cut_shuff_de[allgenes_tmp_cut_shuff_de$Ensembl_ID %in% SARS_COV_2_list_DE_Ens,]
    counts = c(counts,nrow(allgenes_tmp_cut_shuff_de_in_list))
    if (nrow(allgenes_tmp_cut_shuff_de_in_list) >= DE_count_in_list){
      counter = counter + 1
    }
  }
  print(i)
  print(counter/100001)
  print(nrow(allgenes_tmp_DE_in_list))
  print(mean(counts))
  
  enrich_i = c(counter/1001,nrow(allgenes_tmp_DE_in_list),mean(counts))
  enrich_i = as.data.frame(enrich_i)
  enrich_i = t(as.data.frame(enrich_i))
  colnames(enrich_i) = c("pvalue","Obs","Exp")
  rownames(enrich_i) = i
  enrichment = rbind(enrichment,enrich_i)
  print(enrichment)
}


SARS_COV_2_list_Entrez = allgenes2[(allgenes2$Gene %in% SARS_COV_2_list),][,"Entrez_ID"]

newList = c()
for (i in SARS_COV_2_list_Entrez){
  newString = strsplit(as.character(i),', ')
  newList = c(newList,newString[[1]][1])
}

SARS_COV_2_list_Entrez = newList


for (i in comps){
  FDR_col = paste(i,"FDR",sep=".")
  logFC_col = paste(i,"log2FC",sep=".")
  allgenes_tmp = allgenes2[(allgenes2[,FDR_col] <= FDR_cut),]
  allgenes_tmp = allgenes_tmp[!(is.na(allgenes_tmp$Ensembl_ID)),]
  allgenes_tmp = allgenes_tmp[!(is.na(allgenes_tmp$Entrez_ID)),]
  for (j in allgenes_tmp$Ensembl_ID){
    newRow = allgenes_tmp[(allgenes_tmp$Ensembl_ID == j),]
    EntID = newRow$Entrez_ID
    EntString = strsplit(as.character(EntID),', ')
    if (length(EntString[[1]]) > 1){
      newRowFinal = newRow
      for (k in seq(1,length(EntString[[1]])-1)){
        newRowFinal = rbind(newRowFinal,newRow)
        print(k)
        print(newRow)
      }
      newRow = newRowFinal
      allgenes_tmp = allgenes_tmp[!(allgenes_tmp$Ensembl_ID == j),]
      newRow$Entrez_ID = EntString[[1]]
      allgenes_tmp = rbind(allgenes_tmp,newRow)
    }
  }
  allgenes_tmp = allgenes_tmp[order(allgenes_tmp[,logFC_col]),]
  
}



allgenes_tmp = allgenes2[!(is.na(allgenes2$Ensembl_ID)),]
allgenes_tmp = allgenes_tmp[!(is.na(allgenes_tmp$Entrez_ID)),]
for (j in allgenes_tmp$Ensembl_ID){
  newRow = allgenes_tmp[(allgenes_tmp$Ensembl_ID == j),]
  EntID = newRow$Entrez_ID
  EntString = strsplit(as.character(EntID),', ')
  if (length(EntString[[1]]) > 1){
    allgenes_tmp = allgenes_tmp[!(allgenes_tmp$Ensembl_ID == j),]
    newRow$Entrez_ID = EntString[[1]][1]
    allgenes_tmp = rbind(allgenes_tmp,newRow)
  }
}

allgenes_tmp$abs = abs(allgenes_tmp$`condition-IFNavsUT.log2FC`)
allgenes_tmp = allgenes_tmp[order(allgenes_tmp$abs,decreasing = T),]
allgenes_tmp = allgenes_tmp[!(duplicated(allgenes_tmp$Entrez_ID)),]
allgenes_tmp = allgenes_tmp[order(allgenes_tmp$IFNa_Y_vs_N.log2FC),]

allgenes_tmp_entrez = allgenes_tmp

#allgenes_tmp = allgenes_tmp[!(is.na(allgenes_tmp$IFNa_Y_vs_N.FDR)),]

pathwaysList = c()
pathwaysList$`00000_SARS_CoV2_interactions` = SARS_COV_2_list_Entrez

ranks = allgenes_tmp$IFNa_Y_vs_N.log2FC
names(ranks) = allgenes_tmp$Entrez_ID
fgseaRes2 <- fgsea(pathways = pathwaysList, 
                   stats = ranks,
                   minSize=10,
                   maxSize=500,
                   nperm=10000)

plotEnrichment(pathwaysList[["00000_SARS_CoV2_interactions"]],
               ranks) + labs(title="SARS_CoV2_interactions")



for (i in comps){
  FDR_col = paste(i,"FDR",sep=".")
  logFC_col = paste(i,"log2FC",sep=".")
  allgenes_tmp2 = allgenes_tmp
  
  allgenes_tmp3 = allgenes_tmp2[allgenes_tmp$Entrez_ID %in% SARS_COV_2_list_Entrez,]
  allgenes_tmp3_DE = allgenes_tmp3[(allgenes_tmp3[,FDR_col] <= 0.1),]
  allgenes_tmp3_DE = allgenes_tmp3_DE[!(is.na(allgenes_tmp3_DE$Entrez_ID)),]
  minCount = min(allgenes_tmp3_DE[,"Mean_Normalized_Counts"])
  print(minCount)
  
  allgenes_tmp2$abs = abs(allgenes_tmp2[,logFC_col])
  allgenes_tmp2 = allgenes_tmp2[allgenes_tmp2$Mean_Normalized_Counts >= 42,]
  allgenes_tmp2 = allgenes_tmp2[order(allgenes_tmp2$abs,decreasing = T),]
  allgenes_tmp2 = allgenes_tmp2[!(duplicated(allgenes_tmp2$Entrez_ID)),]
  allgenes_tmp2 = allgenes_tmp2[order(allgenes_tmp2[,logFC_col]),]
  allgenes_tmp2 = allgenes_tmp2[!(is.na(allgenes_tmp2[,FDR_col])),]
  
  
  ranks = allgenes_tmp2[,logFC_col]
  names(ranks) = allgenes_tmp2$Entrez_ID
  fgseaRes <- fgsea(pathways = pathwaysList, 
                    stats = ranks,
                    minSize=10,
                    maxSize=500,
                    nperm=10000)
  print(i)
  print(fgseaRes)
  
  print(plotEnrichment(pathwaysList[["00000_SARS_CoV2_interactions"]], ranks) + 
          labs(title=paste("SARS_CoV2_interactions : ",i,sep="")))
  
}


SARS_COV_2_list_Ensembl = allgenes2[(allgenes2$Gene %in% SARS_COV_2_list),][,"Ensembl_ID"]

pathwaysList = c()
pathwaysList$`00000_SARS_CoV2_interactions` = SARS_COV_2_list_Ensembl


for (i in comps){
  FDR_col = paste(i,"FDR",sep=".")
  logFC_col = paste(i,"log2FC",sep=".")
  allgenes_tmp2 = allgenes2
  
  allgenes_tmp3 = allgenes_tmp2[allgenes_tmp$Ensembl_ID %in% SARS_COV_2_list_Ensembl,]
  allgenes_tmp3_DE = allgenes_tmp3[(allgenes_tmp3[,FDR_col] <= 0.1),]
  allgenes_tmp3_DE = allgenes_tmp3_DE[!(is.na(allgenes_tmp3_DE$Ensembl_ID)),]
  minCount = min(allgenes_tmp3_DE[,"Mean_Normalized_Counts"])
  print(minCount)
  
  allgenes_tmp2$abs = abs(allgenes_tmp2[,logFC_col])
  allgenes_tmp2 = allgenes_tmp2[allgenes_tmp2$Mean_Normalized_Counts >= 30,]
  allgenes_tmp2 = allgenes_tmp2[order(allgenes_tmp2$abs,decreasing = T),]
  #allgenes_tmp2 = allgenes_tmp2[!(duplicated(allgenes_tmp2$Entrez_ID)),]
  allgenes_tmp2 = allgenes_tmp2[order(allgenes_tmp2[,FDR_col]),]
  allgenes_tmp2 = allgenes_tmp2[!(is.na(allgenes_tmp2[,FDR_col])),]
  
  
  ranks = allgenes_tmp2[,FDR_col]
  names(ranks) = allgenes_tmp2$Ensembl_ID
  fgseaRes <- fgsea(pathways = pathwaysList, 
                    stats = ranks,
                    minSize=10,
                    maxSize=500,
                    nperm=10000)
  print(i)
  print(fgseaRes)
  
  print(plotEnrichment(pathwaysList[["00000_SARS_CoV2_interactions"]], ranks) + 
          labs(title=paste("SARS_CoV2_interactions : ",i,sep="")))
  
}



#### SARS heatmaps ####
#### OLD ####


####

#syms <- unique(syms)
genes <- annot[SARS_COV_2_list_DE_Ens,]$Gene
M <- stabilizedBR[SARS_COV_2_list_DE_Ens,]
#UT, IL17, IFNa, IFNg, IL13, CO
M = M[,c(31:36,25:30,7:12,13:18,19:24,1:6)]
#M <- M[,c(1:12,19:24,31:36)]
#M <-subset(M,,names)
#M <- M - rowMeans(M, na.rm=T)

#DE.heatmap(M, dendrogram="both", fact=conds[,c("condition"), drop=F], zlim =c(-2,2), lhei=c(1,.25,5,1), height = 1200, width=1200, main=paste("DE SARS-CoV-2 interacting genes - FDR < ", FDR_cut,sep=""), file.prefix="DE/Heatmaps/DE_heatmap_no_symbols.SARS", breaks=.2, Symbols=NA, min.col = "blue", mid.col = "white", max.col = "red", margins = c(12,12,12,12))


M = merge(annot,M, by.x="Ensembl_ID",by.y='row.names')
rownames(M) = M$Gene
M = M[,seq(13,48)]

M = M[order(rownames(M)),]

M_UT = M[,c(1:6)]

M <- M - rowMeans(M_UT, na.rm=T)

targets <- as.data.frame(as.matrix(colnames(M)))

conds <- cbind(targets, do.call(rbind, strsplit(as.character(targets$V1),'_')))
rownames(conds) <- conds[,c(1)]
conds <- conds[,-c(1)]
condsSmall = conds
condsSmall = condsSmall[,-c(2)]
condsSmall = as.data.frame(condsSmall)
colnames(condsSmall) = "condition"


ha_c = HeatmapAnnotation(df = condsSmall, 
                         col = list(condition = cols),
                         show_legend = T)

VP = as.data.frame(allgenes_SARS2_tmp[,c("Type","Viral_Protein")])
rownames(VP) = allgenes_SARS2_tmp$Gene
colnames(VP) = c("Protein_Class","Viral_Protein")

#VP = VP[order(rownames(VP)),]

VP$Viral_Protein = factor(VP$Viral_Protein)

cols_VP <- gg_color_hue(length(unique(VP$Viral_Protein)))
names(cols_VP) = levels(VP$Viral_Protein)

cols_VP_2 <- gg_color_hue(length(unique(VP$Protein_Class)))
VP$Protein_Class = factor(VP$Protein_Class,levels = c("Envelope","Membrane","Nucleocapsid","NSP","ORF","Spike"))
names(cols_VP_2) = levels(VP$Protein_Class)

VP = VP[rownames(VP) %in% rownames(M),]

VP$Protein_Class = as.vector(VP$Protein_Class)
VP$Protein_Class[rownames(VP) == "CTSB"] = "Proteases"
VP$Protein_Class[rownames(VP) == "PCSK7"] = "Proteases"
VP$Protein_Class[rownames(VP) == "TMPRSS2"] = "Proteases"

ha_r = rowAnnotation(df = VP, 
                     col = list(Viral_Protein = cols_VP,
                                Protein_Class = cols_VP_2),
                     show_legend = T)


ht = Heatmap(M, name = "Normalized \nExpression", 
             top_annotation = ha_c, right_annotation = ha_r,
             col = colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red")), 
             cluster_rows = T, cluster_columns = F, show_heatmap_legend = T,
             show_column_names = T, show_row_names = T)

pdf("DE/Heatmaps/DE_heatmap_nosymbols_complex.COVID.bothAnnotations.normToUT.pdf",height = 16, width = 11)
print(ht)
dev.off()

#VP$Protein_Large_Class[VP$Viral_Protein == "Spike"] = "Structural"
#VP$Protein_Class[VP$Viral_Protein == "Spike"] = "Structural"
VP2 = VP[,-c(2)]
VP2 = as.data.frame(VP2)
colnames(VP2) = "Protein_Class"
rownames(VP2) = rownames(VP)

VP2$Protein_Large_Class = as.vector(VP2$Protein_Class)
VP2$Protein_Large_Class[VP2$Protein_Large_Class == "Envelope"] = "Structural"
VP2$Protein_Large_Class[VP2$Protein_Large_Class == "Membrane"] = "Structural"
VP2$Protein_Large_Class[VP2$Protein_Large_Class == "Nucleocapsid"] = "Structural"
VP2$Protein_Large_Class[VP2$Protein_Large_Class == "Spike"] = "Structural"
VP2$Protein_Large_Class[VP2$Protein_Large_Class == "ORF"] = "Accessory"

VP2$Protein_Large_Class[rownames(VP2) == "ACE2"] = "ACE2"
VP2$Protein_Large_Class[rownames(VP2) == "TMPRSS2"] = "TMPRSS2"
VP2$Protein_Large_Class[rownames(VP2) == "PCSK7"] = "Proteases"
VP2$Protein_Large_Class[rownames(VP2) == "CTSB"] = "Proteases"

VP2$Protein_Large_Class = factor(VP2$Protein_Large_Class, levels=c("ACE2","TMPRSS2","Proteases","Structural","NSP","Accessory"))

cols_VP_2 <- gg_color_hue(length(unique(VP2$Protein_Large_Class))-2)
cols_VP_2 = c("#000000","gray",cols_VP_2)
VP$Protein_Large_Class = factor(VP2$Protein_Large_Class,levels = c("ACE2","TMPRSS2","Structural","NSP","Accessory"))
names(cols_VP_2) = levels(VP2$Protein_Large_Class)

VP3 = VP2$Protein_Large_Class
VP3 = as.data.frame(VP3)
colnames(VP3) = "Protein_Class"
rownames(VP3) = rownames(VP2)
names(VP3$Protein_Class) = rownames(VP3)

ha_r = rowAnnotation(df = VP3, 
                     col = list(
                       Protein_Class = cols_VP_2),
                     show_legend = T)


ht = Heatmap(M, name = "Normalized \nExpression", 
             #top_annotation = ha_c, 
             right_annotation = ha_r,
             col = colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red")), 
             cluster_rows = T, cluster_columns = F, show_heatmap_legend = T,
             show_column_names = T, show_row_names = T,row_dend_width = unit(20, "mm"))

pdf("~/Box/Erle group presentations/LRB/SARS-CoV-2 paper/Walter_COV_2_figures/DE_heatmap_nosymbols_complex.COVID.oneAnnotations.normToUT.v4.pdf",height = 16, width = 11)
print(ht)
dev.off()

####
diffed_genes = c("CLCA1","POSTN","IFI6","USP18","CCL20","CSF3","CXCL9","IRF8")
diffed_genes_ens = annot[annot$Gene %in% diffed_genes,]
diffed_genes  = diffed_genes_ens$Gene
diffed_genes_ens = diffed_genes_ens$Ensembl_ID


#syms <- unique(syms)
genes <- annot[diffed_genes_ens,]$Gene
M <- stabilizedBR[diffed_genes_ens,]
#UT, IL17, IFNa, IFNg, IL13, CO
M = M[,c(31:36,25:30,7:12,13:18,19:24,1:6)]
#M <- M[,c(1:12,19:24,31:36)]
#M <-subset(M,,names)
#M <- M - rowMeans(M, na.rm=T)

#DE.heatmap(M, dendrogram="both", fact=conds[,c("condition"), drop=F], zlim =c(-2,2), lhei=c(1,.25,5,1), height = 1200, width=1200, main=paste("DE SARS-CoV-2 interacting genes - FDR < ", FDR_cut,sep=""), file.prefix="DE/Heatmaps/DE_heatmap_no_symbols.SARS", breaks=.2, Symbols=NA, min.col = "blue", mid.col = "white", max.col = "red", margins = c(12,12,12,12))

genes_Cat = c("IL13","IL17","IL17","IFNa","IL13","IFNg","IFNg","IFNa")
names(genes_Cat) = genes
genes_Cat = as.data.frame(genes_Cat)

M = merge(annot,M, by.x="Ensembl_ID",by.y='row.names')
rownames(M) = M$Gene
M = M[,seq(13,48)]

M = M[order(rownames(M)),]

M_UT = M[,c(1:6)]

M <- M - rowMeans(M_UT, na.rm=T)

M = M[match(diffed_genes, rownames(M)),]

targets <- as.data.frame(as.matrix(colnames(M)))

conds <- cbind(targets, do.call(rbind, strsplit(as.character(targets$V1),'_')))
rownames(conds) <- conds[,c(1)]
conds <- conds[,-c(1)]
condsSmall = conds
condsSmall = condsSmall[,-c(2)]
condsSmall = as.data.frame(condsSmall)
colnames(condsSmall) = "condition"

diffed_genes = c("CLCA1","POSTN","IFI6","USP18","CCL20","CSF3","CXCL9","IRF8")
genes_Cat = c("IL13","IL13","IFNa","IFNa","IL17","IL17","IFNg","IFNg")
names(genes_Cat) = diffed_genes
genes_Cat = as.data.frame(genes_Cat)

colnames(genes_Cat) = "Cytokine_Response"
genes_Cat$Cytokine_Response = factor(genes_Cat$Cytokine_Response,levels = c("IL17","IL13","IFNa","IFNg"))



ha_c = HeatmapAnnotation(df = condsSmall, 
                         col = list(condition = cols),
                         show_legend = T)

cols_VP_2 <- gg_color_hue(length(unique(genes_Cat$Cytokine_Response)))
names(cols_VP_2) = levels(genes_Cat$Cytokine_Response)

genes_Cat = genes_Cat[match(rownames(M), rownames(genes_Cat)),]
genes_Cat = as.data.frame(genes_Cat)
colnames(genes_Cat) = "Cytokine_Response"
rownames(genes_Cat) = rownames(M)

ha_r = rowAnnotation(df = genes_Cat, 
                     col = list(Cytokine_Response = cols_VP_2),
                     show_legend = T)


ht = Heatmap(M, name = "Normalized \nExpression", 
             #top_annotation = ha_c, 
             right_annotation = ha_r,
             col = colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red")), 
             cluster_rows = F, cluster_columns = F, show_heatmap_legend = T,
             show_column_names = T, show_row_names = T,row_dend_width = unit(20, "mm"))

pdf("~/Box/Erle group presentations/LRB/SARS-CoV-2 paper/Walter_COV_2_figures/DE_heatmap_nosymbols_complex.Treatment_Markers.oneAnnotations.normToUT.pdf",height = 3, width = 11)
print(ht)
dev.off()


VP = as.data.frame(allgenes_SARS2_tmp[,c("Type","Viral_Protein")])
rownames(VP) = allgenes_SARS2_tmp$Gene
colnames(VP) = c("Protein_Class","Viral_Protein")

#VP = VP[order(rownames(VP)),]

VP$Viral_Protein = factor(VP$Viral_Protein)

cols_VP <- gg_color_hue(length(unique(VP$Viral_Protein)))
names(cols_VP) = levels(VP$Viral_Protein)

VP = VP[rownames(VP) %in% rownames(M),]

cols_VP_2 <- gg_color_hue(length(unique(VP$Protein_Class)))
VP$Protein_Class = factor(VP$Protein_Class,levels = c("Proteases","Envelope","Membrane","Nucleocapsid","NSP","ORF","Spike"))
names(cols_VP_2) = levels(VP$Protein_Class)


VP2 = VP[,-c(2)]
VP2 = as.data.frame(VP2)
colnames(VP2) = "Protein_Class"
rownames(VP2) = rownames(VP)

VP2$Protein_Large_Class = as.vector(VP2$Protein_Class)
VP2$Protein_Large_Class[VP2$Protein_Large_Class == "Envelope"] = "Structural"
VP2$Protein_Large_Class[VP2$Protein_Large_Class == "Membrane"] = "Structural"
VP2$Protein_Large_Class[VP2$Protein_Large_Class == "Nucleocapsid"] = "Structural"
VP2$Protein_Large_Class[VP2$Protein_Large_Class == "Spike"] = "Structural"
VP2$Protein_Large_Class[VP2$Protein_Large_Class == "ORF"] = "Accessory"

VP2$Protein_Large_Class[rownames(VP2) == "ACE2"] = "ACE2"
VP2$Protein_Large_Class[rownames(VP2) == "TMPRSS2"] = "TMPRSS2"
VP2$Protein_Large_Class = factor(VP2$Protein_Large_Class, levels=c("ACE2","TMPRSS2","Proteases","Structural","NSP","Accessory"))

cols_VP_2 <- gg_color_hue(length(unique(VP2$Protein_Large_Class))-2)
cols_VP_2 = c("#000000","#000000",cols_VP_2)
VP$Protein_Large_Class = factor(VP2$Protein_Large_Class,levels = c("ACE2","TMPRSS2","Proteases","Structural","NSP","Accessory"))
names(cols_VP_2) = levels(VP2$Protein_Large_Class)

VP3 = VP2$Protein_Large_Class
VP3 = as.data.frame(VP3)
colnames(VP3) = "Protein_Class"

ha_r = rowAnnotation(df = VP3, 
                     col = list(
                       Protein_Class = cols_VP_2),
                     show_legend = T)


ht = Heatmap(M, name = "Normalized \nExpression", 
             #top_annotation = ha_c, 
             right_annotation = ha_r,
             col = colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red")), 
             cluster_rows = T, cluster_columns = F, show_heatmap_legend = T,
             show_column_names = T, show_row_names = T,row_dend_width = unit(20, "mm"))

pdf("~/Box/Erle group presentations/LRB/SARS-CoV-2 paper/Walter_COV_2_figures/DE_heatmap_nosymbols_complex.COVID.oneAnnotations.normToUT.TMPRSS_ACE2.v3.pdf",height = 16, width = 11)
print(ht)
dev.off()



####

condsSmall <- as.data.frame(conds[,c("condition")])
colnames(condsSmall) <- "condition"
rownames(condsSmall) <- rownames(conds)

cols <- gg_color_hue(6)
names(cols) = levels(condsSmall$condition)

ha = HeatmapAnnotation(df = condsSmall, 
                       col = list(condition = cols),
                       show_legend = T)

ht = Heatmap(M, name = "Normalized \nExpression", 
             top_annotation = ha, 
             col = colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red")), 
             cluster_rows = T, cluster_columns = T, show_heatmap_legend = T,
             show_column_names = T, show_row_names = F)

pdf("~/Box/Erle group presentations/LRB/SARS-CoV-2 paper/Walter_COV_2_figures/DE_heatmap_nosymbols_complex.pdf",height = 12, width = 11)
print(ht)
dev.off()



###############################################
#### Generate heatmaps for all comparisons ####
###############################################




sigcounts <- sigcounts[nrow(sigcounts):1,]
c = 1
for (i in grep(".FDR", colnames(allgenes2))) {
	#j <- colnames(allgenes)[i]
	name <- gsub(".FDR","", colnames(allgenes2)[i])
	list <- allgenes2[which(allgenes2[,i] < fdrcutoff),]$Ensembl_ID
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
		genes <- annot[list,]$Gene
		if (nrow(M) > 0) {
			DE.heatmap(M, dendrogram="both", fact=condsSelect[,c("condition"), drop=F], zlim =c(-2,2), lhei=c(1,.5,6,1), height = 1200, width=1200, main=paste("DE genes - from ", name, " - FDR < ", fdrcutoff), file.prefix=paste("DE/Heatmaps/",name,"/DE_heatmap_nosymbols.",name,".only",sep=""), breaks=.2, Symbols=NULL, min.col = "blue", mid.col = "white", max.col = "red")
			DE.heatmap(M, dendrogram="both", fact=condsSelect[,c("condition"), drop=F], zlim =c(-2,2), lhei=c(1,.5,y,1),  width=1200,  height=x, main=paste("DE genes - from ", name, " - FDR < ", fdrcutoff), file.prefix=paste("DE/Heatmaps/",name,"/DE_heatmap_symbols",name,".only.withGenes",sep=""), breaks=.2, Symbols=genes, min.col = "blue", mid.col = "white", max.col = "red")
		}
	}
	if (as.list(strsplit(name, "-|_|vs"))[[1]][c(1)] != "condition" && grepl("_vs_",name) == TRUE){
		cond <- as.list(strsplit(name, "_"))[[1]][c(1)]
		group <- as.list(strsplit(name, "_"))[[1]][c(2,4)]
		condsSelect <- conds[which(conds[,cond] %in% group),]
		M <- stabilizedBR[,colnames(stabilizedBR) %in% rownames(condsSelect)]
		M <- M[list,]
		M <- M - rowMeans(M, na.rm=T)
		genes <- annot[list,]$Gene
		if (nrow(M) > 0) {
			DE.heatmap(M, dendrogram="both", fact=condsSelect[,c("condition"), drop=F], zlim =c(-2,2), lhei=c(1,.5,6,1), height = 1200, width=1200, main=paste("DE genes - from ", name, " - FDR < ", fdrcutoff), file.prefix=paste("DE/Heatmaps/",name,"/DE_heatmap_nosymbols.",name,".only",sep=""), breaks=.2, Symbols=NULL, min.col = "blue", mid.col = "white", max.col = "red")
			DE.heatmap(M, dendrogram="both", fact=condsSelect[,c("condition"), drop=F], zlim =c(-2,2), lhei=c(1,.5,y,1),  width=1200,  height=x,  main=paste("DE genes - from ", name, " - FDR < ", fdrcutoff), file.prefix=paste("DE/Heatmaps/",name,"/DE_heatmap_symbols",name,".only.withGenes",sep=""), breaks=.2, Symbols=genes, min.col = "blue", mid.col = "white", max.col = "red")
		}
	}
	M <- stabilizedBR[list,]
	M <- M - rowMeans(M, na.rm=T)
	genes <- annot[list,]$Gene
	if (nrow(M) > 0) {
		DE.heatmap(M, dendrogram="both", fact=conds[,c("condition"), drop=F], zlim =c(-2,2), lhei=c(1,.5,6,1), height = 1200, width=1200, main=paste("DE genes - from ", name, " - FDR < ", fdrcutoff), file.prefix=paste("DE/Heatmaps/",name,"/DE_heatmap_nosymbols.",name,".allSamples",sep=""), breaks=.2, Symbols=NULL, min.col = "blue", mid.col = "white", max.col = "red")
		DE.heatmap(M, dendrogram="both", fact=conds[,c("condition"), drop=F], zlim =c(-2,2), lhei=c(1,.5,y,1),  width=1200,  height=x,  main=paste("DE genes - from ", name, " - FDR < ", fdrcutoff), file.prefix=paste("DE/Heatmaps/",name,"/DE_heatmap_symbols",name,".allSamples.withGenes",sep=""), breaks=.2, Symbols=genes, min.col = "blue", mid.col = "white", max.col = "red")
	}
	c = c + 1
}
	


stabilizedForBox <- t(stabilizedBR)
stabilizedForBoxAnnotated <- merge(stabilizedForBox, conds, by="row.names")


rownames(stabilizedForBoxAnnotated) <- stabilizedForBoxAnnotated$Row.names

stabilizedForBoxAnnotatedIL13 = stabilizedForBoxAnnotated[c(19:24,31:36),]

png("SOCS1.box.png",width = 3,height = 4,units = "in",res = 300)
ggplot(data=stabilizedForBoxAnnotatedIL13) + geom_boxplot(aes(x=condition,y=ENSG00000185338,fill=condition)) + 
  labs(y="Log2 Stabilized Counts") +
  scale_fill_manual(breaks = c("UT", "IL13"), values=c("blue", "darkgreen"))+
  theme_bw()
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

#### Peaks for KD ####

up = read.table("~/data/Bonser_ChIP_cleanedPeaks/Analysis_with_old_peakFix/RNA_results/correct_version/up.txt",header=F,sep="\t")
down = read.table("~/data/Bonser_ChIP_cleanedPeaks/Analysis_with_old_peakFix/RNA_results/correct_version/down.txt",header=F,sep="\t")
  
peakFile = read.table("~/data/Bonser_ChIP_cleanedPeaks/Analysis_with_old_peakFix/RNA_results/correct_version/condition-D24_IL13vsD24_UT.DESeq2_results.sort.vsAll.bulkRNA.bed",
                      header=F,sep="\t")
allPeaks = read.table("~/data/Bonser_ChIP_cleanedPeaks/Analysis_with_old_peakFix/DE/All_genes.Bonser_ChIP_cleanedPeaks.txt",
                      header=T,sep="\t")

allGenes = read.table("~/data/Erle_Bonser_0529/Analysis/DE/All_genes.Erle_Bonser_0529.txt",header=T,sep="\t",row.names = 1)
allGenesIL13Sig = allGenes[allGenes$condition.IL13vsUT.FDR <= 0.1,]
allGenesIL13Sig = allGenesIL13Sig[!(is.na(allGenesIL13Sig$Gene)),]
cutBaseMean_Gene = min(allGenesIL13Sig$Mean_Normalized_Counts)

allGenesOverMin = allGenes[allGenes$Mean_Normalized_Counts >= cutBaseMean_Gene,]

allPeakIL13Sig = allPeaks[allPeaks$condition.D24_IL13vsD24_UT.FDR <= 0.1,]
allPeakIL13Sig = allPeakIL13Sig[!(is.na(allPeakIL13Sig$chr)),]
cutBaseMean = min(allPeakIL13Sig$Mean_Normalized_Counts)

allPeaksOverMin = allPeaks[allPeaks$Mean_Normalized_Counts >= cutBaseMean,]

colnames(peakFile) = c("chr","start","stop","peak_id","log2FC","padj","tss","ensembl_id","gene_name","gene_log2FC","gene_padj","distance")

peakFileFilter = peakFile[!(is.na(peakFile$padj)),]
peakFileFilter = peakFileFilter[!(is.na(peakFileFilter$gene_padj)),]

#peakFileFilter = peakFile[peakFile$peak_id %in% allPeaksOverMin$peak_ID,]
#peakFileFilter = peakFileFilter[peakFileFilter$ensembl_id %in% rownames(allGenesOverMin),]

colnames(peakFileFilter) = c("chr","start","stop","peak_id","log2FC","padj","tss","ensembl_id","gene_name","gene_log2FC","gene_padj","distance")
peakFileFilter$DE_cat = "No_Change"
peakFileFilter$DE_cat[peakFileFilter$log2FC > 0 & peakFileFilter$padj <= 0.1] = "Increase"
peakFileFilter$DE_cat[peakFileFilter$log2FC < 0 & peakFileFilter$padj <= 0.1] = "Decrease"

peakFileFilter$DE_cat = factor(peakFileFilter$DE_cat,levels = c("No_Change","Increase","Decrease"))

peakFileFilter$pair = paste(peakFileFilter$peak_id,peakFileFilter$ensembl_id,sep="_")
peakFileFilterU = peakFileFilter[!duplicated(peakFileFilter$pair),]

prot_genes = allgenes2[allgenes2$Biotype == "protein_coding",]

peakFileFilterU = peakFileFilterU[peakFileFilterU$ensembl_id %in% prot_genes$Ensembl_ID,drop=F,]

ggplot(peakFileFilterU, aes(y=gene_log2FC, x=DE_cat)) + geom_boxplot(outlier.shape = NA) + 
  theme_classic() + geom_hline(yintercept=0) + ylim(-5,5)
ggsave("peak_vs_nearby_gene.boxplot.pdf", height=5, width=5, units=c("in"))

ggplot(peakFileFilter, aes(y=gene_log2FC, x=DE_cat)) + geom_boxplot() + 
  theme_classic() + geom_hline(yintercept=0)
ggsave("peak_vs_nearby_gene.boxplot.withOutliers.pdf", height=5, width=5, units=c("in"))

bulk = read.table("~/data/Bonser_ChIP_cleanedPeaks/Analysis_with_old_peakFix/RNA_results/correct_version/bulk_results.tss.name.strand.sort.ChIP.txt",sep='\t',header=F)

colnames(bulk) = c("Gene_chr","Gene_start","Gene_end","Ensembl_ID","Score","Strand","Ensembl_trans_ID","Gene","Gene_log2FC","Gene_FC","Gene_pval","Gene_FDR",
                   "Peak_chr","Peak_start","Peak_end","Peak_ID","Peak_log2FC","Peak_FC","Pead_pval","Peak_FDR","Distance")

bulk$DE_status = "Not_DE"       
bulk$DE_status[bulk$Gene_FDR <= 0.05] = "DE"
bulkFilter = bulk[abs(bulk$Distance) <= 5e5,]
bulkFilter$DE_status = factor(bulkFilter$DE_status, levels = c("Not_DE","DE"))
ggplot(bulkFilter, aes(x=Distance, color=DE_status,fill=DE_status)) +
  geom_histogram(aes(y=..density..),alpha=0.5, position="identity", binwidth = 20000)+
  geom_density(alpha=.2) + theme_classic() + xlim(-300000,200000)
ggsave("Gene_V_vs_nearby_gene.boxplot.pdf", height=5, width=5, units=c("in"))
