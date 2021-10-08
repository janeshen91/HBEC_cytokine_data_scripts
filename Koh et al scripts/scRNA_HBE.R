##### Full scRNA-seq  workflow ####

library(Seurat)
options(future.globals.maxSize= 8891289600)
library(dplyr)
library(Matrix)
library(ComplexHeatmap)
library(circlize)
library(plyr)
library(ggplot2)
library("devtools")
library(clusterExperiment)
library(scone)
library(zinbwave)
library(circlize)
library(reshape)
library(MAST)
library(data.table)

project <- c("Erle_Bonser_scRNA_HG38_seurat3_05_15_2019")
dir <- paste(paste("/data/projects/",project, sep = "/"),"Analysis_Batch_Regress_filter_Seurat3_SCT_cellTypes", sep = "/")
dir.create(dir)
setwd(dir)
library(gplots)
library(cowplot)
library(reshape2)
library(future)
plan("multiprocess", workers = 8)

library("calibrate")


dir.create("DE")
dir.create("QC")
dir.create("DE/MAplots/")
dir.create("DE/VolcanoPlots/")
dir.create("DE/MMplots/")
dir.create("DE/pathways/")
dir.create("DE/Heatmaps/")
dir.create("DE/DEfiles/")
dir.create("QC/clusterCalls/")
dir.create("QC/round1/")


data <- Read10X("/Volumes/Pegasus/Projects/Erle_Bonser_scRNA_HG38_cellranger_2.1_04_03_2018//HBE_round2/outs/raw_gene_bc_matrices_mex/GRCh38/")


colnames(data) <- gsub("-10", "-HBE10", colnames(data))
colnames(data) <- gsub("-1", "-HBE01", colnames(data))
colnames(data) <- gsub("-2", "-HBE02", colnames(data))
colnames(data) <- gsub("-3", "-HBE03", colnames(data))
colnames(data) <- gsub("-4", "-HBE04", colnames(data))
colnames(data) <- gsub("-5", "-HBE05", colnames(data))
colnames(data) <- gsub("-6", "-HBE06", colnames(data))
colnames(data) <- gsub("-7", "-HBE07", colnames(data))
colnames(data) <- gsub("-8", "-HBE08", colnames(data))
colnames(data) <- gsub("-9", "-HBE09", colnames(data))

cellNames <- colnames(data)


Inds <- c("Rep1","Rep2","Rep3","Rep6")

c1 <- read.table("/Volumes/Pegasus/Projects/Erle_Bonser_scRNA_HG38_cellranger_2.1_04_03_2018//HBE1_demuxlet.default_all_GT_noMissing.ref.alphaTune2.top50k.best", header=T, row.names = 1)
rownames(c1) <- gsub("-1", "-HBE01", rownames(c1))
c1$BC <- do.call(rbind, strsplit(as.character(rownames(c1)),"-"))[,c(1)]
c1$batch <- do.call(rbind, strsplit(as.character(rownames(c1)),"-"))[,c(2)]
conds1 <- cbind(c1, do.call(rbind, strsplit(as.character(c1$BEST),'-')))
colnames(conds1)[24] <- "dropStatus"

conds1$Ind <- paste(conds1$'2', conds1$'3', sep=":")
conds1$Ind[conds1$dropStatus == 'SNG' & conds1$'2' == "Rep1"] <- "Rep1"
conds1$Ind[conds1$dropStatus == 'SNG' & conds1$'2' == "Rep2"] <- "Rep2"
conds1$Ind[conds1$dropStatus == 'SNG' & conds1$'2' == "Rep3"] <- "Rep3"
conds1$Ind[conds1$dropStatus == 'SNG' & conds1$'2' == "Rep6"] <- "Rep6"
conds1$Day <- conds1$Ind
conds1$Day <- gsub("Rep1", "D00", conds1$Day)
conds1$Day <- gsub("Rep2", "D24", conds1$Day)
conds1$Day <- gsub("Rep6", "D24", conds1$Day)
conds1$Day <- gsub("Rep3", "D24", conds1$Day)
conds1$Treatment <- conds1$Ind
conds1$Treatment <- gsub("Rep1", "UT", conds1$Treatment)
conds1$Treatment <- gsub("Rep2", "IL17", conds1$Treatment)
conds1$Treatment <- gsub("Rep6", "IFN", conds1$Treatment)
conds1$Treatment <- gsub("Rep3", "CO", conds1$Treatment)
conds1SNG <- conds1[which(conds1$'2' %in% Inds & conds1$dropStatus == "SNG"),]

condsRed1 <- conds1[,c(1:5,7,14,20:30)]
condsRed1SNG <- conds1SNG[,c(1:5,7,14,20:30)]


c2 <- read.table("/Volumes/Pegasus/Projects/Erle_Bonser_scRNA_HG38_cellranger_2.1_04_03_2018/HBE2_demuxlet.default_all_GT_noMissing.ref.alphaTune2.top50k.best", header=T, row.names = 1)
rownames(c2) <- gsub("1", "HBE02", rownames(c2))
c2$BC <- do.call(rbind, strsplit(as.character(rownames(c2)),"-"))[,c(1)]
c2$batch <- do.call(rbind, strsplit(as.character(rownames(c2)),"-"))[,c(2)]
conds2 <- cbind(c2, do.call(rbind, strsplit(as.character(c2$BEST),'-')))
colnames(conds2)[24] <- "dropStatus"

conds2$Ind <- paste(conds2$'2', conds2$'3', sep=":")
conds2$Ind[conds2$dropStatus == 'SNG' & conds2$'2' == "Rep1"] <- "Rep1"
conds2$Ind[conds2$dropStatus == 'SNG' & conds2$'2' == "Rep2"] <- "Rep2"
conds2$Ind[conds2$dropStatus == 'SNG' & conds2$'2' == "Rep3"] <- "Rep3"
conds2$Ind[conds2$dropStatus == 'SNG' & conds2$'2' == "Rep6"] <- "Rep6"
conds2$Day <- conds2$Ind
conds2$Day <- gsub("Rep1", "D03", conds2$Day)
conds2$Day <- gsub("Rep2", "D00", conds2$Day)
conds2$Day <- gsub("Rep6", "D24", conds2$Day)
conds2$Day <- gsub("Rep3", "D24", conds2$Day)
conds2$Treatment <- conds2$Ind
conds2$Treatment <- gsub("Rep1", "UT", conds2$Treatment)
conds2$Treatment <- gsub("Rep2", "UT", conds2$Treatment)
conds2$Treatment <- gsub("Rep6", "IL17", conds2$Treatment)
conds2$Treatment <- gsub("Rep3", "IFN", conds2$Treatment)
conds2SNG <- conds2[which(conds2$'2' %in% Inds & conds2$dropStatus == "SNG"),]

condsRed2 <- conds2[,c(1:5,7,14,20:30)]
condsRed2SNG <- conds2SNG[,c(1:5,7,14,20:30)]


c3 <- read.table("/Volumes/Pegasus/Projects/Erle_Bonser_scRNA_HG38_cellranger_2.1_04_03_2018/HBE3_demuxlet.default_all_GT_noMissing.ref.alphaTune2.top50k.best", header=T, row.names = 1)
rownames(c3) <- gsub("1", "HBE03", rownames(c3))
c3$BC <- do.call(rbind, strsplit(as.character(rownames(c3)),"-"))[,c(1)]
c3$batch <- do.call(rbind, strsplit(as.character(rownames(c3)),"-"))[,c(2)]
conds3 <- cbind(c3, do.call(rbind, strsplit(as.character(c3$BEST),'-')))
colnames(conds3)[24] <- "dropStatus"

conds3$Ind <- paste(conds3$'2', conds3$'3', sep=":")
conds3$Ind[conds3$dropStatus == 'SNG' & conds3$'2' == "Rep1"] <- "Rep1"
conds3$Ind[conds3$dropStatus == 'SNG' & conds3$'2' == "Rep2"] <- "Rep2"
conds3$Ind[conds3$dropStatus == 'SNG' & conds3$'2' == "Rep3"] <- "Rep3"
conds3$Ind[conds3$dropStatus == 'SNG' & conds3$'2' == "Rep6"] <- "Rep6"
conds3$Day <- conds3$Ind
conds3$Day <- gsub("Rep1", "D09", conds3$Day)
conds3$Day <- gsub("Rep2", "D03", conds3$Day)
conds3$Day <- gsub("Rep6", "D00", conds3$Day)
conds3$Day <- gsub("Rep3", "D24", conds3$Day)
conds3$Treatment <- conds3$Ind
conds3$Treatment <- gsub("Rep1", "UT", conds3$Treatment)
conds3$Treatment <- gsub("Rep2", "UT", conds3$Treatment)
conds3$Treatment <- gsub("Rep6", "UT", conds3$Treatment)
conds3$Treatment <- gsub("Rep3", "IL17", conds3$Treatment)

conds3SNG <- conds3[which(conds3$'2' %in% Inds & conds3$dropStatus == "SNG"),]

condsRed3 <- conds3[,c(1:5,7,14,20:30)]
condsRed3SNG <- conds3SNG[,c(1:5,7,14,20:30)]



c4 <- read.table("/Volumes/Pegasus/Projects/Erle_Bonser_scRNA_HG38_cellranger_2.1_04_03_2018/HBE4_demuxlet.default_all_GT_noMissing.ref.alphaTune2.top50k.best", header=T, row.names = 1)
rownames(c4) <- gsub("1", "HBE04", rownames(c4))
c4$BC <- do.call(rbind, strsplit(as.character(rownames(c4)),"-"))[,c(1)]
c4$batch <- do.call(rbind, strsplit(as.character(rownames(c4)),"-"))[,c(2)]
conds4 <- cbind(c4, do.call(rbind, strsplit(as.character(c4$BEST),'-')))
colnames(conds4)[24] <- "dropStatus"

conds4$Ind <- paste(conds4$'2', conds4$'3', sep=":")
conds4$Ind[conds4$dropStatus == 'SNG' & conds4$'2' == "Rep1"] <- "Rep1"
conds4$Ind[conds4$dropStatus == 'SNG' & conds4$'2' == "Rep2"] <- "Rep2"
conds4$Ind[conds4$dropStatus == 'SNG' & conds4$'2' == "Rep3"] <- "Rep3"
conds4$Ind[conds4$dropStatus == 'SNG' & conds4$'2' == "Rep6"] <- "Rep6"
conds4$Day <- conds4$Ind
conds4$Day <- gsub("Rep1", "D24", conds4$Day)
conds4$Day <- gsub("Rep2", "D09", conds4$Day)
conds4$Day <- gsub("Rep6", "D03", conds4$Day)
conds4$Day <- gsub("Rep3", "D00", conds4$Day)
conds4$Treatment <- conds4$Ind
conds4$Treatment <- gsub("Rep1", "UT", conds4$Treatment)
conds4$Treatment <- gsub("Rep2", "UT", conds4$Treatment)
conds4$Treatment <- gsub("Rep6", "UT", conds4$Treatment)
conds4$Treatment <- gsub("Rep3", "UT", conds4$Treatment)

conds4SNG <- conds4[which(conds4$'2' %in% Inds & conds4$dropStatus == "SNG"),]

condsRed4 <- conds4[,c(1:5,7,14,20:30)]
condsRed4SNG <- conds4SNG[,c(1:5,7,14,20:30)]


c5 <- read.table("/Volumes/Pegasus/Projects/Erle_Bonser_scRNA_HG38_cellranger_2.1_04_03_2018/HBE5_demuxlet.default_all_GT_noMissing.ref.alphaTune2.top50k.best", header=T, row.names = 1)
rownames(c5) <- gsub("1", "HBE05", rownames(c5))
c5$BC <- do.call(rbind, strsplit(as.character(rownames(c5)),"-"))[,c(1)]
c5$batch <- do.call(rbind, strsplit(as.character(rownames(c5)),"-"))[,c(2)]
conds5 <- cbind(c5, do.call(rbind, strsplit(as.character(c5$BEST),'-')))
colnames(conds5)[24] <- "dropStatus"

conds5$Ind <- paste(conds5$'2', conds5$'3', sep=":")
conds5$Ind[conds5$dropStatus == 'SNG' & conds5$'2' == "Rep1"] <- "Rep1"
conds5$Ind[conds5$dropStatus == 'SNG' & conds5$'2' == "Rep2"] <- "Rep2"
conds5$Ind[conds5$dropStatus == 'SNG' & conds5$'2' == "Rep3"] <- "Rep3"
conds5$Ind[conds5$dropStatus == 'SNG' & conds5$'2' == "Rep6"] <- "Rep6"
conds5$Day <- conds5$Ind
conds5$Day <- gsub("Rep1", "D24", conds5$Day)
conds5$Day <- gsub("Rep2", "D24", conds5$Day)
conds5$Day <- gsub("Rep6", "D09", conds5$Day)
conds5$Day <- gsub("Rep3", "D03", conds5$Day)
conds5$Treatment <- conds5$Ind
conds5$Treatment <- gsub("Rep1", "IL13", conds5$Treatment)
conds5$Treatment <- gsub("Rep2", "UT", conds5$Treatment)
conds5$Treatment <- gsub("Rep6", "UT", conds5$Treatment)
conds5$Treatment <- gsub("Rep3", "UT", conds5$Treatment)
conds5SNG <- conds5[which(conds5$'2' %in% Inds & conds5$dropStatus == "SNG"),]

condsRed5 <- conds5[,c(1:5,7,14,20:30)]
condsRed5SNG <- conds5SNG[,c(1:5,7,14,20:30)]


c6 <- read.table("/Volumes/Pegasus/Projects/Erle_Bonser_scRNA_HG38_cellranger_2.1_04_03_2018/HBE6_demuxlet.default_all_GT_noMissing.ref.alphaTune2.top50k.best", header=T, row.names = 1)
rownames(c6) <- gsub("1", "HBE06", rownames(c6))
c6$BC <- do.call(rbind, strsplit(as.character(rownames(c6)),"-"))[,c(1)]
c6$batch <- do.call(rbind, strsplit(as.character(rownames(c6)),"-"))[,c(2)]
conds6 <- cbind(c6, do.call(rbind, strsplit(as.character(c6$BEST),'-')))
colnames(conds6)[24] <- "dropStatus"

conds6$Ind <- paste(conds6$'2', conds6$'3', sep=":")
conds6$Ind[conds6$dropStatus == 'SNG' & conds6$'2' == "Rep1"] <- "Rep1"
conds6$Ind[conds6$dropStatus == 'SNG' & conds6$'2' == "Rep2"] <- "Rep2"
conds6$Ind[conds6$dropStatus == 'SNG' & conds6$'2' == "Rep3"] <- "Rep3"
conds6$Ind[conds6$dropStatus == 'SNG' & conds6$'2' == "Rep6"] <- "Rep6"
conds6$Day <- conds6$Ind
conds6$Day <- gsub("Rep1", "D24", conds6$Day)
conds6$Day <- gsub("Rep2", "D24", conds6$Day)
conds6$Day <- gsub("Rep6", "D24", conds6$Day)
conds6$Day <- gsub("Rep3", "D09", conds6$Day)
conds6$Treatment <- conds6$Ind
conds6$Treatment <- gsub("Rep1", "CO", conds6$Treatment)
conds6$Treatment <- gsub("Rep2", "IL13", conds6$Treatment)
conds6$Treatment <- gsub("Rep6", "UT", conds6$Treatment)
conds6$Treatment <- gsub("Rep3", "UT", conds6$Treatment)
conds6SNG <- conds6[which(conds6$'2' %in% Inds & conds6$dropStatus == "SNG"),]

condsRed6 <- conds6[,c(1:5,7,14,20:30)]
condsRed6SNG <- conds6SNG[,c(1:5,7,14,20:30)]


c7 <- read.table("/Volumes/Pegasus/Projects/Erle_Bonser_scRNA_HG38_cellranger_2.1_04_03_2018/HBE7_demuxlet.default_all_GT_noMissing.ref.alphaTune2.top50k.best", header=T, row.names = 1)
rownames(c7) <- gsub("1", "HBE07", rownames(c7))
c7$BC <- do.call(rbind, strsplit(as.character(rownames(c7)),"-"))[,c(1)]
c7$batch <- do.call(rbind, strsplit(as.character(rownames(c7)),"-"))[,c(2)]
conds7 <- cbind(c7, do.call(rbind, strsplit(as.character(c7$BEST),'-')))
colnames(conds7)[24] <- "dropStatus"

conds7$Ind <- paste(conds7$'2', conds7$'3', sep=":")
conds7$Ind[conds7$dropStatus == 'SNG' & conds7$'2' == "Rep1"] <- "Rep1"
conds7$Ind[conds7$dropStatus == 'SNG' & conds7$'2' == "Rep2"] <- "Rep2"
conds7$Ind[conds7$dropStatus == 'SNG' & conds7$'2' == "Rep3"] <- "Rep3"
conds7$Ind[conds7$dropStatus == 'SNG' & conds7$'2' == "Rep6"] <- "Rep6"
conds7$Day <- conds7$Ind
conds7$Day <- gsub("Rep1", "D24", conds7$Day)
conds7$Day <- gsub("Rep2", "D24", conds7$Day)
conds7$Day <- gsub("Rep6", "D24", conds7$Day)
conds7$Day <- gsub("Rep3", "D24", conds7$Day)
conds7$Treatment <- conds7$Ind
conds7$Treatment <- gsub("Rep1", "IFN", conds7$Treatment)
conds7$Treatment <- gsub("Rep2", "CO", conds7$Treatment)
conds7$Treatment <- gsub("Rep6", "IL13", conds7$Treatment)
conds7$Treatment <- gsub("Rep3", "UT", conds7$Treatment)
conds7SNG <- conds7[which(conds7$'2' %in% Inds & conds7$dropStatus == "SNG"),]

condsRed7 <- conds7[,c(1:5,7,14,20:30)]
condsRed7SNG <- conds7SNG[,c(1:5,7,14,20:30)]



c8 <- read.table("/Volumes/Pegasus/Projects/Erle_Bonser_scRNA_HG38_cellranger_2.1_04_03_2018/HBE8_demuxlet.default_all_GT_noMissing.ref.alphaTune2.top50k.best", header=T, row.names = 1)
rownames(c8) <- gsub("1", "HBE08", rownames(c8))
c8$BC <- do.call(rbind, strsplit(as.character(rownames(c8)),"-"))[,c(1)]
c8$batch <- do.call(rbind, strsplit(as.character(rownames(c8)),"-"))[,c(2)]
conds8 <- cbind(c8, do.call(rbind, strsplit(as.character(c8$BEST),'-')))
colnames(conds8)[24] <- "dropStatus"

conds8$Ind <- paste(conds8$'2', conds8$'3', sep=":")
conds8$Ind[conds8$dropStatus == 'SNG' & conds8$'2' == "Rep1"] <- "Rep1"
conds8$Ind[conds8$dropStatus == 'SNG' & conds8$'2' == "Rep2"] <- "Rep2"
conds8$Ind[conds8$dropStatus == 'SNG' & conds8$'2' == "Rep3"] <- "Rep3"
conds8$Ind[conds8$dropStatus == 'SNG' & conds8$'2' == "Rep6"] <- "Rep6"
conds8$Day <- conds8$Ind
conds8$Day <- gsub("Rep1", "D24", conds8$Day)
conds8$Day <- gsub("Rep2", "D24", conds8$Day)
conds8$Day <- gsub("Rep6", "D24", conds8$Day)
conds8$Day <- gsub("Rep3", "D24", conds8$Day)
conds8$Treatment <- conds8$Ind
conds8$Treatment <- gsub("Rep1", "IL17", conds8$Treatment)
conds8$Treatment <- gsub("Rep2", "IFN", conds8$Treatment)
conds8$Treatment <- gsub("Rep6", "CO", conds8$Treatment)
conds8$Treatment <- gsub("Rep3", "IL13", conds8$Treatment)
conds8SNG <- conds8[which(conds8$'2' %in% Inds & conds8$dropStatus == "SNG"),]

condsRed8 <- conds8[,c(1:5,7,14,20:30)]
condsRed8SNG <- conds8SNG[,c(1:5,7,14,20:30)]


c9 <- read.table("/Volumes/Pegasus/Projects/Erle_Bonser_scRNA_HG38_cellranger_2.1_04_03_2018/HBE9_demuxlet.default_all_GT_noMissing.ref.alphaTune2.top50k.best", header=T, row.names = 1)
rownames(c9) <- gsub("1", "HBE09", rownames(c9))
c9$BC <- do.call(rbind, strsplit(as.character(rownames(c9)),"-"))[,c(1)]
c9$batch <- do.call(rbind, strsplit(as.character(rownames(c9)),"-"))[,c(2)]
conds9 <- cbind(c9, do.call(rbind, strsplit(as.character(c9$BEST),'-')))
colnames(conds9)[24] <- "dropStatus"

conds9$Ind <- paste(conds9$'2', conds9$'3', sep=":")
conds9$Ind[conds9$dropStatus == 'SNG' & conds9$'2' == "Rep1"] <- "Rep1"
conds9$Ind[conds9$dropStatus == 'SNG' & conds9$'2' == "Rep2"] <- "Rep2"
conds9$Ind[conds9$dropStatus == 'SNG' & conds9$'2' == "Rep3"] <- "Rep3"
conds9$Ind[conds9$dropStatus == 'SNG' & conds9$'2' == "Rep6"] <- "Rep6"
conds9$Day <- conds9$Ind
conds9$Day <- gsub("Rep1", "D09", conds9$Day)
conds9$Day <- gsub("Rep2", "D24", conds9$Day)
conds9$Day <- gsub("Rep6", "D03", conds9$Day)
conds9$Day <- gsub("Rep3", "D24", conds9$Day)
conds9$Treatment <- conds9$Ind
conds9$Treatment <- gsub("Rep1", "UT", conds9$Treatment)
conds9$Treatment <- gsub("Rep2", "IFN", conds9$Treatment)
conds9$Treatment <- gsub("Rep6", "UT", conds9$Treatment)
conds9$Treatment <- gsub("Rep3", "IL13", conds9$Treatment)
conds9SNG <- conds9[which(conds9$'2' %in% Inds & conds9$dropStatus == "SNG"),]

condsRed9 <- conds9[,c(1:5,7,14,20:30)]
condsRed9SNG <- conds9SNG[,c(1:5,7,14,20:30)]

c10 <- read.table("/Volumes/Pegasus/Projects/Erle_Bonser_scRNA_HG38_cellranger_2.1_04_03_2018/HBE10_demuxlet.default_all_GT_noMissing.ref.alphaTune2.top50k.best", header=T, row.names = 1)
rownames(c10) <- gsub("1", "HBE10", rownames(c10))
c10$BC <- do.call(rbind, strsplit(as.character(rownames(c10)),"-"))[,c(1)]
c10$batch <- do.call(rbind, strsplit(as.character(rownames(c10)),"-"))[,c(2)]
conds10 <- cbind(c10, do.call(rbind, strsplit(as.character(c10$BEST),'-')))
colnames(conds10)[24] <- "dropStatus"

conds10$Ind <- paste(conds10$'2', conds10$'3', sep=":")
conds10$Ind[conds10$dropStatus == 'SNG' & conds10$'2' == "Rep1"] <- "Rep1"
conds10$Ind[conds10$dropStatus == 'SNG' & conds10$'2' == "Rep2"] <- "Rep2"
conds10$Ind[conds10$dropStatus == 'SNG' & conds10$'2' == "Rep3"] <- "Rep3"
conds10$Ind[conds10$dropStatus == 'SNG' & conds10$'2' == "Rep6"] <- "Rep6"
conds10$Day <- conds10$Ind
conds10$Day <- gsub("Rep1", "D24", conds10$Day)
conds10$Day <- gsub("Rep2", "D03", conds10$Day)
conds10$Day <- gsub("Rep6", "D00", conds10$Day)
conds10$Day <- gsub("Rep3", "D24", conds10$Day)
conds10$Treatment <- conds10$Ind
conds10$Treatment <- gsub("Rep1", "IL13", conds10$Treatment)
conds10$Treatment <- gsub("Rep2", "UT", conds10$Treatment)
conds10$Treatment <- gsub("Rep6", "UT", conds10$Treatment)
conds10$Treatment <- gsub("Rep3", "IL17", conds10$Treatment)
conds10SNG <- conds10[which(conds10$'2' %in% Inds & conds10$dropStatus == "SNG"),]

condsRed10 <- conds10[,c(1:5,7,14,20:30)]
condsRed10SNG <- conds10SNG[,c(1:5,7,14,20:30)]



condsSNG <- rbind(condsRed1SNG, condsRed2SNG, condsRed3SNG, condsRed4SNG, condsRed5SNG, condsRed6SNG, condsRed7SNG, condsRed8SNG, condsRed9SNG, condsRed10SNG)
conds <- rbind(condsRed1, condsRed2, condsRed3, condsRed4, condsRed5, condsRed6, condsRed7, condsRed8, condsRed9, condsRed10)

condsSNG <- condsSNG[which(rownames(condsSNG) %in% cellNames),]
conds <- conds[which(rownames(conds) %in% cellNames),]

nrow(conds)
ncol(data)

dataSNG <- data[,which(colnames(data) %in% rownames(condsSNG))]
dataSNGColSums <-  Matrix::colSums(dataSNG)
condsSNG$nUMI <- dataSNGColSums

conds <- conds[which(conds$batch != "."),]

data <- data[,which( colnames(data) %in% rownames(conds)),]
conds <- conds[which(rownames(conds) %in% colnames(data)),]

dataColSums <- Matrix::colSums(data)
conds$nUMI <-  dataColSums

conds$DBLmSNG <- conds$LLK12 - conds$SNG.LLK1
condsSNG$DBLmSNG <- condsSNG$LLK12 - condsSNG$SNG.LLK1

conds$priorSNGmDBL <- conds$PRB.SNG1- conds$PRB.DBL


#### NUMBER OF SNP AND LLK filter ####


condsFilter <- condsSNG[which(condsSNG$SNG.LLK1 <= -1 & condsSNG$N.SNP >= 10 & condsSNG$nUMI >= 100 & condsSNG$PRB.DBL <= 0.7),]
condsF <- conds[which(conds$N.SNP >= 10 & conds$nUMI >= 1000),]


condsFinal <- condsFilter[which(rownames(condsFilter) %in% colnames(data)),]
dataFinal <- data[,which(colnames(data) %in% rownames(condsFinal))]

condsFinal$longName <- paste(paste(condsFinal$BC, condsFinal$batch, sep="_"), condsFinal$Ind, condsFinal$Day, condsFinal$Treatment, sep="-")
condsFinal$Condition <- paste(condsFinal$Day,condsFinal$Treatment, sep="_")
condsFinal$ConditionWithInd <- paste(condsFinal$Ind, condsFinal$Day,condsFinal$Treatment, sep="_")


condsFinal$Ind <- gsub("Rep1","Sub1", condsFinal$Ind)
condsFinal$Ind <- gsub("Rep2","Sub2", condsFinal$Ind)
condsFinal$Ind <- gsub("Rep3","Sub3", condsFinal$Ind)
condsFinal$Ind <- gsub("Rep6","Sub6", condsFinal$Ind)

condsFinal$ConditionWithInd <- gsub("Rep1","Sub1", condsFinal$ConditionWithInd)
condsFinal$ConditionWithInd <- gsub("Rep2","Sub2", condsFinal$ConditionWithInd)
condsFinal$ConditionWithInd <- gsub("Rep3","Sub3", condsFinal$ConditionWithInd)
condsFinal$ConditionWithInd <- gsub("Rep6","Sub6", condsFinal$ConditionWithInd)



condsFinalTime <- condsFinal[which(condsFinal$Treatment == "UT"),]
condsFinalTreatment <- condsFinal[which(condsFinal$Day == "D24"),]


dataTime <- dataFinal[,which(colnames(dataFinal) %in% rownames(condsFinalTime))]
dataTreatment <- dataFinal[,which(colnames(dataFinal) %in% rownames(condsFinalTreatment))]

dataTreatment <- dataTreatment[,Matrix::colSums(dataTreatment > 0) > 200]
dataTreatment <- dataTreatment[Matrix::rowSums(dataTreatment > 0) > 10,]

#dataTreatment = readRDS(file= "../raw_data.treatment.rds")
#condsFinalTreatment  = readRDS(file = "../conditions.treatment.rds")

soTreatment <- CreateSeuratObject(counts  = dataTreatment, min.cells = 10, 
                                  min.features = 100, 
                                  project = "Treatment", names.field = 1, 
                                  names.delim = "_", meta.data=condsFinalTreatment)




markerGenes <- c("CDC6","COL17A1","MCM6","SNAI2","TWIST1","KRT17","KRT5","BPIFB1",
                 "SLC5A5","TCN1","ZBP1","IFI44L","ISG15","IFIT2","AGR3","CCDC176",
                 "DNAH7","FOXJ1","CCDC153","CCDC170","SNTN","MUC5B","S100A8","SCGB1A1",
                 "SCGB3A1","FCGBP","IFIT1","IFI6","FOXA3","CST1","MUC5AC","NOS2",
                 "SLC26A4","TFF3","ALOX15","SERPINB10","SOCS1","CCL26","CLDN3","AQP3",
                 "CCND1","KRT7","UPK1B","TSPAN8","TSPAN19","GSTA1","GSTA2","CCDC146",
                 "CCDC173","KLF4","KRT23","SERPINB2","CCDC39","DNAAF1","DNAH10","DNAH11",
                 "DNAH5","DNAH6","POSTN","SPDEF","CFTR")


ccGene <- read.table(file = "../cell_cycle_genes.txt", header = F,sep = "\t")
ccGeneList <- unique(as.vector(ccGene$V2))
ccGeneList <- ccGeneList[which(ccGeneList %in% rownames(soTreatment@assays$RNA@data))]

ribo <- read.table(file="../ribo_list.txt",sep="\t", header=F)

ribo.genes <- unique(as.vector(ribo$V1))
ribo.genes <- ribo.genes[which(ribo.genes %in% rownames(soTreatment@assays$RNA@data))]


markerGenes <- unique(markerGenes)



marker <- soTreatment@assays$RNA@var.features

MTlist <- marker[grep("MT-",marker)]
RPSlist <- marker[grepl("RPS",marker)]
RPLlist <- marker[grepl("RPL",marker)]

#removeList <- c(MTlist, RPSlist, RPLlist)

#marker <- marker[!marker %in% removeList]

marker <- unique(c(marker, markerGenes))

IL13table <- read.table("../IL13_Y_vs_N.xls",header=T, row.names = 1, sep='\t')
IFNatable <- read.table("../IFNa_Y_vs_N.xls",header=T, row.names = 1, sep='\t')
IFNgtable <- read.table("../IFNg_Y_vs_N.xls",header=T, row.names = 1, sep='\t')

IL13table <- IL13table[which(IL13table$IL13_Y_vs_N.log2FC > 0),]
IFNatable <- IFNatable[which(IFNatable$IFNa_Y_vs_N.log2FC > 0),]
IFNgtable <- IFNgtable[which(IFNgtable$IFNg_Y_vs_N.log2FC > 0),]

IL13genes <- IL13table$Gene
IFNagenes <- as.vector(IFNatable$Gene)
IFNggenes <- as.vector(IFNgtable$Gene)
IFNgenes <- unique(c(IFNggenes,IFNagenes))

IL13genesShort <- markerGenes[markerGenes %in% IL13genes]
IFNgenesShort <- markerGenes[markerGenes %in% IFNgenes]

IL13genesUniq <- IL13genesShort[!IL13genesShort %in% IFNgenesShort]
IFNgenesUniq <- IFNgenesShort[!IFNgenesShort %in% IL13genesShort]
COgenesUniq <- IL13genesShort[IL13genesShort %in% IFNgenesShort]


basalMarkers <- unique(c("CDC6","COL17A1","MCM6","SNAI2","TWIST1","KRT17","KRT5","AQP3","CCND1","S100A2","CLDN1","TP63","ITGA5","ITGA6","EGFR"))
secretoryMarkers <- unique(c("BPIFB1","TCN1","SCGB1A1","SCGB3A1","TSPAN8","WFDC2"))
ciliatedMarkers <- unique(c("ARG3","CCDC176","ERICH3","DNAH7","FOXJ1","CCDC153","CCDC170","SNTN","TSPAN19","CCDC146","CCDC173","CFAP53","CCDC39","GSTA1","GSTA2"))

bigList <- c("IL13genesUniq","IFNgenesUniq","COgenesUniq","basalMarkers","secretoryMarkers","ciliatedMarkers","ccGeneList","ribo.genes","mito.genes")
bigListNoCC <- c("IL13genesUniq","IFNgenesUniq","COgenesUniq","basalMarkers","secretoryMarkers","ciliatedMarkers")

#list <- c(IL13genesUniq,IFNgenesUniq,COgenesUniq,basalMarkers,secretoryMarkers,ciliatedMarkers,ccGeneList,ribo.genes)

#### End Lists ####

dim(soTreatment@assays$RNA@data)

mito.genes <- grep("^MT-", rownames(soTreatment@assays$RNA@data), value = T)
percent.mito <- colSums(as.array(expm1(soTreatment@assays$RNA@data[mito.genes, ])))/colSums(as.array(expm1(soTreatment@assays$RNA@data)))
#total.UMI <- colSums(as.array(soTreatment@raw.data))
soTreatment <- AddMetaData(soTreatment, percent.mito, "percent.mito")

ribo <- read.table(file="../Analysis_Batch_Regress_filter_Seurat3/ribo_list.txt",sep="\t", header=F)

ribo.genes <- unique(as.vector(ribo$V1))
ribo.genes <- ribo.genes[which(ribo.genes %in% rownames(soTreatment@assays$RNA@data))]

percent.ribo <- colSums(as.array(expm1(soTreatment@assays$RNA@data[ribo.genes, ])))/colSums(as.array(expm1(soTreatment@assays$RNA@data)))
soTreatment <- AddMetaData(soTreatment, percent.ribo, "percent.ribo")

condsFinal <- soTreatment@meta.data
condsFinal$Ind <- gsub("Rep1","Sub1", condsFinal$Ind)
condsFinal$Ind <- gsub("Rep2","Sub2", condsFinal$Ind)
condsFinal$Ind <- gsub("Rep3","Sub3", condsFinal$Ind)
condsFinal$Ind <- gsub("Rep6","Sub6", condsFinal$Ind)

soTreatment@meta.data <- condsFinal

#VlnPlot(object = soTime, features.plot = c("nGene","RD.TOTL"), group.by = "ConditionWithInd", nCol = 2, do.return = TRUE, x.lab.rot = T)
VlnPlot(object = soTreatment, features = c("nCount_RNA","nFeature_RNA"), ncol = 1,pt.size = 0, group.by = "ConditionWithInd",  x.lab.rot = T,log = T)
VlnPlot(object = soTreatment, features = c("percent.mito","percent.ribo"),ncol = 1, pt.size = 0, group.by = "batch", nCol = 3, do.return = TRUE, x.lab.rot = T)


#soTime <- ScaleData(soTime, vars.to.regress = c("batch"))

soTreatment@meta.data$batch_day <- soTreatment@meta.data$batch

soTreatment@meta.data$batch_day[soTreatment@meta.data$batch_day == "HBE01"] <- "B2"
soTreatment@meta.data$batch_day[soTreatment@meta.data$batch_day == "HBE02"] <- "B2"
soTreatment@meta.data$batch_day[soTreatment@meta.data$batch_day == "HBE03"] <- "B2"
soTreatment@meta.data$batch_day[soTreatment@meta.data$batch_day == "HBE04"] <- "B2"
soTreatment@meta.data$batch_day[soTreatment@meta.data$batch_day == "HBE05"] <- "B1"
soTreatment@meta.data$batch_day[soTreatment@meta.data$batch_day == "HBE06"] <- "B2"
soTreatment@meta.data$batch_day[soTreatment@meta.data$batch_day == "HBE07"] <- "B1"
soTreatment@meta.data$batch_day[soTreatment@meta.data$batch_day == "HBE08"] <- "B2"
soTreatment@meta.data$batch_day[soTreatment@meta.data$batch_day == "HBE09"] <- "B3"
soTreatment@meta.data$batch_day[soTreatment@meta.data$batch_day == "HBE10"] <- "B3"



soTreatment <- subset(soTreatment, subset = nFeature_RNA >= 500 & nFeature_RNA < 7500 & nCount_RNA > 1200 & percent.mito < .5 & percent.ribo < .35)


for ( i in unique(soTreatment@meta.data$Treatment)){
  print(i)
  cells <- rownames(soTreatment@meta.data[(soTreatment@meta.data$Treatment == i),])
  
  soTMP <- subset(soTreatment,cells = cells)
  soTMP <- subset(soTMP)
  
  
  
  #soTMP <- NormalizeData(soTMP,normalization.method = "LogNormalize",scale.factor = 1000)
  
  soTMP <- FindVariableFeatures(object = soTMP, mean.function = ExpMean, dispersion.function = LogVMR, 
                                selection.method = "vst",nfeatures = 2000,verbose = F)
  soTMP = ScaleData(object = soTMP,vars.to.regress = c("nCount_RNA","percent.mito","percent.ribo"),verbose = F)
  
  soTMP <- RunPCA(object  = soTMP, ndims.print = 1:6, 
                  nfeatures.print = 5, npcs = 100,features = VariableFeatures(object = soTMP))
  
  
  soTMP <- FindNeighbors(soTMP, dims = 1:100,verbose = F)
  soTMP<- RunUMAP(soTMP, dims = 1:100,verbose = F)
  soTMP <- FindClusters(soTMP, resolution = 0.6,verbose = F)
  
  VlnPlot(object = soTMP, features = c("nCount_RNA","nFeature_RNA","percent.mito","percent.ribo"), pt.size = 0, do.return = TRUE, x.lab.rot = T)
  ggsave(paste("QC/QC.violins.cluster.D24_",i,"_v3.pdf",sep=""), height=7, width=7, units=c("in"))
  
  DimPlot(soTMP, reduction = "umap",label = T)
  ggsave(paste("QC/UMAP.named.cluster.D24_",i,"_v3.pdf",sep=""), height=7, width=7, units=c("in"))
  DimPlot(soTMP, reduction = "umap", label = T,group.by ="Ind")
  ggsave(paste("QC/UMAP.named.cluster.D24_",i,"_v3.byInd.pdf", sep=""), height=7, width=7, units=c("in"))
  
  
  
  umap <- soTMP@reductions$umap@cell.embeddings
  for (j in bigList) {
    k <- get(j) 
    m <- as.matrix(soTMP@assays$RNA@data)
    m <- m[which(rownames(m) %in% k),]
    if (class(m) == "matrix") {
      n <- colMeans(m)
    } else { n <- m }
    umap <- merge(umap, as.matrix(n), by='row.names')
    rownames(umap) <- umap[,c(1)]
    umap <- umap[,-c(1)]
    colnames(umap) <- gsub("V1",j,colnames(umap))
    umap <- umap[order(umap[,j]),]
    p <- ggplot2::ggplot(umap, aes_string(x = "UMAP_1", y = "UMAP_2")) +
      geom_point(aes(color=eval(parse(text = j))), size=1) + scale_colour_gradient(low = "lightgray",high = "red") +
      xlab("Dimension 1") + labs(colour = c("Normalized \nExpression")) +
      ylab("Dimension 2") 
    png(paste("QC/round1/markers/UMAP.plot.marker.", j, "_",i,".png", sep=""), height=6, width=7,units = "in",res = 150)
    print(p)
    dev.off()
  }
  
  soTMP <- SCTransform(soTMP, verbose = FALSE)
  
  conds = soTMP@meta.data
  conds = merge(conds,umap,by="row.names")
  
  allAves = as.data.frame(setNames(replicate(4,numeric(0), simplify = F),letters[0:4]))
  for (m in levels(soTMP$seurat_clusters)){
    condsTMP = conds[(conds$seurat_clusters == m),]
    aves = colMeans(condsTMP[,c("basalMarkers","secretoryMarkers","ciliatedMarkers","ccGeneList")])
    allAves = rbind(allAves, aves)
  }
  colnames(allAves) = c("basalMarkers","secretoryMarkers","ciliatedMarkers","ccGeneList")
  rownames(allAves) = paste("cluster_",levels(soTMP$seurat_clusters),sep="")
  allAves$secObasal = allAves$secretoryMarkers/allAves$basalMarkers
  allAves$cellType = "Intermediate"
  #allAves$cellType[(allAves$secretoryMarkers >= 1*allAves$basalMarkers)] = "Secretory_low"
  #allAves$cellType[(allAves$secretoryMarkers >= 1.225*allAves$basalMarkers)] = "Secretory_l2"
  allAves$cellType[(allAves$secretoryMarkers >= 1.4*allAves$basalMarkers)] = "Secretory"
  allAves$cellType[(allAves$basalMarkers >= 2*allAves$secretoryMarkers)] = "Basal"
  allAves$cellType[(allAves$ciliatedMarkers >= .25)] = "Ciliated"
  allAves$cellType[(allAves$ccGeneList >= .4)] = "Proliferating_Basal"
  
  soTMP@meta.data$cellType = soTMP@meta.data$seurat_clusters
  for (n in levels(soTMP$seurat_clusters)){
    cellTypeN = allAves[(as.numeric(n)+1),"cellType"]
    soTMP@meta.data$cellType = factor(soTMP@meta.data$cellType, 
                                      levels = c("Proliferating_Basal","Basal","Intermediate","Secretory","Ciliated"))
    soTMP@meta.data$cellType[(soTMP@meta.data$seurat_clusters == as.numeric(n))] = cellTypeN
  }
  assign(paste('soD24_',i,"_v3",sep=""),soTMP)
  saveRDS(soTMP, file = paste("soD24_",i,"_v3.rds",sep=""))
  DimPlot(soTMP, reduction = "umap",label = T, group.by = "cellType")
  ggsave(paste("QC/UMAP.named.cluster.D24_",i,"_v3.byCellType.pdf", sep=""), height=7, width=7, units=c("in"))
  print(i)
  print(allAves)
}


soD24_UT_v3 <- FindNeighbors(soD24_UT_v3, dims = 1:100,verbose = F)
soD24_UT_v3<- RunUMAP(soD24_UT_v3, dims = 1:100,verbose = F)
soD24_UT_v3 <- FindClusters(soD24_UT_v3, resolution = 2,verbose = F)
DimPlot(soD24_UT_v3, reduction = "umap",label = T)
DimPlot(soD24_UT_v3, reduction = "umap", label = T,group.by ="Ind")

VlnPlot(object = soTMP, features = c("nCount_RNA","nFeature_RNA","percent.mito","percent.ribo"), pt.size = 0, do.return = TRUE, x.lab.rot = T)
ggsave(paste("QC/QC.violins.cluster.D24_",i,"_v3.pdf",sep=""), height=7, width=7, units=c("in"))

####### split by condition and individual #######

soTreatment@meta.data$ConditionWithInd_Batch = paste(soTreatment@meta.data$ConditionWithInd,soTreatment@meta.data$batch_day,sep="_")

cells <- rownames(soTreatment@meta.data[!(soTreatment@meta.data$ConditionWithInd_Batch == "Sub3_D24_IL17_B2"),])
soTreatmentClean <- subset(soTreatment,cells = cells)

cells <- rownames(soTreatmentClean@meta.data[(soTreatmentClean@meta.data$Treatment == "UT"),])
soTreatmentClean <- subset(soTreatmentClean,cells = cells)

for ( i in unique(soTreatmentClean@meta.data$ConditionWithInd_Batch)){
  print(i)
  cells <- rownames(soTreatmentClean@meta.data[(soTreatmentClean@meta.data$ConditionWithInd_Batch == i),])
  
  soTMP <- subset(soTreatmentClean,cells = cells)
  soTMP <- subset(soTMP)
  
  
  
  #soTMP <- NormalizeData(soTMP,normalization.method = "LogNormalize",scale.factor = 1000)
  
  soTMP <- FindVariableFeatures(object = soTMP, mean.function = ExpMean, dispersion.function = LogVMR, 
                                selection.method = "vst",nfeatures = 2000,verbose = F)
  soTMP = ScaleData(object = soTMP,vars.to.regress = c("nCount_RNA","percent.mito","percent.ribo"),verbose = F)
  
  soTMP <- RunPCA(object  = soTMP, ndims.print = 1:6, 
                  nfeatures.print = 5, npcs = 100,features = VariableFeatures(object = soTMP))
  
  
  soTMP <- FindNeighbors(soTMP, dims = 1:100,verbose = F)
  soTMP<- RunUMAP(soTMP, dims = 1:100,verbose = F)
  soTMP <- FindClusters(soTMP, resolution = 0.6,verbose = F)
  
  VlnPlot(object = soTMP, features = c("nCount_RNA","nFeature_RNA","percent.mito","percent.ribo"), pt.size = 0, do.return = TRUE, x.lab.rot = T)
  ggsave(paste("QC/QC.violins.cluster.D24_",i,"_v3.pdf",sep=""), height=7, width=7, units=c("in"))
  
  DimPlot(soTMP, reduction = "umap",label = T)
  ggsave(paste("QC/UMAP.named.cluster.D24_",i,"_v3.pdf",sep=""), height=7, width=7, units=c("in"))
  DimPlot(soTMP, reduction = "umap", label = T,group.by ="Ind")
  ggsave(paste("QC/UMAP.named.cluster.D24_",i,"_v3.byInd.pdf", sep=""), height=7, width=7, units=c("in"))
  
  
  
  umap <- soTMP@reductions$umap@cell.embeddings
  for (j in bigList) {
    k <- get(j) 
    m <- as.matrix(soTMP@assays$RNA@data)
    m <- m[which(rownames(m) %in% k),]
    if (class(m) == "matrix") {
      n <- colMeans(m)
    } else { n <- m }
    umap <- merge(umap, as.matrix(n), by='row.names')
    rownames(umap) <- umap[,c(1)]
    umap <- umap[,-c(1)]
    colnames(umap) <- gsub("V1",j,colnames(umap))
    umap <- umap[order(umap[,j]),]
    p <- ggplot2::ggplot(umap, aes_string(x = "UMAP_1", y = "UMAP_2")) +
      geom_point(aes(color=eval(parse(text = j))), size=1) + scale_colour_gradient(low = "lightgray",high = "red") +
      xlab("Dimension 1") + labs(colour = c("Normalized \nExpression")) +
      ylab("Dimension 2") 
    png(paste("QC/round1/markers/UMAP.plot.marker.", j, "_",i,".png", sep=""), height=6, width=7,units = "in",res = 150)
    print(p)
    dev.off()
  }
  
  soTMP <- SCTransform(soTMP, verbose = FALSE)
  
  conds = soTMP@meta.data
  conds = merge(conds,umap,by="row.names")
  
  allAves = as.data.frame(setNames(replicate(4,numeric(0), simplify = F),letters[0:4]))
  for (m in levels(soTMP$seurat_clusters)){
    condsTMP = conds[(conds$seurat_clusters == m),]
    aves = colMeans(condsTMP[,c("basalMarkers","secretoryMarkers","ciliatedMarkers","ccGeneList")])
    allAves = rbind(allAves, aves)
  }
  colnames(allAves) = c("basalMarkers","secretoryMarkers","ciliatedMarkers","ccGeneList")
  rownames(allAves) = paste("cluster_",levels(soTMP$seurat_clusters),sep="")
  allAves$secObasal = allAves$secretoryMarkers/allAves$basalMarkers
  allAves$cellType = "Intermediate"
  #allAves$cellType[(allAves$secretoryMarkers >= 1*allAves$basalMarkers)] = "Secretory_low"
  #allAves$cellType[(allAves$secretoryMarkers >= 1.225*allAves$basalMarkers)] = "Secretory_l2"
  allAves$cellType[(allAves$secretoryMarkers >= 1.4*allAves$basalMarkers)] = "Secretory"
  allAves$cellType[(allAves$basalMarkers >= 2*allAves$secretoryMarkers)] = "Basal"
  allAves$cellType[(allAves$ciliatedMarkers >= .25)] = "Ciliated"
  allAves$cellType[(allAves$ccGeneList >= .4)] = "Proliferating_Basal"
  
  soTMP@meta.data$cellType = soTMP@meta.data$seurat_clusters
  for (n in levels(soTMP$seurat_clusters)){
    cellTypeN = allAves[(as.numeric(n)+1),"cellType"]
    soTMP@meta.data$cellType = factor(soTMP@meta.data$cellType, 
                                      levels = c("Proliferating_Basal","Basal","Intermediate","Secretory","Ciliated"))
    soTMP@meta.data$cellType[(soTMP@meta.data$seurat_clusters == as.numeric(n))] = cellTypeN
  }
  assign(paste('soD24_',i,"_v3",sep=""),soTMP)
  saveRDS(soTMP, file = paste("soD24_",i,"_v3.rds",sep=""))
  DimPlot(soTMP, reduction = "umap",label = T, group.by = "cellType")
  ggsave(paste("QC/UMAP.named.cluster.D24_",i,"_v3.byCellType.pdf", sep=""), height=7, width=7, units=c("in"))
  print(i)
  print(allAves)
}


soList2 = c(soD24_Sub1_D24_UT_B2_v3,soD24_Sub2_D24_UT_B1_v3,soD24_Sub6_D24_UT_B2_v3,soD24_Sub3_D24_UT_B1_v3,
            soD24_Sub6_D24_IFN_B2_v3,soD24_Sub2_D24_IL17_B2_v3,soD24_Sub3_D24_CO_B2_v3,soD24_Sub3_D24_IFN_B2_v3,soD24_Sub6_D24_IL17_B2_v3,
            soD24_Sub1_D24_IL13_B1_v3,soD24_Sub2_D24_IL13_B2_v3,soD24_Sub1_D24_CO_B2_v3,
            soD24_Sub1_D24_IFN_B1_v3,soD24_Sub6_D24_IL13_B1_v3,soD24_Sub2_D24_CO_B1_v3,  
            soD24_Sub3_D24_IL13_B2_v3,soD24_Sub2_D24_IFN_B2_v3,soD24_Sub6_D24_CO_B2_v3,soD24_Sub1_D24_IL17_B2_v3,soD24_Sub3_D24_IL13_B3_v3,
            soD24_Sub2_D24_IFN_B3_v3,soD24_Sub3_D24_IL17_B3_v3,soD24_Sub1_D24_IL13_B3_v3)

soList3 = c(soD24_Sub1_D24_UT_B2_v3,soD24_Sub2_D24_UT_B1_v3,soD24_Sub6_D24_UT_B2_v3,soD24_Sub3_D24_UT_B1_v3)

soFeatures3 <- SelectIntegrationFeatures(object.list = soList3, nfeatures = 8000)
soListPrep3 <- PrepSCTIntegration(object.list = soList3, anchor.features = soFeatures3)

reference_dataset <- c(1,2,3,4)

#### soTreatmentInt3 is only UT samples for reference cell type assignment ####

soAnchors3 <- FindIntegrationAnchors(object.list = soListPrep3, normalization.method = "SCT", 
                                     anchor.features = soFeatures3)

soTreatmentInt3 <- IntegrateData(anchorset = soAnchors3, normalization.method = "SCT")

soTreatmentInt3 <- RunPCA(object = soTreatmentInt3, verbose = FALSE,npcs = 100)
soTreatmentInt3 <- FindNeighbors(soTreatmentInt3, dims = 1:100,verbose = F)
soTreatmentInt3 <- FindClusters(soTreatmentInt3, resolution = 0.25,verbose = F)
soTreatmentInt3 <- RunUMAP(object = soTreatmentInt3, dims = 1:100)

UT_Int.allMarkers = FindAllMarkers(soTreatmentInt3,logfc.threshold = 0.5,only.pos = T)
UT_Int.allMarkersBackUP = UT_Int.allMarkers

for (i in levels(soTreatmentInt3@active.ident)){
  clusterMarkers = rownames(head(UT_Int.allMarkers[(UT_Int.allMarkers$cluster == i),],n=9))
  FeaturePlot(soTreatmentInt3,features = clusterMarkers,pt.size = 0.1)
  ggsave(paste("QC/UMAP.UT_SCT_merged_cluster_",i,".pdf",sep=""), height=7, width=9, units=c("in"))
}

VlnPlot(object = soTreatmentInt3, features = c("nCount_RNA","nFeature_RNA"), ncol = 1,pt.size = 0,  x.lab.rot = T,log = T)
VlnPlot(object = soTreatmentInt3, features = c("percent.mito","percent.ribo"),ncol = 1, pt.size = 0, do.return = TRUE, x.lab.rot = T)

FeaturePlot(soTreatmentInt3,features = c("nCount_RNA","nFeature_RNA","percent.mito","percent.ribo"))

UT_Int.5v2Markers = FindMarkers(soTreatmentInt3,logfc.threshold = 0.5,only.pos = F,ident.1 = 5,ident.2 = 2)
UT_Int.allMarkers[(UT_Int.allMarkers$cluster == 5),]

dir.create("../Analysis_Batch_Regress_SCT/DE/UT_integrated_markers")

top10List = c()

for (i in unique(allMarkers$cluster)) {
  print(i)
  Tmp <-  allMarkers[(allMarkers$cluster == i),]
  print(head(Tmp))
  write.table(Tmp, sep="\t",quote = F,row.names = T,col.names = NA,
              file = paste("../Analysis_Batch_Regress_SCT/DE/UT_integrated_markers/","UT_Int.Cluster_",i,"Markers.txt"))
  genesList = rownames(head(Tmp,n=10))
  top10List = c(top10List,genesList)
}



soTreatmentInt3@meta.data$cellType = soTreatmentInt3@meta.data$integrated_snn_res.0.25
soTreatmentInt3@meta.data$cellType = as.vector(soTreatmentInt3@meta.data$cellType)
soTreatmentInt3@meta.data$cellType[soTreatmentInt3@meta.data$cellType == 0] = "Intermediate"
soTreatmentInt3@meta.data$cellType[soTreatmentInt3@meta.data$cellType == 1] = "Basal"
soTreatmentInt3@meta.data$cellType[soTreatmentInt3@meta.data$cellType == 2] = "Secretory"
soTreatmentInt3@meta.data$cellType[soTreatmentInt3@meta.data$cellType == 3] = "Ciliated"
soTreatmentInt3@meta.data$cellType[soTreatmentInt3@meta.data$cellType == 4] = "Proliferating_Basal"
soTreatmentInt3@meta.data$cellType[soTreatmentInt3@meta.data$cellType == 5] = "Unknown"
soTreatmentInt3@meta.data$cellType = factor(soTreatmentInt3@meta.data$cellType,
                                            levels = c("Proliferating_Basal","Basal","Intermediate",
                                                       "Secretory","Ciliated","Unknown"))

DimPlot(soTreatmentInt3, reduction = "umap",label = T)
ggsave("QC/UMAP.cluster_Ind_Batch_SCT.UT.pdf", height=7, width=7, units=c("in"))
DimPlot(soTreatmentInt3, reduction = "umap",label = T,group.by = "cellType")
ggsave("QC/UMAP.cellType_Ind_Batch_SCT.UT.pdf", height=7, width=7, units=c("in"))
DimPlot(soTreatmentInt3, reduction = "umap",label = T,group.by = "Treatment")
ggsave("QC/UMAP.Treatment_Ind_Batch_SCT.UT.pdf", height=7, width=7, units=c("in"))
DimPlot(soTreatmentInt3, reduction = "umap",label = T,group.by = "Ind")
ggsave("QC/UMAP.Individual_Ind_Batch_SCT.UT.pdf", height=7, width=7, units=c("in"))
DimPlot(soTreatmentInt3, reduction = "umap",label = T)
ggsave("QC/UMAP.org_cluster_Ind_Batch_SCT.UT.pdf", height=7, width=7, units=c("in"))

saveRDS(file = "soTreatmentInt3.rds",soTreatmentInt3)

DoHeatmap(soTreatmentInt3, features = top10List)
ggsave("../Analysis_Batch_Regress_SCT/DE/UT_integrated_markers/Top10.heatmap.pdf",height = 8,width = 8,units = "in",dpi = 150)


######    Assign cell types in UT individuals #####

#soD24_Sub1_D24_UT_B2_v3 = readRDS("soD24_Sub1_D24_UT_B2_v3.rds")
#soD24_Sub2_D24_UT_B2_v3 = readRDS("soD24_Sub2_D24_UT_B1_v3.rds")
#soD24_Sub3_D24_UT_B2_v3 = readRDS("soD24_Sub3_D24_UT_B1_v3.rds")
#soD24_Sub6_D24_UT_B2_v3 = readRDS("soD24_Sub6_D24_UT_B2_v3.rds")

Sub1Cells = soTreatmentInt3@meta.data[(soTreatmentInt3@meta.data$Ind == "Sub1"),]
Sub1CellsVect = Sub1Cells$cellType
names(Sub1CellsVect) = rownames(Sub1Cells)

Sub2Cells = soTreatmentInt3@meta.data[(soTreatmentInt3@meta.data$Ind == "Sub2"),]
Sub2CellsVect = Sub2Cells$cellType
names(Sub2CellsVect) = rownames(Sub2Cells)

Sub3Cells = soTreatmentInt3@meta.data[(soTreatmentInt3@meta.data$Ind == "Sub3"),]
Sub3CellsVect = Sub3Cells$cellType
names(Sub3CellsVect) = rownames(Sub3Cells)

Sub6Cells = soTreatmentInt3@meta.data[(soTreatmentInt3@meta.data$Ind == "Sub6"),]
Sub6CellsVect = Sub6Cells$cellType
names(Sub6CellsVect) = rownames(Sub6Cells)

soD24_Sub1_D24_UT_B2_v3 = AddMetaData(soD24_Sub1_D24_UT_B2_v3,Sub1CellsVect,col.name="cellType")
soD24_Sub2_D24_UT_B1_v3 = AddMetaData(soD24_Sub2_D24_UT_B2_v3,Sub2CellsVect,col.name="cellType")
soD24_Sub6_D24_UT_B2_v3 = AddMetaData(soD24_Sub6_D24_UT_B2_v3,Sub6CellsVect,col.name="cellType")
soD24_Sub3_D24_UT_B1_v3 = AddMetaData(soD24_Sub3_D24_UT_B2_v3,Sub3CellsVect,col.name="cellType")

condsInt3 = soTreatmentInt3@meta.data

ggplot(condsInt3, aes(x = seurat_clusters, fill=Ind, color=Ind)) + 
  geom_bar(position="dodge", width = .75) + theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8)) +
  scale_x_discrete(name ="cluster") + theme_bw()
ggsave("QC/SCT.UT.clusterCounts.png",width = 6, height = 4,dpi = 150,units = "in")


conds.new <- condsInt3 %>% 
  group_by(Ind,seurat_clusters) %>% 
  dplyr::summarise(n = n()) %>%
  mutate(freq = n / sum(n))

ggplot(conds.new,aes(seurat_clusters,freq,fill=Ind))+
  geom_bar(stat="identity",position='dodge', width = .75) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8)) +
  scale_x_discrete(name ="cluster") + 
  scale_y_continuous(name = "percent of cells") + theme_bw()
ggsave("QC/SCT.UT.clusterPercent.png",width = 6, height = 4,dpi = 150,units = "in")

conds.new <- condsInt3 %>% 
  group_by(Ind,cellType) %>% 
  dplyr::summarise(n = n()) %>%
  mutate(freq = n / sum(n))

ggplot(conds.new,aes(cellType,freq,fill=Ind))+
  geom_bar(stat="identity",position='dodge', width = .75) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8)) +
  scale_x_discrete(name ="cellType") + 
  scale_y_continuous(name = "percent of cells") + theme_bw()
ggsave("QC/SCT.UT.cellTypePercent.png",width = 6, height = 4,dpi = 150,units = "in")


#######   Integrate the rest of the samples ###########


soTreatment@meta.data$ConditionWithInd_Batch = paste(soTreatment@meta.data$ConditionWithInd,soTreatment@meta.data$batch_day,sep="_")

cells <- rownames(soTreatment@meta.data[!(soTreatment@meta.data$ConditionWithInd_Batch == "Sub3_D24_IL17_B2"),])
soTreatmentClean2 <- subset(soTreatment,cells = cells)

cells <- rownames(soTreatmentClean2@meta.data[!(soTreatmentClean2@meta.data$Treatment == "UT"),])
soTreatmentClean2 <- subset(soTreatmentClean2,cells = cells)

for ( i in unique(soTreatmentClean2@meta.data$ConditionWithInd_Batch)){
  print(i)
  cells <- rownames(soTreatmentClean2@meta.data[(soTreatmentClean2@meta.data$ConditionWithInd_Batch == i),])
  
  soTMP <- subset(soTreatmentClean2,cells = cells)
  soTMP <- subset(soTMP)
  
  
  
  #soTMP <- NormalizeData(soTMP,normalization.method = "LogNormalize",scale.factor = 1000)
  
  soTMP <- FindVariableFeatures(object = soTMP, mean.function = ExpMean, dispersion.function = LogVMR, 
                                selection.method = "vst",nfeatures = 2000,verbose = F)
  soTMP = ScaleData(object = soTMP,vars.to.regress = c("nCount_RNA","percent.mito","percent.ribo"),verbose = F)
  
  soTMP <- RunPCA(object  = soTMP, ndims.print = 1:6, 
                  nfeatures.print = 5, npcs = 100,features = VariableFeatures(object = soTMP))
  
  
  soTMP <- FindNeighbors(soTMP, dims = 1:100,verbose = F)
  soTMP<- RunUMAP(soTMP, dims = 1:100,verbose = F)
  soTMP <- FindClusters(soTMP, resolution = 0.6,verbose = F)
  
  VlnPlot(object = soTMP, features = c("nCount_RNA","nFeature_RNA","percent.mito","percent.ribo"), pt.size = 0, do.return = TRUE, x.lab.rot = T)
  ggsave(paste("QC/QC.violins.cluster.D24_",i,"_v3.pdf",sep=""), height=7, width=7, units=c("in"))
  
  DimPlot(soTMP, reduction = "umap",label = T)
  ggsave(paste("QC/UMAP.named.cluster.D24_",i,"_v3.pdf",sep=""), height=7, width=7, units=c("in"))
  DimPlot(soTMP, reduction = "umap", label = T,group.by ="Ind")
  ggsave(paste("QC/UMAP.named.cluster.D24_",i,"_v3.byInd.pdf", sep=""), height=7, width=7, units=c("in"))
  
  
  
  umap <- soTMP@reductions$umap@cell.embeddings
  for (j in bigList) {
    k <- get(j) 
    m <- as.matrix(soTMP@assays$RNA@data)
    m <- m[which(rownames(m) %in% k),]
    if (class(m) == "matrix") {
      n <- colMeans(m)
    } else { n <- m }
    umap <- merge(umap, as.matrix(n), by='row.names')
    rownames(umap) <- umap[,c(1)]
    umap <- umap[,-c(1)]
    colnames(umap) <- gsub("V1",j,colnames(umap))
    umap <- umap[order(umap[,j]),]
    p <- ggplot2::ggplot(umap, aes_string(x = "UMAP_1", y = "UMAP_2")) +
      geom_point(aes(color=eval(parse(text = j))), size=1) + scale_colour_gradient(low = "lightgray",high = "red") +
      xlab("Dimension 1") + labs(colour = c("Normalized \nExpression")) +
      ylab("Dimension 2") 
    png(paste("QC/round1/markers/UMAP.plot.marker.", j, "_",i,".png", sep=""), height=6, width=7,units = "in",res = 150)
    print(p)
    dev.off()
  }
  
  soTMP <- SCTransform(soTMP, verbose = FALSE)
  
  conds = soTMP@meta.data
  conds = merge(conds,umap,by="row.names")
  assign(paste('soD24_',i,"_v3",sep=""),soTMP)
  saveRDS(soTMP, file = paste("soD24_",i,"_v3.rds",sep=""))
  print(i)
}



soList4 = c(soD24_Sub1_D24_UT_B2_v3,soD24_Sub2_D24_UT_B1_v3,soD24_Sub6_D24_UT_B2_v3,soD24_Sub3_D24_UT_B1_v3,
            soD24_Sub6_D24_IFN_B2_v3,soD24_Sub2_D24_IL17_B2_v3,soD24_Sub3_D24_CO_B2_v3,soD24_Sub3_D24_IFN_B2_v3,soD24_Sub6_D24_IL17_B2_v3,
            soD24_Sub1_D24_IL13_B1_v3,soD24_Sub2_D24_IL13_B2_v3,soD24_Sub1_D24_CO_B2_v3,
            soD24_Sub1_D24_IFN_B1_v3,soD24_Sub6_D24_IL13_B1_v3,soD24_Sub2_D24_CO_B1_v3,  
            soD24_Sub3_D24_IL13_B2_v3,soD24_Sub2_D24_IFN_B2_v3,soD24_Sub6_D24_CO_B2_v3,soD24_Sub1_D24_IL17_B2_v3,soD24_Sub3_D24_IL13_B3_v3,
            soD24_Sub2_D24_IFN_B3_v3,soD24_Sub3_D24_IL17_B3_v3,soD24_Sub1_D24_IL13_B3_v3)


soFeatures4 <- SelectIntegrationFeatures(object.list = soList4, nfeatures = 8000)
soListPrep4 <- PrepSCTIntegration(object.list = soList4, anchor.features = soFeatures4)

reference_dataset <- c(1,2,3,4)

soAnchors4 <- FindIntegrationAnchors(object.list = soListPrep4, normalization.method = "SCT", 
                                     anchor.features = soFeatures4,reference = 1)
soTreatmentInt4 <- IntegrateData(anchorset = soAnchors4, normalization.method = "SCT")

soTreatmentInt4 <- RunPCA(object = soTreatmentInt4, verbose = FALSE,npcs = 100,assay = "integrated")
soTreatmentInt4 <- FindNeighbors(soTreatmentInt4, dims = 1:100,verbose = F,assay = "integrated")
soTreatmentInt4 <- FindClusters(soTreatmentInt4, resolution = 0.25,verbose = F)
soTreatmentInt4 <- RunUMAP(object = soTreatmentInt4, dims = 1:100,assay = "integrated")

All_Int.allMarkers = FindAllMarkers(soTreatmentInt4,only.pos = T)
All_Int.6_Markers = FindMarkers(soTreatmentInt4,ident.1 = 6,ident.2 = 3, only.pos = F)

for (i in levels(soTreatmentInt4@active.ident)){
  clusterMarkers = rownames(head(All_Int.allMarkers[(All_Int.allMarkers$cluster == i),],n=9))
  FeaturePlot(soTreatmentInt4,features = clusterMarkers,pt.size = 0.1)
  ggsave(paste("QC/UMAP.ALL_SCT_merged_cluster_",i,".pdf",sep=""), height=7, width=9, units=c("in"))
}


VlnPlot(object = soTreatmentInt4, features = c("nCount_RNA","nFeature_RNA"), pt.size = 0,  ncol = 1,x.lab.rot = T,log = T)
ggsave("QC/Violin.ALL_SCT_merged_features.countsQC.pdf", height=7, width=9, units=c("in"))
VlnPlot(object = soTreatmentInt4, features = c("percent.mito","percent.ribo"),pt.size = 0, do.return = TRUE, x.lab.rot = T,log=T)
ggsave("QC/Violin.ALL_SCT_merged_features.percentQC.pdf", height=7, width=9, units=c("in"))

FeaturePlot(soTreatmentInt4,features = c("nCount_RNA","nFeature_RNA","percent.mito","percent.ribo"))
ggsave("QC/UMAP.ALL_SCT_merged_features.pdf", height=7, width=9, units=c("in"))

soTreatmentInt4@meta.data$cellType = soTreatmentInt4@meta.data$integrated_snn_res.0.25
soTreatmentInt4@meta.data$cellType = as.vector(soTreatmentInt4@meta.data$cellType)
soTreatmentInt4@meta.data$cellType[soTreatmentInt4@meta.data$cellType == 0] = "Basal"
soTreatmentInt4@meta.data$cellType[soTreatmentInt4@meta.data$cellType == 4] = "Proliferating_Basal"
soTreatmentInt4@meta.data$cellType[soTreatmentInt4@meta.data$cellType == 2] = "Secretory"
soTreatmentInt4@meta.data$cellType[soTreatmentInt4@meta.data$cellType == 3] = "Ciliated"
soTreatmentInt4@meta.data$cellType[soTreatmentInt4@meta.data$cellType == 1] = "Intermediate"
soTreatmentInt4@meta.data$cellType[soTreatmentInt4@meta.data$cellType == 5] = "Unknown"
soTreatmentInt4@meta.data$cellType[soTreatmentInt4@meta.data$cellType == 6] = "LowQ_Ciliated"
soTreatmentInt4@meta.data$cellType = factor(soTreatmentInt4@meta.data$cellType,
                                            levels = c("Proliferating_Basal","Basal","Intermediate",
                                                       "Secretory","Ciliated","Unknown","LowQ_Ciliated"))

DimPlot(soTreatmentInt4, reduction = "umap",label = T)
ggsave("QC/UMAP.all_treats.cluster_Ind_Batch_SCT.UT.pdf", height=7, width=7, units=c("in"))
DimPlot(soTreatmentInt4, reduction = "umap",label = T,group.by = "cellType")
ggsave("QC/UMAP.all_treats.cellType_Ind_Batch_SCT.UT.pdf", height=7, width=7, units=c("in"))
DimPlot(soTreatmentInt4, reduction = "umap",label = T,group.by = "Treatment")
ggsave("QC/UMAP.all_treats.Treatment_Ind_Batch_SCT.UT.pdf", height=7, width=7, units=c("in"))
DimPlot(soTreatmentInt4, reduction = "umap",label = T,group.by = "Ind")
ggsave("QC/UMAP.all_treats.Individual_Ind_Batch_SCT.UT.pdf", height=7, width=7, units=c("in"))
DimPlot(soTreatmentInt4, reduction = "umap",label = T)
ggsave("QC/UMAP.all_treats.org_cluster_Ind_Batch_SCT.UT.pdf", height=7, width=7, units=c("in"))

soTreatmentInt4@meta.data$cellType_Treatment =  paste(soTreatmentInt4@meta.data$cellType, soTreatmentInt4@meta.data$Treatment)

DoHeatmap(soTreatmentInt4, features = top10List,group.by = "cellType_Treatment",size = 4) + NoLegend()
ggsave("../Analysis_Batch_Regress_SCT/DE/UT_integrated_markers/Top10.heatmap.allTreatments.groupByCellTypeTreatment.png",height = 8,width = 8,units = "in",dpi = 150)


condsInt4 = soTreatmentInt4@meta.data


dt  <- aggregate(rep(1, length(paste0(condsInt4$seurat_clusters, condsInt4$Treatment))),
                 by=list(condsInt4$seurat_clusters, condsInt4$Treatment), sum)
dtCO <- dt[which(dt$Group.2 == "CO"),]
dtIL13 <- dt[which(dt$Group.2 == "IL13"),]
dtIL17 <- dt[which(dt$Group.2 == "IL17"),]
dtIFN <- dt[which(dt$Group.2 == "IFN"),]
dtUT <- dt[which(dt$Group.2 == "UT"),]

dtM <- merge(dtCO,dtIL13,by="Group.1",all = T)
dtM <- merge(dtM,dtIL17,by="Group.1",all = T)
dtM <- merge(dtM,dtIFN,by="Group.1",all = T)
dtM <- merge(dtM,dtUT,by="Group.1",all = T)
dtM[is.na(dtM)] = 0
dtM = dtM[,-c(2,4,6,8,10)]
colnames(dtM) <- c("seurat_clusters","CO","IL13","IL17","IFN","UT")

write.table(dtM, file="cluster_counts.xls",sep="\t",row.names = F, col.names = T,quote=F)

rownames(dtM) = dtM$seurat_clusters
dtM = dtM[,-c(1)]
dtM$CO = as.numeric(dtM$CO)
dtM$IL13 = as.numeric(dtM$IL13)
dtM$IL17 = as.numeric(dtM$IL17)
dtM$IFN = as.numeric(dtM$IFN)
dtM$UT = as.numeric(dtM$UT)

dtP = dtM/colSums(dtM)
dtP$seurat_clusters = rownames(dtP)
dtPMelt = melt(dtP)
colnames(dtPMelt) = c("seurat_clusters","Treatment","Percent")

dtPMelt$seurat_clusters = factor(dtPMelt$seurat_clusters, levels=c(0,seq(1:27)))

p = ggplot(dtPMelt,aes(seurat_clusters,Percent*100,fill=Treatment))+
  geom_bar(stat="identity",position='dodge', width = .75) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8)) +
  scale_x_discrete(name ="clusters") + 
  scale_y_continuous(name = "percent of cells") + theme_bw()

#p <- ggplot(dtPMelt, aes(x = cluster, y = Percent*100, fill = factor(Day))) +
#  geom_bar(stat="identity", width = 0.7) + 
#  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8)) +
#  scale_x_discrete(name ="cluster") + scale_y_continuous(name = "percent of cells")

png("QC/Barplot.byCluster.png",res = 150, width = 1000,height = 1000)
print(p)
dev.off()




dt  <- aggregate(rep(1, length(paste0(condsInt4$Treatment, condsInt4$cellType))),
                 by=list(condsInt4$Treatment, condsInt4$cellType), sum)
dtPB <- dt[which(dt$Group.2 == "Proliferating_Basal"),]
dtB <- dt[which(dt$Group.2 == "Basal"),]
dtI <- dt[which(dt$Group.2 == "Intermediate"),]
dtS <- dt[which(dt$Group.2 == "Secretory"),]
#dtPC <- dt[which(dt$Group.2 == "Preciliated"),]
dtC <- dt[which(dt$Group.2 == "Ciliated"),]
#dtIon <- dt[which(dt$Group.2 == "Ionocyte"),]
dtU <- dt[which(dt$Group.2 == "Unknown"),]
dtLQC <- dt[which(dt$Group.2 == "LowQ_Ciliated"),]

dtM <- merge(dtPB,dtB,by="Group.1",all = T)
dtM <- merge(dtM,dtI,by="Group.1",all = T)
dtM <- merge(dtM,dtS,by="Group.1",all = T)
#dtM <- merge(dtM,dtPC,by="Group.1",all = T)
dtM <- merge(dtM,dtC,by="Group.1",all = T)
#dtM <- merge(dtM,dtIon,by="Group.1",all = T)
dtM <- merge(dtM,dtU,by="Group.1",all = T)
dtM <- merge(dtM,dtLQC,by="Group.1",all = T)
dtM[is.na(dtM)] = 0
dtM = dtM[,-c(2,4,6,8,10,12,14,16,18)]
colnames(dtM) <- c("Treatment","Proliferating_Basal","Basal","Intermediate","Secretory","Ciliated","Unknown","LowQ_Ciliated")

write.table(dtM, file="treatment_counts.by_cellType.xls",sep="\t",row.names = F, col.names = T,quote=F)

rownames(dtM) = dtM$Treatment
dtM = dtM[,-c(1)]
dtM$Proliferating_Basal = as.numeric(dtM$Proliferating_Basal)
dtM$Basal = as.numeric(dtM$Basal)
dtM$Intermediate = as.numeric(dtM$Intermediate)
dtM$Secretory = as.numeric(dtM$Secretory)
#dtM$Preciliated = as.numeric(dtM$Preciliated)
dtM$Ciliated = as.numeric(dtM$Ciliated)
#dtM$Ionocyte = as.numeric(dtM$Ionocyte)
dtM$Unknown = as.numeric(dtM$Unknown)
dtM$LowQ_Ciliated = as.numeric(dtM$LowQ_Ciliated)

dtP = dtM/rowSums(dtM)
dtP$Treament = rownames(dtP)
dtPMelt = melt(dtP)
colnames(dtPMelt) = c("Treatment","cellType","Percent")

dtPMelt$Treatment = factor(dtPMelt$Treatment, levels=c("UT","IL17","IL13","IFN","CO"))

p = ggplot(dtPMelt,aes(cellType,Percent*100,fill=Treatment))+
  geom_bar(stat="identity",position='dodge', width = .75) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8)) +
  scale_x_discrete(name ="clusters") + 
  scale_y_continuous(name = "percent of cells") + theme_bw()


png("QC/Barplot.cellTypesByTreatment.afterSCT.png",res = 150, width = 1000,height = 1000)
print(p)
dev.off()

pdf("QC/Barplot.cellTypesByTreatment.afterSCT.pdf",width = 6,height = 4)
print(p)
dev.off()

soTreatmentInt <- FindNeighbors(soTreatmentInt) 
soTreatmentInt <- FindClusters(soTreatmentInt, resolution = 2)

DimPlot(soTreatmentInt4, reduction = "umap",label = T)
ggsave("QC/UMAP.new_cluster_Treatment_SCT.pdf", height=7, width=7, units=c("in"))

saveRDS(file = "soTreatmentInt3.rds",soTreatmentInt3)
saveRDS(file = "soTreatmentInt4.rds",soTreatmentInt4)

#### The above is for cell type asignment ####
#### The below is redoing pipeline in less supervised manner carrying over cell IDs from the above####

soTreatmentSCT = readRDS("../Analysis_Batch_Regress_SCT/soTreatmentInt4.rds")
cellTypes = soTreatmentSCT@meta.data$cellType
names(cellTypes) = rownames(soTreatmentSCT@meta.data)

soTreatment = AddMetaData(soTreatment,metadata = cellTypes,col.name = "cellType")

soTreatment <- subset(soTreatment, subset = nFeature_RNA >= 500 & nFeature_RNA < 7500 & 
                        nCount_RNA > 1200 & percent.mito < .5 & percent.ribo < .35)
soTreatment$ConditionWithInd_Batch = paste(soTreatment@meta.data$ConditionWithInd,
                                           soTreatment@meta.data$batch_day, sep="_")


soTreatmentBD = soTreatment

soTreatmentBD <- FindVariableFeatures(object = soTreatmentBD, mean.function = ExpMean, 
                                      dispersion.function = LogVMR, 
                                      selection.method = "vst",nfeatures = 10000)

soTreatmentBD <- ScaleData(soTreatmentBD, 
                           vars.to.regress = c("batch_day","percent.mito",
                                               "nUMI","percent.ribo"), 
                           features =rownames(soTreatmentBD_test@assays$RNA@data))

VlnPlot(object = soTreatmentBD, features = c("nCount_RNA"), pt.size = 0, 
        group.by = "ConditionWithInd", nCol = 1, do.return = TRUE, x.lab.rot = T)
VlnPlot(object = soTreatmentBD, features = c("percent.mito"), pt.size = 0, 
        group.by = "batch", nCol = 3, do.return = TRUE, x.lab.rot = T)
VlnPlot(object = soTreatmentBD, features = c("nCount_RNA"), pt.size = 0, 
        group.by = "batch", nCol = 3, do.return = TRUE, x.lab.rot = T)
VlnPlot(object = soTreatmentBD, features = c("RD.TOTL"), pt.size = 0, 
        group.by = "batch", nCol = 3, do.return = TRUE, x.lab.rot = T)
VlnPlot(object = soTreatmentBD, features = c("percent.ribo"), pt.size = 0, 
        group.by = "batch", nCol = 3, do.return = TRUE, x.lab.rot = T)


#mito.genes2 <- readRDS("../../Erle_Bonser_scRNA_01_31_2019_clean_up/Analysis_BatchRegress/mito.genes.rds")
markerGenes <- c("CDC6","COL17A1","MCM6","SNAI2","TWIST1","KRT17","KRT5","BPIFB1",
                 "SLC5A5","TCN1","ZBP1","IFI44L","ISG15","IFIT2","AGR3","CCDC176",
                 "DNAH7","FOXJ1","CCDC153","CCDC170","SNTN","MUC5B","S100A8","SCGB1A1",
                 "SCGB3A1","FCGBP","IFIT1","IFI6","FOXA3","CST1","MUC5AC","NOS2",
                 "SLC26A4","TFF3","ALOX15","SERPINB10","SOCS1","CCL26","CLDN3","AQP3",
                 "CCND1","KRT7","UPK1B","TSPAN8","TSPAN19","GSTA1","GSTA2","CCDC146",
                 "CCDC173","KLF4","KRT23","SERPINB2","CCDC39","DNAAF1","DNAH10","DNAH11",
                 "DNAH5","DNAH6","POSTN","SPDEF","CFTR")


ccGene <- read.table(file = "../cell_cycle_genes.txt", header = F,sep = "\t")
ccGeneList <- unique(as.vector(ccGene$V2))
ccGeneList <- ccGeneList[which(ccGeneList %in% rownames(soTreatmentBD@assays$RNA@data))]

ribo <- read.table(file="../ribo_list.txt",sep="\t", header=F)

ribo.genes <- unique(as.vector(ribo$V1))
ribo.genes <- ribo.genes[which(ribo.genes %in% rownames(soTreatmentBD@assays$RNA@data))]


markerGenes <- unique(markerGenes)



marker <- soTreatmentBD@assays$RNA@var.features

MTlist <- marker[grep("MT-",marker)]
RPSlist <- marker[grepl("RPS",marker)]
RPLlist <- marker[grepl("RPL",marker)]

#removeList <- c(MTlist, RPSlist, RPLlist)

#marker <- marker[!marker %in% removeList]

marker <- unique(c(marker, markerGenes))

IL13table <- read.table("../IL13_Y_vs_N.xls",header=T, row.names = 1, sep='\t')
IFNatable <- read.table("../IFNa_Y_vs_N.xls",header=T, row.names = 1, sep='\t')
IFNgtable <- read.table("../IFNg_Y_vs_N.xls",header=T, row.names = 1, sep='\t')
#
IL13table <- IL13table[which(IL13table$IL13_Y_vs_N.log2FC > 0),]
IFNatable <- IFNatable[which(IFNatable$IFNa_Y_vs_N.log2FC > 0),]
IFNgtable <- IFNgtable[which(IFNgtable$IFNg_Y_vs_N.log2FC > 0),]

IL13genes <- IL13table$Gene
IFNagenes <- as.vector(IFNatable$Gene)
IFNggenes <- as.vector(IFNgtable$Gene)
IFNgenes <- unique(c(IFNggenes,IFNagenes))

IL13genesShort <- markerGenes[markerGenes %in% IL13genes]
IFNgenesShort <- markerGenes[markerGenes %in% IFNgenes]

IL13genesUniq <- IL13genesShort[!IL13genesShort %in% IFNgenesShort]
IFNgenesUniq <- IFNgenesShort[!IFNgenesShort %in% IL13genesShort]
COgenesUniq <- IL13genesShort[IL13genesShort %in% IFNgenesShort]


basalMarkers <- unique(c("CDC6","COL17A1","MCM6","SNAI2","TWIST1","KRT17","KRT5","AQP3","CCND1","S100A2","CLDN1","TP63","ITGA5","ITGA6","EGFR"))
secretoryMarkers <- unique(c("BPIFB1","TCN1","SCGB1A1","SCGB3A1","TSPAN8","WFDC2"))
ciliatedMarkers <- unique(c("ARG3","CCDC176","ERICH3","DNAH7","FOXJ1","CCDC153","CCDC170","SNTN","TSPAN19","CCDC146","CCDC173","CFAP53","CCDC39","GSTA1","GSTA2"))

bigList <- c("IL13genesUniq","IFNgenesUniq","COgenesUniq","basalMarkers","secretoryMarkers","ciliatedMarkers","ccGeneList","ribo.genes","mito.genes")
bigListNoCC <- c("IL13genesUniq","IFNgenesUniq","COgenesUniq","basalMarkers","secretoryMarkers","ciliatedMarkers")

list <- c(IL13genesUniq,IFNgenesUniq,COgenesUniq,basalMarkers,secretoryMarkers,ciliatedMarkers,ccGeneList,ribo.genes)

for (i in bigList) {
  j <- get(i)
  j <- j[which(j %in% rownames(soTreatmentBD@assays$RNA@scale.data))]
  assign(x = i,value = j)
}


list <- unique(list)

list <- list[which(list %in% rownames(soTreatmentBD@assays$RNA@scale.data))]
list <- unique(c(soTreatmentBD@assays$RNA@var.features, list))

soTreatmentBD@assays$RNA@var.features <- list


soTreatmentBD <- RunPCA(object = soTreatmentBD, ndims.print = 1:6, 
                        nfeatures.print = 5, npcs = 125,
                        features = VariableFeatures(object = soTreatmentBD))

VizDimLoadings(soTreatmentBD, dims = 1:2, reduction = "pca")

DimPlot(soTreatmentBD, reduction = "pca")

DimHeatmap(soTreatmentBD, dims = 1:9, cells = 500, balanced = TRUE)
DimHeatmap(soTreatmentBD, dims = 10:18, cells = 500, balanced = TRUE)
DimHeatmap(soTreatmentBD, dims = 19:27, cells = 500, balanced = TRUE)
DimHeatmap(soTreatmentBD, dims = 28:36, cells = 500, balanced = TRUE)
DimHeatmap(soTreatmentBD, dims = 37:45, cells = 500, balanced = TRUE)
DimHeatmap(soTreatmentBD, dims = 46:54, cells = 500, balanced = TRUE)
DimHeatmap(soTreatmentBD, dims = 55:63, cells = 500, balanced = TRUE)
DimHeatmap(soTreatmentBD, dims = 64:72, cells = 500, balanced = TRUE)
DimHeatmap(soTreatmentBD, dims = 73:81, cells = 500, balanced = TRUE)
DimHeatmap(soTreatmentBD, dims = 82:90, cells = 500, balanced = TRUE)
DimHeatmap(soTreatmentBD, dims = 91:99, cells = 500, balanced = TRUE)
DimHeatmap(soTreatmentBD, dims = 100:108, cells = 500, balanced = TRUE)
DimHeatmap(soTreatmentBD, dims = 109:117, cells = 500, balanced = TRUE)
DimHeatmap(soTreatmentBD, dims = 118:125, cells = 500, balanced = TRUE)

soTreatmentBD <- JackStraw(soTreatmentBD, num.replicate = 100,dims = 125)
soTreatmentBD <- ScoreJackStraw(soTreatmentBD,dims = 1:125)

JackStrawPlot(soTreatmentBD, dims = 1:15)
JackStrawPlot(soTreatmentBD, dims = 16:30)
JackStrawPlot(soTreatmentBD, dims = 31:45)
JackStrawPlot(soTreatmentBD, dims = 46:60)
JackStrawPlot(soTreatmentBD, dims = 61:75)
JackStrawPlot(soTreatmentBD, dims = 76:90)
JackStrawPlot(soTreatmentBD, dims = 91:99)
JackStrawPlot(soTreatmentBD, dims = 100:108)
JackStrawPlot(soTreatmentBD, dims = 109:117)
JackStrawPlot(soTreatmentBD, dims = 118:125)
#ElbowPlot(soTreatmentBD)

pVals = as.data.frame(soTreatmentBD@reductions$pca@jackstraw$overall.p.values)
pVals$pAdj = p.adjust(pVals$Score)
cutPC = pVals[(pVals$pAdj <= 0.1),]
cutPC = max(cutPC$PC)

soTreatmentBD = soTreatmentBD2

cutPC = 50

soTreatmentBD <- FindNeighbors(soTreatmentBD, dims = 1:cutPC)
soTreatmentBD<- RunUMAP(soTreatmentBD, dims = 1:cutPC)

DimPlot(soTreatmentBD, reduction = "umap",label = T,group.by = "Treatment")
ggsave("QC/UMAP.BD.Treatment.named.pdf", height=7, width=7, units=c("in"))
DimPlot(soTreatmentBD, reduction = "umap",label = T,group.by = "Ind")
ggsave("QC/UMAP.BD.Subject.named.pdf", height=7, width=7, units=c("in"))
DimPlot(soTreatmentBD, reduction = "umap",label = T,group.by = "cellType")
ggsave("QC/UMAP.BD.cellType.named.pdf", height=7, width=7, units=c("in"))
DimPlot(soTreatmentBD2, reduction = "umap",label = T,group.by = "batch")
ggsave("QC/UMAP.BD.batch.named.pdf", height=7, width=7, units=c("in"))

FeaturePlot(soTreatmentBD,"XBP1")




soTreatmentBD <- FindClusters(soTreatmentBD, resolution = 2)
DimPlot(soTreatmentBD, reduction = "umap",group.by = "RNA_snn_res.2")
ggsave("QC/UMAP.BD.unnamed.2.pdf", height=7, width=7, units=c("in"))
DimPlot(soTreatmentBD, reduction = "umap",group.by = "RNA_snn_res.2",label = T)
ggsave("QC/UMAP.BD.named.2.pdf", height=7, width=7, units=c("in"))

soTreatmentBD <- FindClusters(soTreatmentBD, resolution = 1.5)
DimPlot(soTreatmentBD, reduction = "umap",group.by = "RNA_snn_res.1.5")
ggsave("QC/UMAP.BD.unnamed.1.5.pdf", height=7, width=7, units=c("in"))
DimPlot(soTreatmentBD, reduction = "umap",group.by = "RNA_snn_res.1.5",label = T)
ggsave("QC/UMAP.BD.named.1.5.pdf", height=7, width=7, units=c("in"))

soTreatmentBD <- FindClusters(soTreatmentBD, resolution = 1.25)
DimPlot(soTreatmentBD, reduction = "umap",group.by = "RNA_snn_res.1.25")
ggsave("QC/UMAP.BD.unnamed.1.25.pdf", height=7, width=7, units=c("in"))
DimPlot(soTreatmentBD, reduction = "umap",group.by = "RNA_snn_res.1.25",label = T)
ggsave("QC/UMAP.BD.named.1.25.pdf", height=7, width=7, units=c("in"))

soTreatmentBD <- FindClusters(soTreatmentBD, resolution = 1)
DimPlot(soTreatmentBD, reduction = "umap",group.by = "RNA_snn_res.1")
ggsave("QC/UMAP.BD.unnamed.1.pdf", height=7, width=7, units=c("in"))
DimPlot(soTreatmentBD, reduction = "umap",group.by = "RNA_snn_res.1",label = T)
ggsave("QC/UMAP.BD.named.1.pdf", height=7, width=7, units=c("in"))

soTreatmentBD <- FindClusters(soTreatmentBD, resolution = 0.8)
DimPlot(soTreatmentBD, reduction = "umap",group.by = "RNA_snn_res.0.8")
ggsave("QC/UMAP.BD.unnamed.0.8.pdf", height=7, width=7, units=c("in"))
DimPlot(soTreatmentBD, reduction = "umap",group.by = "RNA_snn_res.0.8",label = T)
ggsave("QC/UMAP.BD.named.0.8.pdf", height=7, width=7, units=c("in"))



soTreatmentBD <- RunTSNE(object = soTreatmentBD, dims.use = 1:cutPC, do.fast = TRUE)
DimPlot(object = soTreatmentBD, do.label = T,reduction = "tsne",group.by = "RNA_snn_res.1.25")
ggsave("QC/TSNE.BD.named.pdf", height=7, width=7, units=c("in"))
DimPlot(object = soTreatmentBD, do.label = F,reduction = "tsne")
ggsave("QC/TSNE.BD.unnamed.pdf", height=7, width=7, units=c("in"))
DimPlot(soTreatmentBD2, reduction = "umap",label = T,group.by = "cellType")
ggsave("QC/TSNE.BD.Treatment.named.pdf", height=7, width=7, units=c("in"))



FeaturePlot(soTreatmentBD,
            features = c("nCount_RNA","nFeature_RNA","percent.mito","percent.ribo"))
ggsave("QC/UMAP.BD.QC_features.pdf", height=7, width=7, units=c("in"))
VlnPlot(soTreatmentBD,
        features = c("nCount_RNA","nFeature_RNA"),
        group.by = "RNA_snn_res.1.25",ncol = 1,pt.size = 0)
ggsave("QC/VlnPlt.BD.QC_features.pdf", height=7, width=7, units=c("in"))
VlnPlot(soTreatmentBD,
        features = c("percent.mito","percent.ribo"),
        group.by = "RNA_snn_res.1.25",ncol = 1,pt.size = 0)
ggsave("QC/VlnPlt.BD.QC_MitoRiboPct.pdf", height=7, width=7, units=c("in"))

VlnPlot(soTreatmentBD3,features="XBP1",group.by = "cellType", split.by = "Treatment",pt.size = 0)

active.ident = soTreatmentBD2@meta.data$RNA_snn_res.1.25
names(active.ident)  = rownames(soTreatmentBD2@meta.data)
soTreatmentBD2@active.ident = active.ident

#cluster20 = FindMarkers(soTreatmentBD,ident.1 = 20,logfc.threshold = 0.5,min.cells.feature = .25,only.pos = T)
#cluster25 = FindMarkers(soTreatmentBD,ident.1 = 25,logfc.threshold = 0.5,min.cells.feature = .25,only.pos = T)
#cluster06 = FindMarkers(soTreatmentBD,ident.1 = 6,logfc.threshold = 0.5,min.cells.feature = .25,only.pos = T)
#cluster09 = FindMarkers(soTreatmentBD,ident.1 = 9,logfc.threshold = 0.5,min.cells.feature = .25,only.pos = T)
#cluster19 = FindMarkers(soTreatmentBD,ident.1 = 19,logfc.threshold = 0.5,min.cells.feature = .25,only.pos = T)
cluster23 = FindMarkers(soTreatmentBD,ident.1 = 23,logfc.threshold = 0.5,min.cells.feature = .25,only.pos = T)
cluster22 = FindMarkers(soTreatmentBD,ident.1 = 22,logfc.threshold = 0.5,min.cells.feature = .25,only.pos = T)
cluster25 = FindMarkers(soTreatmentBD,ident.1 = 25,logfc.threshold = 0.5,min.cells.feature = .25,only.pos = T)
cluster26 = FindMarkers(soTreatmentBD,ident.1 = 26,logfc.threshold = 0.5,min.cells.feature = .25,only.pos = T)


cells <- rownames(soTreatmentBD@meta.data[!(soTreatmentBD@meta.data$ConditionWithInd_Batch == "Sub3_D24_IL17_B2"),])
soTreatmentBD2 <- subset(soTreatmentBD,cells = cells)

soTreatmentBD2@meta.data$cellType = factor(soTreatmentBD2@meta.data$cellType,levels = c(levels=c("Proliferating_Basal","Basal","Intermediate","Secretory","Preciliated","Ciliated","Ionocyte","Unknown","LowQ_Ciliated")))

soTreatmentBD2@meta.data$cellType[soTreatmentBD2@meta.data$RNA_snn_res.1.25 == 25] = "Preciliated"
soTreatmentBD2@meta.data$cellType[soTreatmentBD2@meta.data$RNA_snn_res.1.25 == 26] = "Ionocyte"

AllMarkersBD2 = FindAllMarkers(soTreatmentBD2,only.pos = T,logfc.threshold = .5,min.pct = .25)

dir.create("DE/BD2_integrated_markers/")
top10List = c()

sigcounts <- as.data.frame(setNames(replicate(3,numeric(0), simplify = F),letters[0:3]))


for (i in unique(AllMarkersBD2$cluster)) {
  print(i)
  Tmp <-  AllMarkersBD2[(AllMarkersBD2$cluster == i),]
  print(head(Tmp))
  padjCount = nrow(Tmp[Tmp$p_val_adj <= 0.1,])
  pvalCount = nrow(Tmp[Tmp$p_val <= 0.1,])
  rowList = data.frame(row.names = i,cluster=i,fdr=padjCount,raw_p = pvalCount)
  sigcounts = rbind(sigcounts,rowList)
  write.table(Tmp, sep="\t",quote = F,row.names = T,col.names = NA,
              file = paste("DE/BD2_integrated_markers/","Cluster_",i,"Markers.txt"))
  Tmp2 = head(Tmp,n=6)
  geneList = Tmp2$gene
  top10List = c(top10List,geneList)
  
}

DoHeatmap(soTreatmentBD2, features = top10List,size = 4) + NoLegend()
ggsave("DE/BD2_integrated_markers/Top10.heatmap.png",height = 16,width = 16,units = "in",dpi = 150)


DimPlot(soTreatmentBD2, reduction = "umap",label = T,group.by = "cellType")
ggsave("QC/UMAP.BD2.cellType.named.additions.lessbig.pdf", height=10, width=12, units=c("in"))

DimPlot(soTreatmentBD2, reduction = "umap",label = T,group.by = "Treatment")
ggsave("QC/UMAP.BD2.Treatment.named.additions.lessbig.pdf", height=10, width=12, units=c("in"))

DimPlot(soTreatmentBD2, reduction = "umap",label = T,group.by = "RNA_snn_res.1.25")
ggsave("QC/UMAP.BD2.RNA_snn_res.1.25.named.additions.lessbig.pdf", height=10, width=12, units=c("in"))

DotPlot(soTreatmentBD2, features = top10List, dot.scale =  4, split.by = "Treatment", 
        cols = c("gray","blue","green","red","yellow")) + RotatedAxis()
ggsave("DE/BD2_integrated_markers/Top10.DotPlot.pdf",height = 16,width = 20,units = "in",dpi = 300)

DotPlot(soTreatmentBD2, features = top10List, dot.scale =  4) + RotatedAxis()
ggsave("DE/BD2_integrated_markers/Top10.DotPlot.noSpit.pdf",height = 6,width = 20,units = "in",dpi = 300)


data = soTreatmentBD2@assays$RNA@data
data = as.data.frame(data)
condsDB2 = soTreatmentBD2@meta.data

dotMat <- as.data.frame(setNames(replicate(4,numeric(0), simplify = F),letters[0:4]))
dotMatScaled <- as.data.frame(setNames(replicate(4,numeric(0), simplify = F),letters[0:4]))


top10ListFilter = top10List[(top10List %in% rownames(data))]
top10ListFilter = unique(top10ListFilter)

dataAves = rowMeans(data[top10ListFilter,])

for (i in levels(soTreatmentBD2@meta.data$RNA_snn_res.1.25)){
  condsTmp = condsDB2[condsDB2$RNA_snn_res.1.25 == i,]
  cellsTmp = rownames(condsTmp)
  dataTmp = data[top10ListFilter,cellsTmp]
  data10 = dataTmp
  data10[data10 > 0 ] = 1
  pctExp = rowSums(data10)/ncol(data10)
  aveExp = rowMeans(dataTmp)
  scaledAveExp = aveExp - dataAves
  dF = data.frame(row.names = paste(names(pctExp),i,sep="_"),gene = names(pctExp),cluster = i,aveExp = aveExp,pctExp = pctExp)
  dF_scaled = data.frame(row.names = paste(names(pctExp),i,sep="_"),gene = names(pctExp),cluster = i,aveExp = scaledAveExp,pctExp = pctExp)
  dF = dF[match(top10ListFilter, dF$gene),]
  dF_scaled = dF_scaled[match(top10ListFilter, dF_scaled$gene),]
  #print(dF[dF$gene == "CFTR",])
  #rownames(dF) = paste(dF$gene,dF$cluster,sep="_")
  dotMat = rbind(dotMat,dF)
  dotMatScaled = rbind(dotMatScaled,dF_scaled)
}

dotMat$gene = factor(dotMat$gene,ordered = T,levels = rev(top10ListFilter))
dotMatScaled$gene = factor(dotMatScaled$gene,ordered = T,levels = rev(top10ListFilter))

#levels(dotMat$gene) = top10ListFilter

dotMatScaled$aveExp[dotMatScaled$aveExp > 2] = 2
dotMatScaled$aveExp[dotMatScaled$aveExp < -2] = -2

ggplot(dotMatScaled, aes(y = gene, x = cluster)) +
  geom_point(aes(size = pctExp, colour=aveExp)) +
  scale_size_area() +
  scale_colour_gradientn(colours = c("purple","grey","orange"),limits=c(-2,2)) +
  theme_bw()
ggsave("DE/BD2_integrated_markers/test.pdf",height = 20,width = 15,units = "in",dpi = 300)

#### cell counts ####

soTreatmentBD2@meta.data$Treatment_Cluster <- paste(soTreatmentBD2@meta.data$Treatment, soTreatmentBD2@meta.data$RNA_snn_res.1.25, sep="_")

conds <- soTreatmentBD2@meta.data



#conds$genotype <- factor(conds$genotype, levels = c("WT","KO"))
conds$Treatment = factor(conds$Treatment)

p <- ggplot(conds, aes(x = RNA_snn_res.1.25, fill=Treatment, color=Treatment)) + 
  geom_bar(position="dodge", width = .75) + theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8)) +
  scale_x_discrete(name ="cluster")

png("QC/Barplot.TreatmentByCluster.png",res = 150, width = 1000,height = 1000)
print(p)
dev.off()

conds$cellType = factor(conds$cellType,levels=c("Proliferating_Basal","Basal","Intermediate","Secretory","Preciliated","Ciliated","Ionocyte","Unknown","LowQ_Ciliated"))

p <- ggplot(conds, aes(x = cellType, fill=Treatment, color=Treatment)) + 
  geom_bar(position="dodge", width = .75) + theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8)) +
  scale_x_discrete(name ="cluster")

png("QC/Barplot.byCelltype.png",res = 150, width = 1000,height = 1000)
print(p)
dev.off()

conds$cluster = conds$RNA_snn_res.1.25

conds = conds[conds$cellType %in% c("Proliferating_Basal","Basal","Intermediate","Secretory","Preciliated","Ciliated","Ionocyte"),]
conds$cellType = factor(conds$cellType,levels=c("Proliferating_Basal","Basal","Intermediate","Secretory","Preciliated","Ciliated","Ionocyte"))
conds$Treatment = factor(conds$Treatment, levels=c("UT","IL17","IL13","IFN","CO"))

conds.new <- conds %>% 
  group_by(Treatment,cellType,Ind) %>% 
  dplyr::summarise(n = n()) %>%
  mutate(freq = n / sum(n))

cellCountsBySubAndTreat <- as.data.frame(setNames(replicate(3,numeric(0), simplify = F),letters[0:3]))
cellCountsBySubAndTreatNum <- as.data.frame(setNames(replicate(2,numeric(0), simplify = F),letters[0:2]))
counter = 0
for (i in c("UT","IL17","IL13","IFN","CO")){
  for (j in c("Sub1","Sub2","Sub3","Sub6")){
    conds.new.tmp = conds.new[conds.new$Treatment == i,]
    conds.new.tmp = conds.new.tmp[conds.new.tmp$Ind == j,]
    sumN = sum(conds.new.tmp$n)
    for (k in levels(conds$cellType)){
      counter = counter + 1
      kn = as.numeric(conds.new.tmp[conds.new.tmp$cellType == k,"n"])
      freqn = kn/sumN
      m = t(as.matrix(c(i,j,k)))
      rownames(m) = counter
      colnames(m) = c("Treatment","Subject","Celltype")
      cellCountsBySubAndTreat = rbind(cellCountsBySubAndTreat,m)
      cellCountsBySubAndTreatNum = rbind(cellCountsBySubAndTreatNum,c(kn,freqn))
    }
  }
}

cellTypeV = soTreatmentBD2@meta.data$cellType
names(cellTypeV) = rownames(soTreatmentBD2@meta.data)
soTreatmentBD2@active.ident = cellTypeV




cellCountsBySubAndTreatNum[is.na(cellCountsBySubAndTreatNum)] <- 0

cellCountsBySubAndTreat = cbind(cellCountsBySubAndTreat,cellCountsBySubAndTreatNum)

colnames(cellCountsBySubAndTreat) = c("Treatment","Subject","Celltype","N","Percent")
cellCountsBySubAndTreat$Percent = cellCountsBySubAndTreat$Percent * 100

ggplot(cellCountsBySubAndTreat,aes(y=Percent,x=Celltype,colour=Treatment),aes_string(shape=Subject)) +
  #geom_point() + 
  geom_boxplot(aes(color=Treatment),outlier.shape = NA, coef = 0) +
  geom_point(position = position_jitterdodge(),size=.7,alpha=0.6) +
  theme_bw()
ggsave("QC/SCT.cellTypePercentBySample.pdf",width = 8, height = 6,dpi = 150,units = "in")
#geom_jitter(aes(shape=Subject),position = position_jitterdodge(dodge.width = .75, jitter.width = 0.1,jitter.height = 0))

ggplot(cellCountsBySubAndTreat,aes(y=Percent,x=Celltype)) +
  #geom_point() + 
  geom_boxplot(aes(fill=Treatment)) 
geom_jitter(position = position_jitterdodge(dodge.width = 1.2, jitter.width = 0.5,jitter.height = 0))

ggplot(conds.new,aes(cellType,freq,fill=Ind))+
  geom_bar(stat="identity",position='dodge', width = .75) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8)) +
  scale_x_discrete(name ="cellType") + 
  scale_y_continuous(name = "percent of cells") + theme_bw()
ggsave("QC/SCT.UT.cellTypePercent.png",width = 6, height = 4,dpi = 150,units = "in")


conds.new <- conds %>% 
  group_by(Treatment,cluster) %>% 
  dplyr::summarise(n = n()) %>%
  mutate(freq = n / sum(n))

p <- ggplot(conds.new,aes(cluster,freq,fill=Treatment))+
  geom_bar(stat="identity",position='dodge', width = .75) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8)) +
  scale_x_discrete(name ="cluster") + scale_y_continuous(name = "percent of cells")


conds.new <- conds %>% 
  group_by(Treatment,cellType) %>% 
  dplyr::summarise(n = n()) %>%
  mutate(freq = n / sum(n))

p <- ggplot(conds.new,aes(cellType,freq,fill=Treatment))+
  geom_bar(stat="identity",position='dodge', width = .75) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8)) +
  scale_x_discrete(name ="cluster") + scale_y_continuous(name = "percent of cells")




dt  <- aggregate(rep(1, length(paste0(conds$cluster, conds$Treatment))),
                 by=list(conds$cluster, conds$Treatment), sum)
dtCO <- dt[which(dt$Group.2 == "CO"),]
dtIL13 <- dt[which(dt$Group.2 == "IL13"),]
dtIL17 <- dt[which(dt$Group.2 == "IL17"),]
dtIFN <- dt[which(dt$Group.2 == "IFN"),]
dtUT <- dt[which(dt$Group.2 == "UT"),]

dtM <- merge(dtCO,dtIL13,by="Group.1",all = T)
dtM <- merge(dtM,dtIL17,by="Group.1",all = T)
dtM <- merge(dtM,dtIFN,by="Group.1",all = T)
dtM <- merge(dtM,dtUT,by="Group.1",all = T)
dtM[is.na(dtM)] = 0
dtM = dtM[,-c(2,4,6,8,10)]
colnames(dtM) <- c("cluster","CO","IL13","IL17","IFN","UT")

write.table(dtM, file="cluster_counts.xls",sep="\t",row.names = F, col.names = T,quote=F)

rownames(dtM) = dtM$cluster
dtM = dtM[,-c(1)]
dtM$CO = as.numeric(dtM$CO)
dtM$IL13 = as.numeric(dtM$IL13)
dtM$IL17 = as.numeric(dtM$IL17)
dtM$IFN = as.numeric(dtM$IFN)
dtM$UT = as.numeric(dtM$UT)

dtP = dtM/rowSums(dtM)
dtP$cluster = rownames(dtP)
dtPMelt = melt(dtP)
colnames(dtPMelt) = c("cluster","Treatment","Percent")

dtPMelt$cluster = factor(dtPMelt$cluster, levels=c(0,seq(1:27)))

p = ggplot(dtPMelt, aes(x = cluster, y = Percent*100, fill = factor(Treatment))) +
  geom_bar(stat="identity", width = 0.7) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8)) +
  scale_x_discrete(name ="cluster") + scale_y_continuous(name = "percent of cells")

#p <- ggplot(dtPMelt, aes(x = cluster, y = Percent*100, fill = factor(Day))) +
#  geom_bar(stat="identity", width = 0.7) + 
#  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8)) +
#  scale_x_discrete(name ="cluster") + scale_y_continuous(name = "percent of cells")

png("QC/Barplot.byCluster.png",res = 150, width = 1000,height = 1000)
print(p)
dev.off()




dt  <- aggregate(rep(1, length(paste0(conds$Treatment, conds$cellType))),
                 by=list(conds$Treatment, conds$cellType), sum)
dtPB <- dt[which(dt$Group.2 == "Proliferating_Basal"),]
dtB <- dt[which(dt$Group.2 == "Basal"),]
dtI <- dt[which(dt$Group.2 == "Intermediate"),]
dtS <- dt[which(dt$Group.2 == "Secretory"),]
dtPC <- dt[which(dt$Group.2 == "Preciliated"),]
dtC <- dt[which(dt$Group.2 == "Ciliated"),]
dtIon <- dt[which(dt$Group.2 == "Ionocyte"),]
dtU <- dt[which(dt$Group.2 == "Unknown"),]
dtLQC <- dt[which(dt$Group.2 == "LowQ_Ciliated"),]

dtM <- merge(dtPB,dtB,by="Group.1",all = T)
dtM <- merge(dtM,dtI,by="Group.1",all = T)
dtM <- merge(dtM,dtS,by="Group.1",all = T)
dtM <- merge(dtM,dtPC,by="Group.1",all = T)
dtM <- merge(dtM,dtC,by="Group.1",all = T)
dtM <- merge(dtM,dtIon,by="Group.1",all = T)
dtM <- merge(dtM,dtU,by="Group.1",all = T)
dtM <- merge(dtM,dtLQC,by="Group.1",all = T)
dtM[is.na(dtM)] = 0
dtM = dtM[,-c(2,4,6,8,10,12,14,16,18)]
colnames(dtM) <- c("Treatment","Proliferating_Basal","Basal","Intermediate","Secretory","Preciliated","Ciliated","Ionocyte","Unknown","LowQ_Ciliated")

write.table(dtM, file="treatment_counts.by_cellType.xls",sep="\t",row.names = F, col.names = T,quote=F)

rownames(dtM) = dtM$Treatment
dtM = dtM[,-c(1)]
dtM$Proliferating_Basal = as.numeric(dtM$Proliferating_Basal)
dtM$Basal = as.numeric(dtM$Basal)
dtM$Intermediate = as.numeric(dtM$Intermediate)
dtM$Secretory = as.numeric(dtM$Secretory)
dtM$Preciliated = as.numeric(dtM$Preciliated)
dtM$Ciliated = as.numeric(dtM$Ciliated)
dtM$Ionocyte = as.numeric(dtM$Ionocyte)
dtM$Unknown = as.numeric(dtM$Unknown)
dtM$LowQ_Ciliated = as.numeric(dtM$LowQ_Ciliated)

dtP = dtM/rowSums(dtM)
dtP$Treament = rownames(dtP)
dtPMelt = melt(dtP)
colnames(dtPMelt) = c("Treatment","cellType","Percent")

dtPMelt$Treatment = factor(dtPMelt$Treatment, levels=c("UT","IL17","IL13","IFN","CO"))

p = ggplot(dtPMelt,aes(cellType,Percent*100,fill=Treatment))+
  geom_bar(stat="identity",position='dodge', width = .75) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8)) +
  scale_x_discrete(name ="cellType") + 
  scale_y_continuous(name = "percent of cells by treatment") 


#p <- ggplot(dtPMelt, aes(x = cluster, y = Percent*100, fill = factor(Day))) +
#  geom_bar(stat="identity", width = 0.7) + 
#  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8)) +
#  scale_x_discrete(name ="cluster") + scale_y_continuous(name = "percent of cells")

pdf("QC/BD2.Barplot.cellTypesByTreatment.pdf")
print(p)
dev.off()



####### Subject specific plots #######

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


#df_to_plot <- soTreatmentBD2@dr$tsne@cell.embeddings
df_to_plot <- soTreatmentBD2@reductions$pca@cell.embeddings
soTreatmentBD2@meta.data$umap_1 <- soTreatmentBD2@reductions$umap@cell.embeddings[,1]
soTreatmentBD2@meta.data$umap_2 <- soTreatmentBD2@reductions$umap@cell.embeddings[,2]
df_to_plot2 <- merge(df_to_plot, soTreatmentBD2@meta.data, by='row.names')
rownames(df_to_plot2) <- df_to_plot2$Row.names 
df_to_plot2 <- df_to_plot2[,-c(1)]

df_to_plot2 <- merge(df_to_plot2, as.matrix(t(soTreatmentBD2@assays$RNA@counts)), by='row.names')
rownames(df_to_plot2) <- df_to_plot2[,c(1)]
df_to_plot2 <- df_to_plot2[,-c(1)]

df_to_plot2$shape <- df_to_plot2$Ind

df_to_plot2$shape <- gsub("Sub1","0", df_to_plot2$shape)
df_to_plot2$shape <- gsub("Sub2","1", df_to_plot2$shape)
df_to_plot2$shape <- gsub("Sub3","2", df_to_plot2$shape)
df_to_plot2$shape <- gsub("Sub6","4", df_to_plot2$shape)

df_to_plot2$shape <- as.integer(df_to_plot2$shape)

#condOrder <- as.vector(unique(df_to_plot2$ConditionWithInd))[c(2,5,3,4,8,1,7,6)]

condsOrder <- c("D24_UT","D24_IFN","D24_IL13","D24_IL17","D24_CO")
subOrder <- c("Sub1","Sub2","Sub3","Sub6")
shapes <- c(0,1,2,4)


df_to_plot2$Condition <- factor(df_to_plot2$Condition, levels=condsOrder)
df_to_plot2$Ind <- factor(df_to_plot2$Ind, levels=subOrder)


scatCols <- gg_color_hue(5)
greyCols <- c(rep("grey90", 5))


ggplot2::ggplot(df_to_plot2, aes_string(x = "umap_1", y = "umap_2")) +
  geom_point(aes(color=Condition, shape=Ind), size=2, alpha=0.4) +
  xlab("Dimension 1") + scale_color_manual(values=scatCols, name="Condition") +
  scale_shape_manual(name = "Subject", values = c(0, 1, 2, 4)) +
  ylab("Dimension 2") + ggtitle("All Samples") +  
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text = element_text(size=6), axis.title=element_text(size=7))
ggsave("QC/UMAP.byCondition.pdf", height=7, width=7, units=c("in"))

ggplot2::ggplot(df_to_plot2, aes_string(x = "umap_1", y = "umap_2")) +
  geom_point(aes(color=Condition), size=2, alpha=0.4) +
  xlab("Dimension 1") + scale_color_manual(values=scatCols, name="Condition") +
  #scale_shape_manual(name = "Subject", values = c(0, 1, 2, 4)) +
  ylab("Dimension 2") + ggtitle("All Samples") +  
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text = element_text(size=6), axis.title=element_text(size=7)) +
  theme_bw()
ggsave("QC/TSNE.byCondition.bw.pdf", height=7, width=7, units=c("in"))

p0 <- 	ggplot2::ggplot(df_to_plot2, aes_string(x = "umap_1", y = "umap_2")) +
  geom_point(aes(color=Condition, shape=Ind), size=2, alpha=0.4) +
  xlab("Dimension 1") + scale_color_manual(values=scatCols, name="Condition") +
  scale_shape_manual(name = "Subject", values = c(0, 1, 2, 4)) +
  ylab("Dimension 2") + #ggtitle("All Samples") +  
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text = element_text(size=6), axis.title=element_text(size=7), legend.position="none")

for (i in subOrder) {
  text <- i
  pTemp <- ggplot() + annotate("text", x = 4, y = 25, size=8, label = i) + 
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y =element_blank(),
          axis.line = element_blank())
  assign(paste("p",i,"_label",sep=""), pTemp)
}

for (i in condsOrder) {
  text <- i
  pTemp <- ggplot() + annotate("text", x = 4, y = 25, size=8, label = i) + 
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y =element_blank(),
          axis.line = element_blank())
  assign(paste("p",i,"_label",sep=""), pTemp)
}

i = 0
for (j in condsOrder) {
  i = i + 1
  #scatColsTemp <- scatCols
  #scatColsTemp[i] <- scatCols[i]
  cellCount <- nrow(df_to_plot2[which(df_to_plot2$Condition == j),])
  df_to_plot4 <- df_to_plot2[which(df_to_plot2[,"Condition"] == j),]
  #tmpDF2 <- df_to_plot2[which(df_to_plot2[,"cond2"] != j),]
  #df_to_plot4  <- rbind(tmpDF2, tmpDF1)
  p <- ggplot(df_to_plot4, aes_string(x = "umap_1", y = "umap_2")) +
    geom_point(aes(color=Condition, shape=Ind), size=2, alpha=0.6) +
    #geom_point(data = subset(df_to_plot4, cond2 == j), aes(color=cond2), size=2, alpha=0.6) +
    xlab("Dimension 1") + scale_color_manual(values=scatCols[i], name="Condition") +
    scale_shape_manual(name = "Subject", values = c(0, 1, 2, 4)) +
    ylab("Dimension 2") + #ggtitle(paste(j, ", cell count = ", cellCount, sep="" )) + 
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text = element_text(size=6), axis.title=element_text(size=7), legend.position="none") 
  print(p)
  assign(paste("p",j,sep=""), p)
  print(paste("p",j,sep=""))
}


i = 0
for (j in subOrder) {
  i = i + 1
  #scatColsTemp <- scatCols
  #scatColsTemp[i] <- scatCols[i]
  cellCount <- nrow(df_to_plot2[which(df_to_plot2$Ind == j),])
  df_to_plot4 <- df_to_plot2[which(df_to_plot2[,"Ind"] == j),]
  #tmpDF2 <- df_to_plot2[which(df_to_plot2[,"cond2"] != j),]
  #df_to_plot4  <- rbind(tmpDF2, tmpDF1)
  p <- ggplot(df_to_plot4, aes_string(x = "umap_1", y = "umap_2")) +
    geom_point(aes(color=Condition, shape=Ind), size=2, alpha=0.6) +
    #geom_point(data = subset(df_to_plot4, cond2 == j), aes(color=cond2), size=2, alpha=0.6) +
    xlab("Dimension 1") + scale_color_manual(values=scatCols, name="Condition") +
    scale_shape_manual(name = "Subject", values = shapes[i]) +
    ylab("Dimension 2") + #ggtitle(paste(j, ", cell count = ", cellCount, sep="" )) + 
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text = element_text(size=6), axis.title=element_text(size=7), legend.position="none") 
  print(p)
  assign(paste("p",j,sep=""), p)
  print(paste("p",j,sep=""))
}

i = 0
for (j in condsOrder) {
  l = 0
  i = i + 1
  for (k in subOrder) {
    l = l + 1
    scatColsTemp <- greyCols
    scatColsTemp[i] <- scatCols[i]
    cellCount <- nrow(df_to_plot2[which(df_to_plot2$Condition == j),])
    df_to_plot4 <- df_to_plot2[which(df_to_plot2[,"Condition"] == j & df_to_plot2[,"Ind"] == k),]
    
    #df_to_plot4 <- rbind(tmpDF2, tmpDF1)
    p <- ggplot(df_to_plot4, aes_string(x = "umap_1", y = "umap_2")) +
      geom_point(aes(color=Condition, shape = Ind), size=2, alpha=0.6) +
      #geom_point(data = subset(df_to_plot4, cond2 == j), aes(color=cond2), size=2, alpha=0.6) +
      scale_shape_manual(name = "Subject", values = shapes[l]) +
      xlab("Dimension 1") + scale_color_manual(values=scatCols[i], name="Condition") +
      ylab("Dimension 2") + #ggtitle(paste(condsOrder[i], ", cell count = ", cellCount, sep="" )) + 
      theme(plot.title = element_text(hjust = 0.5)) +
      theme(axis.text = element_text(size=6), axis.title=element_text(size=7), legend.position="none") 
    print(p)
    assign(paste("p",j,k,sep="_"), p)
    print(paste("p",j,k,sep="_"))
  }
}


plot_grid(p0, pD24_UT, pD24_IFN,pD24_IL17,pD24_IL13,pD24_CO,
          pSub1, p_D24_UT_Sub1, p_D24_IFN_Sub1, p_D24_IL17_Sub1, p_D24_IL13_Sub1, p_D24_CO_Sub1,
          pSub2, p_D24_UT_Sub2, p_D24_IFN_Sub2, p_D24_IL17_Sub2, p_D24_IL13_Sub2, p_D24_CO_Sub2,
          pSub3, p_D24_UT_Sub3, p_D24_IFN_Sub3, p_D24_IL17_Sub3, p_D24_IL13_Sub3, p_D24_CO_Sub3,
          pSub6, p_D24_UT_Sub6, p_D24_IFN_Sub6, p_D24_IL17_Sub6, p_D24_IL13_Sub6, p_D24_CO_Sub6,
          ncol = 6, nrow = 5)

ggsave("QC/UMAP.Subject_Tissue_groupings.pdf", height=13, width=15, units=c("in"))


condsOrder <- c("D24_UT","D24_IFN","D24_IL13","D24_IL17","D24_CO")
subOrder <- c("Sub1","Sub2","Sub3","Sub6")
shapes <- c(0,1,2,4,7)

scatCols <- gg_color_hue(4)
greyCols <- c(rep("grey90", 4))


ggplot2::ggplot(df_to_plot2, aes_string(x = "umap_1", y = "umap_2")) +
  geom_point(aes(color=Ind), size=2, alpha=0.4) +
  xlab("Dimension 1") + scale_color_manual(values=scatCols, name="Condition") +
  #scale_shape_manual(name = "Subject", values = c(0, 1, 2, 4,7)) +
  ylab("Dimension 2") + ggtitle("All Samples") +  
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text = element_text(size=6), axis.title=element_text(size=7))
ggsave("QC/TSNE.bySubject.pdf", height=7, width=7, units=c("in"))

p0 <- 	ggplot2::ggplot(df_to_plot2, aes_string(x = "umap_1", y = "umap_2")) +
  geom_point(aes(color=Ind), size=2, alpha=0.4) +
  xlab("Dimension 1") + scale_color_manual(values=scatCols, name="Condition") +
  #scale_shape_manual(name = "Subject", values = c(0, 1, 2, 4,7)) +
  ylab("Dimension 2") + #ggtitle("All Individuals") +  
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme(axis.text = element_text(size=6), axis.title=element_text(size=7), legend.position="none")

for (i in subOrder) {
  text <- i
  pTemp <- ggplot() + annotate("text", x = 4, y = 25, size=8, label = i) + 
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y =element_blank(),
          axis.line = element_blank())
  assign(paste("p",i,"_label",sep=""), pTemp)
}

for (i in condsOrder) {
  text <- i
  pTemp <- ggplot() + annotate("text", x = 4, y = 25, size=8, label = i) + 
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y =element_blank(),
          axis.line = element_blank())
  assign(paste("p",i,"_label",sep=""), pTemp)
}

i = 0
for (j in condsOrder) {
  i = i + 1
  #scatColsTemp <- scatCols
  #scatColsTemp[i] <- scatCols[i]
  cellCount <- nrow(df_to_plot2[which(df_to_plot2$Condition == j),])
  df_to_plot4 <- df_to_plot2[which(df_to_plot2[,"Condition"] == j),]
  titleName <- gsub("D24_","",j)
  #tmpDF2 <- df_to_plot2[which(df_to_plot2[,"cond2"] != j),]
  #df_to_plot4  <- rbind(tmpDF2, tmpDF1)
  p <- ggplot(df_to_plot4, aes_string(x = "umap_1", y = "umap_2")) +
    geom_point(aes(color=Ind), size=2, alpha=0.6) +
    #geom_point(data = subset(df_to_plot4, cond2 == j), aes(color=cond2), size=2, alpha=0.6) +
    xlab("Dimension 1") + scale_color_manual(values=scatCols, name="Condition") +
    #scale_shape_manual(name = "Subject", values = c(0, 1, 2, 4)) +
    ylab("Dimension 2") + ggtitle(titleName) + 
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text = element_text(size=6), axis.title=element_text(size=7), legend.position="none") 
  print(p)
  assign(paste("p",j,sep=""), p)
  print(paste("p",j,sep=""))
}


i = 0
for (j in subOrder) {
  i = i + 1
  #scatColsTemp <- scatCols
  #scatColsTemp[i] <- scatCols[i]
  cellCount <- nrow(df_to_plot2[which(df_to_plot2$Ind == j),])
  df_to_plot4 <- df_to_plot2[which(df_to_plot2[,"Ind"] == j),]
  #tmpDF2 <- df_to_plot2[which(df_to_plot2[,"cond2"] != j),]
  #df_to_plot4  <- rbind(tmpDF2, tmpDF1)
  p <- ggplot(df_to_plot4, aes_string(x = "umap_1", y = "umap_2")) +
    geom_point(aes(color=Ind), size=2, alpha=0.6) +
    #geom_point(data = subset(df_to_plot4, cond2 == j), aes(color=cond2), size=2, alpha=0.6) +
    xlab("Dimension 1") + scale_color_manual(values=scatCols[i], name="Condition") +
    #scale_shape_manual(name = "Subject", values = shapes[i]) +
    ylab("Dimension 2") + ggtitle(j) + 
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text = element_text(size=6), axis.title=element_text(size=7), legend.position="none") 
  print(p)
  assign(paste("p",j,sep=""), p)
  print(paste("p",j,sep=""))
}

i = 0
for (j in condsOrder) {
  l = 0
  i = i + 1
  for (k in subOrder) {
    l = l + 1
    scatColsTemp <- greyCols
    scatColsTemp[i] <- scatCols[i]
    cellCount <- nrow(df_to_plot2[which(df_to_plot2$Condition == j),])
    df_to_plot4 <- df_to_plot2[which(df_to_plot2[,"Condition"] == j & df_to_plot2[,"Ind"] == k),]
    
    #df_to_plot4 <- rbind(tmpDF2, tmpDF1)
    p <- ggplot(df_to_plot4, aes_string(x = "umap_1", y = "umap_2")) +
      geom_point(aes(color=Ind), size=2, alpha=0.6) +
      #geom_point(data = subset(df_to_plot4, cond2 == j), aes(color=cond2), size=2, alpha=0.6) +
      #scale_shape_manual(name = "Subject", values = shapes[l]) +
      xlab("Dimension 1") + scale_color_manual(values=scatCols[l], name="Condition") +
      ylab("Dimension 2") + #ggtitle(paste(condsOrder[i], ", cell count = ", cellCount, sep="" )) + 
      theme(plot.title = element_text(hjust = 0.5)) +
      theme(axis.text = element_text(size=6), axis.title=element_text(size=7), legend.position="none") 
    print(p)
    assign(paste("p",j,k,sep="_"), p)
    print(paste("p",j,k,sep="_"))
  }
}


plot_grid(p0,pSub1,pSub2,pSub3,pSub6,
          pD24_UT, p_D24_UT_Sub1, p_D24_UT_Sub2, p_D24_UT_Sub3, p_D24_UT_Sub6,
          pD24_IFN, p_D24_IFN_Sub1, p_D24_IFN_Sub2, p_D24_IFN_Sub3, p_D24_IFN_Sub6,
          pD24_IL17, p_D24_IL17_Sub1, p_D24_IL17_Sub2, p_D24_IL17_Sub3, p_D24_IL17_Sub6,
          pD24_IL13, p_D24_IL13_Sub1, p_D24_IL13_Sub2, p_D24_IL13_Sub3, p_D24_IL13_Sub6,
          pD24_CO, p_D24_CO_Sub1, p_D24_CO_Sub2, p_D24_CO_Sub3, p_D24_CO_Sub6,
          ncol = 5, nrow = 6)

ggsave("QC/UMAP.Subject_Treatment_groupings.png", height=15, width=15, units=c("in"))



#### MAST DE TESTS ####


geneNames = rownames(as.matrix(soTreatmentBD2@assays$RNA@data))
names(geneNames) = geneNames
geneNames = as.data.frame(geneNames)

soTreatmentBD2@meta.data$IL13 = "No"
soTreatmentBD2@meta.data$IL13[soTreatmentBD2@meta.data$Treatment == "IL13"] = "IL13"
soTreatmentBD2@meta.data$IL13[soTreatmentBD2@meta.data$Treatment == "CO"] = "IL13"

soTreatmentBD2@meta.data$IL17 = "No"
soTreatmentBD2@meta.data$IL17[soTreatmentBD2@meta.data$Treatment == "IL17"] = "IL17"

soTreatmentBD2@meta.data$IFN = "No"
soTreatmentBD2@meta.data$IFN[soTreatmentBD2@meta.data$Treatment == "IFN"] = "IFN"
soTreatmentBD2@meta.data$IFN[soTreatmentBD2@meta.data$Treatment == "CO"] = "IFN"

dataMat = as.matrix(soTreatmentBD2@assays$RNA@counts)

dataMat = (dataMat/colSums((dataMat)))*1e6
dataMat = log2(dataMat + 1)

cellToTestList = levels(soTreatmentBD2@meta.data$cellType)[c(1,2,3,4,5,6)]

cellList = rownames(soTreatmentBD2@meta.data[(soTreatmentBD2@meta.data$cellType %in% cellToTestList),])

conds = soTreatmentBD2@meta.data[cellList,]
conds$cellType = factor(conds$cellType,levels=cellToTestList)

dataMat = dataMat[,cellList]

scaRaw <- FromMatrix(exprsArray =  dataMat, 
                     cData = conds, 
                     fData = geneNames)

sca = scaRaw

scaSample <- sca[sample(which(freq(sca)>.1), 20),]
flat <- as(scaSample, 'data.table')
ggplot(flat, aes(x=value))+geom_density() +facet_wrap(~geneNames, scale='free_y')

thres <- thresholdSCRNACountMatrix(assay(sca), cutbins=c(.01,.1,.2,.5,1,2,3,4,5,6,7,10), min_per_bin = 30)
par(mfrow=c(5,4))
plot(thres)

#assays(sca) <- list(thresh=thres$counts_threshold, tpm=assay(sca))

freq_expressed <- 0.01
FCTHRESHOLD <- log2(1.25)

freqSca = freq(sca)

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
dir.create("DE/DEfiles/mastAllTerms")
dir.create("DE/DEfiles/mastAllTerms/All")
dir.create("DE/DEfiles/mastAllTerms/DEonly")

sigcounts <- as.data.frame(setNames(replicate(2,numeric(0), simplify = F),letters[0:2]))

condsNames = colnames(zlmMod@coefD)
#condsNames = condsNames[c(2:12,25:56)]
condsNames = condsNames[c(26:37)]
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
  
  #fdrPass = nrow(fcHurdleSig)
  #rawPass = nrow(fcHurdle[(fcHurdle$`Pr(>Chisq)` < 0.1),])
  #sigcounts <- rbind(sigcounts,c(fdrPass,rawPass))
  
  print(i)
  print(c(fdrPass,rawPass))
  print(head(fcHurdleSig))
  
  write.table(fcHurdle,file = paste("DE/DEfiles/mastAllTerms/All/",i,".all.full.txt",sep=""),row.names = F,col.names = T,sep = "\t",quote=F)
  write.table(fcHurdleSig,file = paste("DE/DEfiles/mastAllTerms/DEonly/",i,".DE.full.txt",sep=""),row.names = F,col.names = T,sep = "\t",quote=F)
}


#### MAST as conditions ####
#### This is for final celltype-specific DE tables ####

condsList = colnames(zlmCond@coefD)
summaryCondTable <- summary(zlmCond)


summaryCondTableDt <- summaryCondTable$datatable
fcHurdle <- merge(summaryTableDt[contrast==i & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                  summaryTableDt[contrast==i & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid')
fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
setorder(fcHurdle, fdr)


baseLevel = levels(zlmCond@sca@colData$celType_Treatment)[1]
baseLevel = gsub("Proliferating_Basal","ProliferatingBasal",baseLevel)
allLevels = levels(zlmCond@sca@colData$celType_Treatment)
allLevels = gsub("Proliferating_Basal","ProliferatingBasal",allLevels)
denominators = allLevels[grep("UT",allLevels)]
numerators = allLevels[-(grep("UT",allLevels))]
numeratorsCell = unlist(strsplit( numerators, "_"))[seq(1,length(unlist(strsplit( numerators, "_"))),2)]

dir.create("DE/DEfiles/mastAllPairs")
sigCountsPairs = c()
sigCountsNames = c()
pairWiseTests = c()
for (denom in denominators){
  print(denom)
  cellType = gsub("_UT","",denom)
  
  num = numerators[(numeratorsCell == cellType)]
  num = gsub("ProliferatingBasal","Proliferating_Basal",num)
  denom = gsub("ProliferatingBasal","Proliferating_Basal",denom)
  if (denom == baseLevel){
    for (testNum in num){
      print(testNum)
      waldRes = waldTest(zlmCond, Hypothesis(paste("celType_Treatment",testNum,sep="")))
      waldRes2 = waldRes[,,'Pr(>Chisq)']
      hurdleVec = waldRes2[,"hurdle"]
      hurdleVec = as.data.frame(hurdleVec)
      hurdleVec$primerid = rownames(hurdleVec)
      colnames(hurdleVec) = c("hurdle","primerid")
      logfcTable = summaryCondTableDt[(contrast==paste("celType_Treatment",testNum,sep="") & component=='logFC'),]
      fullTable = merge(logfcTable,hurdleVec, by= "primerid")
      fullTable$hurdle.padj = p.adjust(fullTable$hurdle,method = "fdr")
      fullTable = fullTable[order(fullTable$hurdle.padj),]
      reducTable = as.data.frame(fullTable[,c(1,6,8,9)])
      reducTable$FC = 2^reducTable$coef
      reducTable = reducTable[,c(1,2,5,3,4)]
      colnames(reducTable) = c("Gene","log2FC","FC","hurdle_p","hurdle_fdr")
      pairWiseTests[[paste(testNum,".vs.",denom,sep="")]] = reducTable
      write.table(reducTable,paste("DE/DEfiles/mastAllPairs/comp_",testNum,".vs.",denom,".txt",sep=""),quote = F,sep = "\t",row.names = F,col.names = T)
      
      numDE = nrow(reducTable[reducTable$hurdle_fdr <= 0.1,])
      sigCountsPairs = c(sigCountsPairs,numDE)
      sigCountsNames = c(sigCountsNames,paste(testNum,".vs.",denom,sep=""))
      
    }
  } else {
    logfcTableDenom = summaryCondTableDt[(contrast==paste("celType_Treatment",denom,sep="") & component=='logFC'),]
    for (testNum in num){
      print(testNum)
      waldRes = waldTest(zlmCond, Hypothesis(paste("celType_Treatment",testNum,"-",
                                                   "celType_Treatment",denom,sep="")))
      waldRes2 = waldRes[,,'Pr(>Chisq)']
      hurdleVec = waldRes2[,"hurdle"]
      hurdleVec = as.data.frame(hurdleVec)
      hurdleVec$primerid = rownames(hurdleVec)
      colnames(hurdleVec) = c("hurdle","primerid")
      logfcTableNum = summaryCondTableDt[(contrast==paste("celType_Treatment",testNum,sep="") & component=='logFC'),]
      fcMerge = as.data.frame(merge(logfcTableNum, logfcTableDenom,by="primerid",suffixes=paste(".",c(testNum,denom),sep="")))
      fcMerge$coef = fcMerge[,paste("coef.",testNum,sep="")] - fcMerge[,paste("coef.",denom,sep="")]
      fullTable = merge(fcMerge,hurdleVec, by= "primerid")
      fullTable$hurdle.padj = p.adjust(fullTable$hurdle,method = "fdr")
      fullTable = fullTable[order(fullTable$hurdle.padj),]
      lenFullTable = ncol(fullTable)
      reducTable = as.data.frame(fullTable[,c(1,seq(lenFullTable-2,lenFullTable))])
      reducTable$FC = 2^reducTable$coef
      reducTable = reducTable[,c(1,2,5,3,4)]
      colnames(reducTable) = c("Gene","log2FC","FC","hurdle_p","hurdle_fdr")
      pairWiseTests[[paste(testNum,".vs.",denom,sep="")]] = reducTable
      
      write.table(reducTable,paste("DE/DEfiles/mastAllPairs/comp_",testNum,".vs.",denom,".txt",sep=""),quote = F,sep = "\t",row.names = F,col.names = T)
      
      numDE = nrow(reducTable[reducTable$hurdle_fdr <= 0.1,])
      sigCountsPairs = c(sigCountsPairs,numDE)
      sigCountsNames = c(sigCountsNames,paste(testNum,".vs.",denom,sep=""))
    }
  }
}



names(sigCountsPairs) = sigCountsNames
sigCountsPairs = as.data.frame(sigCountsPairs)

IL13conds = pairWiseTests[["Secretory_IL13.vs.Secretory_UT"]]
IL13model = mastTables[["IL13IL13+cellTypeSecretory:IL13IL13"]]

mergeIL13 = merge(IL13conds,IL13model,by.x="Gene",by.y="primerid")

ggplot(mergeIL13,aes(x=log2FC,y=coef.IL13IL13)) + geom_point(alpha=0.2)


#### all pairs sigcounts ####

sigcountsPairs <- as.data.frame(setNames(replicate(2,numeric(0), simplify = F),letters[0:2]))
nameVec = c()
testNames = names(pairWiseTests)
testNames = testNames[order(testNames)]
for (i in testNames) {
  fcHurdleReduc = pairWiseTests[[i]]
  fcHurdleReducSig = fcHurdleReduc[(fcHurdleReduc$hurdle_fdr <= 0.1),]
  fcHurdleReducSig = fcHurdleReducSig[!(is.na(fcHurdleReducSig$Gene)),]
  fdrPass = nrow(fcHurdleReducSig)
  fcHurdleReducSig = fcHurdleReduc[(fcHurdleReduc$hurdle_p <= 0.1),]
  fcHurdleReducSig = fcHurdleReducSig[!(is.na(fcHurdleReducSig$Gene)),]
  rawPass = nrow(fcHurdleReducSig)
  sigcountsPairs <- rbind(sigcountsPairs,c(fdrPass,rawPass))
  nameVec = c(nameVec,i)
}
rownames(sigcountsPairs) = nameVec

sigCountsPairsReOrder = sigcountsPairs[c(20,4,12,24,16,8,19,3,11,23,15,7,18,2,10,22,14,6,17,1,9,21,13,5),]
sigCountsPairsReOrder = as.data.frame(sigCountsPairsReOrder)
colnames(sigCountsPairsReOrder) = c("FDR < 0.1","pVal < 0.1")
write.table(sigCountsPairsReOrder, file="DE/DEfiles/mastAllPairs/sigcounts.txt",
            sep="\t",row.names = T,col.names = F,quote=F)

#### all term sigcounts ####

sigcountsMast <- as.data.frame(setNames(replicate(2,numeric(0), simplify = F),letters[0:2]))
nameVec = c()
testNames = names(mastTables)
testNames = testNames[order(testNames)]
for (i in testNames) {
  fcHurdleReduc = mastTables[[i]]
  if (ncol(fcHurdleReduc) == 11) {
    fdrPass = nrow(fcHurdleReduc[(fcHurdleReduc$hurdleFDR <= 0.1),])
    rawPass = nrow(fcHurdleReduc[(fcHurdleReduc$hurdle <= 0.1),])
    sigcountsMast <- rbind(sigcountsMast,c(fdrPass,rawPass))
    nameVec = c(nameVec,i)
  }
  else {
    fdrPass = nrow(fcHurdleReduc[(fcHurdleReduc$fdr <= 0.1),])
    rawPass = nrow(fcHurdleReduc[(fcHurdleReduc$fdr <= 0.1),])
    sigcountsMast <- rbind(sigcountsMast,c(fdrPass,rawPass))
    nameVec = c(nameVec,i)
  }
}
rownames(sigcountsMast) = nameVec


#### celltype markers counts ####

sigcountsMastCellTypes <- as.data.frame(setNames(replicate(2,numeric(0), simplify = F),letters[0:2]))
nameVec = c()
testNames = names(mastTables)[1:7]
testNames = c("cellTypeProliferating_Basal",testNames)

for (i in testNames) {
  fcHurdleReduc = mastTables[[i]]
  fdrPass = nrow(fcHurdleReduc[(fcHurdleReduc$fdr <= 0.1),])
  rawPass = nrow(fcHurdleReduc[(fcHurdleReduc$`Pr(>Chisq)` <= 0.1),])
  print(c(i,fdrPass,rawPass))
  sigcountsMastCellTypes <- rbind(sigcountsMastCellTypes,c(fdrPass,rawPass))
  nameVec = c(nameVec,i)
}

rownames(sigcountsMastCellTypes) = nameVec
colnames(sigcountsMastCellTypes) = c("fdr < 0.1","pval < 0.1")

#sigcountsMastCellTypes = sigcountsMastCellTypes[c(2,6,5,1,3,4,7)]

write.table(sigcountsMastCellTypes, file="DE/DEfiles/mastAllTerms/sigcounts.CellTypes.txt",
            sep="\t",row.names = T,col.names = NA,quote=F)


#### cell type markers ####


cellTypeTests = names(mastTablesDE)[1:5]
mastCellTypeMarkers = c()

cellTypeTests = c("cellTypeProliferating_Basal",cellTypeTests)

for (i in cellTypeTests){
  dataTab = mastTablesDE[[i]]
  dataTab = dataTab[dataTab$coef > 0,]
  rownames(dataTab) = dataTab$primerid
  markerList = dataTab$primerid
  remaining = cellTypeTests[!(cellTypeTests %in% i)]
  for (j in remaining){
    dataNext = mastTablesDE[[j]]
    rownames(dataNext) = dataNext$primerid
    for (gene in markerList){
      fc = dataTab[(dataTab$primerid == gene),coef]
      fcOther = dataNext[(dataNext$primerid == gene),coef]
      if (length(fcOther) > 0){
        if( fcOther >= fc*.5) {
          markerList = markerList[!(markerList %in% gene)]
        }
      }
    }
  }
  dataTabNew = dataTab[(dataTab$primerid %in% markerList),]
  mastCellTypeMarkers[[i]] = dataTabNew
  if ("MCIDAS" %in% dataTabNew$primerid){
    print(dataTabNew[(dataTabNew$primerid == "MCIDAS"),])
  }
  print(i)
  print(nrow(dataTabNew))
}


