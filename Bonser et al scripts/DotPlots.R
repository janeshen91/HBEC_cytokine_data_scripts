library(ggplot2)
library(dplyr)
library(Seurat)

##read in the RDS files that have the scRNAseq data
COV2<-readRDS("../../iCloud Drive (Archive)//Documents/scRNAseq/Walter's scripts and data/COV2_cellType_and_cellTypeXtreatment_results.rds")
#read in the names of the genes we want to plot
cellTypeSpefici<-read.csv(file="/Users/janeshen/Documents/covid19cellatlas/ITGBspecificGenes.txt",header=F)


#get only the genes we're interested in
rows<-COV2$Gene %in% cellTypeSpefici$V1
COV2genes<-unique(COV2$Gene)
CellTypeSpecific<-COV2[rows,]

#get only the cell types we're interested in
CellTypeSpecific3<-subset(CellTypeSpecific,CellTypeSpecific$cellType %in% c("Basal","Secretory","Ciliated"))
CellTypeSpecific3<-subset(CellTypeSpecific3,CellTypeSpecific3$cellType_Treatment %in% c("Basal_UT","Ciliated_UT","Secretory_UT"))


##re-organize the genes by average gene expression
CellTypeSpecific4<-CellTypeSpecific3 %>% group_by(cellType) %>% 
  mutate(mx = max(cellType.avg_logFC)) %>% 
  arrange(desc(mx), desc(cellType.avg_logFC))

orderGene<-as.character(CellTypeSpecific4$Gene)


CellTypeSpecific5<-aggregate(.~Gene, data=CellTypeSpecific4[,c("Gene","aveExp")],FUN=max)
colnames(CellTypeSpecific5)<-c("Gene","maxAveExp")
CellTypeSpecific6<-merge(CellTypeSpecific4,CellTypeSpecific5,by="Gene")
CellTypeSpecific6$cellType<-as.character(CellTypeSpecific6$cellType)
CellTypeSpecific6$CellTypeMaxExpres<-ifelse(CellTypeSpecific6$maxAveExp==CellTypeSpecific6$aveExp,CellTypeSpecific6$cellType,"NA")

####subset out the basal, ciliated and secretory genes and then organize the order of the genes plotted by the average logfc in these celltypes
basal<-subset(CellTypeSpecific6,CellTypeSpecific6$CellTypeMaxExpres=="Basal")
dim(basal)
basalOrder<-basal[order(basal$cellType.avg_logFC,decreasing = TRUE),]
orderGeneBasal<-as.character(basalOrder$Gene)
###
ciliated<-subset(CellTypeSpecific6,CellTypeSpecific6$CellTypeMaxExpres=="Ciliated")
CiliatedOrder<-ciliated[order(ciliated$cellType.avg_logFC,decreasing = TRUE),]
orderGeneCiliated<-as.character(CiliatedOrder$Gene)

###
Secretory<-subset(CellTypeSpecific6,CellTypeSpecific6$CellTypeMaxExpres=="Secretory")
SecretoryOrder<-Secretory[order(Secretory$cellType.avg_logFC,decreasing = TRUE),]
orderGeneSecretory<-as.character(SecretoryOrder$Gene)

##overall
orderOverall<-c(orderGeneBasal,orderGeneSecretory,orderGeneCiliated)



##basal only

Basal<-subset(CellTypeSpecific3,CellTypeSpecific3$cellType=="Basal")
head(Basal)
Basal<-Basal[order(Basal$cellType.avg_logFC,decreasing = TRUE),] 
orderGeneBasal<-as.character(Basal$Gene)


Secretory<-subset(CellTypeSpecific3,CellTypeSpecific3$cellType=="Secretory")
Secretory<-Secretory[order(Secretory$cellType.avg_logFC,decreasing = TRUE),] 
head(Secretory)
orderGeneSecretory<-as.character(Secretory$Gene)

Ciliated<-subset(CellTypeSpecific3,CellTypeSpecific3$cellType=="Ciliated")
Ciliated<-Ciliated[order(Ciliated$cellType.avg_logFC,decreasing = TRUE),] 
head(Ciliated)
orderGeneCiliated<-as.character(Ciliated$Gene)


overallOrder<-c(orderGeneBasal[c(1:6)],orderGeneSecretory[c(1:9)],orderGeneCiliated[c(1:16)])

####plot
p = ggplot(CellTypeSpecific4, aes(y = Gene, x = cellType)) +
  geom_point(aes(size = pctExp, colour=cellType.avg_logFC),) +  scale_size_area() +
  geom_point(data=CellTypeSpecific3[(CellTypeSpecific3$cellType.p_val_adj <= 0.1),],
             size=6,pch=21, fill=NA, colour="black", stroke=1.25)+
  scale_colour_gradient2(low = "blue", mid = "grey", high = "red", 
                         na.value = "grey") +
  theme_bw() +theme(axis.text.x = element_text(angle = 60, hjust = 1,color="black",size=12),
                    axis.title=element_text(size=16), 
                    axis.text.y = element_text(color="black",size=12,face="italic"),
                    legend.text=element_text(size=10),legend.title=element_text(size=14),
                    panel.grid = element_blank(),plot.title = element_text(hjust = 0.5))+labs(size="Percent\nExpressed", colour="Log2fc") + ylab("Gene Name") +
  guides(color = guide_colorbar(order = 1),size = guide_legend(order=2))+
  scale_x_discrete(limits=orderCellTypes)+scale_y_discrete(limits=rev(orderOverall))+
  geom_hline(aes(yintercept=15.5),color="gray") +
  geom_hline(aes(yintercept=23.5),color="gray")
p
ggsave(paste("/Users/janeshen/Box/Erle group presentations/Jane Shen/COV_2_figures/OriginalData.SCV2.dotPlot.pdf",sep=""),height = 8,width = 8,units = "in",dpi = 300)

