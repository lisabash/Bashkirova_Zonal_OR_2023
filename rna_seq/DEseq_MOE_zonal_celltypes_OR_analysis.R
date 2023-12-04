library(GenomicFeatures)
library(DESeq2)
library(Rsamtools) #needed for BamFileList
library(rtracklayer) #needed to import GTF as GRanges object
library(GenomicAlignments) #used for summarizeOverlaps
library( "RColorBrewer" ) #for making heatmaps
library(gridExtra)
library(ggrepel) # for labeling points on plot with gene name
library("pheatmap") #used for heatmaps
library(dplyr)
library(tibble)
library(tidyverse)

# zonal color palette
zonal_fine_colors_fish_classIfirst<-c("#990000","#FF0000","#FFCC00","#33CC33","#0099FF","#0000FF")
zonal_fine_colors_dark_fish_classIfirst<-c("#660000","#990000","#999900","#003300","#0066CC","#000066")

#load in annotation
genes <-import("/media/storageA/kevin/annotation/genes+soria+pcdh.gtf")
genes_txdb <- makeTxDbFromGRanges(genes) # creates database
exonsByGene <- exonsBy(genes_txdb, by="gene") #GRangesList object, run exonsByGene$Olfr1507 to see the annotated exons
gene_width <- data.frame(sum(width(exonsByGene)))

#read in my data
sampleTable <- read.delim("/media/storageA/lisa/2020_zonal_paper/analysis/RNAseq_zonal_dev/200424/DEseq2_sampleTable_zonal_dev_200423.txt",header =T)
sampleTable

bamfiles <- filenames <- file.path(sampleTable$data) #creates character vector with paths to bam files
bamfiles <- BamFileList(bamfiles, yieldSize=1000000)


se <- summarizeOverlaps(features=exonsByGene, reads=bamfiles,
                        mode="Union",
                        singleEnd=FALSE,
                        ignore.strand=FALSE) #for truseq data include preprocess.reads=invertStrand


colData(se) <-DataFrame(sampleTable) #is this additional metadata necessary?
colnames(se)<-sampleTable$library

#one factor design formula is simple: ~ condition
dds <- DESeqDataSet(se, design = ~ condition)
dds$condition <- relevel(dds$condition, ref="Z1Mash") #put control condition first, otherwise it will do put the levels alphabetically
### remove genes with no mapped reads in any condition
dds <- dds[ rowSums(counts(dds)) > 1, ]

#differential expression analysis
dds <- DESeq(dds) #estimates the size factors (to control for library size), dispersion estimation for each gene, and fitting to a linear model. Returns a DESeqDataSet.
Olfr <- grep( "^Olfr", rownames(dds), value = TRUE)
nonOlfr <- grep( "^Olfr", rownames(dds), invert= TRUE, value = TRUE)

#############################
# get FPKM AND TPM from DDS #
#############################
### Calculate FPKM
fpkm_scaled <- as.data.frame(fpkm(dds, robust = TRUE)) # robus normalizes on DEseq "scaling factor"
fpkm_scaled2<-fpkm_scaled %>% rownames_to_column(var="GeneID") %>% gather(condition,fpkm,-GeneID) %>% separate(condition,sep="_",into=c("V1","V2","rep")) %>% 
  unite("condition", V1, V2)
fpkm_scaled.mean <- fpkm_scaled2 %>% group_by(condition,GeneID) %>% summarize(mean_fpkm=mean(fpkm))
fpkm_scaled.mean.table <- spread(fpkm_scaled.mean,condition,mean_fpkm)

### Calculate unnormalized FPKM to generate TPM
fpkm_unscaled <- as.data.frame(fpkm(dds, robust = FALSE)) %>% rownames_to_column(var="GeneID")
fpkm_unscaled <- fpkm_unscaled %>% gather(condition,fpkm,-GeneID) %>% separate(condition,sep="_",into=c("V1","V2","rep")) %>% 
  unite("condition", V1, V2)
fpkm_unscaled.mean<- fpkm_unscaled  %>% group_by(condition,GeneID) %>% summarize(mean_fpkm=mean(fpkm))
fpkm.mean.table <- spread(fpkm_unscaled.mean,condition,mean_fpkm)

### Calculate TPM
tpm <- fpkm_unscaled %>% group_by(condition,rep) %>% summarize(total_fpkm = sum(fpkm)) %>% full_join(fpkm_unscaled) %>% mutate(tpm = (fpkm / total_fpkm) *1000000) %>% dplyr::select(GeneID,condition,rep,tpm)
tpm.mean <- tpm %>% group_by(condition,GeneID) %>% summarize(mean_tpm = mean(tpm))
tpm.sd <- tpm %>% group_by(condition,GeneID) %>% summarize(sd = sd(tpm))
tpm.table <- tpm %>% mutate(lib = paste0(condition,"_",rep)) %>% group_by(lib,GeneID) %>% dplyr::select(GeneID,condition = lib,tpm) %>% spread(condition,tpm)
tpm.mean.table <- spread(tpm.mean,condition,mean_tpm)
tpm.sd.table <- spread(tpm.sd,condition,sd)

###########################################
# LINEPLOTS OF NFI EXPRESSION  (Figure 4B)#
###########################################

nfi.tpm.table<-tpm.mean.table %>% filter(GeneID %in% c("Nfia","Nfib","Nfix")) %>% gather(sample,tpm,-GeneID) %>% separate(sample,sep="_",into=c("zone","celltype"), remove = FALSE) %>% unite("group",zone,GeneID, remove = FALSE)
nfi.tpm.table
nfi.sd.table<-tpm.sd.table %>% filter(GeneID %in% c("Nfia","Nfib","Nfix")) %>% gather(sample,sd,-GeneID) %>% separate(sample,sep="_",into=c("zone","celltype"), remove = FALSE) %>% unite("group",zone,GeneID, remove = FALSE)
nfi.sd.table
nfi.mean.sd.table<-left_join(nfi.tpm.table,nfi.sd.table)

p.nfi.mean.split<-ggplot(nfi.mean.sd.table,mapping=aes(celltype,tpm,group=group))+geom_line(mapping=aes(color=zone), size=1)+
  geom_point(aes(col=zone, shape=zone), size=4)+ggtitle("NFI family expression")+
  theme(axis.text.y = element_text(colour= "black", size = 12),axis.text.x = element_text(colour= "black", size = 14), axis.title.x = element_blank(),axis.title.y = element_text(color="black", size=14),plot.title = element_text(hjust = 0.5, size=20),strip.text.x = element_text(size = 14, colour = "black"))+
  ylab("TPM")+scale_color_manual(values=c("#EC001D","#23BD7D","#0505B7"))+ 
  facet_grid(cols = vars(GeneID))
p.nfi.mean.split
#ggsave(file="lineplot.NFI_TPM.split.mean_wider.pdf",p.nfi.mean.split, width = 10, height = 4, units = "in" ,useDingbats=FALSE)

#############################
# TMP OLFR PLOTS (Figure 1E)#
#############################

# load zonal annotation
zonal_anno_fine_tb<-read_tsv("/media/storageA/lisa/2020_zonal_paper/annotation/zonal_anno_TanLocus_soriaNames/zonal_anno_2020_fine_TanLocus_SoriaName.bed", col_names = c("chr","start","stop","name","qscore","strand","zone"))
zonal_anno_fine_tb<-zonal_anno_fine_tb %>% select(name,zone)
tpm.mean.table.zone<-right_join(tpm.mean.table,zonal_anno_fine_tb,by=c("GeneID" = "name"))
tpm.mean.table.zone$zone<-factor(tpm.mean.table.zone$zone, levels=c("classI","1","2","3","4","5"))

# INP cells
p.z1.mashngn_classI<-ggplot(tpm.mean.table.zone, mapping = aes(zone,Z1_mashNgn, color=zone, fill = zone)) + ylim(0,3)+
  geom_boxplot(alpha = 0.7,  lwd =1.2 , outlier.shape=NA)+ #geom_boxplot(alpha = 0.1, outlier.shape = NA , lwd =1.2 ) +
  scale_color_manual(values=zonal_fine_colors_dark_fish_classIfirst) + 
  scale_fill_manual(values = zonal_fine_colors_fish_classIfirst) +#+ geom_jitter()
  theme(legend.position="none",axis.text.y = element_text(colour= "black", size = 12),axis.text.x = element_text(colour= "black", size = 14), axis.title.x = element_blank(),axis.title.y = element_text(colour= "black", size = 14),plot.title = element_text(hjust = 0.5, size=20))+
  ggtitle("Zone1 MashNgn")
p.z1.mashngn_classI

p.z2.mashngn_classI<-ggplot(tpm.mean.table.zone, mapping = aes(zone,Z2_mashNgn, color=zone, fill = zone)) + ylim(0,3)+
  geom_boxplot(alpha = 0.7,  lwd =1.2 , outlier.shape=NA)+ #geom_boxplot(alpha = 0.1, outlier.shape = NA , lwd =1.2 ) +
  scale_color_manual(values=zonal_fine_colors_dark_fish_classIfirst) + 
  scale_fill_manual(values = zonal_fine_colors_fish_classIfirst) +#+ geom_jitter()
  theme(legend.position="none",axis.text.y = element_text(colour= "black", size = 12),axis.text.x = element_text(colour= "black", size = 14), axis.title.x = element_blank(),axis.title.y = element_text(colour= "black", size = 14),plot.title = element_text(hjust = 0.5, size=20))+
  ggtitle("Zone2 MashNgn")
p.z2.mashngn_classI

p.z5.mashngn_classI<-ggplot(tpm.mean.table.zone, mapping = aes(zone,Z5_mashNgn, color=zone, fill = zone)) + ylim(0,3)+
  geom_boxplot(alpha = 0.7,  lwd =1.2 , outlier.shape=NA)+ #geom_boxplot(alpha = 0.1, outlier.shape = NA , lwd =1.2 ) +
  scale_color_manual(values=zonal_fine_colors_dark_fish_classIfirst) + 
  scale_fill_manual(values = zonal_fine_colors_fish_classIfirst) +#+ geom_jitter()
  theme(legend.position="none",axis.text.y = element_text(colour= "black", size = 12),axis.text.x = element_text(colour= "black", size = 14), axis.title.x = element_blank(),axis.title.y = element_text(colour= "black", size = 14),plot.title = element_text(hjust = 0.5, size=20))+
  ggtitle("Zone5 MashNgn")
p.z5.mashngn_classI

grid.arrange(p.z1.mashngn_classI, p.z2.mashngn_classI, p.z5.mashngn_classI, ncol = 3,nrow=1)
#g<-arrangeGrob(p.z1.mashngn_classI, p.z2.mashngn_classI, p.z5.mashngn_classI, ncol = 3,nrow=1)
#ggsave(file="boxplot.ORs_mashNgn_only_classI_colors2022.pdf",g, width = 12, height = 4, units = "in" )

# iOSN cells
p.z1.ngn_classI<-ggplot(tpm.mean.table.zone, mapping = aes(zone,Z1_NgnBright, color=zone, fill = zone)) + ylim(0,30)+
  geom_boxplot(alpha = 0.7,  lwd =1.2 , outlier.shape=NA)+ #geom_boxplot(alpha = 0.1, outlier.shape = NA , lwd =1.2 ) +
  scale_color_manual(values=zonal_fine_colors_dark_fish_classIfirst) + 
  scale_fill_manual(values = zonal_fine_colors_fish_classIfirst) +#+ geom_jitter()
  theme(legend.position="none",axis.text.y = element_text(colour= "black", size = 12),axis.text.x = element_text(colour= "black", size = 14), axis.title.x = element_blank(),axis.title.y = element_text(colour= "black", size = 14),plot.title = element_text(hjust = 0.5, size=20))+
  ggtitle("Zone1 Ngn")
p.z1.ngn_classI

p.z2.ngn_classI<-ggplot(tpm.mean.table.zone, mapping = aes(zone,Z2_NgnBright, color=zone, fill = zone)) + ylim(0,15)+
  geom_boxplot(alpha = 0.7,  lwd =1.2 , outlier.shape=NA)+ #geom_boxplot(alpha = 0.1, outlier.shape = NA , lwd =1.2 ) +
  scale_color_manual(values=zonal_fine_colors_dark_fish_classIfirst) + 
  scale_fill_manual(values = zonal_fine_colors_fish_classIfirst) +#+ geom_jitter()
  theme(legend.position="none",axis.text.y = element_text(colour= "black", size = 12),axis.text.x = element_text(colour= "black", size = 14), axis.title.x = element_blank(),axis.title.y = element_text(colour= "black", size = 14),plot.title = element_text(hjust = 0.5, size=20))+
  ggtitle("Zone2 Ngn")
p.z2.ngn_classI

p.z5.ngn_classI<-ggplot(tpm.mean.table.zone, mapping = aes(zone,Z5_NgnBright, color=zone, fill = zone)) + ylim(0,30)+
  geom_boxplot(alpha = 0.7,  lwd =1.2 , outlier.shape=NA)+ #geom_boxplot(alpha = 0.1, outlier.shape = NA , lwd =1.2 ) +
  scale_color_manual(values=zonal_fine_colors_dark_fish_classIfirst) + 
  scale_fill_manual(values = zonal_fine_colors_fish_classIfirst) +#+ geom_jitter()
  theme(legend.position="none",axis.text.y = element_text(colour= "black", size = 12),axis.text.x = element_text(colour= "black", size = 14), axis.title.x = element_blank(),axis.title.y = element_text(colour= "black", size = 14),plot.title = element_text(hjust = 0.5, size=20))+
  ggtitle("Zone5 Ngn")
p.z5.ngn_classI

grid.arrange(p.z1.ngn_classI, p.z2.ngn_classI, p.z5.ngn_classI, ncol = 3,nrow=1)
#g<-arrangeGrob(p.z1.ngn_classI, p.z2.ngn_classI, p.z5.ngn_classI, ncol = 3,nrow=1)
#ggsave(file="boxplot.ORs_ngn_only_classI_colors2022_2.pdf",g, width = 12, height = 4, units = "in" )

# mOSN cells
p.z1.ngndim_classI<-ggplot(tpm.mean.table.zone, mapping = aes(zone,Z1_NgnDim, color=zone, fill = zone)) + ylim(0,40)+
  geom_boxplot(alpha = 0.7,  lwd =1.2 , outlier.shape=NA)+ #geom_boxplot(alpha = 0.1, outlier.shape = NA , lwd =1.2 ) +
  scale_color_manual(values=zonal_fine_colors_dark_fish_classIfirst) + 
  scale_fill_manual(values = zonal_fine_colors_fish_classIfirst) +#+ geom_jitter()
  theme(legend.position="none",axis.text.y = element_text(colour= "black", size = 12),axis.text.x = element_text(colour= "black", size = 14), axis.title.x = element_blank(),axis.title.y = element_text(colour= "black", size = 14),plot.title = element_text(hjust = 0.5, size=20))+
  ggtitle("Zone1 NgnDim")
p.z1.ngndim_classI

p.z2.ngndim_classI<-ggplot(tpm.mean.table.zone, mapping = aes(zone,Z2_NgnDim, color=zone, fill = zone)) + ylim(0,40)+
  geom_boxplot(alpha = 0.7,  lwd =1.2 , outlier.shape=NA)+ #geom_boxplot(alpha = 0.1, outlier.shape = NA , lwd =1.2 ) +
  scale_color_manual(values=zonal_fine_colors_dark_fish_classIfirst) + 
  scale_fill_manual(values = zonal_fine_colors_fish_classIfirst) +#+ geom_jitter()
  theme(legend.position="none",axis.text.y = element_text(colour= "black", size = 12),axis.text.x = element_text(colour= "black", size = 14), axis.title.x = element_blank(),axis.title.y = element_text(colour= "black", size = 14),plot.title = element_text(hjust = 0.5, size=20))+
  ggtitle("Zone1 NgnDim")
p.z2.ngndim_classI

p.z5.ngndim_classI<-ggplot(tpm.mean.table.zone, mapping = aes(zone,Z5_NgnDim, color=zone, fill = zone)) + ylim(0,40)+
  geom_boxplot(alpha = 0.7,  lwd =1.2 , outlier.shape=NA)+ #geom_boxplot(alpha = 0.1, outlier.shape = NA , lwd =1.2 ) +
  scale_color_manual(values=zonal_fine_colors_dark_fish_classIfirst) + 
  scale_fill_manual(values = zonal_fine_colors_fish_classIfirst) +#+ geom_jitter()
  theme(legend.position="none",axis.text.y = element_text(colour= "black", size = 12),axis.text.x = element_text(colour= "black", size = 14), axis.title.x = element_blank(),axis.title.y = element_text(colour= "black", size = 14),plot.title = element_text(hjust = 0.5, size=20))+
  ggtitle("Zone1 NgnDim")
p.z5.ngndim_classI

grid.arrange(p.z1.ngndim_classI, p.z2.ngndim_classI, p.z5.ngndim_classI, ncol = 3,nrow=1)
#g<-arrangeGrob(p.z1.ngndim_classI, p.z2.ngndim_classI, p.z5.ngndim_classI, ncol = 3,nrow=1)
#ggsave(file="boxplot.ORs_ngndim_only_classI_colors2022_ylim40.pdf",g, width = 12, height = 4, units = "in" )

####################################
# DIFFERENTIAL EXPRESSION ANALYSIS #
####################################

#  FIND DIFFERENTIALLY EXPRESSED TRANSCRIPTION FACTORS

# TFs from Gene Ontology database annotation “DNA binding transcription factor activity”
TFs<-read.table("/media/storageA/lisa/2020_zonal_paper/annotation/GO/GO_DNA-binding_transcription_factor_activity_MF.txt")
TF.vect<-TFs$V1

########################################################
# padj<0.05, 2 fold change between zone1 and zone5 cells
########################################################

# MAKE TABLE OF DIFFERENTIALLY EXPRESSED TFS (Supplementary File 1)

# Find diff expressed TFs
# GBC cells
res05 <- results(dds, contrast=c("condition","Z5Mash","Z1Mash"),alpha=0.05) #more stringent significance cutoff, by default alpha=0.1
res.mash5v1.alpha05<-rownames_to_column(as.data.frame(res05))
res.mash5v1.alpha05.up<-filter(res.mash5v1.alpha05,rowname %in% TF.vect) %>% filter(padj < 0.05) %>% filter(log2FoldChange > 1)
res.mash5v1.alpha05.down<-filter(res.mash5v1.alpha05,rowname %in% TF.vect) %>% filter(padj < 0.05) %>% filter(log2FoldChange < -1) 
res.mash5v1.alpha05.UPandDOWN<-filter(res.mash5v1.alpha05,rowname %in% TF.vect) %>% filter(padj < 0.05) %>% filter(log2FoldChange < -1 |log2FoldChange > 1) %>% mutate('celltype'="GBC")
mash.5v1.names<-unique(c(res.mash5v1.alpha05.up$rowname,res.mash5v1.alpha05.down$rowname))

# INP cells
res05 <- results(dds, contrast=c("condition","Z5MashNgn","Z1MashNgn"),alpha=0.05) #more stringent significance cutoff, by default alpha=0.1
res.mashNgn5v1.alpha05<-rownames_to_column(as.data.frame(res05))
res.mashNgn5v1.alpha05.up<-filter(res.mashNgn5v1.alpha05,rowname %in% TF.vect) %>% filter(padj < 0.05) %>% filter(log2FoldChange > 1) 
res.mashNgn5v1.alpha05.down<-filter(res.mashNgn5v1.alpha05,rowname %in% TF.vect) %>% filter(padj < 0.05) %>% filter(log2FoldChange < -1) 
res.mashNgn5v1.alpha05.UPandDOWN<-filter(res.mashNgn5v1.alpha05,rowname %in% TF.vect) %>% filter(padj < 0.05) %>% filter(log2FoldChange < -1 |log2FoldChange > 1) %>% mutate('celltype'="INP")
mashNgn.5v1.names<-unique(c(res.mashNgn5v1.alpha05.up$rowname,res.mashNgn5v1.alpha05.down$rowname))

# iOSN cells
res05 <- results(dds, contrast=c("condition","Z5NgnBright","Z1NgnBright"),alpha=0.05) #more stringent significance cutoff, by default alpha=0.1
res.NgnBright5v1.alpha05<-rownames_to_column(as.data.frame(res05))
res.NgnBright5v1.alpha05.up<-filter(res.NgnBright5v1.alpha05,rowname %in% TF.vect) %>% filter(padj < 0.05) %>% filter(log2FoldChange > 1) 
res.NgnBright5v1.alpha05.down<-filter(res.NgnBright5v1.alpha05,rowname %in% TF.vect) %>% filter(padj < 0.05) %>% filter(log2FoldChange < -1) 
res.NgnBright5v1.alpha05.UPandDOWN<-filter(res.NgnBright5v1.alpha05,rowname %in% TF.vect) %>% filter(padj < 0.05) %>% filter(log2FoldChange < -1 |log2FoldChange > 1) %>% mutate('celltype'="iOSN")
NgnBright.5v1.names<-unique(c(res.NgnBright5v1.alpha05.up$rowname,res.NgnBright5v1.alpha05.down$rowname))

# mOSN cells
res05 <- results(dds, contrast=c("condition","Z5NgnDim","Z1NgnDim"),alpha=0.05) #more stringent significance cutoff, by default alpha=0.1
summary(res05)
res.NgnDim5v1.alpha05<-rownames_to_column(as.data.frame(res05))
res.NgnDim5v1.alpha05.up<-filter(res.NgnDim5v1.alpha05,rowname %in% TF.vect) %>% filter(padj < 0.05) %>% filter(log2FoldChange > 1) 
res.NgnDim5v1.alpha05.down<-filter(res.NgnDim5v1.alpha05,rowname %in% TF.vect) %>% filter(padj < 0.05) %>% filter(log2FoldChange < -1) 
res.NgnDim5v1.alpha05.UPandDOWN<-filter(res.NgnDim5v1.alpha05,rowname %in% TF.vect) %>% filter(padj < 0.05) %>% filter(log2FoldChange < -1 |log2FoldChange > 1) %>% mutate('celltype'="OSN")
NgnDim.5v1.names<-unique(c(res.NgnDim5v1.alpha05.up$rowname,res.NgnDim5v1.alpha05.down$rowname))

##################################
#make table (Supplementary File 1)
##################################
tpm.mash10.5v1_table<-tpm.mean.table %>% filter(GeneID %in% mash.5v1.names) %>% select(GeneID,Z1_mash,Z2_mash,Z5_mash) %>% arrange(desc(Z5_mash))
log2mash_table<-filter(res.mash5v1.alpha05,rowname %in% mash.5v1.names) %>% select(rowname,log2FoldChange,padj) 
tpm.mash10.5v1_table<-left_join(tpm.mash10.5v1_table,log2mash_table, by=c("GeneID" = "rowname" )) %>% 
  mutate("cell_type" = "GBC") %>% rename("zone1_TPM"=Z1_mash,"zone2_TPM"=Z2_mash,"zone5_TPM"=Z5_mash)
tpm.mashNgn10.5v1_table<-tpm.mean.table %>% filter(GeneID %in% mashNgn.5v1.names) %>% select(GeneID,Z1_mashNgn,Z2_mashNgn,Z5_mashNgn) %>% arrange(desc(Z5_mashNgn))
log2mashNgn_table<-filter(res.mashNgn5v1.alpha05,rowname %in% mashNgn.5v1.names) %>% select(rowname,log2FoldChange,padj )
tpm.mashNgn10.5v1_table<-left_join(tpm.mashNgn10.5v1_table,log2mashNgn_table, by=c("GeneID" = "rowname" )) %>% 
  mutate("cell_type" = "INP") %>% rename("zone1_TPM"=Z1_mashNgn,"zone2_TPM"=Z2_mashNgn,"zone5_TPM"=Z5_mashNgn)
tpm.NgnBright10.5v1_table<-tpm.mean.table %>% filter(GeneID %in% NgnBright.5v1.names) %>% select(GeneID,Z1_NgnBright,Z2_NgnBright,Z5_NgnBright) %>% arrange(desc(Z5_NgnBright))
log2NgnBright_table<-filter(res.NgnBright5v1.alpha05,rowname %in% NgnBright.5v1.names) %>% select(rowname,log2FoldChange,padj)
tpm.NgnBright10.5v1_table<-left_join(tpm.NgnBright10.5v1_table,log2NgnBright_table, by=c("GeneID" = "rowname" ))%>% 
  mutate("cell_type" = "iOSN") %>% rename("zone1_TPM"=Z1_NgnBright,"zone2_TPM"=Z2_NgnBright,"zone5_TPM"=Z5_NgnBright)
tpm.NgnDim10.5v1_table<-tpm.mean.table %>% filter(GeneID %in% NgnDim.5v1.names) %>% select(GeneID,Z1_NgnDim,Z2_NgnDim,Z5_NgnDim) %>% arrange(desc(Z5_NgnDim))
log2NgnDim_table<-filter(res.NgnDim5v1.alpha05,rowname %in% NgnDim.5v1.names) %>% select(rowname,log2FoldChange,padj)
tpm.NgnDim10.5v1_table<-left_join(tpm.NgnDim10.5v1_table,log2NgnDim_table, by=c("GeneID" = "rowname" ))%>% 
  mutate("cell_type" = "mOSN") %>% rename("zone1_TPM"=Z1_NgnDim,"zone2_TPM"=Z2_NgnDim,"zone5_TPM"=Z5_NgnDim)

DE.TFs_GO_DNAbinding_TFactivity_2fold<-rbind(tpm.mash10.5v1_table,tpm.mashNgn10.5v1_table,tpm.NgnBright10.5v1_table,tpm.NgnDim10.5v1_table)


########################################################
# padj<0.05, 3 fold change between zone1 and zone5 cells
########################################################
# USE TO MAKE HEATMAP OF A MORE SELECTIVE GROUP OF DIFFERENTIALLY EXPRESSED TFS (FIGURE 4A)

# Find diff expressed TFs
# GBC cells
res05 <- results(dds, contrast=c("condition","Z5Mash","Z1Mash"),alpha=0.05) #more stringent significance cutoff, by default alpha=0.1
summary(res05)
res.mash5v1.alpha05<-rownames_to_column(as.data.frame(res05))
res.mash5v1.alpha05.up<-filter(res.mash5v1.alpha05,rowname %in% TF.vect) %>% filter(padj < 0.05) %>% filter(log2FoldChange > 1.5)
dim(res.mash5v1.alpha05.up)
#23 TFs up in zone5
res.mash5v1.alpha05.down<-filter(res.mash5v1.alpha05,rowname %in% TF.vect) %>% filter(padj < 0.05) %>% filter(log2FoldChange < -1.5) 
dim(res.mash5v1.alpha05.down)
#17 TFs up in zone5
res.mash5v1.alpha05.UPandDOWN<-filter(res.mash5v1.alpha05,rowname %in% TF.vect) %>% filter(padj < 0.05) %>% filter(log2FoldChange < -1.5 |log2FoldChange > 1.5) %>% mutate('celltype'="GBC")
dim(res.mash5v1.alpha05.UPandDOWN)
#40 TFs differentially expresesed
mash.5v1.names<-unique(c(res.mash5v1.alpha05.up$rowname,res.mash5v1.alpha05.down$rowname))


# INP cells
res05 <- results(dds, contrast=c("condition","Z5MashNgn","Z1MashNgn"),alpha=0.05) #more stringent significance cutoff, by default alpha=0.1
summary(res05)
res.mashNgn5v1.alpha05<-rownames_to_column(as.data.frame(res05))
res.mashNgn5v1.alpha05.up<-filter(res.mashNgn5v1.alpha05,rowname %in% TF.vect) %>% filter(padj < 0.05) %>% filter(log2FoldChange > 1.5) 
dim(res.mashNgn5v1.alpha05.up)
#20 TFs up zone5
res.mashNgn5v1.alpha05.down<-filter(res.mashNgn5v1.alpha05,rowname %in% TF.vect) %>% filter(padj < 0.05) %>% filter(log2FoldChange < -1.5) 
dim(res.mashNgn5v1.alpha05.down)
#25 TFs down zone5
res.mashNgn5v1.alpha05.UPandDOWN<-filter(res.mashNgn5v1.alpha05,rowname %in% TF.vect) %>% filter(padj < 0.05) %>% filter(log2FoldChange < -1.5 |log2FoldChange > 1.5) %>% mutate('celltype'="INP")
dim(res.mashNgn5v1.alpha05.UPandDOWN)
#45 TFs diffrentially expressed
mashNgn.5v1.names<-unique(c(res.mashNgn5v1.alpha05.up$rowname,res.mashNgn5v1.alpha05.down$rowname))

# iOSN cells
res05 <- results(dds, contrast=c("condition","Z5NgnBright","Z1NgnBright"),alpha=0.05) #more stringent significance cutoff, by default alpha=0.1
summary(res05)
res.NgnBright5v1.alpha05<-rownames_to_column(as.data.frame(res05))
res.NgnBright5v1.alpha05.up<-filter(res.NgnBright5v1.alpha05,rowname %in% TF.vect) %>% filter(padj < 0.05) %>% filter(log2FoldChange > 1.5) 
dim(res.NgnBright5v1.alpha05.up)
#38 TFs up zone5
res.NgnBright5v1.alpha05.down<-filter(res.NgnBright5v1.alpha05,rowname %in% TF.vect) %>% filter(padj < 0.05) %>% filter(log2FoldChange < -1.5) 
dim(res.NgnBright5v1.alpha05.down)
#24 TFs down zone5
res.NgnBright5v1.alpha05.UPandDOWN<-filter(res.NgnBright5v1.alpha05,rowname %in% TF.vect) %>% filter(padj < 0.05) %>% filter(log2FoldChange < -1.5 |log2FoldChange > 1.5) %>% mutate('celltype'="iOSN")
dim(res.NgnBright5v1.alpha05.UPandDOWN)
#62 TFs differentially expressed
NgnBright.5v1.names<-unique(c(res.NgnBright5v1.alpha05.up$rowname,res.NgnBright5v1.alpha05.down$rowname))

# mOSN cells
res05 <- results(dds, contrast=c("condition","Z5NgnDim","Z1NgnDim"),alpha=0.05) #more stringent significance cutoff, by default alpha=0.1
summary(res05)
res.NgnDim5v1.alpha05<-rownames_to_column(as.data.frame(res05))
res.NgnDim5v1.alpha05.up<-filter(res.NgnDim5v1.alpha05,rowname %in% TF.vect) %>% filter(padj < 0.05) %>% filter(log2FoldChange > 1.5) 
dim(res.NgnDim5v1.alpha05.up)
#14 TFs up zone5
res.NgnDim5v1.alpha05.down<-filter(res.NgnDim5v1.alpha05,rowname %in% TF.vect) %>% filter(padj < 0.05) %>% filter(log2FoldChange < -1.5) 
dim(res.NgnDim5v1.alpha05.down)
#5 TFs down zone5
res.NgnDim5v1.alpha05.UPandDOWN<-filter(res.NgnDim5v1.alpha05,rowname %in% TF.vect) %>% filter(padj < 0.05) %>% filter(log2FoldChange < -1.5 |log2FoldChange > 1.5) %>% mutate('celltype'="OSN")
dim(res.NgnDim5v1.alpha05.UPandDOWN)
#19 TFs differentially expressed
NgnDim.5v1.names<-unique(c(res.NgnDim5v1.alpha05.up$rowname,res.NgnDim5v1.alpha05.down$rowname))

###########
#make table
###########
tpm.mash10.5v1_table<-tpm.mean.table %>% filter(GeneID %in% mash.5v1.names) %>% select(GeneID,Z1_mash,Z2_mash,Z5_mash) %>% arrange(desc(Z5_mash))
log2mash_table<-filter(res.mash5v1.alpha05,rowname %in% mash.5v1.names) %>% select(rowname,log2FoldChange,padj) 
tpm.mash10.5v1_table<-left_join(tpm.mash10.5v1_table,log2mash_table, by=c("GeneID" = "rowname" )) %>% 
  mutate("cell_type" = "GBC") %>% rename("zone1_TPM"=Z1_mash,"zone2_TPM"=Z2_mash,"zone5_TPM"=Z5_mash)
tpm.mashNgn10.5v1_table<-tpm.mean.table %>% filter(GeneID %in% mashNgn.5v1.names) %>% select(GeneID,Z1_mashNgn,Z2_mashNgn,Z5_mashNgn) %>% arrange(desc(Z5_mashNgn))
log2mashNgn_table<-filter(res.mashNgn5v1.alpha05,rowname %in% mashNgn.5v1.names) %>% select(rowname,log2FoldChange,padj )
tpm.mashNgn10.5v1_table<-left_join(tpm.mashNgn10.5v1_table,log2mashNgn_table, by=c("GeneID" = "rowname" )) %>% 
  mutate("cell_type" = "INP") %>% rename("zone1_TPM"=Z1_mashNgn,"zone2_TPM"=Z2_mashNgn,"zone5_TPM"=Z5_mashNgn)
tpm.NgnBright10.5v1_table<-tpm.mean.table %>% filter(GeneID %in% NgnBright.5v1.names) %>% select(GeneID,Z1_NgnBright,Z2_NgnBright,Z5_NgnBright) %>% arrange(desc(Z5_NgnBright))
log2NgnBright_table<-filter(res.NgnBright5v1.alpha05,rowname %in% NgnBright.5v1.names) %>% select(rowname,log2FoldChange,padj)
tpm.NgnBright10.5v1_table<-left_join(tpm.NgnBright10.5v1_table,log2NgnBright_table, by=c("GeneID" = "rowname" ))%>% 
  mutate("cell_type" = "iOSN") %>% rename("zone1_TPM"=Z1_NgnBright,"zone2_TPM"=Z2_NgnBright,"zone5_TPM"=Z5_NgnBright)
tpm.NgnDim10.5v1_table<-tpm.mean.table %>% filter(GeneID %in% NgnDim.5v1.names) %>% select(GeneID,Z1_NgnDim,Z2_NgnDim,Z5_NgnDim) %>% arrange(desc(Z5_NgnDim))
log2NgnDim_table<-filter(res.NgnDim5v1.alpha05,rowname %in% NgnDim.5v1.names) %>% select(rowname,log2FoldChange,padj)
tpm.NgnDim10.5v1_table<-left_join(tpm.NgnDim10.5v1_table,log2NgnDim_table, by=c("GeneID" = "rowname" ))%>% 
  mutate("cell_type" = "mOSN") %>% rename("zone1_TPM"=Z1_NgnDim,"zone2_TPM"=Z2_NgnDim,"zone5_TPM"=Z5_NgnDim)

DE.TFs_GO_DNAbinding_TFactivity_3fold<-rbind(tpm.mash10.5v1_table,tpm.mashNgn10.5v1_table,tpm.NgnBright10.5v1_table,tpm.NgnDim10.5v1_table)

################################
# TFs TPM heatmaps  (Figure 4A)#
################################
# use only TFs that are expressed with at least 15 TPM 

#set up color scale
#add "dummy" gene to set padj and log2FoldChange for color scale color scale, to be later removed from the figure
dummy_gene_mash<-tibble(GeneID=c("dummy_l","dummy_h"),Z1_mash=c(50,50),Z2_mash=c(50,50),Z5_mash=c(50,50),log2FoldChange=c(-6.04,4.87),padj=c(4.44e-36,0.05))
dummy_gene_mashngn<-tibble(GeneID=c("dummy_l","dummy_h"),Z1_mashNgn=c(50,50),Z2_mashNgn=c(50,50),Z5_mashNgn=c(50,50),log2FoldChange=c(-6.04,4.87),padj=c(4.44e-36,0.05))
dummy_gene_ngn<-tibble(GeneID=c("dummy_l","dummy_h"),Z1_NgnBright=c(50,50),Z2_NgnBright=c(50,50),Z5_NgnBright=c(50,50),log2FoldChange=c(-6.04,4.87),padj=c(4.44e-36,0.05))
dummy_gene_ngndim<-tibble(GeneID=c("dummy_l","dummy_h"),Z1_NgnDim=c(50,50),Z2_NgnDim=c(50,50),Z5_NgnDim=c(50,50),log2FoldChange=c(-6.04,4.87),padj=c(4.44e-36,0.05))

#breaksList = c(seq(0, 200, by = 2),seq(202, 400, by = 6))
breaksList = c(seq(0, 200, by = 2),seq(202, 400, by = 6),seq(406, 556, by = 30))
#hmcol <- colorRampPalette(brewer.pal(9, "Greens"))(length(breaksList))
hmcol <- colorRampPalette(c("white","black"))(length(breaksList))

#GBC heatmap
tpm.mash10.5v1<-tpm.mean.table %>% filter(GeneID %in% mash.5v1.names) %>% select(GeneID,Z1_mash,Z2_mash,Z5_mash) %>% filter(Z1_mash > 15 | Z2_mash > 15 | Z5_mash > 15 ) %>% arrange(desc(Z5_mash))
log2mash<-filter(res.mash5v1.alpha05,rowname %in% mash.5v1.names) %>% select(rowname,log2FoldChange,padj) 
tpm.mash10.5v1<-left_join(tpm.mash10.5v1,log2mash, by=c("GeneID" = "rowname" ))
filter(tpm.mash10.5v1,log2FoldChange > 0) %>% dim()
filter(tpm.mash10.5v1,log2FoldChange < 0) %>% dim()
#20
tpm.mash10.5v1<-rbind(tpm.mash10.5v1,dummy_gene_mash)
tpm.mash10.5v1.mat<-as.matrix(tpm.mash10.5v1[,2:6])
rownames(tpm.mash10.5v1.mat)<-tpm.mash10.5v1$GeneID
#pal <- list(log2FoldChange = setNames(brewer.pal(nlevels(tpm.mash10.5v1$log2FoldChange),"RdBu"),levels(tpm.mash10.5v1$log2FoldChange)),padj = setNames(brewer.pal(9, "Purples"),levels(tpm.mash10.5v1$padj)))
pal <- list(log2FoldChange = setNames(brewer.pal(9, "RdBu"),levels(tpm.mash10.5v1$log2FoldChange)),padj = setNames(rev(brewer.pal(9, "Purples")),levels(tpm.mash10.5v1$padj)))
anno<-tpm.mash10.5v1 %>% select(GeneID,log2FoldChange,padj)
anno<-data.frame(anno, row.names = 1)
#p.mash.tpm10.5v1<-pheatmap(t(tpm.mash10.5v1.mat[,1:3]), color = hmcol, annotation_col = anno, annotation_colors = pal, border_color = "black",cluster_rows = FALSE, cluster_cols = FALSE, main = "Mash (TPM 10 cutoff, Z5 vs Z1 sig diff)")
p.mash.tpm10.5v1_dummy<-pheatmap(t(tpm.mash10.5v1.mat[,1:3]), color = hmcol, breaks = breaksList,annotation_col = anno, annotation_colors = pal, border_color = "black",cluster_rows = FALSE, cluster_cols = FALSE, main = "Mash (TPM 15 cutoff, Z5 vs Z1 sig diff)")

#INP heatmap
tpm.mashNgn10.5v1<-tpm.mean.table %>% filter(GeneID %in% mashNgn.5v1.names) %>% select(GeneID,Z1_mashNgn,Z2_mashNgn,Z5_mashNgn) %>% filter(Z1_mashNgn > 15 | Z2_mashNgn > 15 | Z5_mashNgn > 15 ) %>% arrange(desc(Z5_mashNgn))
log2mashNgn<-filter(res.mashNgn5v1.alpha05,rowname %in% mashNgn.5v1.names) %>% select(rowname,log2FoldChange,padj )
tpm.mashNgn10.5v1<-left_join(tpm.mashNgn10.5v1,log2mashNgn, by=c("GeneID" = "rowname" ))
filter(tpm.mashNgn10.5v1,log2FoldChange > 0) %>% dim()
#13
filter(tpm.mashNgn10.5v1,log2FoldChange < 0) %>% dim()
#15
dim(tpm.mashNgn10.5v1)
#28
tpm.mashNgn10.5v1<-rbind(tpm.mashNgn10.5v1,dummy_gene_mashngn)
tpm.mashNgn10.5v1.mat<-as.matrix(tpm.mashNgn10.5v1[,2:5])
rownames(tpm.mashNgn10.5v1.mat)<-tpm.mashNgn10.5v1$GeneID
pal <- list(log2FoldChange = setNames(brewer.pal(9, "RdBu"),levels(tpm.mashNgn10.5v1$log2FoldChange)),padj = setNames(rev(brewer.pal(9, "Purples")),levels(tpm.mashNgn10.5v1$padj)))
anno<-tpm.mashNgn10.5v1 %>% select(GeneID,log2FoldChange,padj)
anno<-data.frame(anno, row.names = 1)
#p.mashNgn.tpm10.5v1<-pheatmap(t(tpm.mashNgn10.5v1.mat[,1:3]), color = hmcol, annotation_col = anno, annotation_colors = pal, border_color = "black",cluster_rows = FALSE, cluster_cols = FALSE, main = "MashNgn (TPM 10 cutoff, Z5 vs Z1 sig diff)")
p.mashNgn.tpm10.5v1_dummy<-pheatmap(t(tpm.mashNgn10.5v1.mat[,1:3]), color = hmcol, breaks = breaksList,annotation_col = anno, annotation_colors = pal, border_color = "black",cluster_rows = FALSE, cluster_cols = FALSE, main = "MashNgn (TPM 15 cutoff, Z5 vs Z1 sig diff)")

#iOSN heatmap
tpm.NgnBright10.5v1<-tpm.mean.table %>% filter(GeneID %in% NgnBright.5v1.names) %>% select(GeneID,Z1_NgnBright,Z2_NgnBright,Z5_NgnBright) %>% filter(Z1_NgnBright > 15 | Z2_NgnBright > 15 | Z5_NgnBright > 15 ) %>% arrange(desc(Z5_NgnBright))
log2NgnBright<-filter(res.NgnBright5v1.alpha05,rowname %in% NgnBright.5v1.names) %>% select(rowname,log2FoldChange,padj)
tpm.NgnBright10.5v1<-left_join(tpm.NgnBright10.5v1,log2NgnBright, by=c("GeneID" = "rowname" ))
filter(tpm.NgnBright10.5v1,log2FoldChange > 0) %>% dim()
#12
filter(tpm.NgnBright10.5v1,log2FoldChange < 0) %>% dim()
#13
dim(tpm.NgnBright10.5v1)
#25
tpm.NgnBright10.5v1<-rbind(tpm.NgnBright10.5v1,dummy_gene_ngn)
tpm.NgnBright10.5v1.mat<-as.matrix(tpm.NgnBright10.5v1[,2:5])
rownames(tpm.NgnBright10.5v1.mat)<-tpm.NgnBright10.5v1$GeneID
pal <- list(log2FoldChange = setNames(brewer.pal(9, "RdBu"),levels(tpm.NgnBright10.5v1$log2FoldChange)),padj = setNames(rev(brewer.pal(9, "Purples")),levels(tpm.NgnBright10.5v1$padj)))
anno<-tpm.NgnBright10.5v1 %>% select(GeneID,log2FoldChange,padj)
anno<-data.frame(anno, row.names = 1)
#p.NgnBright.tpm10.5v1<-pheatmap(t(tpm.NgnBright10.5v1.mat[,1:3]), color = hmcol, annotation_col = anno, annotation_colors = pal,border_color = "black",cluster_rows = FALSE, cluster_cols = FALSE, main = "NgnBright (TPM 10 cutoff, Z5 vs Z1 sig diff)")
p.NgnBright.tpm10.5v1_dummy<-pheatmap(t(tpm.NgnBright10.5v1.mat[,1:3]), color = hmcol, breaks = breaksList,annotation_col = anno, annotation_colors = pal,border_color = "black",cluster_rows = FALSE, cluster_cols = FALSE, main = "NgnBright (TPM 15 cutoff, Z5 vs Z1 sig diff)")

#mOSN heatmap
tpm.NgnDim10.5v1<-tpm.mean.table %>% filter(GeneID %in% NgnDim.5v1.names) %>% select(GeneID,Z1_NgnDim,Z2_NgnDim,Z5_NgnDim) %>% filter(Z1_NgnDim > 15 | Z2_NgnDim > 15 | Z5_NgnDim > 15 ) %>% arrange(desc(Z5_NgnDim))
log2NgnDim<-filter(res.NgnDim5v1.alpha05,rowname %in% NgnDim.5v1.names) %>% select(rowname,log2FoldChange,padj)
tpm.NgnDim10.5v1<-left_join(tpm.NgnDim10.5v1,log2NgnDim, by=c("GeneID" = "rowname" ))
filter(tpm.NgnDim10.5v1,log2FoldChange > 0) %>% dim()
#12
filter(tpm.NgnDim10.5v1,log2FoldChange < 0) %>% dim()
#13
dim(tpm.NgnDim10.5v1)
#9
tpm.NgnDim10.5v1<-rbind(tpm.NgnDim10.5v1,dummy_gene_ngndim)
tpm.NgnDim10.5v1.mat<-as.matrix(tpm.NgnDim10.5v1[,2:5])
rownames(tpm.NgnDim10.5v1.mat)<-tpm.NgnDim10.5v1$GeneID
pal <- list(log2FoldChange = setNames(brewer.pal(9, "RdBu"),levels(tpm.NgnDim10.5v1$log2FoldChange)),padj = setNames(rev(brewer.pal(9, "Purples")),levels(tpm.NgnDim10.5v1$padj)))
anno<-tpm.NgnDim10.5v1 %>% select(GeneID,log2FoldChange,padj)
anno<-data.frame(anno, row.names = 1)
p.NgnDim.tpm10.5v1_dummy<-pheatmap(t(tpm.NgnDim10.5v1.mat[,1:3]), color = hmcol, breaks = breaksList,annotation_col = anno, annotation_colors = pal, border_color = "black",cluster_rows = FALSE, cluster_cols = FALSE, main = "NgnDim (TPM 15 cutoff, Z5 vs Z1 sig diff)")

plot_list<-list(p.mash.tpm10.5v1_dummy[[4]],p.mashNgn.tpm10.5v1_dummy[[4]],p.NgnBright.tpm10.5v1_dummy[[4]],p.NgnDim.tpm10.5v1_dummy[[4]])
grid.arrange(arrangeGrob(grobs= plot_list,nrow=4))

#write to pdf
#g<-arrangeGrob(grobs= plot_list,nrow=4)
#ggsave(file="heatmap.Z5vZ1updown.log2_scale_dummy_GO_TF_3fold_15TPM_black_white.pdf",g, width = 14, height = 10, units = "in" )

