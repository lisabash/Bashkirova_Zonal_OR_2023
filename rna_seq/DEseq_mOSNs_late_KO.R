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

#zonal OR color scale
zonal_fine_colors_fish2<-c("#990000","#FF0000","#FFCC00","#33CC33","#00CCCC","#0000FF")
zonal_fine_colors_dark_fish2<-c("#660000","#990000","#999900","#003300","#003366","#000066")

#load in annotation
genes <-import("/media/storageA/kevin/annotation/genes+soria+pcdh.gtf")
genes_txdb <- makeTxDbFromGRanges(genes) # creates database
exonsByGene <- exonsBy(genes_txdb, by="gene") #GRangesList object
gene_width <- data.frame(sum(width(exonsByGene)))

#read in my data
sampleTable <- read.delim("/media/storageA/lisa/2020_zonal_paper/analysis/RNAseq_nfiKO/200910_lateKO/single_end/DEseq.lateNFIko.200910.sample.table_figure.txt",header =T)
sampleTable

bamfiles <- filenames <- file.path(sampleTable$data) #creates character vector with paths to bam files
bamfiles <- BamFileList(bamfiles, yieldSize=1000000)

se <- summarizeOverlaps(features=exonsByGene, reads=bamfiles,
                        mode="Union",
                        singleEnd=TRUE) #for truseq data include preprocess.reads=invertStrand


colData(se) <-DataFrame(sampleTable) #is this additional metadata necessary?
colnames(se)<-sampleTable$library

#one factor design formula is simple: ~ condition
#to do multifactor design (for example to control for a batch), then use design formula ~ batch + condition
#se$date<-factor(se$date)
dds <- DESeqDataSet(se, design = ~ geno)
dds$geno <- relevel(dds$geno, ref="wt") #put control condition first, otherwise it will do put the levels alphabetically
### remove genes with no mapped reads in any condition
dds <- dds[ rowSums(counts(dds)) > 1, ] #should I use O reads instead?

#differential expression analysis
dds <- DESeq(dds) #estimates the size factors (to control for library size), dispersion estimation for each gene, and fitting to a linear model. Returns a DESeqDataSet.
Olfr <- grep( "^Olfr", rownames(dds), value = TRUE)
nonOlfr <- grep( "^Olfr", rownames(dds), invert= TRUE, value = TRUE)

#############################
### Calculate FPKM
#############################
fpkm_scaled <- as.data.frame(fpkm(dds, robust = TRUE)) # robus normalizes on DEseq "scaling factor"

### Calculate normalized mean FPKM 
fpkm_scaled2<-fpkm_scaled %>% rownames_to_column(var="GeneID") %>% 
  gather(condition,fpkm,-GeneID) %>% separate(condition,sep="_",into=c("V1","V2","rep","date")) %>% 
  unite("condition", V1, V2) 
fpkm_scaled.mean <- fpkm_scaled2 %>% group_by(condition,GeneID) %>% summarize(mean_fpkm=mean(fpkm))
fpkm_scaled.mean.table <- spread(fpkm_scaled.mean,condition,mean_fpkm)

### Calculate unnormalized FPKM to generate TPM
fpkm_unscaled <- as.data.frame(fpkm(dds, robust = FALSE)) %>% rownames_to_column(var="GeneID")
fpkm_unscaled <- fpkm_unscaled %>% gather(condition,fpkm,-GeneID) %>% 
  separate(condition,sep="_",into=c("V1","V2","rep","date")) %>% 
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

zonal_anno_fine_tb<-read_tsv("/media/storageA/lisa/2020_zonal_paper/annotation/zonal_anno_TanLocus_soriaNames/zonal_anno_2020_fine_TanLocus_SoriaName.bed", col_names = c("chr","start","stop","name","score","strand","zone"))
zonal_anno_fine_tb<-zonal_anno_fine_tb %>% select(name,zone)
tpm.mean.table.zone<-right_join(tpm.mean.table,zonal_anno_fine_tb,by=c("GeneID" = "name"))
tpm.mean.table.zone_log2<-tpm.mean.table.zone %>% mutate("log2"= log2((RNA.NfiABXko_OMPCre+0.01)/(RNA.control_OMPCre+0.01)))
tpm.mean.table.zone_log2$zone<-factor(tpm.mean.table.zone_log2$zone,levels = c("classI","1","2","3","4","5"))

########## BOXPLOT OR EXPRESSION (Figure 4D)

p.log<-ggplot(tpm.mean.table.zone_log2, mapping = aes(zone,log2, color=zone, fill = zone)) + ylim(-10,10)+
  geom_boxplot(alpha = 0.7,  lwd =1.2, outlier.shape = NA)+ #geom_boxplot(alpha = 0.1, outlier.shape = NA , lwd =1.2 ) +
  scale_color_manual(values=zonal_fine_colors_dark_fish2) + 
  scale_fill_manual(values = zonal_fine_colors_fish2)+
  theme(legend.position="none",axis.text.y = element_text(colour= "black", size = 12),axis.text.x = element_text(colour= "black", size = 14), axis.title.x = element_blank(),axis.title.y = element_text(colour= "black", size = 14),plot.title = element_text(hjust = 0.5, size=20))+
  ggtitle("NfiABX lateKO/wt log2")+ ylab("TPM")
p.log
#ggsave(file="boxplot_late_nfiKO_log2_MOE_scale10_SEalign.pdf",p.log, width = 4, height = 4, units = "in" ,useDingbats = FALSE)

##### WILCOXON TEST
lk_wilcox5<-tpm.mean.table.zone_log2 %>% 
  filter(zone == 5) %>% 
  #select(-zone, -log2) %>%
  gather("sample","value", -GeneID,-zone, -log2) %>% drop_na()
wilcox.test(value ~ sample,lk_wilcox5) 
#p-value = 0.1652

lk_wilcox4<-tpm.mean.table.zone_log2 %>% 
  filter(zone == 4) %>% 
  gather("sample","value", -GeneID,-zone, -log2) %>% drop_na()
wilcox.test(value ~ sample,lk_wilcox4) 
#p-value = 0.9862

lk_wilcox3<-tpm.mean.table.zone_log2 %>% 
  filter(zone == 3) %>% 
  gather("sample","value", -GeneID,-zone, -log2) %>% drop_na()
wilcox.test(value ~ sample,lk_wilcox3) 
#p-value = 0.6585

lk_wilcox2<-tpm.mean.table.zone_log2 %>% 
  filter(zone == 2) %>% 
  gather("sample","value", -GeneID,-zone, -log2) %>% drop_na()
wilcox.test(value ~ sample,lk_wilcox2) 
#p-value = 0.3849

lk_wilcox1<-tpm.mean.table.zone_log2 %>% 
  filter(zone == 1) %>% 
  gather("sample","value", -GeneID,-zone, -log2) %>% drop_na()
wilcox.test(value ~ sample,lk_wilcox1) 
#p-value = 0.2611

lk_wilcoxclassi<-tpm.mean.table.zone_log2 %>% 
  filter(zone == "classI") %>% 
  gather("sample","value", -GeneID,-zone, -log2) %>% drop_na()
wilcox.test(value ~ sample,lk_wilcoxclassi) 
#p-value = 0.3577
