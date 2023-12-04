library(GenomicFeatures)
library(DESeq2)
library(Rsamtools) #needed for BamFileList
library(rtracklayer) #needed to import GTF as GRanges object
library(GenomicAlignments) #used for summarizeOverlaps
library( "RColorBrewer" ) #for making heatmaps
library(gridExtra)
library("pheatmap") #used for heatmaps
library(ggrepel) #for genes points on plot
library(dplyr)
library(tibble)
library(tidyverse)

#color palette used for ORs by zone
zonal_fine_colors_fish<-c("#FF0000","#FFCC00","#33CC33","#00CCCC","#0000FF","#999999")
zonal_fine_colors_dark_fish<-c("#990000","#999900","#003300","#003366","#000066","black")

#load in annotation
genes <-import("/media/storageA/kevin/annotation/genes+soria+pcdh.gtf")
genes_txdb <- makeTxDbFromGRanges(genes) # creates database
exonsByGene <- exonsBy(genes_txdb, by="gene") #GRangesList object
gene_width <- data.frame(sum(width(exonsByGene)))

#read in my data
sampleTable <- read.delim("/media/storageA/lisa/2020_zonal_paper/analysis/RNAseq_nfiKO/210128_deseq_bonly_het/DEseq.NfiKO_zone5.sample_table_210125.txt",header =T)
bamfiles <- filenames <- file.path(sampleTable$data) #creates character vector with paths to bam files
bamfiles <- BamFileList(bamfiles, yieldSize=1000000)
sampleTable

se <- summarizeOverlaps(features=exonsByGene, reads=bamfiles,
                        mode="Union",
                        singleEnd=FALSE,
                        ignore.strand=FALSE) #for truseq data include preprocess.reads=invertStrand

colData(se) <-DataFrame(sampleTable) #is this additional metadata necessary?
colnames(se)<-sampleTable$library

dds <- DESeqDataSet(se, design = ~ condition)
dds$condition <- relevel(dds$condition, ref="wt_Z5") #put control condition first, otherwise it will put the levels alphabetically
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
write.table(fpkm_scaled, file="z5_nfi_FPKM.txt", sep="\t",quote=FALSE)
fpkm_scaled2<-fpkm_scaled %>% rownames_to_column(var="GeneID") %>% gather(condition,fpkm,-GeneID) %>% separate(condition,sep="_",into=c("condition","rep"))
fpkm_scaled.mean <- fpkm_scaled2 %>% group_by(condition,GeneID) %>% summarize(mean_fpkm=mean(fpkm))
fpkm_scaled.mean.table <- spread(fpkm_scaled.mean,condition,mean_fpkm)
write.table(fpkm_scaled.mean.table, file="z5_nfi_FPKM_mean.txt",quote=FALSE)

### Calculate unnormalized FPKM to generate TPM
fpkm_unscaled <- as.data.frame(fpkm(dds, robust = FALSE)) %>% rownames_to_column(var="GeneID")
fpkm_unscaled <- fpkm_unscaled %>% gather(condition,fpkm,-GeneID) %>% separate(condition,sep="_",into=c("condition","rep"))
fpkm_unscaled.mean<- fpkm_unscaled  %>% group_by(condition,GeneID) %>% summarize(mean_fpkm=mean(fpkm))
fpkm.mean.table <- spread(fpkm_unscaled.mean,condition,mean_fpkm)

### Calculate TPM
tpm <- fpkm_unscaled %>% group_by(condition,rep) %>% summarize(total_fpkm = sum(fpkm)) %>% full_join(fpkm_unscaled) %>% mutate(tpm = (fpkm / total_fpkm) *1000000) %>% dplyr::select(GeneID,condition,rep,tpm)
tpm.mean <- tpm %>% group_by(condition,GeneID) %>% summarize(mean_tpm = mean(tpm))
tpm.sd <- tpm %>% group_by(condition,GeneID) %>% summarize(sd = sd(tpm))
tpm.table <- tpm %>% mutate(lib = paste0(condition,"_",rep)) %>% group_by(lib,GeneID) %>% dplyr::select(GeneID,condition = lib,tpm) %>% spread(condition,tpm)
write.table(tpm.table, file="z5_nfi_TPM.txt",quote=FALSE)
tpm.mean.table <- spread(tpm.mean,condition,mean_tpm)
tpm.sd.table <- spread(tpm.sd,condition,sd)
write.table(as.data.frame(tpm.mean.table), file="z5_nfi_TPM_mean_210128.txt")

# load OR zonal annotation
zonal_anno_fine_tb<-read_tsv("/media/storageA/lisa/2020_zonal_paper/annotation/zonal_anno_TanLocus_soriaNames/zonal_anno_2020_fine_TanLocus_SoriaName_sorted.bed", col_names = c("chr","start","stop","name","num","strand","zone"))
zonal_anno_fine_tb<-zonal_anno_fine_tb %>% select(name,zone)

tpm.table.zone<-right_join(tpm.table,zonal_anno_fine_tb,by=c("GeneID" = "name"))
tpm.mean.table.zone<-right_join(tpm.mean.table,zonal_anno_fine_tb,by=c("GeneID" = "name"))
tpm.mean.table.zone$zone<-factor(tpm.mean.table.zone$zone, levels = c("classI","1","2","3","4","5"))


#################################################################
# Boxplots OR expression in NFI cKO genotypes by zone (Figure 4E)
#################################################################

p.abx<-ggplot(tpm.mean.table.zone, mapping = aes(zone,NfiABXko, color=zone, fill = zone)) + #ylim(0,300)+
  geom_boxplot(alpha = 0.7,  lwd =1.2, outlier.shape = NA)+ #geom_boxplot(alpha = 0.1, outlier.shape = NA , lwd =1.2 ) +
  scale_color_manual(values=zonal_fine_colors_dark_fish) + scale_fill_manual(values = zonal_fine_colors_fish) + geom_jitter(alpha =0.7)+
  theme(legend.position="none",axis.text.y = element_text(colour= "black", size = 12),axis.text.x = element_text(colour= "black", size = 14), axis.title.x = element_blank(),axis.title.y = element_text(colour= "black", size = 14),plot.title = element_text(hjust = 0.5, size=20))+
  ggtitle("zone5 NfiABXko")+ ylab("TPM")+geom_text_repel(aes(label=if_else(NfiABXko > 80, as.character(GeneID),"")), size=3)
p.abx

p.ab<-ggplot(tpm.mean.table.zone, mapping = aes(zone,NfiABko, color=zone, fill = zone)) + #ylim(0,300)+
  geom_boxplot(alpha = 0.7,  lwd =1.2 ,outlier.shape = NA )+ #geom_boxplot(alpha = 0.1, outlier.shape = NA , lwd =1.2 ) +
  scale_color_manual(values=zonal_fine_colors_dark_fish) + scale_fill_manual(values = zonal_fine_colors_fish) + geom_jitter(alpha =0.7)+
  theme(legend.position="none",axis.text.y = element_text(colour= "black", size = 12),axis.text.x = element_text(colour= "black", size = 14), axis.title.x = element_blank(),axis.title.y = element_text(colour= "black", size = 14),plot.title = element_text(hjust = 0.5, size=20))+
  ggtitle("zone5 NfiABko")+ ylab("TPM")+geom_text_repel(aes(label=if_else(NfiABko > 80, as.character(GeneID),"")), size=3)
p.ab

p.x<-ggplot(tpm.mean.table.zone, mapping = aes(zone,NfiXko, color=zone, fill = zone)) + #ylim(0,1000)+
  geom_boxplot(alpha = 0.7,  lwd =1.2,outlier.shape = NA  )+ #geom_boxplot(alpha = 0.1, outlier.shape = NA , lwd =1.2 ) +
  scale_color_manual(values=zonal_fine_colors_dark_fish) + scale_fill_manual(values = zonal_fine_colors_fish) + geom_jitter(alpha =0.7)+
  theme(legend.position="none",axis.text.y = element_text(colour= "black", size = 12),axis.text.x = element_text(colour= "black", size = 14), axis.title.x = element_blank(),axis.title.y = element_text(colour= "black", size = 14),plot.title = element_text(hjust = 0.5, size=20))+
  ggtitle("zone5 NfiXko")+ ylab("TPM")+geom_text_repel(aes(label=if_else(NfiXko > 80, as.character(GeneID),"")), size=3)+
  scale_y_continuous(breaks=c(100,200, 300, 400,500,600,700,800,900),limits=c(0, 1000))
p.x

p.b<-ggplot(tpm.mean.table.zone, mapping = aes(zone,NfiBko, color=zone, fill = zone)) + #ylim(0,1000)+
  geom_boxplot(alpha = 0.7,  lwd =1.2 ,outlier.shape = NA )+ #geom_boxplot(alpha = 0.1, outlier.shape = NA , lwd =1.2 ) +
  scale_color_manual(values=zonal_fine_colors_dark_fish) + scale_fill_manual(values = zonal_fine_colors_fish) + geom_jitter(alpha =0.7)+
  theme(legend.position="none",axis.text.y = element_text(colour= "black", size = 12),axis.text.x = element_text(colour= "black", size = 14), axis.title.x = element_blank(),axis.title.y = element_text(colour= "black", size = 14),plot.title = element_text(hjust = 0.5, size=20))+
  ggtitle("zone5 NfiBko")+ ylab("TPM")+geom_text_repel(aes(label=if_else(NfiBko > 80, as.character(GeneID),"")), size=3)+
  scale_y_continuous(breaks=c(100,200, 300, 400,500,600,700,800,900),limits=c(0, 1000))
p.b

p.abxhet<-ggplot(tpm.mean.table.zone, mapping = aes(zone,NfiABXhet, color=zone, fill = zone)) + #ylim(0,1250)+
  geom_boxplot(alpha = 0.7,  lwd =1.2 ,outlier.shape = NA )+ #geom_boxplot(alpha = 0.1, outlier.shape = NA , lwd =1.2 ) +
  scale_color_manual(values=zonal_fine_colors_dark_fish) + scale_fill_manual(values = zonal_fine_colors_fish) + geom_jitter(alpha =0.7)+
  theme(legend.position="none",axis.text.y = element_text(colour= "black", size = 12),axis.text.x = element_text(colour= "black", size = 14), axis.title.x = element_blank(),axis.title.y = element_text(colour= "black", size = 14),plot.title = element_text(hjust = 0.5, size=20))+
  ggtitle("zone5 NfiABXhet")+ ylab("TPM")+geom_text_repel(aes(label=if_else(NfiABXhet > 80, as.character(GeneID),"")), size=3)+
  scale_y_continuous(breaks=c(200,400,600,800,1000,1200),limits=c(0, 1250))
p.abxhet

p.wt<-ggplot(tpm.mean.table.zone, mapping = aes(zone,wt, color=zone, fill = zone)) + ylim(0,1000)+
  geom_boxplot(alpha = 0.7,  lwd =1.2 , outlier.shape = NA)+ #geom_boxplot(alpha = 0.1, outlier.shape = NA , lwd =1.2 ) +
  scale_color_manual(values=zonal_fine_colors_dark_fish) + scale_fill_manual(values = zonal_fine_colors_fish) + geom_jitter(alpha =0.7)+
  theme(legend.position="none",axis.text.y = element_text(colour= "black", size = 12),axis.text.x = element_text(colour= "black", size = 14), axis.title.x = element_blank(),axis.title.y = element_text(colour= "black", size = 14),plot.title = element_text(hjust = 0.5, size=20))+
  ggtitle("zone5 WT")+ ylab("TPM")+geom_text_repel(aes(label=if_else(wt > 80, as.character(GeneID),"")), size=3)+ 
  scale_y_continuous(breaks=c(100,200, 300, 400,500,600,700,800,900),limits=c(0, 1000))
p.wt

grid.arrange(p.abx, p.ab, p.abxhet, p.b, p.x, p.wt ,ncol = 6,nrow=1)
#write to pdf
#g<-arrangeGrob(p.abx, p.ab, p.abxhet, p.x, p.wt ,ncol = 6,nrow=1)
#ggsave(file="z5_nfiKO_labeled_points_large_classIfirst.pdf",g, width = 18, height = 8, units = "in" ,useDingbats = FALSE)

#####################################################################
### get OR significance count for Pie charts (Supplemental Figure 4C)
#####################################################################

#### ABX
res <- results(dds, contrast=c("condition","NFIABXko_Z5","wt_Z5"),cooksCutoff=FALSE)
summary(res)
res <- res[order(res$padj),] #order results by pvalue
df.res<-rownames_to_column(as.data.frame(res))
df.res.Olfr<- df.res%>% filter(rowname %in% Olfr)

df.res.Olfr.zone<-inner_join(df.res.Olfr,zonal_anno_fine_tb,by=c("rowname" = "name"))
df.res.Olfr.zone.sig<-df.res.Olfr.zone %>% mutate(fold=if_else(log2FoldChange < 0,"down","up")) %>%
  mutate(sig= if_else(!(is.na(padj)), if_else(padj < 0.05,"sig","not"),"not")) #%>% #View()
group_by(zone,fold,sig) %>% summarise(n())

abx_summary<-df.res.Olfr.zone.sig %>% select(rowname,sig,zone,fold) %>% 
  group_by(zone,sig,fold)%>% summarise(n())

abx_summary %>% mutate("sel"= if_else(sig == "not", "not_sig",if_else(fold == "up", "sig_up","sig_down"))) %>%
  group_by(zone,sel) 

df<-data.frame("zone"= c("zone1","zone1","zone1","zone2","zone2","zone2",
                         "zone3","zone3","zone3","zone4","zone4","zone4",
                         "zone5","zone5","zone5","classI","classI","classI"), 
               "sel" = c("not_sig","sig_up","sig_down","not_sig","sig_up","sig_down",
                         "not_sig","sig_up","sig_down","not_sig","sig_up","sig_down",
                         "not_sig","sig_up","sig_down","not_sig","sig_up","sig_down"),
               "num" = c(228,5,0,  167,107,0,  74,88,0,  51,9,82,  1,0,42,  125,4,0))

df_old<-data.frame("zone"= c("zone1","zone1","zone1","zone2","zone2","zone2",
                             "zone3","zone3","zone3","zone4","zone4","zone4",
                             "zone5","zone5","zone5","classI","classI","classI"), 
                   "sel" = c("not_sig","sig_up","sig_down","not_sig","sig_up","sig_down",
                             "not_sig","sig_up","sig_down","not_sig","sig_up","sig_down",
                             "not_sig","sig_up","sig_down","not_sig","sig_up","sig_down"),
                   "num" = c(228,6,0,  168,105,0,  70,88,0,  51,9,82,  1,0,41,  125,4,0))

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold"),
    legend.position = "none"
  )

abx.1<-ggplot(filter(df,zone == "zone1"), aes(x="", y=num, fill=sel))+
  geom_bar(width = 1, stat = "identity")+coord_polar("y", start=0)+ 
  scale_fill_manual(values=c("gray","blue3","red3"))+blank_theme
abx.1
abx.2<-ggplot(filter(df,zone == "zone2"), aes(x="", y=num, fill=sel))+
  geom_bar(width = 1, stat = "identity")+coord_polar("y", start=0)+ 
  scale_fill_manual(values=c("gray","blue3","red3"))+blank_theme
abx.3<-ggplot(filter(df,zone == "zone3"), aes(x="", y=num, fill=sel))+
  geom_bar(width = 1, stat = "identity")+coord_polar("y", start=0)+ 
  scale_fill_manual(values=c("gray","blue3","red3"))+blank_theme
abx.4<-ggplot(filter(df,zone == "zone4"), aes(x="", y=num, fill=sel))+
  geom_bar(width = 1, stat = "identity")+coord_polar("y", start=0)+ 
  scale_fill_manual(values=c("gray","blue3","red3"))+blank_theme
abx.5<-ggplot(filter(df,zone == "zone5"), aes(x="", y=num, fill=sel))+
  geom_bar(width = 1, stat = "identity")+coord_polar("y", start=0)+ 
  scale_fill_manual(values=c("gray","blue3","red3"))+blank_theme
abx.classI<-ggplot(filter(df,zone == "classI"), aes(x="", y=num, fill=sel))+
  geom_bar(width = 1, stat = "identity")+coord_polar("y", start=0)+ 
  scale_fill_manual(values=c("gray","blue3","red3"))+blank_theme
grid.arrange(abx.1,abx.2,abx.3,abx.4,abx.5,abx.classI, ncol = 6)

#### AB
res <- results(dds, contrast=c("condition","NFIABko_Z5","wt_Z5"),cooksCutoff=FALSE)
summary(res)
res <- res[order(res$padj),] #order results by pvalue
df.res<-rownames_to_column(as.data.frame(res))
df.res.Olfr<- df.res%>% filter(rowname %in% Olfr)

df.res.Olfr.zone<-inner_join(df.res.Olfr,zonal_anno_fine_tb,by=c("rowname" = "name"))
df.res.Olfr.zone.sig<-df.res.Olfr.zone %>% mutate(fold=if_else(log2FoldChange < 0,"down","up")) %>%
  mutate(sig= if_else(!(is.na(padj)), if_else(padj < 0.05,"sig","not"),"not")) #%>%
group_by(zone,fold,sig) %>% summarise(n())

ab_summary<-df.res.Olfr.zone.sig %>% select(rowname,sig,zone,fold) %>% 
  group_by(zone,sig,fold)%>% summarise(n())

ab_summary %>% mutate("sel"= if_else(sig == "not", "not_sig",if_else(fold == "up", "sig_up","sig_down"))) %>%
  group_by(zone,sel) 

df<-data.frame("zone"= c("zone1","zone1","zone1","zone2","zone2","zone2",
                         "zone3","zone3","zone3","zone4","zone4","zone4",
                         "zone5","zone5","zone5","classI","classI","classI"), 
               "sel" = c("not_sig","sig_up","sig_down","not_sig","sig_up","sig_down",
                         "not_sig","sig_up","sig_down","not_sig","sig_up","sig_down",
                         "not_sig","sig_up","sig_down","not_sig","sig_up","sig_down"),
               "num" = c(230,3,0,  243,31,0,  89,72,0,  96,28,18,  2,0,41,  126,3,0))
class(df$num)


ab.1<-ggplot(filter(df,zone == "zone1"), aes(x="", y=num, fill=sel))+
  geom_bar(width = 1, stat = "identity")+coord_polar("y", start=0)+ 
  scale_fill_manual(values=c("gray","blue3","red3"))+blank_theme
ab.1
ab.2<-ggplot(filter(df,zone == "zone2"), aes(x="", y=num, fill=sel))+
  geom_bar(width = 1, stat = "identity")+coord_polar("y", start=0)+ 
  scale_fill_manual(values=c("gray","blue3","red3"))+blank_theme
ab.3<-ggplot(filter(df,zone == "zone3"), aes(x="", y=num, fill=sel))+
  geom_bar(width = 1, stat = "identity")+coord_polar("y", start=0)+ 
  scale_fill_manual(values=c("gray","blue3","red3"))+blank_theme
ab.4<-ggplot(filter(df,zone == "zone4"), aes(x="", y=num, fill=sel))+
  geom_bar(width = 1, stat = "identity")+coord_polar("y", start=0)+ 
  scale_fill_manual(values=c("gray","blue3","red3"))+blank_theme
ab.5<-ggplot(filter(df,zone == "zone5"), aes(x="", y=num, fill=sel))+
  geom_bar(width = 1, stat = "identity")+coord_polar("y", start=0)+ 
  scale_fill_manual(values=c("gray","blue3","red3"))+blank_theme
ab.classI<-ggplot(filter(df,zone == "classI"), aes(x="", y=num, fill=sel))+
  geom_bar(width = 1, stat = "identity")+coord_polar("y", start=0)+ 
  scale_fill_manual(values=c("gray","blue3","red3"))+blank_theme
grid.arrange(ab.1,ab.2,ab.3,ab.4,ab.5,ab.classI, ncol = 6)

#### B
res <- results(dds, contrast=c("condition","NFIBko_Z5","wt_Z5"),cooksCutoff=FALSE)
summary(res)
res <- res[order(res$padj),] #order results by pvalue
df.res<-rownames_to_column(as.data.frame(res))
df.res.Olfr<- df.res%>% filter(rowname %in% Olfr)

df.res.Olfr.zone<-inner_join(df.res.Olfr,zonal_anno_fine_tb,by=c("rowname" = "name"))
df.res.Olfr.zone.sig<-df.res.Olfr.zone %>% mutate(fold=if_else(log2FoldChange < 0,"down","up")) %>%
  mutate(sig= if_else(!(is.na(padj)), if_else(padj < 0.05,"sig","not"),"not")) #%>%
group_by(zone,fold,sig) %>% summarise(n())

df.res.Olfr.zone  %>% summarise(n())

b_summary<-df.res.Olfr.zone.sig %>% select(rowname,sig,zone,fold) %>% 
  group_by(zone,sig,fold)%>% summarise(n())

b_summary %>% mutate("sel"= if_else(sig == "not", "not_sig",if_else(fold == "up", "sig_up","sig_down"))) %>%
  group_by(zone,sel) 

df<-data.frame("zone"= c("zone1","zone1","zone1","zone2","zone2","zone2",
                         "zone3","zone3","zone3","zone4","zone4","zone4",
                         "zone5","zone5","zone5","classI","classI","classI"), 
               "sel" = c("not_sig","sig_up","sig_down","not_sig","sig_up","sig_down",
                         "not_sig","sig_up","sig_down","not_sig","sig_up","sig_down",
                         "not_sig","sig_up","sig_down","not_sig","sig_up","sig_down"),
               "num" = c(233,0,0,  270,4,0,  147,13,1,  136,5,1,  41,0,2,  129,0,0))
class(df$num)

b.1<-ggplot(filter(df,zone == "zone1"), aes(x="", y=num, fill=sel))+
  geom_bar(width = 1, stat = "identity")+coord_polar("y", start=0)+ 
  scale_fill_manual(values=c("gray","blue3","red3"))+blank_theme
b.1
b.2<-ggplot(filter(df,zone == "zone2"), aes(x="", y=num, fill=sel))+
  geom_bar(width = 1, stat = "identity")+coord_polar("y", start=0)+ 
  scale_fill_manual(values=c("gray","blue3","red3"))+blank_theme
b.3<-ggplot(filter(df,zone == "zone3"), aes(x="", y=num, fill=sel))+
  geom_bar(width = 1, stat = "identity")+coord_polar("y", start=0)+ 
  scale_fill_manual(values=c("gray","blue3","red3"))+blank_theme
b.4<-ggplot(filter(df,zone == "zone4"), aes(x="", y=num, fill=sel))+
  geom_bar(width = 1, stat = "identity")+coord_polar("y", start=0)+ 
  scale_fill_manual(values=c("gray","blue3","red3"))+blank_theme
b.5<-ggplot(filter(df,zone == "zone5"), aes(x="", y=num, fill=sel))+
  geom_bar(width = 1, stat = "identity")+coord_polar("y", start=0)+ 
  scale_fill_manual(values=c("gray","blue3","red3"))+blank_theme
b.classI<-ggplot(filter(df,zone == "classI"), aes(x="", y=num, fill=sel))+
  geom_bar(width = 1, stat = "identity")+coord_polar("y", start=0)+ 
  scale_fill_manual(values=c("gray","blue3","red3"))+blank_theme
grid.arrange(b.1,b.2,b.3,b.4,b.5,b.classI, ncol = 6)

#### X
res <- results(dds, contrast=c("condition","NFIXko_Z5","wt_Z5"),cooksCutoff=FALSE)
summary(res)
res <- res[order(res$padj),] #order results by pvalue
df.res<-rownames_to_column(as.data.frame(res))
df.res.Olfr<- df.res%>% filter(rowname %in% Olfr)

df.res.Olfr.zone<-inner_join(df.res.Olfr,zonal_anno_fine_tb,by=c("rowname" = "name"))
df.res.Olfr.zone.sig<-df.res.Olfr.zone %>% mutate(fold=if_else(log2FoldChange < 0,"down","up")) %>%
  mutate(sig= if_else(!(is.na(padj)), if_else(padj < 0.05,"sig","not"),"not")) #%>%group_by(zone,fold,sig) %>% summarise(n())

df.res.Olfr.zone  %>% summarise(n())

x_summary<-df.res.Olfr.zone.sig %>% select(rowname,sig,zone,fold) %>% 
  group_by(zone,sig,fold)%>% summarise(n())

x_summary %>% mutate("sel"= if_else(sig == "not", "not_sig",if_else(fold == "up", "sig_up","sig_down"))) %>%
  group_by(zone,sel) 

df<-data.frame("zone"= c("zone1","zone1","zone1","zone2","zone2","zone2",
                         "zone3","zone3","zone3","zone4","zone4","zone4",
                         "zone5","zone5","zone5","classI","classI","classI"), 
               "sel" = c("not_sig","sig_up","sig_down","not_sig","sig_up","sig_down",
                         "not_sig","sig_up","sig_down","not_sig","sig_up","sig_down",
                         "not_sig","sig_up","sig_down","not_sig","sig_up","sig_down"),
               "num" = c(231,0,2,  271,0,3,  145,0,16,  137,0,5,  43,0,0,  129,0,0))
class(df$num)

x.1<-ggplot(filter(df,zone == "zone1"), aes(x="", y=num, fill=sel))+
  geom_bar(width = 1, stat = "identity")+coord_polar("y", start=0)+ 
  scale_fill_manual(values=c("gray","blue3","red3"))+blank_theme
x.1
x.2<-ggplot(filter(df,zone == "zone2"), aes(x="", y=num, fill=sel))+
  geom_bar(width = 1, stat = "identity")+coord_polar("y", start=0)+ 
  scale_fill_manual(values=c("gray","blue3","red3"))+blank_theme
x.3<-ggplot(filter(df,zone == "zone3"), aes(x="", y=num, fill=sel))+
  geom_bar(width = 1, stat = "identity")+coord_polar("y", start=0)+ 
  scale_fill_manual(values=c("gray","blue3","red3"))+blank_theme
x.4<-ggplot(filter(df,zone == "zone4"), aes(x="", y=num, fill=sel))+
  geom_bar(width = 1, stat = "identity")+coord_polar("y", start=0)+ 
  scale_fill_manual(values=c("gray","blue3","red3"))+blank_theme
x.5<-ggplot(filter(df,zone == "zone5"), aes(x="", y=num, fill=sel))+
  geom_bar(width = 1, stat = "identity")+coord_polar("y", start=0)+ 
  scale_fill_manual(values=c("gray","blue3","red3"))+blank_theme
x.classI<-ggplot(filter(df,zone == "classI"), aes(x="", y=num, fill=sel))+
  geom_bar(width = 1, stat = "identity")+coord_polar("y", start=0)+ 
  scale_fill_manual(values=c("gray","blue3","red3"))+blank_theme
grid.arrange(x.1,x.2,x.3,x.4,x.5,x.classI, ncol = 6)

#### ABX het
res <- results(dds, contrast=c("condition","NFIABXhet_Z5","wt_Z5"),cooksCutoff=FALSE)
summary(res)
res <- res[order(res$padj),] #order results by pvalue
df.res<-rownames_to_column(as.data.frame(res))
df.res.Olfr<- df.res%>% filter(rowname %in% Olfr)

df.res.Olfr.zone<-inner_join(df.res.Olfr,zonal_anno_fine_tb,by=c("rowname" = "name"))
df.res.Olfr.zone.sig<-df.res.Olfr.zone %>% mutate(fold=if_else(log2FoldChange < 0,"down","up")) %>%
  mutate(sig= if_else(!(is.na(padj)), if_else(padj < 0.05,"sig","not"),"not")) #%>%
group_by(zone,fold,sig) %>% summarise(n())

df.res.Olfr.zone  %>% summarise(n())

abx_het_summary<-df.res.Olfr.zone.sig %>% select(rowname,sig,zone,fold) %>% 
  group_by(zone,sig,fold)%>% summarise(n())

abx_het_summary %>% mutate("sel"= if_else(sig == "not", "not_sig",if_else(fold == "up", "sig_up","sig_down"))) %>%
  group_by(zone,sel) 

df<-data.frame("zone"= c("zone1","zone1","zone1","zone2","zone2","zone2",
                         "zone3","zone3","zone3","zone4","zone4","zone4",
                         "zone5","zone5","zone5","classI","classI","classI"), 
               "sel" = c("not_sig","sig_up","sig_down","not_sig","sig_up","sig_down",
                         "not_sig","sig_up","sig_down","not_sig","sig_up","sig_down",
                         "not_sig","sig_up","sig_down","not_sig","sig_up","sig_down"),
               "num" = c(233,0,0,  267,6,1,  146,15,0,  125,16,1,  37,0,6,  129,0,0))
class(df$num)

abx_het.1<-ggplot(filter(df,zone == "zone1"), aes(x="", y=num, fill=sel))+
  geom_bar(width = 1, stat = "identity")+coord_polar("y", start=0)+ 
  scale_fill_manual(values=c("gray","blue3","red3"))+blank_theme
abx_het.1
abx_het.2<-ggplot(filter(df,zone == "zone2"), aes(x="", y=num, fill=sel))+
  geom_bar(width = 1, stat = "identity")+coord_polar("y", start=0)+ 
  scale_fill_manual(values=c("gray","blue3","red3"))+blank_theme
abx_het.3<-ggplot(filter(df,zone == "zone3"), aes(x="", y=num, fill=sel))+
  geom_bar(width = 1, stat = "identity")+coord_polar("y", start=0)+ 
  scale_fill_manual(values=c("gray","blue3","red3"))+blank_theme
abx_het.4<-ggplot(filter(df,zone == "zone4"), aes(x="", y=num, fill=sel))+
  geom_bar(width = 1, stat = "identity")+coord_polar("y", start=0)+ 
  scale_fill_manual(values=c("gray","blue3","red3"))+blank_theme
abx_het.5<-ggplot(filter(df,zone == "zone5"), aes(x="", y=num, fill=sel))+
  geom_bar(width = 1, stat = "identity")+coord_polar("y", start=0)+ 
  scale_fill_manual(values=c("gray","blue3","red3"))+blank_theme
abx_het.classI<-ggplot(filter(df,zone == "classI"), aes(x="", y=num, fill=sel))+
  geom_bar(width = 1, stat = "identity")+coord_polar("y", start=0)+ 
  scale_fill_manual(values=c("gray","blue3","red3"))+blank_theme
grid.arrange(abx_het.1,abx_het.2,abx_het.3,abx_het.4,abx_het.5,abx_het.classI, ncol = 6)

grid.arrange(abx.1,abx.2,abx.3,abx.4,abx.5,abx.classI, 
             ab.1,ab.2,ab.3,ab.4,ab.5,ab.classI,
             abx_het.1,abx_het.2,abx_het.3,abx_het.4,abx_het.5,abx_het.classI,
             b.1,b.2,b.3,b.4,b.5,b.classI,
             x.1,x.2,x.3,x.4,x.5,x.classI,
             ncol = 6, nrow = 5)
g<-arrangeGrob(abx.1,abx.2,abx.3,abx.4,abx.5,abx.classI, 
               ab.1,ab.2,ab.3,ab.4,ab.5,ab.classI,
               abx_het.1,abx_het.2,abx_het.3,abx_het.4,abx_het.5,abx_het.classI,
               b.1,b.2,b.3,b.4,b.5,b.classI,
               x.1,x.2,x.3,x.4,x.5,x.classI,
               ncol = 6, nrow = 5)
ggsave(file="pie.z5_nfiKO_genos_abxko_abko_abxhet_bko_xko.pdf",g, width = 8, height = 8, units = "in" ,useDingbats = FALSE)


#################################################################
# Volcano plots NFI ABX cKO non-OR genes (Supplemental Figure 4F)
#################################################################

res <- results(dds, contrast=c("condition","NFIABXko_Z5","wt_Z5"))
summary(res)
res <- res[order(res$padj),] #order results by pvalue
df.res<-rownames_to_column(as.data.frame(res))

# to look at changes in mOSN markers we want to exclude genes that have a zonal expression pattern
### load in zonal blacklist: significantly differentially expressed genes with p-adj 0.05 between z5 mOSNs vs z1 mOSNs or z5 mOSNs vs z2 mOSNs
z5v1_blacklist<-read_delim(file = "/media/storageA/lisa/2020_zonal_paper/analysis/RNAseq_zonal_dev/200528_dev_markers/results.Z5vZ1.fromDIMonlyDE2_sig_nonOlfr.tsv" ,delim = "\t")
z5v1_blacklist_markers<-z5v1_blacklist$rowname
z5v2_blacklist<-read_delim(file = "/media/storageA/lisa/2020_zonal_paper/analysis/RNAseq_zonal_dev/200528_dev_markers/results.Z5vZ2.fromDIMonlyDE2_sig_nonOlfr.tsv" ,delim = "\t")
z5v2_blacklist_markers<-z5v2_blacklist$rowname
blacklist<-unique(c(z5v1_blacklist_markers,z5v2_blacklist_markers))
#484 blacklisted non Olfr genes

#load in mOSN markers and exclude blacklist genes
# use the 200 genes most enriched in mOSNs
OMP_markers<-read_delim(file = "/media/storageA/lisa/2020_zonal_paper/analysis/RNAseq_zonal_dev/200519/NgnDim_p05_2fold_TPM10_notZonal_markers.tsv" ,delim = "\t")
OMP_markers<-OMP_markers %>% filter(!(rowname %in% blacklist)) %>% head(n=200)
OMP_markers<-OMP_markers$rowname

df.res.OMP<-df.res %>% filter(rowname %in% OMP_markers) %>% 
  mutate("log"=if_else(log2FoldChange > -1 & log2FoldChange < 1,"not2fold","2fold"))%>% mutate("sig"=if_else(padj < 0.05 & log == "2fold","sigfold",if_else(padj >= 0.05, "notsig","signotfold" ))) #%>%

# load in zone 5 upregulated DE genes
Z5_v1or2markers<-read_delim(file = "/media/storageA/lisa/2020_zonal_paper/analysis/RNAseq_zonal_dev/200519/Z5ngnDim_vs1or2_markers.tsv" ,delim = "\t")
Z5_v1or2markers<-Z5_v1or2markers$rowname

df.res.Z5v1or2<-df.res %>% filter(rowname %in% Z5_v1or2markers) %>% 
  mutate("log"=if_else(log2FoldChange > -1 & log2FoldChange < 1,"not2fold","2fold"))%>% mutate("sig"=if_else(padj < 0.05 & log == "2fold","sigfold",if_else(padj >= 0.05, "notsig","signotfold" )))

# load in z1 up DE genes
Z1or2_v5markers<-read_delim(file = "/media/storageA/lisa/2020_zonal_paper/analysis/RNAseq_zonal_dev/200519/Z1or2ngnDim_vs5_markers.tsv" ,delim = "\t")
Z1or2_v5markers<-Z1or2_v5markers$rowname

df.res.Z1or2v5<-df.res %>% filter(rowname %in% Z1or2_v5markers) %>% 
  mutate("log"=if_else(log2FoldChange > -1 & log2FoldChange < 1,"not2fold","2fold"))%>% mutate("sig"=if_else(padj < 0.05 & log == "2fold","sigfold",if_else(padj >= 0.05, "notsig","signotfold" )))


#volcano plots
p.omp<-ggplot(df.res.OMP, aes(x=log2FoldChange,y=-log10(padj),col = sig)) + 
  geom_point(size=3, alpha=0.5) + 
  geom_text_repel(aes(label=if_else(-log10(padj) >25  , as.character(rowname),"")), size=4)+
  xlim(-10,10)+ ylim (0,100)+ 
  theme(plot.title = element_text(hjust = 0.5, color = "black", size =20),legend.position = "none", axis.title.x = element_text(color="black", size=12),axis.title.y = element_text(color="black", size=12), axis.text.x = element_text(color="black",size=10),axis.text.y = element_text(color="black",size=10))+ 
  scale_color_manual(values = c("black", "red","blue")) + 
  geom_vline(xintercept = c(-1,1), color = "gray2",linetype="dashed")+ggtitle("OSN markers")
p.omp

p.zone5<-ggplot(df.res.Z5v1or2, aes(x=log2FoldChange,y=-log10(padj),col = sig)) + 
  geom_point(size=3, alpha=0.5) + 
  geom_text_repel(aes(label=if_else(-log10(padj) > 25 , as.character(rowname),"")), size=4)+
  xlim(-10,10)+ ylim (0,100)+ 
  theme(plot.title = element_text(hjust = 0.5, color = "black", size =20),legend.position = "none", axis.title.x = element_text(color="black", size=12),axis.title.y = element_text(color="black", size=12), axis.text.x = element_text(color="black",size=10),axis.text.y = element_text(color="black",size=10))+ 
  scale_color_manual(values = c("black", "red", "blue"))+ 
  geom_vline(xintercept = c(-1,1), color = "gray2",linetype="dashed")+ggtitle("Zone5 markers")
p.zone5

p.zone1<-ggplot(df.res.Z1or2v5, aes(x=log2FoldChange,y=-log10(padj),col = sig)) + 
  geom_point(size=3, alpha=0.5) + 
  geom_text_repel(aes(label=if_else(-log10(padj) >25 , as.character(rowname),"")), size=4)+
  xlim(-10,10)+ ylim (0,100)+ 
  theme(plot.title = element_text(hjust = 0.5, color = "black", size =20),legend.position = "none", axis.title.x = element_text(color="black", size=12),axis.title.y = element_text(color="black", size=12), axis.text.x = element_text(color="black",size=10),axis.text.y = element_text(color="black",size=10))+ 
  scale_color_manual(values = c("black", "red", "blue"))+ 
  geom_vline(xintercept = c(-1,1), color = "gray2",linetype="dashed")+ggtitle("Zone1-2 markers")
p.zone1

grid.arrange(p.omp,  p.zone5,p.zone1,ncol = 3,nrow=1)
#g<-arrangeGrob(p.omp,  p.zone5,p.zone1,ncol = 3,nrow=1)
#ggsave(file="2023_nfiKO_OMPmarkers_zonalmarkers_volcano_wide.pdf",g, width = 9, height = 5, units = "in" )

