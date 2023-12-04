library(GenomicFeatures)
library(DESeq2)
library(Rsamtools) #needed for BamFileList
library(rtracklayer) #needed to import GTF as GRanges object
library(GenomicAlignments) #used for summarizeOverlaps
library( "RColorBrewer" ) #for making heatmaps
library(gridExtra)
library("pheatmap") #used for heatmaps
library(dplyr)
library(tibble)
library(tidyverse)

#load in annotation
genes <-import("/media/storageA/kevin/annotation/genes+soria+pcdh.gtf")
genes_txdb <- makeTxDbFromGRanges(genes) # creates database
exonsByGene <- exonsBy(genes_txdb, by="gene") #GRangesList object, run exonsByGene$Olfr1507 to see the annotated exons
gene_width <- data.frame(sum(width(exonsByGene)))

#read in my data
sampleTable <- read.delim("/media/storageA/lisa/2020_zonal_paper/analysis/RNAseq_zonal_dev/200424/DEseq2_sampleTable_zonal_dev_200423.txt",header =T)

bamfiles <- filenames <- file.path(sampleTable$data) #creates character vector with paths to bam files
bamfiles <- BamFileList(bamfiles, yieldSize=1000000)
sampleTable

se <- summarizeOverlaps(features=exonsByGene, reads=bamfiles,
                        mode="Union",
                        singleEnd=FALSE,
                        ignore.strand=FALSE) #for truseq data include preprocess.reads=invertStrand


colData(se) <-DataFrame(sampleTable) #is this additional metadata necessary?
colnames(se)<-sampleTable$library

######################
# GET ZONE5 MARKERS ##
######################
##############################################
# sig diff z5NgnDim (mOSN) vs z1NgnDim (mOSN)#
##############################################

#one factor design formula is simple: ~ condition

dds <- DESeqDataSet(se, design = ~ condition)
dds$condition <- relevel(dds$condition, ref="Z1NgnDim") #put control condition first, otherwise it will do put the levels alphabetically
dds <- dds[ rowSums(counts(dds)) > 1, ]

#differential expression analysis
dds <- DESeq(dds) #estimates the size factors (to control for library size), dispersion estimation for each gene, and fitting to a linear model. Returns a DESeqDataSet.
str(dds)


Olfr <- grep( "^Olfr", rownames(dds), value = TRUE)
nonOlfr <- grep( "^Olfr", rownames(dds), invert= TRUE, value = TRUE)

res <-results(dds, name="condition_Z5NgnDim_vs_Z1NgnDim")
resultsNames(dds) #lists all the coefficients


df<-as.data.frame(rownames_to_column(as.data.frame(res)))
write_delim(rownames_to_column(as.data.frame(res)),path = "results.Z5vZ1.tsv",col_names = TRUE,delim = "\t")

#filter out Olfr genes
res.NgnDimv5v1.up5<-rownames_to_column(as.data.frame(res)) %>% filter(log2FoldChange > 1) %>% filter(padj < 0.05) %>% filter(rowname %in% nonOlfr)
res.NgnDimv5v1.down5<-rownames_to_column(as.data.frame(res)) %>% filter(log2FoldChange < -1) %>% filter(padj < 0.05) %>% filter(rowname %in% nonOlfr)

write_delim(data.frame(res.NgnDimv5v1.up5$rowname),path = "Z5vZ1_NgnDim_sig.tsv",col_names = TRUE,delim = "\t")
write_delim(data.frame(res.NgnDimv5v1.down5$rowname),path = "Z1vZ5_NgnDim_sig.tsv",col_names = TRUE,delim = "\t")


TPM.tb<-read.table(file = "/media/storageA/lisa/2020_zonal_paper/analysis/RNAseq_zonal_dev/200424/zonal_dev_TPM_mean.txt", header = TRUE, row.names = 1)

res.NgnDimv5v1.up5<-left_join(res.NgnDimv5v1.up5,TPM.tb, by=c("rowname"="GeneID")) %>% 
  select(rowname,baseMean,log2FoldChange,padj,Z1_NgnDim,Z2_NgnDim,Z5_NgnDim) %>%
  drop_na() %>% filter(Z5_NgnDim > 5) %>% 
  rename("baseMean5v1" = baseMean, "log2FoldChange5v1" = log2FoldChange, "padj5v1" = padj)

res.NgnDimv5v1.down5<-left_join(res.NgnDimv5v1.down5,TPM.tb, by=c("rowname"="GeneID")) %>% 
  select(rowname,baseMean,log2FoldChange,padj,Z1_NgnDim,Z2_NgnDim,Z5_NgnDim) %>%
  drop_na() %>% filter(Z1_NgnDim > 5) %>% 
  rename("baseMean5v1" = baseMean, "log2FoldChange5v1" = log2FoldChange, "padj5v1" = padj)

##############################################
# sig diff z5NgnDim (mOSN) vs z2NgnDim (mOSN)#
##############################################
dds <- DESeqDataSet(se, design = ~ condition)
dds$condition <- relevel(dds$condition, ref="Z2NgnDim") #put control condition first, otherwise it will do put the levels alphabetically
dds <- dds[ rowSums(counts(dds)) > 1, ] #should I use O reads instead?

#differential expression analysis
dds <- DESeq(dds) #estimates the size factors (to control for library size), dispersion estimation for each gene, and fitting to a linear model. Returns a DESeqDataSet.
Olfr <- grep( "^Olfr", rownames(dds), value = TRUE)
nonOlfr <- grep( "^Olfr", rownames(dds), invert= TRUE, value = TRUE)

res <-results(dds, name="condition_Z5NgnDim_vs_Z2NgnDim")
resultsNames(dds) #lists all the coefficients
summary(res)
res <- res[order(res$padj),] #order results by pvalue

#filter out Olfr genes
res.NgnDimv5v2.up5<-rownames_to_column(as.data.frame(res)) %>% filter(log2FoldChange > 1) %>% filter(padj < 0.05)%>% filter(rowname %in% nonOlfr)
res.NgnDimv5v2.down5<-rownames_to_column(as.data.frame(res)) %>% filter(log2FoldChange < -1) %>% filter(padj < 0.05)%>% filter(rowname %in% nonOlfr)

TPM.tb<-read.table(file = "/media/storageA/lisa/2020_zonal_paper/analysis/RNAseq_zonal_dev/200424/zonal_dev_TPM_mean.txt", header = TRUE, row.names = 1)

res.NgnDimv5v2.up5<-left_join(res.NgnDimv5v2.up5,TPM.tb, by=c("rowname"="GeneID")) %>% 
  select(rowname,baseMean,log2FoldChange,padj,Z1_NgnDim,Z2_NgnDim,Z5_NgnDim) %>%
  drop_na() %>% filter(Z5_NgnDim > 5) %>% rename("baseMean5v2" = baseMean, "log2FoldChange5v2" = log2FoldChange, "padj5v2" = padj)

res.NgnDimv5v2.down5<-left_join(res.NgnDimv5v2.down5,TPM.tb, by=c("rowname"="GeneID")) %>% 
  select(rowname,baseMean,log2FoldChange,padj,Z1_NgnDim,Z2_NgnDim,Z5_NgnDim) %>%
  drop_na() %>% filter(Z2_NgnDim > 5) %>% rename("baseMean5v2" = baseMean, "log2FoldChange5v2" = log2FoldChange, "padj5v2" = padj)


########################################
# write Z5 NgnDim (mOSN) Markers files #
########################################
res.NgnDimv5v1or2_up<-full_join(res.NgnDimv5v1.up5,res.NgnDimv5v2.up5)
NgnDimv5v1or2_up_markers<-res.NgnDimv5v1or2_up$rowname
write_delim(res.NgnDimv5v1or2_up,path = "Z5ngnDim_vs1or2_markers.tsv",col_names = TRUE,delim = "\t")

res.NgnDimv5v1and2_up<-inner_join(res.NgnDimv5v1.up5,res.NgnDimv5v2.up5)
NgnDimv5v1and2_up_markers<-res.NgnDimv5v1and2_up$rowname
write_delim(res.NgnDimv5v1and2_up,path = "Z5ngnDim_vs1and2_markers.tsv",col_names = TRUE,delim = "\t")

#############################################
# write Z1 or 2 NgnDim (mOSN) Markers files #
#############################################
res.NgnDimv5v1or2_down<-full_join(res.NgnDimv5v1.down5,res.NgnDimv5v2.down5)
NgnDimv5v1or2_down_markers<-res.NgnDimv5v1or2_down$rowname
write_delim(res.NgnDimv5v1or2_down,path = "Z1or2ngnDim_vs5_markers.tsv",col_names = TRUE,delim = "\t")

res.NgnDimv5v1and2_down<-inner_join(res.NgnDimv5v1.down5,res.NgnDimv5v2.down5)
NgnDimv5v1and2_down_markers<-res.NgnDimv5v1and2_down$rowname
write_delim(res.NgnDimv5v1and2_down,path = "Z1and2ngnDim_vs5_markers.tsv",col_names = TRUE,delim = "\t")


################################################
# LOOK AT CELL TYPE SPECIFIC MARKERS NOT ZONALLY
################################################

# MASH (GBC)
dds <- DESeqDataSet(se, design = ~ zone + tissue)
dds$tissue <- relevel(dds$tissue, ref="Mash") #put control condition first, otherwise it will do put the levels alphabetically
### remove genes with no mapped reads in any condition
dds <- dds[ rowSums(counts(dds)) > 1, ] 

#differential expression analysis
dds <- DESeq(dds) #estimates the size factors (to control for library size), dispersion estimation for each gene, and fitting to a linear model. Returns a DESeqDataSet.
Olfr <- grep( "^Olfr", rownames(dds), value = TRUE)
nonOlfr <- grep( "^Olfr", rownames(dds), invert= TRUE, value = TRUE)

res <-results(dds)
resultsNames(dds) #lists all the coefficients

res <- results(dds, contrast=c("tissue","MashNgn","Mash"))
summary(res)
res <- res[order(res$padj),] #order results by pvalue

res05 <- results(dds, contrast=c("tissue","MashNgn","Mash"),alpha=0.05) #more stringent significance cutoff, by default alpha=0.1
summary(res05)
res05 <- res05[order(res05$padj),] #order results by pvalue

#filter out Olfr genes
res.MashNgnvMash.alpha05_up<-rownames_to_column(as.data.frame(res05)) %>% 
  filter(log2FoldChange > 1) %>% 
  filter(rowname %in% nonOlfr)#%>% filter(baseMean > 30)
#2364 up genes

res.MashNgnvMash.alpha05_down<-rownames_to_column(as.data.frame(res05)) %>% 
  filter(log2FoldChange < -1) %>% 
  filter(rowname %in% nonOlfr)#%>% filter(baseMean > 30)
#2155 down(up in Mash) genes

TPM.tb<-read.table(file = "/media/storageA/lisa/2020_zonal_paper/analysis/RNAseq_zonal_dev/200424/zonal_dev_TPM_mean.txt", header = TRUE, row.names = 1)
res.MashNgnvMash.alpha05_down<-left_join(res.MashNgnvMash.alpha05_down,TPM.tb, by=c("rowname"="GeneID")) %>% 
  select(rowname,baseMean,log2FoldChange,padj,Z1_mash,Z2_mash,Z5_mash,Z1_mashNgn,Z2_mashNgn,Z5_mashNgn) %>%
  drop_na() %>% #filter(baseMean > 30)
  filter(Z1_mash > 10 & Z2_mash > 10 & Z5_mash > 10)

Mash_markers<-res.MashNgnvMash.alpha05_down$rowname
write_delim(res.MashNgnvMash.alpha05_down_v2,path = "Mash_p05_2fold_TPM10_notZonal_markers.tsv",col_names = TRUE,delim = "\t")


# mashNgn (INP)
TPM.tb<-read.table(file = "/media/storageA/lisa/2020_zonal_paper/analysis/RNAseq_zonal_dev/200424/zonal_dev_TPM_mean.txt", header = TRUE, row.names = 1)
res.MashNgnvMash.alpha05_up<-left_join(res.MashNgnvMash.alpha05_up,TPM.tb, by=c("rowname"="GeneID")) %>% 
  select(rowname,baseMean,log2FoldChange,padj,Z1_mash,Z2_mash,Z5_mash,Z1_mashNgn,Z2_mashNgn,Z5_mashNgn) %>%
  drop_na() %>% #filter(baseMean > 30)
  filter(Z1_mashNgn > 10 & Z2_mashNgn > 10 & Z5_mashNgn > 10)

MashNgn_markers<-res.MashNgnvMash.alpha05_up$rowname
write_delim(res.MashNgnvMash.alpha05_up,path = "MashNgn_p05_2fold_TPM10_notZonal_markers.tsv",col_names = TRUE,delim = "\t")

### mashNgn (INP) vsNgnBright (iOSN)
dds <- DESeqDataSet(se, design = ~ zone + tissue)
dds$tissue <- relevel(dds$tissue, ref="MashNgn") #put control condition first, otherwise it will do put the levels alphabetically
### remove genes with no mapped reads in any condition
dds <- dds[ rowSums(counts(dds)) > 1, ] #should I use O reads instead?

#differential expression analysis
dds <- DESeq(dds) #estimates the size factors (to control for library size), dispersion estimation for each gene, and fitting to a linear model. Returns a DESeqDataSet.
Olfr <- grep( "^Olfr", rownames(dds), value = TRUE)
nonOlfr <- grep( "^Olfr", rownames(dds), invert= TRUE, value = TRUE)

res <-results(dds)
resultsNames(dds) #lists all the coefficients

#two commands to specify to specify the contrast we want for the results table in multi-factor design.
#res <- results(dds, name="condition_treated_vs_untreated")
#res <- results(dds, contrast=c("condition","treated","untreated"))
res <- results(dds, contrast=c("tissue","NgnBright","MashNgn"))
summary(res)
res <- res[order(res$padj),] #order results by pvalue

res05 <- results(dds, contrast=c("tissue","NgnBright","MashNgn"),alpha=0.05) #more stringent significance cutoff, by default alpha=0.1
summary(res05)
res05 <- res05[order(res05$padj),] #order results by pvalue

#filter out Olfr genes
res.NgnBrightvMashNgn.alpha05_up<-rownames_to_column(as.data.frame(res05)) %>% 
  filter(log2FoldChange > 1) %>% 
  filter(rowname %in% nonOlfr)%>% #filter(baseMean > 30)
  
  
  res.NgnBrightvMashNgn.alpha05_down<-rownames_to_column(as.data.frame(res05)) %>% 
  filter(log2FoldChange < -1) %>% 
  filter(rowname %in% nonOlfr)%>% #filter(baseMean > 30)
  
  
#up in mashNgn (INP) relative to ngnbright (iOSN)
TPM.tb<-read.table(file = "/media/storageA/lisa/2020_zonal_paper/analysis/RNAseq_zonal_dev/200424/zonal_dev_TPM_mean.txt", header = TRUE, row.names = 1)
res.NgnBrightvMashNgn.alpha05_down<-left_join(res.NgnBrightvMashNgn.alpha05_down,TPM.tb, by=c("rowname"="GeneID")) %>% 
  select(rowname,baseMean,log2FoldChange,padj,Z1_mashNgn,Z2_mashNgn,Z5_mashNgn,Z1_NgnBright,Z2_NgnBright,Z5_NgnBright) %>%
  drop_na() %>% #filter(baseMean > 30)
  filter(Z1_mashNgn > 10 & Z2_mashNgn > 10 & Z5_mashNgn > 10)

res.MashNgn.alpha05_peak<-inner_join(res.MashNgnvMash.alpha05_up,res.NgnBrightvMashNgn.alpha05_down, by="rowname") 
#only 19 genes peak
MashNgn_markers_peak<-res.MashNgn.alpha05_peak$rowname
write_delim(res.MashNgn.alpha05_peak,path = "MashNgn_p05_2fold_TPM10_notZonal_markers_peak.tsv",col_names = TRUE,delim = "\t")

#ngnBright (iOSN) markers
TPM.tb<-read.table(file = "/media/storageA/lisa/2020_zonal_paper/analysis/RNAseq_zonal_dev/200424/zonal_dev_TPM_mean.txt", header = TRUE, row.names = 1)
res.NgnBrightvMashNgn.alpha05_up<-left_join(res.NgnBrightvMashNgn.alpha05_up,TPM.tb, by=c("rowname"="GeneID")) %>% 
  select(rowname,baseMean,log2FoldChange,padj,Z1_mashNgn,Z2_mashNgn,Z5_mashNgn,Z1_NgnBright,Z2_NgnBright,Z5_NgnBright) %>%
  drop_na() %>% #filter(baseMean > 30)
  filter(Z1_NgnBright > 10 & Z2_NgnBright > 10 & Z5_NgnBright > 10)

NgnBright_marker<-res.NgnBrightvMashNgn.alpha05_up$rowname
write_delim(res.NgnBrightvMashNgn.alpha05_up,path = "NgnBright_p05_2fold_TPM10_notZonal_markers.tsv",col_names = TRUE,delim = "\t")

### ngnDim (mOSN)
dds <- DESeqDataSet(se, design = ~ zone + tissue)
dds$tissue <- relevel(dds$tissue, ref="NgnBright") #put control condition first, otherwise it will do put the levels alphabetically
### remove genes with no mapped reads in any condition
dds <- dds[ rowSums(counts(dds)) > 1, ] #should I use O reads instead?

#differential expression analysis
dds <- DESeq(dds) #estimates the size factors (to control for library size), dispersion estimation for each gene, and fitting to a linear model. Returns a DESeqDataSet.
Olfr <- grep( "^Olfr", rownames(dds), value = TRUE)
nonOlfr <- grep( "^Olfr", rownames(dds), invert= TRUE, value = TRUE)

res <-results(dds)
resultsNames(dds) #lists all the coefficients

res <- results(dds, contrast=c("tissue","NgnDim","NgnBright"))
summary(res)
res <- res[order(res$padj),] #order results by pvalue

res05 <- results(dds, contrast=c("tissue","NgnDim","NgnBright"),alpha=0.05) #more stringent significance cutoff, by default alpha=0.1
summary(res05)
res05 <- res05[order(res05$padj),] #order results by pvalue

#filter out Olfr genes
res.NgnDimvNgnBright.alpha05_up<-rownames_to_column(as.data.frame(res05)) %>% 
  filter(log2FoldChange > 1) %>% 
  filter(rowname %in% nonOlfr)#%>% filter(baseMean > 30)

res.NgnDimvNgnBright.alpha05_down<-rownames_to_column(as.data.frame(res05)) %>% 
  filter(log2FoldChange < -1) %>% 
  filter(rowname %in% nonOlfr)#%>% filter(baseMean > 30)

#up in Ngnbright (iOSN) relative to ngndim (mOSN)
TPM.tb<-read.table(file = "/media/storageA/lisa/2020_zonal_paper/analysis/RNAseq_zonal_dev/200424/zonal_dev_TPM_mean.txt", header = TRUE, row.names = 1)
res.NgnDimvNgnBright.alpha05_down<-left_join(res.NgnDimvNgnBright.alpha05_down,TPM.tb, by=c("rowname"="GeneID")) %>% 
  select(rowname,baseMean,log2FoldChange,padj,Z1_NgnBright,Z2_NgnBright,Z5_NgnBright,Z1_NgnDim,Z2_NgnDim,Z5_NgnDim) %>%
  drop_na() %>% #filter(baseMean > 30)
  filter(Z1_NgnBright > 10 & Z2_NgnBright > 10 & Z5_NgnBright > 10)

res.NgnBright.alpha05_peak<-inner_join(res.NgnBrightvMashNgn.alpha05_up,res.NgnDimvNgnBright.alpha05_down, by="rowname") 
#only 62 genes peak
NgnBright_markers_peak<-res.NgnBright.alpha05_peak$rowname
write_delim(res.NgnBright.alpha05_peak,path = "NgnBright_p05_2fold_TPM10_notZonal_markers_peak.tsv",col_names = TRUE,delim = "\t")

#ngnDim (mOSN) markers
TPM.tb<-read.table(file = "/media/storageA/lisa/2020_zonal_paper/analysis/RNAseq_zonal_dev/200424/zonal_dev_TPM_mean.txt", header = TRUE, row.names = 1)
res.NgnDimvNgnBright.alpha05_up<-left_join(res.NgnDimvNgnBright.alpha05_up,TPM.tb, by=c("rowname"="GeneID")) %>% 
  select(rowname,baseMean,log2FoldChange,padj,Z1_NgnBright,Z2_NgnBright,Z5_NgnBright,Z1_NgnDim,Z2_NgnDim,Z5_NgnDim) %>%
  drop_na() %>% #filter(baseMean > 30)
  filter(Z1_NgnDim > 10 & Z2_NgnDim > 10 & Z5_NgnDim > 10)

NgnDim_marker<-res.NgnDimvNgnBright.alpha05_up$rowname
write_delim(res.NgnDimvNgnBright.alpha05_up,path = "NgnDim_p05_2fold_TPM10_notZonal_markers.tsv",col_names = TRUE,delim = "\t")
