setwd("/media/storageA/lisa/2020_zonal_paper/analysis/hic/200925_repeat_heatmaps/")
library(data.table)
library(dplyr)
library(ggplot2)
library(ggpointdensity)
library(tidyr)
library(stringr)
library(readr) # read_delim
library(tibble)
library(cowplot)
library(RColorBrewer)

# resolution of analysis
bin_size <- 50000
dump_dir <- "/media/storageD/distiller/dump.balanced/"


## FUNCTIONS ##
format_trans_path <- function(name){
  path <- str_c(dump_dir,"trans.",name,".",bin_size,".balanced.txt.gz")
  return(path)
}
format_cis_path <- function(name){
  path <- str_c(dump_dir,"cis.",name,".",bin_size,".balanced.txt.gz")
  return(path)
}

bedToBait <- function(fname) {
  y<-read_tsv(fname,col_names = c("chrom", "start", "stop","name"),comment="#") %>%
    mutate(chrNo = as.numeric(str_replace(chrom,"chr",""))) %>% 
    arrange(chrNo,start) %>% 
    mutate(bin = floor(start/bin_size)*bin_size) %>% 
    mutate(bait = str_c(chrNo,"_",bin)) %>% dplyr::select(bait) %>% distinct()
  return(y)
}

bedToBait <- function(x) {
  y <- read.delim(x,header=FALSE,comment.char="#") %>%
    mutate(chrNo = as.numeric(str_replace(V1,"chr",""))) %>%
    filter(!is.na(chrNo)) %>%
    arrange(chrNo,V2) %>%
    mutate(bin = floor(as.numeric(V2)/bin_size)*bin_size) %>% 
    mutate(bait = make_bait_str(chrNo, bin)) %>%
    dplyr::select(bait) %>%
    distinct()
  return(y)
}


bedToBait2 <- function(fname) {
  y<-read_tsv(fname,col_names = c("chrom", "start", "stop","name","score","strand"),comment="#") %>%
    mutate(chrNo = as.numeric(str_replace(chrom,"chr",""))) %>% 
    arrange(chrNo,start) %>% 
    mutate(bin = floor(start/bin_size)*bin_size) %>% 
    mutate(bait = str_c(chrNo,"_",bin)) %>% dplyr::select(bait) %>% distinct()
  return(y)
}
#### LOAD DATA #####
extract_contacts_Trans <- function(name,loci){
  write.table(dplyr::select(loci,bait),file="loci.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
  system(command="sort -k1,1 loci.txt > loci.sort.txt",wait=TRUE)
  join.command <- str_c("echo \"join -1 1 -2 1 loci.sort.txt <(zcat ",dump_dir,"trans.",name,".",bin_size,".balanced.txt.gz) > loci.contacts.txt\" | bash")
  system(command=join.command,wait=TRUE)
  loci_data <- fread("loci.contacts.txt",sep=" ",header=FALSE,col.names=c("bait","prey","contacts"))
  system(command="rm loci.txt loci.sort.txt loci.contacts.txt")
  return(loci_data)
}

extract_contacts_Cis <- function(name,loci){
  write.table(dplyr::select(loci,bait),file="loci.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
  system(command="sort -k1,1 loci.txt > loci.sort.txt",wait=TRUE)
  join.command <- str_c("echo \"join -1 1 -2 1 loci.sort.txt <(zcat ",dump_dir,"cis.",name,".",bin_size,".balanced.txt.gz) > loci.contacts.txt\" | bash")
  system(command=join.command,wait=TRUE)
  loci_data <- fread("loci.contacts.txt",sep=" ",header=FALSE,col.names=c("bait","prey","contacts"))
  system(command="rm loci.txt loci.sort.txt loci.contacts.txt")
  return(loci_data)
}

label_cis_trans <- function(combinations) {
  labeled.combinations <- combinations %>% mutate(bait_sp = bait) %>% mutate(prey_sp = prey) %>% separate(prey_sp,c("p.chr","p.coord"),"_",convert=TRUE) %>% separate(bait_sp,c("b.chr","b.coord"),"_",convert=TRUE) %>% mutate(type= ifelse(p.chr == b.chr,"cis","trans")) %>% 
    mutate(distance = ifelse(type == "cis",abs(p.coord - b.coord),"NA")) %>% mutate(distance = as.numeric(distance))
  return(labeled.combinations)
}

process_data <- function(hic.data){
  hic.counts_per_bait <- group_by(hic.data,bait) %>% summarize(total_contacts_per_bait = sum(contacts))
  hic.counts <- left_join(hic.data,hic.counts_per_bait)
  return(hic.counts)
}
############################
# DEFINING INPUT BED FILES #
############################

#at 50kb
Islands <- bedToBait("/media/storageA/kevin/annotation/Greek_Islands.nonX.50kb.mm10.bed")%>% mutate(type="Islands")
#Islands_slop50kb<-bedToBait("/media/storageA/kevin/annotation/Greek_Islands.slop50kb.mm10.merged.50kb.bed") %>% mutate(type="Islands") %>% drop_na()
#Islands <- Islands_slop50kb

OR_Clusters <- bedToBait("/media/storageA/kevin/annotation/ORClusters.ordered.no-chrX.mm10.50kb.mm10.bed") %>% mutate(type="OR_Clusters")
All_Bins <- bedToBait("/media/storageA/kevin/annotation/mm10_assembled.50kb.bed")

####
#old zonal annotation
#zone1_bins <- bedToBait2("/media/storageA/lisa/annotations/zone1_genes_mm10_nofishOR.bed") %>% mutate(type="zone1") %>% drop_na()
#zone2to3_bins <- bedToBait2("/media/storageA/lisa/annotations/zone2to3_genes_mm10.bed") %>% mutate(type="zone2to3") %>% drop_na()
#zone4to5_bins <- bedToBait2("/media/storageA/lisa/annotations/zone4to5_genes_mm10.bed") %>% mutate(type="zone4to5") %>% drop_na()

# zone annotation
#zone1_bins <- bedToBait("/media/storageA/lisa/2020_zonal_paper/annotation/zonal_anno_2020_zone1NOclassi.bed") %>% mutate(type="zone1")%>% drop_na()
#zone2to3_bins <- bedToBait("/media/storageA/lisa/2020_zonal_paper/annotation/zonal_anno_2020_zone23.bed") %>% mutate(type="zone2to3")%>% drop_na()
#zone4to5_bins <- bedToBait("/media/storageA/lisa/2020_zonal_paper/annotation/zonal_anno_2020_zone45.bed") %>% mutate(type="zone4to5")%>% drop_na()

## newest zonal annotation
zone1_bins <- bedToBait2("/media/storageA/lisa/2020_zonal_paper/annotation/zonal_anno_TanLocus_soriaNames/zonal_anno_2020_fine_TanLocus_SoriaName_zone1.bed") %>% mutate(type="zone1")
zone2to3_bins <- bedToBait2("/media/storageA/lisa/2020_zonal_paper/annotation/zonal_anno_TanLocus_soriaNames/zonal_anno_2020_fine_TanLocus_SoriaName_zone2_3.sorted.bed") %>% mutate(type="zone2to3")
zone4to5_bins <- bedToBait2("/media/storageA/lisa/2020_zonal_paper/annotation/zonal_anno_TanLocus_soriaNames/zonal_anno_2020_fine_TanLocus_SoriaName_zone45.bed") %>% mutate(type="zone4to5")

zonal_bins <-rbind(zone1_bins,zone2to3_bins,zone4to5_bins) %>% group_by(bait) %>% summarize(bin.zones = paste(type, collapse = ",")) %>% inner_join(x=All_Bins) %>% 
  mutate(bin.zones = factor(bin.zones,levels=c("zone1","zone1,zone2to3","zone2to3","zone2to3,zone4to5","zone4to5","zone1,zone4to5","zone1,zone2to3,zone4to5")))
###### Analyze interactions between sets of loci
## define bin combinations
Islands.combinations <- CJ(bait = Islands$bait, prey=Islands$bait) %>% as.data.frame
#OR_Clusters.combinations <- CJ(bait = OR_Clusters$bait, prey=OR_Clusters$bait) %>% as.data.frame %>% anti_join(Islands.combinations)

names_list = c(
  "Z5-NgnBright_replicates_merged",
  "Z1-NgnBright_replicates_merged",
  "Z5.NfiWT.merged",
  "Z5.NfiKO.merged",
  "Z5.OMP_merged",
  "Z1.OMP_merged",
  "P2_merged")

plot_output_filenames = unlist(lapply(names_list, function (name) { return(str_c(name, ".new_anno.mean_contacts.hm.balanced.vmax0.00055.pdf")) }))

for (i in 1:length(names_list)) {
  
  cell_type_name <- names_list[i]
  plot_output_filename = plot_output_filenames[i]
  print(cell_type_name)
  print("running ...")
  
  ### Libraries
  Z5.OMP.name <- cell_type_name
  Z5.OMP.trans_hic_path <-format_trans_path(Z5.OMP.name)
  Z5.OMP.cis_hic_path <-format_cis_path(Z5.OMP.name)
  #####################################
  # Specify bins for which to load all data
  #####################################
  Selected_Bins <- rbind(Islands,OR_Clusters) %>% 
    dplyr::select(bait,type) %>% 
    group_by(bait) %>% 
    summarize(anno = paste(type, collapse = ",")) %>% 
    as.data.frame()
  
  
  ### Load ALL contacts for Selected_Bins
  Z5.OMP.data.trans <- extract_contacts_Trans(Z5.OMP.name,Selected_Bins)
  Z5.OMP.data.cis <- extract_contacts_Cis(Z5.OMP.name,Selected_Bins)
  Z5.OMP.data <- rbind(Z5.OMP.data.trans,Z5.OMP.data.cis) 
  Z5.OMP.data <- process_data(Z5.OMP.data)
  group_by(Z5.OMP.data, total_contacts_per_bait) %>% summarize(n())
  
  # filter out OR clusters that aren't int he data
  # and form combinsations of pairs of remaining OR clusters
  OR_Clusters.good <- left_join(OR_Clusters, Z5.OMP.data) %>% filter(!is.na(contacts)) %>% select(bait, type) %>% distinct()
  OR_Clusters.combinations <- CJ(bait = OR_Clusters.good$bait, prey=OR_Clusters.good$bait) %>% as.data.frame %>% anti_join(Islands.combinations)
  
  ## Extract data
  Z5.OMP.ORClusters_to_ORClusters <- left_join(OR_Clusters.combinations,Z5.OMP.data) %>%  
    dplyr::select(bait,prey,contacts) %>% 
    replace_na(list(contacts = 0)) %>% 
    mutate(anno = "ORClusters_to_ORClusters") %>% label_cis_trans() 
  Z5.OMP.ORClusters_to_ORClusters.zonal <- Z5.OMP.ORClusters_to_ORClusters %>% 
    inner_join(zonal_bins %>% 
                 dplyr::select(prey=bait,prey.zones = bin.zones)) %>% 
    inner_join(zonal_bins %>% 
                 dplyr::select(bait,bait.zones = bin.zones)) %>% 
    arrange(b.chr,b.coord,p.chr,p.coord) %>% 
    mutate(bait = factor(bait), prey = factor(prey))
  
  ##########################
  # aggregate by zone
  #############################
  my_palette <- colorRampPalette(brewer.pal(n = 9, name = 'Reds'))(99)
  col_breaks = c(seq(0,0.00055,length=99))
  col_fun = colorRamp2(col_breaks,my_palette)
  
  Z5.OMP.ORClusters_to_ORClusters.zonalSummary <- Z5.OMP.ORClusters_to_ORClusters.zonal %>% filter(type == "trans") %>% 
    filter(!is.na(prey.zones)) %>% filter(!is.na(bait.zones))%>% 
    group_by(bait.zones,prey.zones) %>% 
    summarize(number_of_pairs = n(),mean_contacts = mean(contacts)) %>% 
    filter(bait.zones %in% c("zone1", "zone2to3", "zone4to5")) %>%  
    filter(prey.zones %in% c("zone1", "zone2to3", "zone4to5"))
  
  
  Z5.OMP.ORClusters_to_ORClusters.zonalSummary.number_of_pairs.mat <- Z5.OMP.ORClusters_to_ORClusters.zonalSummary %>% 
    dplyr::select(bait.zones,prey.zones,number_of_pairs) %>% 
    spread(prey.zones,number_of_pairs) %>% 
    column_to_rownames(var="bait.zones") %>% as.matrix()
  
  Z5.OMP.ORClusters_to_ORClusters.zonalSummary.mean_contacts.mat <- Z5.OMP.ORClusters_to_ORClusters.zonalSummary %>% 
    dplyr::select(bait.zones,prey.zones,mean_contacts) %>% 
    spread(prey.zones,mean_contacts) %>% 
    column_to_rownames(var="bait.zones") %>% as.matrix()
  
  
  heatmap.plot<-Heatmap(Z5.OMP.ORClusters_to_ORClusters.zonalSummary.mean_contacts.mat,
                        column_order=colnames(Z5.OMP.ORClusters_to_ORClusters.zonalSummary.number_of_pairs.mat),
                        row_order=rownames(Z5.OMP.ORClusters_to_ORClusters.zonalSummary.number_of_pairs.mat),
                        show_row_names = TRUE,
                        show_column_names = TRUE,
                        col = col_fun,
                        width = unit(7, "in"), height = unit(7, "in"),
                        border=TRUE
  )
  
  pdf(file=plot_output_filename,width=10,height=10,useDingbats=FALSE)
  draw(heatmap.plot)	
  dev.off()
  
  
}
