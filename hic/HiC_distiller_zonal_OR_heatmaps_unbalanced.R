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
library("ComplexHeatmap")
library(circlize)

######################################################################
# SCRIPT TO GENERATE AGGREGATE HEATMAPS OF OR CLUSTER CONTACTS BY ZONE 
######################################################################

# used in Figure 3B, Figure 6B, and Supp Figure 6B
# input is dumped Hi-C contacts (aligned with distiller), bed files of OR clusters and OR genes annotated by zone


# resolution of analysis
bin_size <- 50000
dump_dir <- "/media/storageD/distiller/dump/"

## FUNCTIONS ##
format_trans_path <- function(name){
  path <- str_c(dump_dir,"trans.",name,".",bin_size,".txt.gz")
  return(path)
}
format_cis_path <- function(name){
  path <- str_c(dump_dir,"cis.",name,".",bin_size,".txt.gz")
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
  join.command <- str_c("echo \"join -1 1 -2 1 loci.sort.txt <(zcat ",dump_dir,"trans.",name,".",bin_size,".txt.gz) > loci.contacts.txt\" | bash")
  system(command=join.command,wait=TRUE)
  loci_data <- fread("loci.contacts.txt",sep=" ",header=FALSE,col.names=c("bait","prey","contacts"))
  system(command="rm loci.txt loci.sort.txt loci.contacts.txt", wait = TRUE)
  return(loci_data)
}


extract_contacts_Cis <- function(name,loci){
  write.table(dplyr::select(loci,bait),file="loci.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
  system(command="sort -k1,1 loci.txt > loci.sort.txt",wait=TRUE)
  join.command <- str_c("echo \"join -1 1 -2 1 loci.sort.txt <(zcat ",dump_dir,"cis.",name,".",bin_size,".txt.gz) > loci.contacts.txt\" | bash")
  system(command=join.command,wait=TRUE)
  loci_data <- fread("loci.contacts.txt",sep=" ",header=FALSE,col.names=c("bait","prey","contacts"))
  system(command="rm loci.txt loci.sort.txt loci.contacts.txt", wait = TRUE)
  return(loci_data)
}

# functions for working with contact data
process_data <- function(hic.data,hic.contacts){
  hic.counts_per_bait <- group_by(hic.data,bait) %>% summarize(total_contacts_per_bait = sum(contacts))
  hic.counts <- left_join(hic.data,hic.counts_per_bait) %>% mutate(fraction_contacts_per_bait = contacts/total_contacts_per_bait) %>% 
    mutate(norm_contacts = 1000000000*contacts/hic.contacts)
  return(hic.counts)
}
label_cis_trans <- function(combinations) {
  labeled.combinations <- combinations %>% mutate(bait_sp = bait) %>% mutate(prey_sp = prey) %>% separate(prey_sp,c("p.chr","p.coord"),"_",convert=TRUE) %>% separate(bait_sp,c("b.chr","b.coord"),"_",convert=TRUE) %>% mutate(type= ifelse(p.chr == b.chr,"cis","trans")) %>% 
    mutate(distance = ifelse(type == "cis",abs(p.coord - b.coord),"NA")) %>% mutate(distance = as.numeric(distance))
  return(labeled.combinations)
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

OR_Clusters_sel <- select(OR_Clusters,bait)
zonal_bins <-rbind(zone1_bins,zone2to3_bins,zone4to5_bins) %>% group_by(bait) %>% summarize(bin.zones = paste(type, collapse = ",")) %>% inner_join(x=OR_Clusters_sel) %>% 
  mutate(bin.zones = factor(bin.zones,levels=c("zone1","zone1,zone2to3","zone2to3","zone2to3,zone4to5","zone4to5","zone1,zone4to5","zone1,zone2to3,zone4to5")))

###### Analyze interactions between sets of loci
## define bin combinations
Islands.combinations <- CJ(bait = Islands$bait, prey=Islands$bait) %>% as.data.frame
OR_Clusters.combinations <- CJ(bait = OR_Clusters$bait, prey=OR_Clusters$bait) %>% as.data.frame %>% anti_join(Islands.combinations)

names_list = c(
  "Z5-NgnBright_replicates_merged",
  "Z1-NgnBright_replicates_merged",
  "Z5.NfiWT.merged",
  "Z5.NfiKO.merged",
  "Z5.OMP_merged",
  "Z1.OMP_merged",
  "P2_merged")
total_contacts_list = c(
  836047057,
  671386897,
  311783592,
  351710826,
  380353788,
  313312520,
  495587086)
plot_output_filenames = unlist(lapply(names_list, function (name) { return(str_c(name, ".new_anno.mean_contacts.hm.vmax10_filterORcluster.pdf")) }))

for (i in 1:length(names_list)) {
  
  cell_type_name <- names_list[i]
  cell_type_total_contacts <- total_contacts_list[i]
  plot_output_filename = plot_output_filenames[i]
  print(cell_type_name)
  print(cell_type_total_contacts)
  print("running ...")
  
  ### Libraries
  Z5.OMP.name <- cell_type_name
  Z5.OMP.total_hic_contacts <- cell_type_total_contacts
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
  Z5.OMP.data <- process_data(Z5.OMP.data,Z5.OMP.total_hic_contacts)
  
  ## Extract data
  Z5.OMP.ORClusters_to_ORClusters <- left_join(OR_Clusters.combinations,Z5.OMP.data) %>%  
    dplyr::select(bait,prey,norm_contacts,fraction_contacts_per_bait) %>% 
    replace_na(list(norm_contacts = 0,fraction_contacts_per_bait=0)) %>% 
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
  col_breaks = c(seq(0,10,length=99))
  col_fun = colorRamp2(col_breaks,my_palette)
  
  Z5.OMP.ORClusters_to_ORClusters.zonalSummary <- Z5.OMP.ORClusters_to_ORClusters.zonal %>% filter(type == "trans") %>% 
    filter(!is.na(prey.zones)) %>% filter(!is.na(bait.zones))%>% 
    group_by(bait.zones,prey.zones) %>% 
    summarize(number_of_pairs = n(),cumulative_contacts=sum(norm_contacts),mean_fraction_contacts=mean(fraction_contacts_per_bait),mean_contacts2 = mean(norm_contacts)) %>% 
    mutate(mean_contacts = cumulative_contacts /number_of_pairs) %>%  
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

