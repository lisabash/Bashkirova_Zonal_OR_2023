library(tidyverse)


# used to determine correlation between interchromosomal Hi-C contacts between OR genes and deposition of histone maks
# Supp Figure 2D
# input is dumped Hi-C contacts (50kb resolution) and file with quantification of ChIP deposition over OR cluster bins at 50kb resolution (using HOMER)
# each OR cluster bin is further annotated by the zone of expression of the resident OR genes
# output is matrix file (used as input to generate heatmap with python script)

# zone annotation
make_bait_str <- function(chrNo, bin) {
  return(sprintf("%d_%d", chrNo, bin))
  # return(str_c(sprintf(chrNo,scientific = FALSE),"_",sprintf("%d", bin)))
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

bedToBait_keepAll_fields <- function(fname) {
  y<-read_tsv(fname,col_names = c("chrom", "start", "stop","name"),comment="#") %>%
    mutate(chrNo = as.numeric(str_replace(chrom,"chr",""))) %>% 
    arrange(chrNo,start) %>% 
    mutate(bin = floor(start/bin_size)*bin_size) %>% 
    mutate(bait = make_bait_str(chrNo, bin)) %>% 
    dplyr::select(bait,name) %>% 
    distinct()
  return(y)
}

baitToBinFomat<-function(zonal_bins){ 
  file<-zonal_bins %>% separate(bait, into = c("chr", "start"), sep = "_", remove = FALSE) %>% 
    mutate("X"="chr") %>% unite("chrom",c(X,chr), sep = "", remove = FALSE)
  file$start <-as.integer(file$start)
  file<-data.frame("chr"=file$chrom, "start"= file$start, "stop"= file$start+50000,"name"=".","score"="0","score"="+","zone" =  file$type, "bait"=file$bait)
  return(file)
}

# resolution of analysis
#bin_size <- 25000
bin_size <- 50000
dump_dir <- "/media/storageD/distiller/dump/"

zone1_bins <- bedToBait("zone1.bed") %>% mutate(type="zone1")
zone2to3_bins <- bedToBait("zone2_3.bed") %>% mutate(type="zone2to3")
zone4to5_bins <- bedToBait("zone4_5.bed") %>% mutate(type="zone4to5")

zonal_bins<-rbind(zone1_bins,zone2to3_bins,zone4to5_bins)
zonal_bins_count<-zonal_bins%>% group_by(bait) %>% summarise(bin_count=n())
zonal_bins<-left_join(zonal_bins, zonal_bins_count) %>% 
  filter(bin_count == 1) %>% dplyr::select(-bin_count)

Islands_slop50kb <- bedToBait("Greek_Islands.slop50kb.mm10.merged.50kb.bed")%>% mutate(type="Islands_slop50kb")
zonal_bins.nonIsland  <- anti_join(zonal_bins,Islands_slop50kb,by="bait")

write_delim(x = zonal_bins.nonIsland, path = "zonal_baits_anno.txt", delim = "\t")
#now there is only 334

zone1_bins_bed<-baitToBinFomat(zone1_bins)
zone2to3_bins_bed<-baitToBinFomat(zone2to3_bins)
zone4to5_bins_bed <- baitToBinFomat(zone4to5_bins)

joined_zonal_bins_bed<-rbind(zone1_bins_bed,zone2to3_bins_bed,zone4to5_bins_bed) %>% 
  filter(bait %in% c(zonal_bins.nonIsland$bait))

write_delim(joined_zonal_bins_bed, "zonal_OR_50kb_bins.bed", delim="\t",col_names = FALSE)

######################
# EXIT R AND RUN HOMER 
#annotatePeaks.pl zonal_OR_50kb_bins.bed  mm10 -d \
#homer.mergedReps.nchip.OE_tk79 \
#homer.mergedReps2021.Z1_OMP_tk79 \
#homer.mergedReps2021.Z2_OMP_tk79 \
#homer.mergedReps2021.Z5_OMP_tk79 \
#homer.nchip.Z5_NgnBright_tk79_r2 \
#homer.mergedReps.nchip.Z5_NgnBright_tk79 > annotatePeaks_zonal_OR_50kb_bin_nchips_tagdesnsity_2021.txt 2> std_error_annotatePeaks_2021.txt

####################################
#### ORDER BAITS BY TK79 TAGDENSITY
####################################
# all ORs intermixed
read_tagdensity_file<-function(tagdensity_file){
  td<-read_delim(file=tagdensity_file,delim="\t",
                 skip = 1, col_names=FALSE) %>% 
    dplyr::select(X2,X3,X20,X21,X22,X23,X24,X25) %>% mutate("start" = X3-1)%>% dplyr::select(-X3)%>%
    unite("bait",c(X2,start), sep = "_") %>% separate(bait, into=c("drop","bait"), sep = "r") %>% 
    rename("OE_tk79"=X20,"Z1_OSN_tk79"=X21,"Z2_OSN_tk79"=X22,"Z5_OSN_tk79"=X23, "Z5_Ngn_tk79r1" = X24,"Z5_Ngn_tk79" = X25) %>% 
    arrange(desc(Z5_Ngn_tk79))
  return(td)
}

zonal_baits_tk79_tagdensity<-read_tagdensity_file("annotatePeaks_zonal_OR_50kb_bin_nchips_tagdesnsity_2021.txt")

Z5_Ngn_tk79_order_bait<-zonal_baits_tk79_tagdensity %>% 
  left_join(zonal_bins) %>% rename("zone"= type) %>%
  dplyr::select(bait, zone,Z5_Ngn_tk79) %>% 
  arrange(desc(Z5_Ngn_tk79))%>% 
  arrange(Z5_Ngn_tk79)
write_delim(Z5_Ngn_tk79_order_bait, path= "baits.Z5Ngntk79_order.zones_intermixed.tsv",delim="\t")


###########################
### read in hic data
###########################

library(data.table)
## FUNCTIONS ##
format_trans_path <- function(name){
  path <- str_c(dump_dir,"trans.",name,".",bin_size,".txt.gz")
  return(path)
}
format_cis_path <- function(name){
  path <- str_c(dump_dir,"cis.",name,".",bin_size,".txt.gz")
  return(path)
}
#### LOAD DATA #####
extract_contacts_Trans <- function(name,loci){
  write.table(dplyr::select(loci,bait),file="loci.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
  system(command="sort -k1,1 loci.txt > loci.sort.txt",wait=TRUE)
  join.command <- str_c("echo \"join -1 1 -2 1 loci.sort.txt <(zcat ",dump_dir,"trans.",name,".",bin_size,".txt.gz) > loci.contacts.txt\" | bash")
  system(command=join.command,wait=TRUE)
  loci_data <- fread("loci.contacts.txt",sep=" ",header=FALSE,col.names=c("bait","prey","contacts"))
  system(command="rm loci.txt loci.sort.txt loci.contacts.txt")
  return(loci_data)
}

extract_contacts_Cis <- function(name,loci){
  write.table(dplyr::select(loci,bait),file="loci.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
  system(command="sort -k1,1 loci.txt > loci.sort.txt",wait=TRUE)
  join.command <- str_c("echo \"join -1 1 -2 1 loci.sort.txt <(zcat ",dump_dir,"cis.",name,".",bin_size,".txt.gz) > loci.contacts.txt\" | bash")
  system(command=join.command,wait=TRUE)
  loci_data <- fread("loci.contacts.txt",sep=" ",header=FALSE,col.names=c("bait","prey","contacts"))
  system(command="rm loci.txt loci.sort.txt loci.contacts.txt")
  return(loci_data)
}

# functions for working with contact data
process_data <- function(hic.data,hic.contacts){
  hic.counts_per_bait <- group_by(hic.data,bait) %>% summarize(total_contacts_per_bait = sum(contacts))
  hic.counts <- left_join(hic.data,hic.counts_per_bait) %>% mutate(fraction_contacts_per_bait = contacts/total_contacts_per_bait) %>% mutate(norm_contacts = 1000000000*contacts/hic.contacts)
  return(hic.counts)
}
label_cis_trans <- function(combinations) {
  labeled.combinations <- combinations %>% mutate(bait_sp = bait) %>% mutate(prey_sp = prey) %>% separate(prey_sp,c("p.chr","p.coord"),"_",convert=TRUE) %>% separate(bait_sp,c("b.chr","b.coord"),"_",convert=TRUE) %>% mutate(type= ifelse(p.chr == b.chr,"cis","trans")) %>% mutate(distance = ifelse(type == "cis",abs(p.coord - b.coord),"NA")) %>% mutate(distance = as.numeric(distance))
  return(labeled.combinations)
}

#at 50kb
Islands <- bedToBait("Greek_Islands.nonX.50kb.mm10.bed")%>% mutate(type="Islands")
OR_Clusters_full <- bedToBait_keepAll_fields("ORClusters.ordered.no-chrX.mm10.50kb.mm10.bed") %>% mutate(type="OR_Clusters")
OR_Clusters<-OR_Clusters_full %>% dplyr::select(-name)
All_Bins <- bedToBait("mm10_assembled.50kb.bed")


Z5.Ngn.name <- "Z5-NgnBright_replicates_merged"
Z5.Ngn.total_hic_contacts <- 836047057 #total nodups
Z5.Ngn.trans_hic_path <-format_trans_path(Z5.Ngn.name)
Z5.Ngn.cis_hic_path <-format_cis_path(Z5.Ngn.name)


Selected_Bins <- rbind(Islands,OR_Clusters) %>% dplyr::select(bait,type) %>% group_by(bait) %>% summarize(anno = paste(type, collapse = ",")) %>% as.data.frame()

Z5.Ngn.data.trans <- extract_contacts_Trans(Z5.Ngn.name,Selected_Bins)
Z5.Ngn.data.cis <- extract_contacts_Cis(Z5.Ngn.name,Selected_Bins)
Z5.Ngn.data <- rbind(Z5.Ngn.data.trans,Z5.Ngn.data.cis) 
Z5.Ngn.data <- process_data(Z5.Ngn.data,Z5.Ngn.total_hic_contacts)

##################################################
###### Analyze interactions between sets of loci
## define bin combinations
######################################################

Zonal_OR.combinations <- CJ(bait = zonal_bins.nonIsland$bait, prey=zonal_bins.nonIsland$bait) %>% as.data.frame 
#remove_intra 
Zonal_OR.combinations.trans<-Zonal_OR.combinations %>% separate(bait, into = c("chr.1","start.1"),sep = "_", remove = FALSE) %>% 
  separate(prey, into = c("chr.2","start.2"),sep = "_", remove = FALSE) %>% 
  mutate("locus" = if_else(chr.1 == chr.2, "cis","trans")) %>% filter(locus == "trans")%>% dplyr::select(bait,prey)
########################
## Extract data
Z5.Ngn.Zonal_OR_to_Zonal_OR_trans <- left_join(Zonal_OR.combinations.trans,Z5.Ngn.data) %>%  dplyr::select(bait,prey,norm_contacts,fraction_contacts_per_bait) %>% replace_na(list(norm_contacts = 0,fraction_contacts_per_bait=0)) %>% mutate(anno = "Zonal_OR_to_Zonal_OR") %>% label_cis_trans() 

##########################
# Heatmaps norm contacts
#############################
Z5.Ngn.Zonal_OR_to_Zonal_OR.trans.mat <- Z5.Ngn.Zonal_OR_to_Zonal_OR_trans %>% arrange(b.chr,b.coord,p.chr,p.coord) %>% dplyr::select(bait,prey,norm_contacts) %>% spread(prey,norm_contacts) %>% column_to_rownames(var="bait") %>% as.matrix()
write.csv(Z5.Ngn.Zonal_OR_to_Zonal_OR.trans.mat, "Z5.Ngn.Zonal_OR_to_Zonal_OR_inter.mat.csv")



