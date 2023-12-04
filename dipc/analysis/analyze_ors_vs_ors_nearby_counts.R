library(tidyverse)
library(gridExtra)

# input is file with quantification of ORs within a specific distance in particle radii for every OR locus (produced by pairwise_distances.sh)

cells = read_tsv(file="../good_cells.txt", col_names=c("cell_name", "sex"))

######## LOAD IN ZONAL ANNOTATION
get_noGI_ORs<- function(fname){
  read_tsv(file = fname, col_names= FALSE) %>% 
    select(X4) %>% rename(name=X4)
}
zone1_noGI<-get_noGI_ORs("/media/storageE/lisa_dipc/bed_files/zone_1_noGIOverlap_50kb.bed")
zone23_noGI<-get_noGI_ORs("/media/storageE/lisa_dipc/bed_files/zone_23_noGIOverlap_50kb.bed")
zone45_noGI<-get_noGI_ORs("/media/storageE/lisa_dipc/bed_files/zone_45_noGIOverlap_50kb.bed")

####################################
# LOAD IN DATA AT 2.5 PARTICLE RADII
####################################

ors_to_ors_counts_list <- vector(mode = "list", length = nrow(cells))
for (row in 1:nrow(cells)) {
  cell_name <- cells[row, "cell_name"]$cell_name
  print(cell_name)
  sex <- cells[row, "sex"]
  counts_fname <- str_c("../compute-nearby-3d/counts/", cell_name, ".", sex, ".orsNoGI_vs_orsNoGI.dist2.5.counts.txt")
  counts_data <- read_tsv(counts_fname, col_names=TRUE)
  counts_data <- rename(counts_data, or_name = "#name")
  counts_data <- counts_data %>% mutate("cell_name" = cell_name)
  ors_to_ors_counts_list[[row]] <- counts_data
}
ors_to_ors_counts <- bind_rows(ors_to_ors_counts_list)
ors_to_ors_counts <- ors_to_ors_counts %>%
  separate(cell_name,c("zone","cell_num"), sep = "OMP.",remove = FALSE) %>%
  separate(or_name,c("chr","pos","name"), sep="_", remove= FALSE) %>% 
  select(-pos,-cell_num,-or_name)

#file format: 3d networks found interogating all ORs (excluding GI slop50) vs ORs from zonal categories 
#(excluding GI slop50kb and ORs from other zones slop50kb)

ors_to_ors_counts_2.5radii <-ors_to_ors_counts %>% mutate(radii="2.5radii")


####################################
# LOAD IN DATA AT 5 PARTICLE RADII
####################################

ors_to_ors_counts_list <- vector(mode = "list", length = nrow(cells))
for (row in 1:nrow(cells)) {
  cell_name <- cells[row, "cell_name"]$cell_name
  print(cell_name)
  sex <- cells[row, "sex"]
  counts_fname <- str_c("../compute-nearby-3d/counts/", cell_name, ".", sex, ".orsNoGI_vs_orsNoGI.dist5.counts.txt")
  counts_data <- read_tsv(counts_fname, col_names=TRUE)
  counts_data <- rename(counts_data, or_name = "#name")
  counts_data <- counts_data %>% mutate("cell_name" = cell_name)
  ors_to_ors_counts_list[[row]] <- counts_data
}
ors_to_ors_counts <- bind_rows(ors_to_ors_counts_list)
ors_to_ors_counts <- ors_to_ors_counts %>%
  separate(cell_name,c("zone","cell_num"), sep = "OMP.",remove = FALSE) %>%
  separate(or_name,c("chr","pos","name"), sep="_", remove= FALSE) %>% 
  select(-pos,-cell_num,-or_name)

#file format: 3d networks found interogating all ORs (excluding GI slop50) vs ORs from zonal categories 
#(excluding GI slop50kb and ORs from other zones slop50kb)

ors_to_ors_counts_5radii <-ors_to_ors_counts %>% mutate(radii="5radii")

####################################
# LOAD IN DATA AT 10 PARTICLE RADII
####################################

ors_to_ors_counts_list <- vector(mode = "list", length = nrow(cells))
for (row in 1:nrow(cells)) {
  cell_name <- cells[row, "cell_name"]$cell_name
  print(cell_name)
  sex <- cells[row, "sex"]
  counts_fname <- str_c("../compute-nearby-3d/counts/", cell_name, ".", sex, ".orsNoGI_vs_orsNoGI.dist10.counts.txt")
  counts_data <- read_tsv(counts_fname, col_names=TRUE)
  counts_data <- rename(counts_data, or_name = "#name")
  counts_data <- counts_data %>% mutate("cell_name" = cell_name)
  ors_to_ors_counts_list[[row]] <- counts_data
}
ors_to_ors_counts <- bind_rows(ors_to_ors_counts_list)
ors_to_ors_counts <- ors_to_ors_counts %>%
  separate(cell_name,c("zone","cell_num"), sep = "OMP.",remove = FALSE) %>%
  separate(or_name,c("chr","pos","name"), sep="_", remove= FALSE) %>% 
  select(-pos,-cell_num,-or_name)

#file format: 3d networks found interogating all ORs (excluding GI slop50) vs ORs from zonal categories 
#(excluding GI slop50kb and ORs from other zones slop50kb)

ors_to_ors_counts_10radii <-ors_to_ors_counts %>% mutate(radii="10radii")


####################################
# LOAD IN DATA AT 20 PARTICLE RADII
####################################

ors_to_ors_counts_list <- vector(mode = "list", length = nrow(cells))
for (row in 1:nrow(cells)) {
  cell_name <- cells[row, "cell_name"]$cell_name
  print(cell_name)
  sex <- cells[row, "sex"]
  counts_fname <- str_c("../compute-nearby-3d/counts/", cell_name, ".", sex, ".orsNoGI_vs_orsNoGI.dist20.counts.txt")
  counts_data <- read_tsv(counts_fname, col_names=TRUE)
  counts_data <- rename(counts_data, or_name = "#name")
  counts_data <- counts_data %>% mutate("cell_name" = cell_name)
  ors_to_ors_counts_list[[row]] <- counts_data
}
ors_to_ors_counts <- bind_rows(ors_to_ors_counts_list)
ors_to_ors_counts <- ors_to_ors_counts %>%
  separate(cell_name,c("zone","cell_num"), sep = "OMP.",remove = FALSE) %>%
  separate(or_name,c("chr","pos","name"), sep="_", remove= FALSE) %>% 
  select(-pos,-cell_num,-or_name)

#file format: 3d networks found interogating all ORs (excluding GI slop50) vs ORs from zonal categories 
#(excluding GI slop50kb and ORs from other zones slop50kb)

ors_to_ors_counts_20radii <-ors_to_ors_counts %>% mutate(radii="20radii")


all.radii<-rbind(ors_to_ors_counts_2.5radii,ors_to_ors_counts_5radii,ors_to_ors_counts_10radii,ors_to_ors_counts_20radii)

#add zonal information
all.radii_zone<-all.radii %>% mutate(OR_zone = if_else(name %in% zone1_noGI$name, "zone1_OR", 
                                               if_else(name %in% zone23_noGI$name, "zone23_OR", 
                                               if_else(name %in% zone45_noGI$name, "zone45_OR", "unknown"))))

#########################
# PLOTS (Supp Figure 3C)
#########################

#largest aggregate
working_fract<-all.radii_zone%>% group_by(cell_name, zone, radii) %>% 
  summarise(max_inter=max(n_other_inter), median_inter = median(n_other_inter), mean_inter = mean(n_other_inter))
working_fract$radii<-factor(working_fract$radii, levels=c("2.5radii","5radii","10radii","20radii"))
p<-ggplot(working_fract, aes(zone,max_inter, color= zone)) + 
  geom_boxplot(outlier.shape = NA ) +
  scale_color_manual(values= c("red4", "blue4")) + 
  ggtitle("max inter")+
  facet_grid(cols = vars(radii))
p
#ggsave(file="boxplot.max_inter_all_radii.pdf",p, width = 6, height = 5, units = "in" )

#largest number of chromosomes in aggregate
working_fract<-all.radii%>% group_by(cell_name, zone, radii) %>% 
  summarise(max_hom=max(n_other_hom))
working_fract$radii<-factor(working_fract$radii, levels=c("2.5radii","5radii","10radii","20radii"))
p<-ggplot(filter(working_fract, radii != "20radii"), aes(zone,max_hom, color= zone)) + 
  geom_boxplot(outlier.shape = NA ) +
  #geom_point(position=position_jitterdodge(jitter.width = 0.3, jitter.height = 0), alpha =0.5)+
  #ylim(0,7) + 
  #scale_fill_manual(values= c("red3", "green3","blue3"))+
  scale_color_manual(values= c("red4", "blue4")) + 
  ggtitle("max inter")+
  facet_grid(cols = vars(radii))
p
#ggsave(file="boxplot.max_hom_all_radii.pdf",p, width = 6, height = 5, units = "in" )

##################################
# PLOTS ZONAL OR (Supp Figure 3D)
##################################
working_fract_zone<-all.radii_zone%>% group_by(cell_name, zone, radii, OR_zone) %>% 
  summarise(max_inter=max(n_other_inter), median_inter = median(n_other_inter), mean_inter = mean(n_other_inter), fract_inter=sum(n_other_inter>0)/n(),total=n())
working_fract_zone$radii<-factor(working_fract_zone$radii, levels=c("2.5radii","5radii","10radii","20radii"))
working_fract_zone$OR_zone<-factor(working_fract_zone$OR_zone, levels=c("zone1_OR","zone23_OR","zone45_OR","unknown"))


p_fract<-ggplot(filter(working_fract_zone, radii == "2.5radii" & OR_zone != "unknown"), aes(OR_zone,fract_inter, color= OR_zone)) + 
  geom_boxplot(outlier.shape = NA ) +
  #geom_point(position=position_jitterdodge(jitter.width = 0.3, jitter.height = 0), alpha =0.5)+
  scale_color_manual(values= c("red4", "green4","blue4")) + 
  ggtitle("fract ORs with inter contacts per cell")+
  facet_grid(cols = vars(zone))
p_fract
#ggsave(file="boxplot.230716.2.5radii.fract_nonzero_by_zone.pdf",p_fract, width = 6, height = 4, units = "in", useDingbats=FALSE )

p_fract<-ggplot(filter(working_fract_zone, radii == "2.5radii" & OR_zone != "unknown" & zone == "z5"), aes(OR_zone,fract_inter, color= OR_zone)) + 
  geom_boxplot(outlier.shape = NA ) +
  #geom_point(position=position_jitterdodge(jitter.width = 0.3, jitter.height = 0), alpha =0.5)+
  scale_color_manual(values= c("red4", "green4","blue4")) + 
  ggtitle("fract ORs with inter contacts per cell")+
  ylim(0,1)+
  facet_grid(cols = vars(zone))
p_fract
#ggsave(file="boxplot.230716.2.5radii.fract_nonzero_by_zone_z5_cells_only.pdf",p_fract, width = 3.5, height = 4, units = "in", useDingbats=FALSE )


df<-all.radii_zone %>% 
  filter(radii == "2.5radii"& zone == "z5" ) %>% 
  group_by(cell_name, OR_zone) %>% 
  summarise(mean_inter=mean(n_other_inter),max_inter=max(n_other_inter),median_inter = median(n_other_inter),fract_inter=sum(n_other_inter>0)/n(),total=n() )

df_z1_z2<-filter(df, OR_zone %in% c("zone1_OR", "zone23_OR"))
wilcox.test(fract_inter ~ OR_zone, df_z1_z2)
#p-value = 0.2829

df_z1_z5<-filter(df, OR_zone %in% c("zone1_OR", "zone45_OR"))
wilcox.test(fract_inter ~ OR_zone, df_z1_z5)
#p-value = 0.0001646

df_z2_z5<-filter(df, OR_zone %in% c("zone23_OR", "zone45_OR"))
wilcox.test(fract_inter ~ OR_zone, df_z2_z5)
#p-value = 0.003812

