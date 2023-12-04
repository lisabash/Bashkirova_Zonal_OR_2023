# load the zonal annotation from Tan 2018 "A Near-Complete Spatial Map of Olfactory Receptors in the Mouse Main Olfactory Epithelium"
# convert continuous zonal annotation into integers

Tan2018_tb<-read_tsv("/Users/elizaveta/Desktop/Zonal_anno/Tan_2018_zonal_anno.txt") %>% rename("name"= 'gene name', "strand" = 'gene strand', "locus" = 'gene locus', "classI" = 'is Class I OR', "zone" = 'inferred zone index')
Tan2018_tb_full<-Tan2018_tb %>% mutate(new_zone = if_else(zone == "low expression","low expression",
                                                  if_else(zone == "unusual", "unusual",
                                                  if_else(zone < 1.5, "1",
                                                  if_else(zone < 2.5, "2",
                                                  if_else(zone < 3.5, "3",
                                                  if_else(zone < 4.5, "4", "5"))))))) %>% arrange(zone)

Tan2018_tb_full %>% group_by(new_zone) %>% summarise(count=n())
# A tibble: 7 x 2
#new_zone       count
#<chr>          <int>
#1 1                373
#2 2                286
#3 3                164
#4 4                144
#5 5                 44
#6 low expression   384
#7 unusual           22

Tan2018_tb_NA<-Tan2018_tb %>% mutate(new_zone = if_else(zone < 1.5, 1,
                                                if_else(zone < 2.5, 2,
                                                if_else(zone < 3.5, 3,
                                                if_else(zone < 4.5, 4, 
                                                if_else(zone <= 5, 5, NULL )))))) %>% 
  
                              mutate(new_zone_coarse = if_else(zone < 1.5, "zone1",
                                                if_else(zone < 3.5, "zone23",
                                                if_else(zone <= 5, "zone45", NULL )))) %>% arrange(zone)

Tan2018_tb_NA %>% group_by(new_zone) %>% summarise(count=n())
# A tibble: 6 x 2
#new_zone count
#<dbl> <int>
#1        1   373
#2        2   286
#3        3   164
#4        4   144
#5        5    44
#6       NA   406

Tan2018_tb_NA %>% group_by(new_zone_coarse) %>% summarise(count=n())
# A tibble: 4 x 2
#new_zone_coarse count
#<chr>           <int>
#1 zone1             373
#2 zone23            450
#3 zone5             188
#4 NA                406

Tan2018_tb_NA %>% filter(is.na(new_zone) == FALSE) %>% summarise(count=n())
# A tibble: 1 x 1
#count
#<int>
#1  1011



Tan2018_tb_NA %>% filter(is.na(new_zone) == FALSE) %>% group_by(classI) %>% summarise(count=n())
# A tibble: 2 x 2
#classI count
#<chr>  <int>
#1 no       896
#2 yes      115

Tan2018_tb_NA %>% filter(is.na(new_zone) == FALSE) %>% filter(classI == "no") %>% group_by(new_zone) %>% summarise(count=n())
# A tibble: 5 x 2
#new_zone count
#<dbl> <int>
#1        1   261
#2        2   283
#3        3   164
#4        4   144
#5        5    44

Tan2018_tb_NA %>% filter(is.na(new_zone) == FALSE) %>% filter(classI == "no") %>% group_by(new_zone_coarse) %>% summarise(count=n())
# A tibble: 3 x 2
#new_zone_coarse count
#<chr>           <int>
#1 zone1             261
#2 zone23            447
#3 zone5             188

Tan2018_tb_NA_classI_fine<-Tan2018_tb_NA %>% mutate(new_zone = ifelse(classI == "yes", "classI", new_zone))
Tan2018_tb_NA_classI_fine %>% group_by(new_zone) %>% summarise(count=n())
# A tibble: 7 x 2
#new_zone count
#<chr>    <int>
#1 1          261
#2 2          283
#3 3          164
#4 4          144
#5 5           44
#6 classI     159
#7 NA         362

Tan2018_tb_NA_classI_coarse<- Tan2018_tb_NA %>% mutate(new_zone_coarse = ifelse(classI == "yes", "classI",new_zone_coarse)) 
Tan2018_tb_NA_classI_coarse %>% group_by(new_zone_coarse) %>% summarise(count=n())
# A tibble: 5 x 2
#new_zone_coarse count
#<chr>           <int>
#1 classI            159
#2 zone1             261
#3 zone23            447
#4 zone5             188
#5 NA                362


# write zonal annotation file #
# will use this zonal assignment for ORs
# but will use a different locus (chr, start, stop) for OR genes (determined by Soria) to have a better mapping of OR promoters

Tan2018_tb_NA_classI_coarse<-Tan2018_tb_NA_classI_coarse  %>%  separate(locus,into = c("chr","start","stop"), sep = ":|-")%>% select(name, chr, start, stop, strand,new_zone_coarse)

write_tsv(filter(Tan2018_tb_NA_classI_coarse,complete.cases(Tan2018_tb_NA_classI_coarse)),"zonal_anno_2020_coarse.tsv",col_names=FALSE)

Tan2018_tb_NA_classI_fine<-Tan2018_tb_NA_classI_fine  %>%  separate(locus,into = c("chr","start","stop"), sep = ":|-")%>% select(name, chr, start, stop, strand,new_zone)

write_tsv(filter(Tan2018_tb_NA_classI_fine,complete.cases(Tan2018_tb_NA_classI_fine)),"zonal_anno_2020_fine.tsv",col_names=FALSE)


