library(tidyverse)
library(ggplot2)
library(gridExtra) 

OR_OR.con<-read_delim(file="/media/storageE/lisa_dipc/contact-densities/ors_ors_unknown_allele_contacts_densities.txt", delim = "\t", col_names = c("cell","sex","OR_OR"))

density_table<-OR_OR.con %>% 
  separate(cell,into = c("zone","num"), sep="OMP", remove=FALSE) %>% 
  select(-num)

density_table_zone1OSN<-density_table %>% filter(zone == "z1")
density_table_zone5OSN<-density_table %>% filter(zone == "z5")

density_table_gather<-density_table %>% gather("contact_type","density",-cell,-zone,-sex)

#Plot
p.all<-ggplot(filter(density_table_gather,contact_type %in% c("OR_OR")), aes(zone,density, color =factor(zone)))+ 
  geom_boxplot(outlier.shape = NA)+ 
  #geom_point(position = position_jitterdodge(jitter.height = 0), alpha = 0.5)+
  scale_color_manual(values = c("red4", "blue4"))+
  ggtitle("zonal contact density")+ xlab("cell zone")+ facet_grid(cols = vars(factor(contact_type, levels = c("OR_OR","GI_GI"))))
p.all
ggsave("boxplot.contact_density_OR_OR.pdf",p.all,width = 4, height = 4.5 , useDingbats = FALSE)



