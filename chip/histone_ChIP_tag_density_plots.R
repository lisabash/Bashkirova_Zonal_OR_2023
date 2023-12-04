library(tidyverse)
library(ggplot2)
library(gridExtra)

#####################################
# boxplot  of ChIP signal (Figure 2B)
#####################################

# zonal OR color palette
zonal_fine_colors_fish_classIfirst<-c("#990000","#FF0000","#FFCC00","#33CC33","#0099FF","#0000FF")
zonal_fine_colors_dark_fish_classIfirst<-c("#660000","#990000","#999900","#003300","#0066CC","#000066")

# read in file counting reads over OR genes in ChIP data (generated with HOMER script "HOMER_histone_ChIP_tag_density.sh")
tagdensity<-read_tsv(file="zones_tagdensity.txt",col_names = c("tk79","tk9","zone"))
tagdensity

# boxplots

p.79<-ggplot(tagdensity, mapping = aes(zone,tk79, color=zone, fill = zone))+ylim(0,0.11)+
  geom_boxplot(alpha = 0.7,  lwd =1.2, outlier.shape = NA)+ 
  scale_color_manual(values=zonal_fine_colors_dark_fish_classIfirst) + 
  scale_fill_manual(values = zonal_fine_colors_fish_classIfirst) + #geom_jitter()+
  theme(legend.position="none",axis.text.y = element_text(colour= "black", size = 12),axis.text.x = element_text(colour= "black", size = 14), axis.title.x = element_blank(),axis.title.y = element_text(colour= "black", size = 14),plot.title = element_text(hjust = 0.5, size=20))+
  ggtitle("tk79")+ ylab("coverage/bp")
p.79

p.9<-ggplot(tagdensity, mapping = aes(zone,tk9, color=zone, fill = zone)) +ylim(0,0.06)+
  geom_boxplot(alpha = 0.7,  lwd =1.2, outlier.shape = NA)+ 
  scale_color_manual(values=zonal_fine_colors_dark_fish_classIfirst) + 
  scale_fill_manual(values = zonal_fine_colors_fish_classIfirst) + #geom_jitter()+
  theme(legend.position="none",axis.text.y = element_text(colour= "black", size = 12),axis.text.x = element_text(colour= "black", size = 14), axis.title.x = element_blank(),axis.title.y = element_text(colour= "black", size = 14),plot.title = element_text(hjust = 0.5, size=20))+
  ggtitle("tk9")+ ylab("coverage/bp")
p.9

grid.arrange(p.9, p.79,ncol = 2,nrow=1)
#g<-arrangeGrob(p.9, p.79,ncol = 2,nrow=1)
#ggsave(file="tagdensity_MOE_tk9_tk79.pdf",g, width = 8, height = 4, units = "in" ,useDingbats = FALSE)
