###this plot calculates beta-diversiy while excluding the original sample and only looks into clustering patterns between the subsamples at 0.25Gb,0.5Gb,0.75Gb,1Gb and 1.25Gb.
#Here the legend label order has also been rectified

####genus level####
library(vegan)
library(ggplot2)
bray_genus_t_NMDS <- read.csv("C:/Users/sneha/OneDrive/Desktop/final_draft/beta_without_raw/genus_without_raw/bray_genus_t_modified_final_without_raw.csv",header=TRUE)

#NMDS
#make community matrix - extract columns with abundance information
com_genus = bray_genus_t_NMDS[,4:ncol(bray_genus_t_NMDS)]

#turn abundance data frame into a matrix
m_com_genus = as.matrix(com_genus)
set.seed(123) #so that same data is reproduced everytime
nmds_genus = metaMDS(m_com_genus, distance = "bray")
nmds_genus
plot(nmds_genus)
#extract NMDS scores (x and y coordinates) for sites from newer versions of vegan package
data.scores.genus = as.data.frame(scores(nmds_genus)$sites)

#add columns to data frame 
data.scores.genus$sample = bray_genus_t_NMDS$sample
data.scores.genus$gigabase = bray_genus_t_NMDS$gigabase
data.scores.genus$group = bray_genus_t_NMDS$group
data.scores.genus$group = factor(data.scores.genus$group, levels = paste0("S", 1:10))

head(data.scores.genus)



xx_genus = ggplot(data.scores.genus, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes( shape = gigabase, colour = group))+ 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "group", y = "NMDS2", shape = "gigabase")  + 
  scale_colour_manual(values = c("red", "yellow","green","grey","pink","blue","darkgreen","lightblue","maroon","purple"))+
  scale_shape_manual(values=c(0,1,2,3,4)) #10 samples, so 10 colours (sample 1 to sample 10); 6 groups, so 5 shapes (0.25Gb,0.5Gb,0.75Gb,1Gb,1.25Gb) 5 shapes as we wont consider raw , we will see how they cluster to 1.25Gb as refernce
xx_genus
#xx_genus + coord_cartesian(xlim = c(0, 2), ylim = c(0, 2))
#to encircle the clustering
xx_genus + stat_ellipse(aes(group = group), geom = "polygon", fill = NA, color = "black")

ggsave("NMDS.svg")

###------done------###

####Species#####

library(vegan)
library(ggplot2)
bray_species_t_NMDS <- read.csv("C:/Users/sneha/OneDrive/Desktop/final_draft/beta_without_raw/species_without_raw/bray_species_t_modified_final_without_raw.csv",header=TRUE)

#NMDS
#make community matrix - extract columns with abundance information
com_species = bray_species_t_NMDS[,4:ncol(bray_species_t_NMDS)]

#turn abundance data frame into a matrix
m_com_species = as.matrix(com_species)
set.seed(123) #so that same data is reproduced everytime
nmds_species = metaMDS(m_com_species, distance = "bray")
nmds_species
plot(nmds_species)
#extract NMDS scores (x and y coordinates) for sites from newer versions of vegan package
data.scores.species = as.data.frame(scores(nmds_species)$sites)

#add columns to data frame 
data.scores.species$sample = bray_species_t_NMDS$sample
data.scores.species$gigabase = bray_species_t_NMDS$gigabase
data.scores.species$group = bray_species_t_NMDS$group
data.scores.species$group = factor(data.scores.species$group, levels = paste0("S", 1:10))

head(data.scores.species)



xx_species = ggplot(data.scores.species, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes( shape = gigabase, colour = group))+ 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "group", y = "NMDS2", shape = "gigabase")  + 
  scale_colour_manual(values = c("red", "yellow","green","grey","pink","blue","darkgreen","lightblue","maroon","purple"))+
  scale_shape_manual(values=c(0,1,2,3,4))
#10 samples, so 10 colours (sample 1 to sample 10); 6 groups, so 5 shapes (0.25Gb,0.5Gb,0.75Gb,1Gb,1.25Gb) we dont take the original/raw samples. Instead we cluster to 1.25Gb
xx_species
#to encircle the clustering
xx_species + stat_ellipse(aes(group = group), geom = "polygon", fill = NA, color = "black")

ggsave("NMDS.svg")
