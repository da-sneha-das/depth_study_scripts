#code to generate taxonomic analysis plots at species level

###--4th draft ---###
###SPECIES LEVEL--###
##TAXONOMY
###----barcodes and sample---
#sample1-s1-000251849  
#sample2-s2-000246807
#sample3-s3-000251080
#sample4-s4-000251757 
#sample5-s5-000251074
#sample6-s6-000246255
#sample7-s7-000246259
#sample8-s8-000371236
#sample9-s9-000371290
#sample10-s10-000360675


#NO of observed genus- wilcox test
##Wilcox test- to compare each depth against the reference depth of 1.25Gb (highest for us)
library(ggplot2)
library(ggpubr)
library(rstatix)
##reading the input file
obs_species <- read.table("C:/Users/sneha/OneDrive/Desktop/final_draft/species/species_count_wilcox.csv", sep = ",", header = T, check.names = F, row.names = 1, comment.char = "", stringsAsFactors = F)
test_t <- t(obs_species)
test_melt <- reshape2::melt(test_t, value.name = "Observed_species")
names(test_melt)[1] <- "Sequencing_Depth"
names(test_melt)[2] <- "Sample"
test_melt$Sequencing_Depth <- factor(test_melt$Sequencing_Depth, levels = c("0.25Gb", "0.5Gb", "0.75Gb", "1Gb", "1.25Gb"))
# Generate boxplot
boxplot <- ggplot(data = test_melt, mapping=
                    aes(x = Sequencing_Depth, y = Observed_species,colour=Sequencing_Depth)) +
  geom_boxplot()+ xlab("Sequencing Depth") + ylab("Number of observed species")+
  geom_point(position = position_jitterdodge(jitter.width = 0.7))
#boxplot
#wilcox test- comparing each group with reference group
stat.test <- stat.test <- test_melt %>%
  wilcox_test(Observed_species ~ Sequencing_Depth, ref.group = "1.25Gb", paired = FALSE) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj") %>%
  add_xy_position()
stat.test

#adding p values to box plot
p <- boxplot + stat_pvalue_manual(
  stat.test, label = "p.adj") # Change label to p.adj.signif is asterisks are preferred
p+ theme_classic()
#-----x-----##

#Kruskal wallis test- shannon, simpson, pielou
library(ggplot2)
library(ggpubr)
library(rstatix)
df_8 <- read.csv("C:/Users/sneha/OneDrive/Desktop/final_draft/species/alpha_diversity_kruskal_species.csv",header=TRUE)
df_8

#shannon_diversity of all samples at all subsampling levels- violin plot
box_plot_shannon_s <- ggplot(data=df_8, mapping=aes(x=Sequencing_Depth,y=Shannon_index,colour=Sequencing_Depth))+ geom_boxplot()+ labs(x="Sequencing Depth",y="Shannon Index")+
  geom_point(position = position_jitterdodge(jitter.width = 0.3))
box_plot_shannon_s+ theme_classic() #no p value in plot
box_plot_shannon_s + stat_compare_means(method = "kruskal.test",label = "p.format") #p value in plot
kruskal.test(Shannon_index ~ Sequencing_Depth, data = df_8)#only prints stats


#simpson_diversity of all samples at all subsampling levels- violin plot
box_plot_simpson_s <- ggplot(data=df_8, mapping=aes(x=Sequencing_Depth,y=Simpson_index,color=Sequencing_Depth))+ geom_boxplot()+labs(x="Sequencing Depth",y="Simpson Index")+
  geom_point(position = position_jitterdodge(jitter.width = 0.3))
box_plot_simpson_s+ theme_classic() #no p value in plot
box_plot_simpson_s + stat_compare_means(method = "kruskal.test",label = "p.format") #p value in plot
kruskal.test(Simpson_index ~ Sequencing_Depth, data = df_8)#only stats

#Pielou_index of all samples at all subsampling levels- violin plot
Box_plot_pielou_s <- ggplot(data=df_8, mapping=aes(x=Sequencing_Depth,y=Pielou_index,colour=Sequencing_Depth))+ geom_boxplot()+labs(x="Sequencing Depth",y="Pielou's Index")+
  geom_point(position = position_jitterdodge(jitter.width = 0.3))
Box_plot_pielou_s+ theme_classic() #no p value in plot
Box_plot_pielou_s + stat_compare_means(method = "kruskal.test",label = "p.format")#p value in plot
kruskal.test(Pielou_index ~ Sequencing_Depth, data = df_8)#only stats


#Distribution of no of genus
box_genus_count_s <- ggplot(data=df_8,mapping=aes(x=Sequencing_Depth,y=Species_count,colour=Sequencing_Depth))+ geom_boxplot()+labs(x="Sequencing Depth",y="Observed Species")
box_genus_count_s #without p value in plot
box_genus_count_s + stat_compare_means(method = "kruskal.test",label = "p.format") #with p value in plot
kruskal.test(Species_count ~ Sequencing_Depth, data = df_8) #just the test, no plot

#Alternate method to do just statistics-kruskal-walis
stat.test <- compare_means(Species_count~ Sequencing_Depth, method= "kruskal.test",data= df_8,p.adjust.method="none") #no plot- just stats
stat.test #displaying the stats


#####-------x----------#####

##-----Post-hoc analysis test---for pairwise comparison--##
#Post-hoc analysis test- DunnTest using Bonferroni for P value correction for Kruskal Wallis result
#Note that- post-hoc analysis test is done, only when we are comparing Kruskal Walis between more than 2 groups examle group A, B and C, and then we get a signficant P value (that is P value less than 0.05)
#When we get p-value which is signficant, then only we do post0hoc test to find out that the difference is exactly between which groups ie is it between A and B or A and C or B and C or are all the 3 groups different from each other
#If the p-value is not significant, that means theres no difference between the groups, hence it makes no sense to do post hoc test
#If we are comparing only 2 groups, only kruskal walis test will suffice and post hoc is notr required bcz we know here the difference is between the 2 groups only. Only if we had more than 2 groups, then the post hoc analysis would help us to find the difference is between which groups
#p-adjustment is done for correction of multiple testing
##------our data post-hoc test----##
#We are doing the post-hoc test only on the no of genus with respect to depth because for no of species we got significant p-value in kruskal walis test.
#For our shannon, simpson and pielou index, the kruskal wallis p-value was not significant,so it meant there's no difference between the sequencing depths, hence we dont do a post-hoc analysis

#library(FSA) #no plot---just the dunn's values for supplementary
#dunn_Test_genus <- dunnTest(Genus_count~Sequencing_Depth,data=df_8,method="none")
#dunn_Test_genus

#Alternate-method-dunn_test---final plot and dunn's test
library(rstatix) ###consider this for draft
obs_dunn_species <- dunn_test(data = df_8, formula = Species_count~Sequencing_Depth,p.adjust.method="none") #this for draft
obs_dunn_species #this for draft
obs_dunn_species <- obs_dunn_species %>% add_xy_position(x = "Sequencing_Depth")
obs_dunn_species

#plot
dunn_plot_no_of_species <- ggboxplot(df_8,y="Species_count",x="Sequencing_Depth",fill="Sequencing_Depth")+ stat_pvalue_manual(obs_dunn_species,hide.ns=FALSE,label="p = {p.adj}")+ labs(x="Sequencing depth",y="Observed Species") #plotting the dunn test p values with the no of genus 
dunn_plot_no_of_species #plot for dunn's test p value comparison between various depths

#####---Line chart------####

#shannon dataframe
df_1 <- read.csv("C:/Users/sneha/OneDrive/Desktop/final_draft/species/shannon_all_species_final.csv",header=TRUE)
#df_1

#simpson dataframe
df_2 <- read.csv("C:/Users/sneha/OneDrive/Desktop/final_draft/species/simpson_all_species_final.csv",header=TRUE)
#df_2

#pielou dataframe
df_3 <- read.csv("C:/Users/sneha/OneDrive/Desktop/final_draft/species/pielou_all_species_final.csv",header=TRUE)
#df_3

#shannon plot
df_1_shannon <-  ggplot(df_1, aes(x=bases_Gb))+ geom_line(aes(y=S1,group=1,colour="S1"))+geom_point(aes(y=S1,group=1,color="pink"))+geom_line(aes(y=S2,group=1,color="S2"))+geom_point(aes(y=S2 ,group=1,color="pink"))+
  geom_line(aes(y=S3,group=1,colour="S3"))+geom_point(aes(y=S3,group=1,color="pink"))+geom_line(aes(y=S4,group=1,colour="S4"))+geom_point(aes(y=S4,group=1,color="pink"))+ geom_line(aes(y=S5,group=1,colour="S5"))+geom_point(aes(y=S5,group=1,color="pink"))+
  geom_line(aes(y=S6,group=1,colour="S6"))+geom_point(aes(y=S6,group=1,color="pink"))+geom_line(aes(y=S7,group=1,colour="S7"))+geom_point(aes(y=S7,group=1,color="pink"))+
  geom_line(aes(y=S8,group=1,colour="S8"))+geom_point(aes(y=S8,group=1,color="pink"))+
  geom_line(aes(y=S9,group=1,colour="S9"))+geom_point(aes(y=S9,group=1,color="pink"))+
  geom_line(aes(y=S10,group=1,colour="S10"))+geom_point(aes(y=S10,group=1,color="pink"))+labs(x="Sequencing Depth (Gb)",y="Shannon Index")
df_1_shannon

#simpson plot
df_2_simpson <-  ggplot(df_2, aes(x=bases_Gb))+ geom_line(aes(y=S1,group=1,colour="S1"))+geom_point(aes(y=S1,group=1,color="pink"))+geom_line(aes(y=S2,group=1,color="S2"))+geom_point(aes(y=S2 ,group=1,color="pink"))+
  geom_line(aes(y=S3,group=1,colour="S3"))+geom_point(aes(y=S3,group=1,color="pink"))+geom_line(aes(y=S4,group=1,colour="S4"))+geom_point(aes(y=S4,group=1,color="pink"))+ geom_line(aes(y=S5,group=1,colour="S5"))+geom_point(aes(y=S5,group=1,color="pink"))+
  geom_line(aes(y=S6,group=1,colour="S6"))+geom_point(aes(y=S6,group=1,color="pink"))+geom_line(aes(y=S7,group=1,colour="S7"))+geom_point(aes(y=S7,group=1,color="pink"))+
  geom_line(aes(y=S8,group=1,colour="S8"))+geom_point(aes(y=S8,group=1,color="pink"))+
  geom_line(aes(y=S9,group=1,colour="S9"))+geom_point(aes(y=S9,group=1,color="pink"))+
  geom_line(aes(y=S10,group=1,colour="S10"))+geom_point(aes(y=S10,group=1,color="pink"))+labs(x="Sequencing Depth (Gb)",y="Simpson Index")
df_2_simpson


#pielou index
df_3_pielou <-  ggplot(df_3, aes(x=bases_Gb))+ geom_line(aes(y=S1,group=1,colour="S1"))+geom_point(aes(y=S1,group=1,color="pink"))+geom_line(aes(y=S2,group=1,color="S2"))+geom_point(aes(y=S2 ,group=1,color="pink"))+
  geom_line(aes(y=S3,group=1,colour="S3"))+geom_point(aes(y=S3,group=1,color="pink"))+geom_line(aes(y=S4,group=1,colour="S4"))+geom_point(aes(y=S4,group=1,color="pink"))+ geom_line(aes(y=S5,group=1,colour="S5"))+geom_point(aes(y=S5,group=1,color="pink"))+
  geom_line(aes(y=S6,group=1,colour="S6"))+geom_point(aes(y=S6,group=1,color="pink"))+geom_line(aes(y=S7,group=1,colour="S7"))+geom_point(aes(y=S7,group=1,color="pink"))+
  geom_line(aes(y=S8,group=1,colour="S8"))+geom_point(aes(y=S8,group=1,color="pink"))+
  geom_line(aes(y=S9,group=1,colour="S9"))+geom_point(aes(y=S9,group=1,color="pink"))+
  geom_line(aes(y=S10,group=1,colour="S10"))+geom_point(aes(y=S10,group=1,color="pink"))+labs(x="Sequencing Depth (Gb)",y="Pielou Index")
df_3_pielou

#####-correlation-####
#SPEARMAN

#correlation between bases and shannon diversity
df_5 <- read.csv("C:/Users/sneha/OneDrive/Desktop/final_draft/species/shannon_correlation_species_final.csv",header=TRUE)
df_5
#correlation between bases and simpson diversity
df_6 <- read.csv("C:/Users/sneha/OneDrive/Desktop/final_draft/species/simpson_correlation_species_final.csv",header=TRUE)
df_6
#correlation between bases and pileou index
df_7 <- read.csv("C:/Users/sneha/OneDrive/Desktop/final_draft/species/pielou_correlation_species_final.csv",header=TRUE)
df_7

##----------sample1-----------##

#correlation- sample1-shannonn
corr_shannon_sample1_s <- cor.test(df_5$bases_Gb,df_5$S1,method="spearman")
corr_shannon_sample1_s
#visualisation using scatter plot-
plot_sample1_shannon_s <- ggscatter(df_5, x = "bases_Gb", y = "S1", 
                                    add = "reg.line", conf.int = TRUE, 
                                    cor.coef = TRUE, cor.method = "spearman",
                                    xlab = "Sequencing Depth (Gb)", ylab = "Shannon Index",color="brown",point=TRUE)
plot_sample1_shannon_s
#correlation- sample1- simpson
corr_simpson_sample1_s <-  cor.test(df_6$bases_Gb,df_6$S1,method="spearman")
corr_simpson_sample1_s 
#visualisation- scatter plot
plot_sample1_simpson_s <-ggscatter(df_6, x = "bases_Gb", y = "S1", 
                                   add = "reg.line", conf.int = TRUE, color="blue", 
                                   cor.coef = TRUE, cor.method = "spearman",
                                   xlab = "Sequencing Depth (Gb)", ylab = "Simpson Index")
plot_sample1_simpson_s

#correlation- sample1- Pielou
corr_pielou_sample1_s <-  cor.test(df_7$bases_Gb,df_7$S1,method="spearman")
corr_pielou_sample1_s
#visualisation- scatter plot
plot_sample1_pielou_s <-ggscatter(df_7, x = "bases_Gb", y = "S1", 
                                  add = "reg.line", conf.int = TRUE, color="darkgreen",
                                  cor.coef = TRUE, cor.method = "spearman",pvalue.thresholds = c(0.01, 0.05, 0.1),
                                  xlab = "Sequencing Depth (Gb)", ylab = "Pielou Index")
plot_sample1_pielou_s

####---sample2------####

#correlation- sample2-shannonn
corr_shannon_sample2_s <- cor.test(df_5$bases_Gb,df_5$S2,method="spearman")
corr_shannon_sample2_s
#visualisation using scatter plot-
plot_sample2_shannon_s <- ggscatter(df_5, x = "bases_Gb", y = "S2", 
                                    add = "reg.line", conf.int = TRUE, color="brown",
                                    cor.coef = TRUE, cor.method = "spearman",
                                    xlab = "Sequencing Depth (Gb)", ylab = "Shannon Index")
plot_sample2_shannon_s
#correlation- sample2- simpson
corr_simpson_sample2_s <-  cor.test(df_6$bases_Gb,df_6$S2,method="spearman")
corr_simpson_sample2_s
#visualisation- scatter plot
plot_sample2_simpson_s <-ggscatter(df_6, x = "bases_Gb", y = "S2", 
                                   add = "reg.line", conf.int = TRUE, color="blue",
                                   cor.coef = TRUE, cor.method = "spearman",
                                   xlab = "Sequencing Depth (Gb)", ylab = "Simpson Index")
plot_sample2_simpson_s

#correlation- sample2- Pielou
corr_pielou_sample2_s <-  cor.test(df_7$bases_Gb,df_7$S2,method="spearman")
corr_pielou_sample2_s
#visualisation- scatter plot
plot_sample2_pielou_s <-ggscatter(df_7, x = "bases_Gb", y = "S2", 
                                  add = "reg.line", conf.int = TRUE,color="darkgreen", 
                                  cor.coef = TRUE, cor.method = "spearman",
                                  xlab = "Sequencing Depth (Gb)", ylab = "Pielou Index")
plot_sample2_pielou_s

######---sample3----------######

#correlation- sample3-shannonn
corr_shannon_sample3_s <- cor.test(df_5$bases_Gb,df_5$S3,method="spearman")
corr_shannon_sample3_s
#visualisation using scatter plot-
plot_sample3_shannon_s <- ggscatter(df_5, x = "bases_Gb", y = "S3", 
                                    add = "reg.line", conf.int = TRUE, color="brown",
                                    cor.coef = TRUE, cor.method = "spearman",
                                    xlab = "Sequencing Depth (Gb)", ylab = "Shannon Index")
plot_sample3_shannon_s
#correlation- sample3-simpson
corr_simpson_sample3_s <-  cor.test(df_6$bases_Gb,df_6$S3,method="spearman")
corr_simpson_sample3_s
#visualisation- scatter plot
plot_sample3_simpson_s <-ggscatter(df_6, x = "bases_Gb", y = "S3", 
                                   add = "reg.line", conf.int = TRUE, color="blue",
                                   cor.coef = TRUE, cor.method = "spearman",
                                   xlab = "Sequencing Depth (Gb)", ylab = "Simpson Index")
plot_sample3_simpson_s

#correlation- sample3- Pielou
corr_pielou_sample3_s <-  cor.test(df_7$bases_Gb,df_7$S3,method="spearman")
corr_pielou_sample3_s
#visualisation- scatter plot
plot_sample3_pielou_s <-ggscatter(df_7, x = "bases_Gb", y = "S3", 
                                  add = "reg.line", conf.int = TRUE,color="darkgreen", 
                                  cor.coef = TRUE, cor.method = "spearman",
                                  xlab = "Sequencing Depth (Gb)", ylab = "Pielou Index")
plot_sample3_pielou_s

########----sample4------------####

#correlation- sample4-shannonn
corr_shannon_sample4_s <- cor.test(df_5$bases_Gb,df_5$S4,method="spearman")
corr_shannon_sample4_s
#visualisation using scatter plot-
plot_sample4_shannon_s <- ggscatter(df_5, x = "bases_Gb", y = "S4", 
                                    add = "reg.line", conf.int = TRUE,color="brown",
                                    cor.coef = TRUE, cor.method = "spearman",
                                    xlab = "Sequencing Depth (Gb)", ylab = "Shannon Index")
plot_sample4_shannon_s
#correlation- sample4- simpson
corr_simpson_sample4_s <-  cor.test(df_6$bases_Gb,df_6$S4,method="spearman")
corr_simpson_sample4_s
#visualisation- scatter plot
plot_sample4_simpson_s <-ggscatter(df_6, x = "bases_Gb", y = "S4", 
                                   add = "reg.line", conf.int = TRUE, color="blue",
                                   cor.coef = TRUE, cor.method = "spearman",
                                   xlab = "Sequencing Depth (Gb)", ylab = "Simpson Index")
plot_sample4_simpson_s

#correlation- sample4- Pielou
corr_pielou_sample4_s <-  cor.test(df_7$bases_Gb,df_7$S4,method="spearman")
corr_pielou_sample4_s
#visualisation- scatter plot
plot_sample4_pielou_s <-ggscatter(df_7, x = "bases_Gb", y = "S4", 
                                  add = "reg.line", conf.int = TRUE,color="darkgreen",
                                  cor.coef = TRUE, cor.method = "spearman",
                                  xlab = "Sequencing Depth (Gb)", ylab = "Pielou Index")
plot_sample4_pielou_s

############-----------sample5----------##############

#correlation- sample5-shannonn
corr_shannon_sample5_s <- cor.test(df_5$bases_Gb,df_5$S5,method="spearman")
corr_shannon_sample5_s
#visualisation using scatter plot-
plot_sample5_shannon_s <- ggscatter(df_5, x = "bases_Gb", y = "S5", 
                                    add = "reg.line", conf.int = TRUE,color="brown",
                                    cor.coef = TRUE, cor.method = "spearman",
                                    xlab = "Sequencing Depth (Gb)", ylab = "Shannon Index")
plot_sample5_shannon_s
#correlation- sample5- simpson
corr_simpson_sample5_s <-  cor.test(df_6$bases_Gb,df_6$S5,method="spearman")
corr_simpson_sample5_s
#visualisation- scatter plot
plot_sample5_simpson_s <-ggscatter(df_6, x = "bases_Gb", y = "S5", 
                                   add = "reg.line", conf.int = TRUE,color="blue",
                                   cor.coef = TRUE, cor.method = "spearman",
                                   xlab = "Sequencing Depth (Gb)", ylab = "Simpson Index")
plot_sample5_simpson_s

#correlation- sample5- Pielou
corr_pielou_sample5_s <-  cor.test(df_7$bases_Gb,df_7$S5,method="spearman")
corr_pielou_sample5_s
#visualisation- scatter plot
plot_sample5_pielou_s <-ggscatter(df_7, x = "bases_Gb", y = "S5", 
                                  add = "reg.line", conf.int = TRUE,color="darkgreen", 
                                  cor.coef = TRUE, cor.method = "spearman",
                                  xlab = "Sequencing Depth (Gb)", ylab = "Pielou Index")
plot_sample5_pielou_s

##########------sample6----------##########

#correlation- sample6-shannon
corr_shannon_sample6_s <- cor.test(df_5$bases_Gb,df_5$S6,method="spearman")
corr_shannon_sample6_s
#visualisation using scatter plot-
plot_sample6_shannon_s <- ggscatter(df_5, x = "bases_Gb", y = "S6", 
                                    add = "reg.line", conf.int = TRUE,color="brown", 
                                    cor.coef = TRUE, cor.method = "spearman",
                                    xlab = "Sequencing Depth (Gb)", ylab = "Shannon Index")
plot_sample6_shannon_s
#correlation- sample6- simpson
corr_simpson_sample6_s <-  cor.test(df_6$bases_Gb,df_6$S6,method="spearman")
corr_simpson_sample6_s
#visualisation- scatter plot
plot_sample6_simpson_s <-ggscatter(df_6, x = "bases_Gb", y = "S6", 
                                   add = "reg.line", conf.int = TRUE,color="blue", 
                                   cor.coef = TRUE, cor.method = "spearman",
                                   xlab = "Sequencing Depth (Gb)", ylab = "Simpson Index")
plot_sample6_simpson_s

#correlation- sample6- Pielou
corr_pielou_sample6_s <-  cor.test(df_7$bases_Gb,df_7$S6,method="spearman")
corr_pielou_sample6_s
#visualisation- scatter plot
plot_sample6_pielou_s <-ggscatter(df_7, x = "bases_Gb", y = "S6", 
                                  add = "reg.line", conf.int = TRUE,color="darkgreen", 
                                  cor.coef = TRUE, cor.method = "spearman",
                                  xlab = "Sequencing Depth (Gb)", ylab = "Pielou Index")
plot_sample6_pielou_s

#############--------sample7----------############

#correlation- sample7-shannonn
corr_shannon_sample7_s <- cor.test(df_5$bases_Gb,df_5$S7,method="spearman")
corr_shannon_sample7_s
#visualisation using scatter plot-
plot_sample7_shannon_s <- ggscatter(df_5, x = "bases_Gb", y = "S7", 
                                    add = "reg.line", conf.int = TRUE,color="brown",
                                    cor.coef = TRUE, cor.method = "spearman",
                                    xlab = "Sequencing Depth (Gb)", ylab = "Shannon Index")
plot_sample7_shannon_s
#correlation- sample7- simpson
corr_simpson_sample7_s <-  cor.test(df_6$bases_Gb,df_6$S7,method="spearman")
corr_simpson_sample7_s
#visualisation- scatter plot
plot_sample7_simpson_s <-ggscatter(df_6, x = "bases_Gb", y = "S7", 
                                   add = "reg.line", conf.int = TRUE,color="blue",
                                   cor.coef = TRUE, cor.method = "spearman",
                                   xlab = "Sequencing Depth (Gb)", ylab = "Simpson Index")
plot_sample7_simpson_s

#correlation- sample7- Pielou
corr_pielou_sample7_s <-  cor.test(df_7$bases_Gb,df_7$S7,method="spearman")
corr_pielou_sample7_s
#visualisation- scatter plot
plot_sample7_pielou_s <-ggscatter(df_7, x = "bases_Gb", y = "S7", 
                                  add = "reg.line", conf.int = TRUE,color="darkgreen",
                                  cor.coef = TRUE, cor.method = "spearman",
                                  xlab = "Sequencing Depth (Gb)", ylab = "Pielou Index")
plot_sample7_pielou_s

#######---sample 8-----###

#correlation- sample8-shannonn
corr_shannon_sample8_s <- cor.test(df_5$bases_Gb,df_5$S8,method="spearman")
corr_shannon_sample8_s
#visualisation using scatter plot-
plot_sample8_shannon_s <- ggscatter(df_5, x = "bases_Gb", y = "S8", 
                                    add = "reg.line", conf.int = TRUE,color="brown",
                                    cor.coef = TRUE, cor.method = "spearman",
                                    xlab = "Sequencing Depth (Gb)", ylab = "Shannon Index")
plot_sample8_shannon_s
#correlation- sample8- simpson
corr_simpson_sample8_s <-  cor.test(df_6$bases_Gb,df_6$S8,method="spearman")
corr_simpson_sample8_s
#visualisation- scatter plot
plot_sample8_simpson_s <-ggscatter(df_6, x = "bases_Gb", y = "S8", 
                                   add = "reg.line", conf.int = TRUE,color="blue",
                                   cor.coef = TRUE, cor.method = "spearman",
                                   xlab = "Sequencing Depth (Gb)", ylab = "Simpson Index")
plot_sample8_simpson_s

#correlation- sample8- Pielou
corr_pielou_sample8_s <-  cor.test(df_7$bases_Gb,df_7$S8,method="spearman")
corr_pielou_sample8_s
#visualisation- scatter plot
plot_sample8_pielou_s <-ggscatter(df_7, x = "bases_Gb", y = "S8", 
                                  add = "reg.line", conf.int = TRUE,color="darkgreen",
                                  cor.coef = TRUE, cor.method = "spearman",
                                  xlab = "Sequencing Depth (Gb)", ylab = "Pielou Index")
plot_sample8_pielou_s


#######---sample 9-----###

#correlation- sample9-shannonn
corr_shannon_sample9_s <- cor.test(df_5$bases_Gb,df_5$S9,method="spearman")
corr_shannon_sample9_s
#visualisation using scatter plot-
plot_sample9_shannon_s <- ggscatter(df_5, x = "bases_Gb", y = "S9", 
                                    add = "reg.line", conf.int = TRUE,color="brown",
                                    cor.coef = TRUE, cor.method = "spearman",
                                    xlab = "Sequencing Depth (Gb)", ylab = "Shannon Index")
plot_sample9_shannon_s
#correlation- sample9- simpson
corr_simpson_sample9_s <-  cor.test(df_6$bases_Gb,df_6$S9,method="spearman")
corr_simpson_sample9_s
#visualisation- scatter plot
plot_sample9_simpson_s <-ggscatter(df_6, x = "bases_Gb", y = "S9", 
                                   add = "reg.line", conf.int = TRUE,color="blue",
                                   cor.coef = TRUE, cor.method = "spearman",
                                   xlab = "Sequencing Depth (Gb)", ylab = "Simpson Index")
plot_sample9_simpson_s

#correlation- sample9- Pielou
corr_pielou_sample9_s <-  cor.test(df_7$bases_Gb,df_7$S9,method="spearman")
corr_pielou_sample9_s
#visualisation- scatter plot
plot_sample9_pielou_s <-ggscatter(df_7, x = "bases_Gb", y = "S9", 
                                  add = "reg.line", conf.int = TRUE,color="darkgreen",
                                  cor.coef = TRUE, cor.method = "spearman",
                                  xlab = "Sequencing Depth (Gb)", ylab = "Pielou Index")
plot_sample9_pielou_s

#######---sample 10-----###

#correlation- sample10-shannonn
corr_shannon_sample10_s <- cor.test(df_5$bases_Gb,df_5$S10,method="spearman")
corr_shannon_sample10_s
#visualisation using scatter plot-
plot_sample10_shannon_s <- ggscatter(df_5, x = "bases_Gb", y = "S10", 
                                     add = "reg.line", conf.int = TRUE,color="brown",
                                     cor.coef = TRUE, cor.method = "spearman",
                                     xlab = "Sequencing Depth (Gb)", ylab = "Shannon Index")
plot_sample10_shannon_s
#correlation- sample10- simpson
corr_simpson_sample10_s <-  cor.test(df_6$bases_Gb,df_6$S10,method="spearman")
corr_simpson_sample10_s
#visualisation- scatter plot
plot_sample10_simpson_s <-ggscatter(df_6, x = "bases_Gb", y = "S10", 
                                    add = "reg.line", conf.int = TRUE,color="blue",
                                    cor.coef = TRUE, cor.method = "spearman",
                                    xlab = "Sequencing Depth (Gb)", ylab = "Simpson Index")
plot_sample10_simpson_s

#correlation- sample10- Pielou
corr_pielou_sample10_s <-  cor.test(df_7$bases_Gb,df_7$S10,method="spearman")
corr_pielou_sample10_s
#visualisation- scatter plot
plot_sample10_pielou_s <-ggscatter(df_7, x = "bases_Gb", y = "S10", 
                                   add = "reg.line", conf.int = TRUE,color="darkgreen",
                                   cor.coef = TRUE, cor.method = "spearman",
                                   xlab = "Sequencing Depth (Gb)", ylab = "Pielou Index")
plot_sample10_pielou_s



###--BETA DIVERSITY----###

###beta-diversity-species

####---These steps are already perfomed---###--dont perform again
#Input file
bray_species <- read.csv("C:/Users/sneha/OneDrive/Desktop/final_draft/species/sample1to10_species_master_final.csv",header=TRUE)
#transpose the file
bray_species_t <- t(bray_species)
#save the transposed file
write.csv(bray_species_t,file = "C:/Users/sneha/OneDrive/Desktop/final_draft/species/bray_species_t_final.csv")
#modify the transposed file(bray_genus_t.csv) manually by adding extra row, column lables as required and save it as a new copy- bray_genus_t_modifed.csv

####----start from here----####
library(vegan)
library(ggplot2)
bray_species_t_NMDS <- read.csv("C:/Users/sneha/OneDrive/Desktop/final_draft/species/bray_species_t_modified_final.csv",header=TRUE)

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
  scale_shape_manual(values=c(0,1,2,3,4,5))
#10 samples, so 10 colours (sample 1 to sample 10); 6 groups, so 6 shapes (0.25Gb,0.5Gb,0.75Gb,1Gb,1.25Gb,raw)
xx_species
#to encircle the clustering
xx_species + stat_ellipse(aes(group = group), geom = "polygon", fill = NA, color = "black")

ggsave("NMDS.svg")


###------done------###






