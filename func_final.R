###code to generate plots for functional analysis
#Main plots- for Wilcox test, used in our project 

#######---function--- final draft---######
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

#NO of observed pathways- wilcox test
##Wilcox test- to compare each depth against the reference depth of 1.25Gb (highest for us)
library(ggplot2)
library(ggpubr)
library(rstatix)
##reading the input file
obs_pathway <- read.table("C:/Users/sneha/OneDrive/Desktop/final_draft/function/pathway_count_wilcox.csv", sep = ",", header = T, check.names = F, row.names = 1, comment.char = "", stringsAsFactors = F)
test_t <- t(obs_pathway)
test_melt <- reshape2::melt(test_t, value.name = "Observed_pathways")
names(test_melt)[1] <- "Sequencing_Depth"
names(test_melt)[2] <- "Sample"
test_melt$Sequencing_Depth <- factor(test_melt$Sequencing_Depth, levels = c("0.25Gb", "0.5Gb", "0.75Gb", "1Gb", "1.25Gb"))
# Generate boxplot
boxplot <- ggplot(data = test_melt, mapping=
                    aes(x = Sequencing_Depth, y = Observed_pathways,colour=Sequencing_Depth)) +
  geom_boxplot()+ xlab("Sequencing Depth") + ylab("Number of Pathways (KEGG)")+
  geom_point(position = position_jitterdodge(jitter.width = 1))
boxplot
#wilcox test- comparing each group with reference group
stat.test <- stat.test <- test_melt %>%
  wilcox_test(Observed_pathways ~ Sequencing_Depth, ref.group = "1.25Gb", paired = FALSE) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance("p.adj") %>%
  add_xy_position()
stat.test

#adding p values to box plot
p <- boxplot + stat_pvalue_manual(
  stat.test, label = "p.adj") # Change label to p.adj.signif is asterisks are preferred
p + theme_classic()
#-----x-----##

##KEGG shannon
library(ggplot2)
df_kegg <- read.csv("C:/Users/sneha/OneDrive/Desktop/final_draft/function/shannon_all_func_final.csv",header=TRUE)
df_kegg_shannon <-  ggplot(df_kegg, aes(x=bases_Gb))+ geom_line(aes(y=S1,group=1,colour="S1"))+geom_point(aes(y=S1,group=1,color="pink"))+geom_line(aes(y=S2,group=1,color="S2"))+geom_point(aes(y=S2 ,group=1,color="pink"))+
  geom_line(aes(y=S3,group=1,colour="S3"))+geom_point(aes(y=S3,group=1,color="pink"))+geom_line(aes(y=S4,group=1,colour="S4"))+geom_point(aes(y=S4,group=1,color="pink"))+ geom_line(aes(y=S5,group=1,colour="S5"))+geom_point(aes(y=S5,group=1,color="pink"))+
  geom_line(aes(y=S6,group=1,colour="S6"))+geom_point(aes(y=S6,group=1,color="pink"))+geom_line(aes(y=S7,group=1,colour="S7"))+geom_point(aes(y=S7,group=1,color="pink"))+
  geom_line(aes(y=S8,group=1,colour="S8"))+geom_point(aes(y=S8,group=1,color="pink"))+
  geom_line(aes(y=S9,group=1,colour="S9"))+geom_point(aes(y=S9,group=1,color="pink"))+
  geom_line(aes(y=S10,group=1,colour="S10"))+geom_point(aes(y=S10,group=1,color="pink"))+labs(x="Sequencing Depth (Gb)",y="Shannon Index")
df_kegg_shannon


#df_1_shannon_ribbon <-  ggplot(df_ribbon, aes(x=bases_Gb))+ geom_ribbon(aes(y=S1,group=1,colour="S1"))+geom_point(aes(y=S1,group=1,color="pink"))+geom_ribbon(aes(y=S2,group=1,color="S2"))+geom_point(aes(y=S2 ,group=1,color="pink"))+
  geom_ribbon(aes(y=S3,group=1,colour="S3"))+geom_point(aes(y=S3,group=1,color="pink"))+geom_ribbon(aes(y=S4,group=1,colour="S4"))+geom_point(aes(y=S4,group=1,color="pink"))+ geom_ribbon(aes(y=S5,group=1,colour="S5"))+geom_point(aes(y=S5,group=1,color="pink"))+
  geom_ribbon(aes(y=S6,group=1,colour="S6"))+geom_point(aes(y=S6,group=1,color="pink"))+geom_ribbon(aes(y=S7,group=1,colour="S7"))+geom_point(aes(y=S7,group=1,color="pink"))+
  geom_ribbon(aes(y=S8,group=1,colour="S8"))+geom_point(aes(y=S8,group=1,color="pink"))+
  geom_ribbon(aes(y=S9,group=1,colour="S9"))+geom_point(aes(y=S9,group=1,color="pink"))+
  geom_ribbon(aes(y=S10,group=1,colour="S10"))+geom_point(aes(y=S10,group=1,color="pink"))+labs(x="Sequencing Depth (Gb)",y="Shannon Index")
#df_1_shannon_ribbon