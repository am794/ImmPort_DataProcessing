##############################################
# ImmPort downstream analysis and statistics #
##############################################

## SDY224: Bcells panel: It has 136 FCS files with 0-10 days time points

#Read the file 
library(openxlsx)
library(ggplot2)
library(reshape2)
#percentages_PP <- read.xlsx("/Volumes/Samsung_T5/Immport_UH2/SDY224/Bcell/Events_percentages_sdy224_bcells.xlsx",sheet=6)[-c(125,126),]
percentages_PP <- read.xlsx("/Users/amandava/Desktop/ImmPort/Events_percentages_sdy224_bcells.xlsx",sheet=6)[-c(125,126),]
perc <- percentages_PP[,c(1,5,6,7,8)]

plasmablasts <- as.data.frame(unclass(percentages_PP[,c(1,5,6,7)]))
plasmablasts$Days <- factor(plasmablasts$Days, levels=unique(plasmablasts$Days))

plasmacells <- percentages_PP[,c(1,5,6,8)]
plasmacells$Days <- factor(plasmacells$Days, levels=unique(plasmacells$Days))

quartz()
ggplot(data=plasmablasts)+geom_boxplot(aes(x=Days,y=Plasmablasts,fill=Days))+ scale_fill_brewer(palette = "Blues")+
  guides(fill = FALSE)+geom_line(aes(x=Days,y=Plasmablasts,group=Subject,color=Subject))+ggtitle("Population percentages of Plasmablast cells (Parent population: CD3- CD19+ CD20)")+
  xlab("Visit")+ylab("Population percentages")

quartz()
ggplot(data=plasmacells)+geom_boxplot(aes(x=Days,y=Plasma.cells,fill=Days))+ scale_fill_brewer(palette = "Blues")+
  geom_line(aes(x=Days,y=Plasma.cells,group=Subject,color=Subject,linetype=Subject))+
  guides(fill = FALSE)+ggtitle("Population percentages of Plasma cells (Parent population: CD3- CD19+ CD20)")+
  xlab("Visit")+ylab("Population percentages")+ylim(0,6.5)
