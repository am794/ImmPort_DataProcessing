#################
# modeling data #
#################

library("lme4")
library("emmeans")

#####Read the files########
sdy144 <- read.xlsx("/Volumes/Samsung_T5/Immport_UH2/Resutls_Analysis/MemoryCD4_CD8/SDY144_MemoryTcells.xlsx",sheet=2)[,-1]
sdy364 <- read.xlsx("/Volumes/Samsung_T5/Immport_UH2/Resutls_Analysis/MemoryCD4_CD8/SDY364_MemoryTcells.xlsx",sheet=2)[,-1]
sdy368 <- read.xlsx("/Volumes/Samsung_T5/Immport_UH2/Resutls_Analysis/MemoryCD4_CD8/SDY368_MemoryTcells.xlsx",sheet=2)[,-1]
sdy387 <- read.xlsx("/Volumes/Samsung_T5/Immport_UH2/Resutls_Analysis/MemoryCD4_CD8/SDY387_MemoryTcells.xlsx",sheet=3)[,-1]
sdy404 <- read.xlsx("/Volumes/Samsung_T5/Immport_UH2/Resutls_Analysis/MemoryCD4_CD8/SDY404_MemoryTcells.xlsx",sheet=2)[,-1]

#####SDY144
head(sdy144)
sdy144$study_time_collected <- paste0("Day_",sdy144$study_time_collected)
melt_sdy144 <- melt(sdy144,id=c("name","study_time_collected","subject_accession","max_subject_age","Study"))
colnames(melt_sdy144)[6:7] <- c("CellType","CellProportion")

sdy144_abtiter <- read.table("/Volumes/Samsung_T5/Immport_UH2/Companion_studies_1/SDY144/Study/Tab/SDY144-DR28_Tab/Tab/hai_result.txt",sep = "\t", header=TRUE)
ab_resp <- function(sdy144_abtiter){
  sdy144_ab_titer_d <- dcast(sdy144_abtiter,SUBJECT_ACCESSION+STUDY_TIME_COLLECTED~VIRUS_STRAIN_PREFERRED,value.var='VALUE_REPORTED',fun.aggregate = mean)
  cnames <- colnames(sdy144_ab_titer_d)
  colnames(sdy144_ab_titer_d)[3:length(cnames)] <- lapply(cnames[3:length(cnames)],function(x) paste0(strsplit(x,"/")[[1]][1],"_",strsplit(x,"/")[[1]][2]))
  sdy144_ab_titer_d$Geometric_mean <- apply(sdy144_ab_titer_d[,3:5],1,function(x) geometric.mean(x,na.rm=TRUE))
  sdy144_ab_titer_d$STUDY_TIME_COLLECTED <- paste0("Day",sdy144_ab_titer_d$STUDY_TIME_COLLECTED)
  resp <- dcast(sdy144_ab_titer_d,SUBJECT_ACCESSION~STUDY_TIME_COLLECTED,value.var='Geometric_mean')
  resp$FoldChange <- resp[,3]/resp[,2]
  resp$responder0[resp[,2] >40 | resp[,2] == 40] = "Responder"
  resp$responder0[resp[,2]<40] = "NonResponder"
  resp$responder30[resp[,3] > 40 | resp[,3] == 40] = "Responder"
  resp$responder30[resp[,3]<40] = "NonResponder"
  resp$responderFC[resp$FoldChange > 4 | resp$FoldChange == 4] ="Responder"
  resp$responderFC[resp$FoldChange < 4] = "NonResponder"
  resp$responder_Day0_FC <- ifelse(resp$responder0 == "Responder" | resp$responderFC == "Responder","Responder","NonResponder")
  return(resp)
}
resp_144 <- ab_resp(sdy144_abtiter)
colnames(resp_144)[1] <- "subject_accession"
merged_sdy144 <- merge(sdy144,resp_144,by ="subject_accession",all.x=TRUE)
merged_sdy144$study_time_collected <- as.factor(merged_sdy144$study_time_collected)
#merged_sdy144$max_subject_age <- as.factor(merged_sdy144$max_subject_age)
merged_sdy144 <- merged_sdy144[order(merged_sdy144$study_time_collected),]
merged_sdy144$responder0 <- as.factor(merged_sdy144$responder0)
merged_sdy144$study_time_collected <- factor(merged_sdy144$study_time_collected,levels = c("Day_0","Day_1","Day_7","Day_28"))
mod1 <- glmer(responder0~study_time_collected+CD4p_Tcells+FoldChange+CD4p_Tcells*study_time_collected+(1|subject_accession),data = merged_sdy144,family = binomial,,control=glmerControl(optimizer="bobyqa",
                                                                                                                                                                              optCtrl=list(maxfun=2e5)))
quartz()
ggplot(merged_sdy144,aes(x=study_time_collected,y=CD4p_Tcells,color=subject_accession))+
  geom_point()+geom_smooth(aes(x=study_time_collected,y=CD4p_Tcells),method="glm")+scale_y_continuous(limits=c(0,1))+
  scale_colour_orig(guide="none")

melted_merged_sdy144 <- na.omit(melt(merged_sdy144,id= c("subject_accession","name","study_time_collected","max_subject_age","Study","Day0","Day30","FoldChange","responder0","responder30","responderFC" ,"responder_Day0_FC")))
melted_merged_sdy144 <- melted_merged_sdy144[order(melted_merged_sdy144$study_time_collected),]
melted_merged_sdy144$study_time_collected <- factor(melted_merged_sdy144$study_time_collected,levels = c("Day_0","Day_1","Day_7","Day_28"))
quartz()
ggplot(melted_merged_sdy144,aes(x=study_time_collected,y=value,color=responder0))+
  geom_line(aes(group=subject_accession))+facet_wrap(~variable,scales = "free")+
  ggtitle("SDY144")+ylab("Cell Proportion")+labs(color='Responder') 

#Visits <- as.character(unique(melted_merged_sdy144$study_time_collected))
#melted_merged_sdy144$study_time_collected <- dplyr::recode(melted_merged_sdy144$study_time_collected,Day_1="Day_01",Day_7="Day_07")

mod3 <- glmer(CD4p_Tcells~study_time_collected+responder0+responder0*study_time_collected+(1|subject_accession),data = merged_sdy144,family=binomial(lin=logit))
em <- emmeans(mod3,~study_time_collected|responder0)
plot(em)


mod11 <- lmer(CD4p_Tcells~study_time_collected+responder0+FoldChange+responder0*study_time_collected+(1|subject_accession),data = merged_sdy144)
em11 <- emmeans(mod11,~study_time_collected|responder0)
pv11 <- contrast(em11,alpha=0.05,interaction = "trt.vs.ctrl1")
em12 <- emmeans(mod11,~responder0|study_time_collected)
pv12 <- contrast(em12,alpha=0.05,interaction = "trt.vs.ctrl1")
em13 <- emmeans(mod11,~responder0*study_time_collected)
pv13 <- contrast(em13,alpha=0.05,interaction = "trt.vs.ctrl1")
em14 <- emmeans(mod11,~study_time_collected|responder0)
p <- plot(em14)
p+coord_flip()

mod21 <- lmer(CD8p_Tcells~study_time_collected+responder0+FoldChange+responder0*study_time_collected+(1|subject_accession),data = merged_sdy144)
em21 <- emmeans(mod21,~study_time_collected|responder0)
pv21 <- contrast(em21,alpha=0.05,interaction = "trt.vs.ctrl1")
em22 <- emmeans(mod21,~responder0|study_time_collected)
pv22 <- contrast(em22,alpha=0.05,interaction = "trt.vs.ctrl1")
em23 <- emmeans(mod21,~responder0*study_time_collected)
pv23 <- contrast(em23,alpha=0.05,interaction = "trt.vs.ctrl1")
em24 <- emmeans(mod21,~study_time_collected|responder0)
p <- plot(em24)
p+coord_flip()


mod11.a <- lmer(CD4p_CD45RAn_CCR7n~study_time_collected+responder0+FoldChange+responder0*study_time_collected+(1|subject_accession),data = merged_sdy144)
em11.a <- emmeans(mod11.a,~study_time_collected|responder0)
pv11.a <- contrast(em11.a,alpha=0.05,interaction = "trt.vs.ctrl1")
em12.a <- emmeans(mod11.a,~responder0|study_time_collected)
pv12.a <- contrast(em12.a,alpha=0.05,interaction = "trt.vs.ctrl1")
em13.a <- emmeans(mod11.a,~responder0*study_time_collected)
pv13.a <- contrast(em13.a,alpha=0.05,interaction = "trt.vs.ctrl1")
em14.a <- emmeans(mod11.a,~study_time_collected|responder0)
p <- plot(em14.a)
p+coord_flip()

mod11.b <- lmer(CD4p_CD45RAp_CCR7n~study_time_collected+responder0+FoldChange+responder0*study_time_collected+(1|subject_accession),data = merged_sdy144)
em11.b <- emmeans(mod11.b,~study_time_collected|responder0)
pv11.b <- contrast(em11.b,alpha=0.05,interaction = "trt.vs.ctrl1")
em12.b <- emmeans(mod11.b,~responder0|study_time_collected)
pv12.b <- contrast(em12.b,alpha=0.05,interaction = "trt.vs.ctrl1")
em13.b <- emmeans(mod11.b,~responder0*study_time_collected)
pv13.b <- contrast(em13.b,alpha=0.05,interaction = "trt.vs.ctrl1")
em14.b <- emmeans(mod11.b,~study_time_collected|responder0)
p <- plot(em14.b)
p+coord_flip()


mod11.c <- lmer(CD4p_CD45RAn_CCR7p~study_time_collected+responder0+FoldChange+responder0*study_time_collected+(1|subject_accession),data = merged_sdy144)
em11.c <- emmeans(mod11.c,~study_time_collected|responder0)
pv11.c <- contrast(em11.c,alpha=0.05,interaction = "trt.vs.ctrl1")
em12.c <- emmeans(mod11.c,~responder0|study_time_collected)
pv12.c <- contrast(em12.c,alpha=0.05,interaction = "trt.vs.ctrl1")
em13.c <- emmeans(mod11.c,~responder0*study_time_collected)
pv13.c <- contrast(em13.c,alpha=0.05,interaction = "trt.vs.ctrl1")
em14.c <- emmeans(mod11.c,~study_time_collected|responder0)
p <- plot(em14.c)
p+coord_flip()


mod11.d <- lmer(CD4p_CD45RAp_CCR7p~study_time_collected+responder0+responder0*study_time_collected+(1|subject_accession),data = merged_sdy144)
em11.d <- emmeans(mod11.d,~study_time_collected|responder0)
pv11.d <- contrast(em11.d,alpha=0.05,interaction = "trt.vs.ctrl1")
em12.d <- emmeans(mod11.d,~responder0|study_time_collected)
pv12.d <- contrast(em12.d,alpha=0.05,interaction = "trt.vs.ctrl1")
em13.d <- emmeans(mod11.d,~responder0*study_time_collected)
pv13.d <- contrast(em13.d,alpha=0.05,interaction = "trt.vs.ctrl1")
em14.d <- emmeans(mod11.d,~study_time_collected|responder0)
p <- plot(em14.d)
p+coord_flip()


mod21.a <- lmer(CD8p_CD45RAn_CCR7n~study_time_collected+responder0+responder0*study_time_collected+(1|subject_accession),data = merged_sdy144)
em21.a <- emmeans(mod21.a,~study_time_collected|responder0)
pv21.a <- contrast(em21.a,alpha=0.05,interaction = "trt.vs.ctrl1")
em22.a <- emmeans(mod21.a,~responder0|study_time_collected)
pv22.a <- contrast(em22.a,alpha=0.05,interaction = "trt.vs.ctrl1")
em23.a <- emmeans(mod21.a,~responder0*study_time_collected)
pv23.a <- contrast(em23.a,alpha=0.05,interaction = "trt.vs.ctrl1")
em24.a <- emmeans(mod21.a,~study_time_collected|responder0)
p <- plot(em24.a)
p+coord_flip()

mod21.b <- lmer(CD8p_CD45RAp_CCR7n~study_time_collected+responder0+responder0*study_time_collected+(1|subject_accession),data = merged_sdy144)
em21.b <- emmeans(mod21.b,~study_time_collected|responder0)
pv21.b <- contrast(em21.b,alpha=0.05,interaction = "trt.vs.ctrl1")
em22.b <- emmeans(mod21.b,~responder0|study_time_collected)
pv22.b <- contrast(em22.b,alpha=0.05,interaction = "trt.vs.ctrl1")
em23.b <- emmeans(mod21.b,~responder0*study_time_collected)
pv23.b <- contrast(em23.b,alpha=0.05,interaction = "trt.vs.ctrl1")
em24.b <- emmeans(mod21.b,~study_time_collected|responder0)
p <- plot(em24.b)
p+coord_flip()

mod21.c <- lmer(CD8p_CD45RAn_CCR7p~study_time_collected+responder0+responder0*study_time_collected+(1|subject_accession),data = merged_sdy144)
em21.c<- emmeans(mod21.c,~study_time_collected|responder0)
pv21.c <- contrast(em21.c,alpha=0.05,interaction = "trt.vs.ctrl1")
em22.c <- emmeans(mod21.c,~responder0|study_time_collected)
pv22.c <- contrast(em22.c,alpha=0.05,interaction = "trt.vs.ctrl1")
em23.c <- emmeans(mod21.c,~responder0*study_time_collected)
pv23.c <- contrast(em23.c,alpha=0.05,interaction = "trt.vs.ctrl1")
em24.c <- emmeans(mod21.c,~study_time_collected|responder0)
p <- plot(em24.c)
p+coord_flip()


mod21.d <- lmer(CD8p_CD45RAp_CCR7p~study_time_collected+responder0+responder0*study_time_collected+(1|subject_accession),data = merged_sdy144)
em21.d<- emmeans(mod21.d,~study_time_collected|responder0)
pv21.d <- contrast(em21.d,alpha=0.05,interaction = "trt.vs.ctrl1")
em22.d <- emmeans(mod21.d,~responder0|study_time_collected)
pv22.d <- contrast(em22.d,alpha=0.05,interaction = "trt.vs.ctrl1")
em23.d <- emmeans(mod21.d,~responder0*study_time_collected)
pv23.d <- contrast(em23.d,alpha=0.05,interaction = "trt.vs.ctrl1")
em24.d <- emmeans(mod21.d,~study_time_collected|responder0)
p <- plot(em24.d)
p+coord_flip()

#####364
head(sdy364)
sdy364$study_time_collected <- paste0("Day_",sdy364$study_time_collected)
melt_sdy364 <- melt(sdy364,id=c("name","study_time_collected","subject_accession","max_subject_age","Study"))
colnames(melt_sdy364)[6:7] <- c("CellType","CellProportion")
sdy364_abtiter <- read.table("/Volumes/Samsung_T5/Immport_UH2/Companion_studies_1/SDY364/Study/Tab/SDY364-DR29_Tab/Tab/hai_result.txt",sep = "\t", header=TRUE)
resp_364 <- ab_resp(sdy364_abtiter)
colnames(resp_364)[1] <- "subject_accession"
merged_sdy364 <- merge(sdy364,resp_364,by ="subject_accession",all.x=TRUE)
merged_sdy364$study_time_collected <- as.factor(merged_sdy364$study_time_collected)
#merged_sdy144$max_subject_age <- as.factor(merged_sdy144$max_subject_age)
merged_sdy364 <- merged_sdy364[order(merged_sdy364$study_time_collected),]
merged_sdy364$responder0 <- as.factor(merged_sdy364$responder0)
merged_sdy364$study_time_collected <- factor(merged_sdy364$study_time_collected,levels = c("Day_0","Day_1","Day_7","Day_28"))
                                                                                                                                                                            
quartz()
ggplot(merged_sdy364,aes(x=study_time_collected,y=CD4p_Tcells,color=subject_accession))+
  geom_point()+geom_smooth(aes(x=study_time_collected,y=CD4p_Tcells),method="glm")+scale_y_continuous(limits=c(0,1))+
  scale_colour_orig(guide="none")

melted_merged_sdy364 <- na.omit(melt(merged_sdy364,id= c("subject_accession","name","study_time_collected","max_subject_age","Study","Day0","Day28","FoldChange","responder0","responder30","responderFC" ,"responder_Day0_FC")))
melted_merged_sdy364 <- melted_merged_sdy364[order(melted_merged_sdy364$study_time_collected),]
melted_merged_sdy364$study_time_collected <- factor(melted_merged_sdy364$study_time_collected,levels = c("Day_0","Day_1","Day_7","Day_28"))
quartz()
ggplot(melted_merged_sdy364,aes(x=study_time_collected,y=value,color=responder0))+
  geom_line(aes(group=subject_accession))+facet_wrap(~variable,scales = "free")+
  ggtitle("SDY364")+ylab("Cell Proportion")+labs(color='Responder') 


mod11_sdy364 <- lmer(CD4p_Tcells~study_time_collected+responder0+responder0*study_time_collected+(1|subject_accession),data = merged_sdy364)
em11_sdy364 <- emmeans(mod11_sdy364,~study_time_collected|responder0)
pv11_sdy364 <- contrast(em11_sdy364,alpha=0.05,interaction = "trt.vs.ctrl1")
em12_sdy364 <- emmeans(mod11_sdy364,~responder0|study_time_collected)
pv12_sdy364 <- contrast(em12_sdy364,alpha=0.05,interaction = "trt.vs.ctrl1")
em13_sdy364 <- emmeans(mod11_sdy364,~responder0*study_time_collected)
pv13_sdy364 <- contrast(em13_sdy364,alpha=0.05,interaction = "trt.vs.ctrl1")
em14_sdy364 <- emmeans(mod11_sdy364,~study_time_collected|responder0)
p <- plot(em14_sdy364)
p+coord_flip()

mod21_sdy364 <- lmer(CD8p_Tcells~study_time_collected+responder0+responder0*study_time_collected+(1|subject_accession),data = merged_sdy364)
em21_sdy364 <- emmeans(mod21_sdy364,~study_time_collected|responder0)
pv21_sdy364 <- contrast(em21_sdy364,alpha=0.05,interaction = "trt.vs.ctrl1")
em22_sdy364 <- emmeans(mod21_sdy364,~responder0|study_time_collected)
pv22_sdy364 <- contrast(em22_sdy364,alpha=0.05,interaction = "trt.vs.ctrl1")
em23_sdy364 <- emmeans(mod21_sdy364,~responder0*study_time_collected)
pv23_sdy364 <- contrast(em23_sdy364,alpha=0.05,interaction = "trt.vs.ctrl1")
em24_sdy364 <- emmeans(mod21_sdy364,~study_time_collected|responder0)
p <- plot(em24_sdy364)
p+coord_flip()

mod11_sdy364.a <- lmer(CD4p_CD45RAn_CCR7n~study_time_collected+responder0+responder0*study_time_collected+(1|subject_accession),data = merged_sdy364)
em11_sdy364.a <- emmeans(mod11_sdy364.a,~study_time_collected|responder0)
pv11_sdy364.a <- contrast(em11_sdy364.a,alpha=0.05,interaction = "trt.vs.ctrl1")
em12_sdy364.a <- emmeans(mod11_sdy364.a,~responder0|study_time_collected)
pv12_sdy364.a <- contrast(em12_sdy364.a,alpha=0.05,interaction = "trt.vs.ctrl1")
em13_sdy364.a <- emmeans(mod11_sdy364.a,~responder0*study_time_collected)
pv13_sdy364.a <- contrast(em13_sdy364.a,alpha=0.05,interaction = "trt.vs.ctrl1")
em14_sdy364.a <- emmeans(mod11_sdy364.a,~study_time_collected|responder0)
p <- plot(em14_sdy364.a)
p+coord_flip()


mod11_sdy364.b <- lmer(CD4p_CD45RAp_CCR7n~study_time_collected+responder0+responder0*study_time_collected+(1|subject_accession),data = merged_sdy364)
em11_sdy364.b <- emmeans(mod11_sdy364.b,~study_time_collected|responder0)
pv11_sdy364.b <- contrast(em11_sdy364.b,alpha=0.05,interaction = "trt.vs.ctrl1")
em12_sdy364.b <- emmeans(mod11_sdy364.b,~responder0|study_time_collected)
pv12_sdy364.b <- contrast(em12_sdy364.b,alpha=0.05,interaction = "trt.vs.ctrl1")
em13_sdy364.b <- emmeans(mod11_sdy364.b,~responder0*study_time_collected)
pv13_sdy364.b <- contrast(em13_sdy364.b,alpha=0.05,interaction = "trt.vs.ctrl1")
em14_sdy364.b <- emmeans(mod11_sdy364.b,~study_time_collected|responder0)
p <- plot(em14_sdy364.b)
p+coord_flip()


mod11_sdy364.c <- lmer(CD4p_CD45RAn_CCR7p~study_time_collected+responder0+responder0*study_time_collected+(1|subject_accession),data = merged_sdy364)
em11_sdy364.c <- emmeans(mod11_sdy364.c,~study_time_collected|responder0)
pv11_sdy364.c <- contrast(em11_sdy364.c,alpha=0.05,interaction = "trt.vs.ctrl1")
em12_sdy364.c <- emmeans(mod11_sdy364.c,~responder0|study_time_collected)
pv12_sdy364.c <- contrast(em12_sdy364.c,alpha=0.05,interaction = "trt.vs.ctrl1")
em13_sdy364.c <- emmeans(mod11_sdy364.c,~responder0*study_time_collected)
pv13_sdy364.c <- contrast(em13_sdy364.c,alpha=0.05,interaction = "trt.vs.ctrl1")
em14_sdy364.c <- emmeans(mod11_sdy364.c,~study_time_collected|responder0)
p <- plot(em14_sdy364.c)
p+coord_flip()


mod11_sdy364.d <- lmer(CD4p_CD45RAp_CCR7p~study_time_collected+responder0+responder0*study_time_collected+(1|subject_accession),data = merged_sdy364)
em11_sdy364.d <- emmeans(mod11_sdy364.d,~study_time_collected|responder0)
pv11_sdy364.d <- contrast(em11_sdy364.d,alpha=0.05,interaction = "trt.vs.ctrl1")
em12_sdy364.d <- emmeans(mod11_sdy364.d,~responder0|study_time_collected)
pv12_sdy364.d <- contrast(em12_sdy364.d,alpha=0.05,interaction = "trt.vs.ctrl1")
em13_sdy364.d <- emmeans(mod11_sdy364.d,~responder0*study_time_collected)
pv13_sdy364.d <- contrast(em13_sdy364.d,alpha=0.05,interaction = "trt.vs.ctrl1")
em14_sdy364.d <- emmeans(mod11_sdy364.d,~study_time_collected|responder0)
p <- plot(em14_sdy364.d)
p+coord_flip()


mod21_sdy364.a <- lmer(CD8p_CD45RAn_CCR7n~study_time_collected+responder0+responder0*study_time_collected+(1|subject_accession),data = merged_sdy364)
em21_sdy364.a <- emmeans(mod21_sdy364.a,~study_time_collected|responder0)
pv21_sdy364.a <- contrast(em21_sdy364.a,alpha=0.05,interaction = "trt.vs.ctrl1")
em22_sdy364.a <- emmeans(mod21_sdy364.a,~responder0|study_time_collected)
pv22_sdy364.a <- contrast(em22_sdy364.a,alpha=0.05,interaction = "trt.vs.ctrl1")
em23_sdy364.a <- emmeans(mod21_sdy364.a,~responder0*study_time_collected)
pv23_sdy364.a <- contrast(em23_sdy364.a,alpha=0.05,interaction = "trt.vs.ctrl1")
em24_sdy364.a <- emmeans(mod21_sdy364.a,~study_time_collected|responder0)
p <- plot(em24_sdy364.a)
p+coord_flip()


mod21_sdy364.b <- lmer(CD8p_CD45RAp_CCR7n~study_time_collected+responder0+responder0*study_time_collected+(1|subject_accession),data = merged_sdy364)
em21_sdy364.b <- emmeans(mod21_sdy364.b,~study_time_collected|responder0)
pv21_sdy364.b <- contrast(em21_sdy364.b,alpha=0.05,interaction = "trt.vs.ctrl1")
em22_sdy364.b <- emmeans(mod21_sdy364.b,~responder0|study_time_collected)
pv22_sdy364.b <- contrast(em22_sdy364.b,alpha=0.05,interaction = "trt.vs.ctrl1")
em23_sdy364.b <- emmeans(mod21_sdy364.b,~responder0*study_time_collected)
pv23_sdy364.b <- contrast(em23_sdy364.b,alpha=0.05,interaction = "trt.vs.ctrl1")
em24_sdy364.b <- emmeans(mod21_sdy364.b,~study_time_collected|responder0)
p <- plot(em24_sdy364.b)
p+coord_flip()


mod21_sdy364.c <- lmer(CD8p_CD45RAn_CCR7p~study_time_collected+responder0+responder0*study_time_collected+(1|subject_accession),data = merged_sdy364)
em21_sdy364.c <- emmeans(mod21_sdy364.c,~study_time_collected|responder0)
pv21_sdy364.c <- contrast(em21_sdy364.c,alpha=0.05,interaction = "trt.vs.ctrl1")
em22_sdy364.c <- emmeans(mod21_sdy364.c,~responder0|study_time_collected)
pv22_sdy364.c <- contrast(em22_sdy364.c,alpha=0.05,interaction = "trt.vs.ctrl1")
em23_sdy364.c <- emmeans(mod21_sdy364.c,~responder0*study_time_collected)
pv23_sdy364.c <- contrast(em23_sdy364.c,alpha=0.05,interaction = "trt.vs.ctrl1")
em24_sdy364.c <- emmeans(mod21_sdy364.c,~study_time_collected|responder0)
p <- plot(em24_sdy364.c)
p+coord_flip()


mod21_sdy364.d <- lmer(CD8p_CD45RAp_CCR7p~study_time_collected+responder0+responder0*study_time_collected+(1|subject_accession),data = merged_sdy364)
em21_sdy364.d <- emmeans(mod21_sdy364.d,~study_time_collected|responder0)
pv21_sdy364.d <- contrast(em21_sdy364.d,alpha=0.05,interaction = "trt.vs.ctrl1")
em22_sdy364.d <- emmeans(mod21_sdy364.d,~responder0|study_time_collected)
pv22_sdy364.d <- contrast(em22_sdy364.d,alpha=0.05,interaction = "trt.vs.ctrl1")
em23_sdy364.d <- emmeans(mod21_sdy364.d,~responder0*study_time_collected)
pv23_sdy364.d <- contrast(em23_sdy364.d,alpha=0.05,interaction = "trt.vs.ctrl1")
em24_sdy364.d <- emmeans(mod21_sdy364.d,~study_time_collected|responder0)
p <- plot(em24_sdy364.d)
p+coord_flip()

#############SDY368
#####364
head(sdy368)
#sdy364$study_time_collected <- paste0("Day_",sdy364$study_time_collected)
melt_sdy368 <- melt(sdy368,id=c("name","study_time_collected","subject_accession","max_subject_age","Study"))
colnames(melt_sdy368)[6:7] <- c("CellType","CellProportion")
sdy368_abtiter <- read.table("/Volumes/Samsung_T5/Immport_UH2/Companion_studies_1/SDY368/Study/Tab/SDY368-DR29_Tab/Tab/hai_result.txt",sep = "\t", header=TRUE)
resp_368 <- ab_resp(sdy368_abtiter)
colnames(resp_368)[1] <- "subject_accession"
merged_sdy368 <- merge(sdy368,resp_368,by ="subject_accession",all.x=TRUE)
merged_sdy368$study_time_collected <- as.factor(merged_sdy368$study_time_collected)
#merged_sdy144$max_subject_age <- as.factor(merged_sdy144$max_subject_age)
merged_sdy368 <- merged_sdy368[order(merged_sdy368$study_time_collected),]
merged_sdy368$responder0 <- as.factor(merged_sdy368$responder0)
merged_sdy368$study_time_collected <- factor(merged_sdy368$study_time_collected,levels = c("Day_0","Day_1","Day_7","Day_28"))

quartz()
ggplot(merged_sdy364,aes(x=study_time_collected,y=CD4p_Tcells,color=subject_accession))+
  geom_point()+geom_smooth(aes(x=study_time_collected,y=CD4p_Tcells),method="glm")+scale_y_continuous(limits=c(0,1))+
  scale_colour_orig(guide="none")

melted_merged_sdy368 <- na.omit(melt(merged_sdy368,id= c("subject_accession","name","study_time_collected","max_subject_age","Study","Day0","Day28","FoldChange","responder0","responder30","responderFC" ,"responder_Day0_FC")))
melted_merged_sdy368 <- melted_merged_sdy368[order(melted_merged_sdy368$study_time_collected),]
melted_merged_sdy368$study_time_collected <- factor(melted_merged_sdy368$study_time_collected,levels = c("Day_0","Day_1","Day_7","Day_28"))
quartz()
ggplot(melted_merged_sdy368,aes(x=study_time_collected,y=value,color=responder0))+
  geom_line(aes(group=subject_accession))+facet_wrap(~variable,scales = "free")+
  ggtitle("SDY368")+ylab("Cell Proportion")+labs(color='Responder') 

mod11_sdy368 <- lmer(CD4p_Tcells~study_time_collected+responder0+responder0*study_time_collected+(1|subject_accession),data = merged_sdy368)
em11_sdy368 <- emmeans(mod11_sdy368,~study_time_collected|responder0)
pv11_sdy368 <- contrast(em11_sdy368,alpha=0.05,interaction = "trt.vs.ctrl1")
em12_sdy368 <- emmeans(mod11_sdy368,~responder0|study_time_collected)
pv12_sdy368 <- contrast(em12_sdy368,alpha=0.05,interaction = "trt.vs.ctrl1")
em13_sdy368 <- emmeans(mod11_sdy368,~responder0*study_time_collected)
pv13_sdy368 <- contrast(em13_sdy368,alpha=0.05,interaction = "trt.vs.ctrl1")
em14_sdy368 <- emmeans(mod11_sdy368,~study_time_collected|responder0)
p <- plot(em14_sdy368)
p+coord_flip()

mod21_sdy368 <- lmer(CD8p_Tcells~study_time_collected+responder0+responder0*study_time_collected+(1|subject_accession),data = merged_sdy368)
em21_sdy368 <- emmeans(mod21_sdy368,~study_time_collected|responder0)
pv21_sdy368 <- contrast(em21_sdy368,alpha=0.05,interaction = "trt.vs.ctrl1")
em22_sdy368 <- emmeans(mod21_sdy368,~responder0|study_time_collected)
pv22_sdy368 <- contrast(em22_sdy368,alpha=0.05,interaction = "trt.vs.ctrl1")
em23_sdy368 <- emmeans(mod21_sdy368,~responder0*study_time_collected)
pv23_sdy368 <- contrast(em23_sdy368,alpha=0.05,interaction = "trt.vs.ctrl1")
em24_sdy368 <- emmeans(mod21_sdy368,~study_time_collected|responder0)
p <- plot(em24_sdy368)
p+coord_flip()

mod11_sdy368.a <- lmer(CD4p_CD45RAn_CCR7n~study_time_collected+responder0+responder0*study_time_collected+(1|subject_accession),data = merged_sdy368)
em11_sdy368.a <- emmeans(mod11_sdy368.a,~study_time_collected|responder0)
pv11_sdy368.a <- contrast(em11_sdy368.a,alpha=0.05,interaction = "trt.vs.ctrl1")
em12_sdy368.a <- emmeans(mod11_sdy368.a,~responder0|study_time_collected)
pv12_sdy368.a <- contrast(em12_sdy368.a,alpha=0.05,interaction = "trt.vs.ctrl1")
em13_sdy368.a <- emmeans(mod11_sdy368.a,~responder0*study_time_collected)
pv13_sdy368.a <- contrast(em13_sdy368.a,alpha=0.05,interaction = "trt.vs.ctrl1")
em14_sdy368.a <- emmeans(mod11_sdy368.a,~study_time_collected|responder0)
p <- plot(em14_sdy368.a)
p+coord_flip()


mod11_sdy368.b <- lmer(CD4p_CD45RAp_CCR7n~study_time_collected+responder0+responder0*study_time_collected+(1|subject_accession),data = merged_sdy368)
em11_sdy368.b <- emmeans(mod11_sdy368.b,~study_time_collected|responder0)
pv11_sdy368.b <- contrast(em11_sdy368.b,alpha=0.05,interaction = "trt.vs.ctrl1")
em12_sdy368.b <- emmeans(mod11_sdy368.b,~responder0|study_time_collected)
pv12_sdy368.b <- contrast(em12_sdy368.b,alpha=0.05,interaction = "trt.vs.ctrl1")
em13_sdy368.b <- emmeans(mod11_sdy368.b,~responder0*study_time_collected)
pv13_sdy368.b <- contrast(em13_sdy368.b,alpha=0.05,interaction = "trt.vs.ctrl1")
em14_sdy368.b <- emmeans(mod11_sdy368.b,~study_time_collected|responder0)
p <- plot(em14_sdy368.b)
p+coord_flip()


mod11_sdy368.c <- lmer(CD4p_CD45RAn_CCR7p~study_time_collected+responder0+responder0*study_time_collected+(1|subject_accession),data = merged_sdy368)
em11_sdy368.c <- emmeans(mod11_sdy368.c,~study_time_collected|responder0)
pv11_sdy368.c <- contrast(em11_sdy368.c,alpha=0.05,interaction = "trt.vs.ctrl1")
em12_sdy368.c <- emmeans(mod11_sdy368.c,~responder0|study_time_collected)
pv12_sdy368.c <- contrast(em12_sdy368.c,alpha=0.05,interaction = "trt.vs.ctrl1")
em13_sdy368.c <- emmeans(mod11_sdy368.c,~responder0*study_time_collected)
pv13_sdy368.c <- contrast(em13_sdy368.c,alpha=0.05,interaction = "trt.vs.ctrl1")
em14_sdy368.c <- emmeans(mod11_sdy368.c,~study_time_collected|responder0)
p <- plot(em14_sdy368.c)
p+coord_flip()


mod11_sdy368.d <- lmer(CD4p_CD45RAp_CCR7p~study_time_collected+responder0+responder0*study_time_collected+(1|subject_accession),data = merged_sdy368)
em11_sdy368.d <- emmeans(mod11_sdy368.d,~study_time_collected|responder0)
pv11_sdy368.d <- contrast(em11_sdy368.d,alpha=0.05,interaction = "trt.vs.ctrl1")
em12_sdy368.d <- emmeans(mod11_sdy368.d,~responder0|study_time_collected)
pv12_sdy368.d <- contrast(em12_sdy368.d,alpha=0.05,interaction = "trt.vs.ctrl1")
em13_sdy368.d <- emmeans(mod11_sdy368.d,~responder0*study_time_collected)
pv13_sdy368.d <- contrast(em13_sdy368.d,alpha=0.05,interaction = "trt.vs.ctrl1")
em14_sdy368.d <- emmeans(mod11_sdy368.d,~study_time_collected|responder0)
p <- plot(em14_sdy368.d)
p+coord_flip()


mod21_sdy368.a <- lmer(CD8p_CD45RAn_CCR7n~study_time_collected+responder0+responder0*study_time_collected+(1|subject_accession),data = merged_sdy368)
em21_sdy368.a <- emmeans(mod21_sdy368.a,~study_time_collected|responder0)
pv21_sdy368.a <- contrast(em21_sdy368.a,alpha=0.05,interaction = "trt.vs.ctrl1")
em22_sdy368.a <- emmeans(mod21_sdy368.a,~responder0|study_time_collected)
pv22_sdy368.a <- contrast(em22_sdy368.a,alpha=0.05,interaction = "trt.vs.ctrl1")
em23_sdy368.a <- emmeans(mod21_sdy368.a,~responder0*study_time_collected)
pv23_sdy368.a <- contrast(em23_sdy368.a,alpha=0.05,interaction = "trt.vs.ctrl1")
em24_sdy368.a <- emmeans(mod21_sdy368.a,~study_time_collected|responder0)
p <- plot(em24_sdy368.a)
p+coord_flip()


mod21_sdy368.b <- lmer(CD8p_CD45RAp_CCR7n~study_time_collected+responder0+responder0*study_time_collected+(1|subject_accession),data = merged_sdy368)
em21_sdy368.b <- emmeans(mod21_sdy368.b,~study_time_collected|responder0)
pv21_sdy368.b <- contrast(em21_sdy368.b,alpha=0.05,interaction = "trt.vs.ctrl1")
em22_sdy368.b <- emmeans(mod21_sdy368.b,~responder0|study_time_collected)
pv22_sdy368.b <- contrast(em22_sdy368.b,alpha=0.05,interaction = "trt.vs.ctrl1")
em23_sdy368.b <- emmeans(mod21_sdy368.b,~responder0*study_time_collected)
pv23_sdy368.b <- contrast(em23_sdy368.b,alpha=0.05,interaction = "trt.vs.ctrl1")
em24_sdy368.b <- emmeans(mod21_sdy368.b,~study_time_collected|responder0)
p <- plot(em24_sdy368.b)
p+coord_flip()


mod21_sdy368.c <- lmer(CD8p_CD45RAn_CCR7p~study_time_collected+responder0+responder0*study_time_collected+(1|subject_accession),data = merged_sdy368)
em21_sdy368.c <- emmeans(mod21_sdy368.c,~study_time_collected|responder0)
pv21_sdy368.c <- contrast(em21_sdy368.c,alpha=0.05,interaction = "trt.vs.ctrl1")
em22_sdy368.c <- emmeans(mod21_sdy368.c,~responder0|study_time_collected)
pv22_sdy368.c <- contrast(em22_sdy368.c,alpha=0.05,interaction = "trt.vs.ctrl1")
em23_sdy368.c <- emmeans(mod21_sdy368.c,~responder0*study_time_collected)
pv23_sdy368.c <- contrast(em23_sdy368.c,alpha=0.05,interaction = "trt.vs.ctrl1")
em24_sdy368.c <- emmeans(mod21_sdy368.c,~study_time_collected|responder0)
p <- plot(em24_sdy368.c)
p+coord_flip()


mod21_sdy368.d <- lmer(CD8p_CD45RAp_CCR7p~study_time_collected+responder0+responder0*study_time_collected+(1|subject_accession),data = merged_sdy368)
em21_sdy368.d <- emmeans(mod21_sdy368.d,~study_time_collected|responder0)
pv21_sdy368.d <- contrast(em21_sdy368.d,alpha=0.05,interaction = "trt.vs.ctrl1")
em22_sdy368.d <- emmeans(mod21_sdy368.d,~responder0|study_time_collected)
pv22_sdy368.d <- contrast(em22_sdy368.d,alpha=0.05,interaction = "trt.vs.ctrl1")
em23_sdy368.d <- emmeans(mod21_sdy368.d,~responder0*study_time_collected)
pv23_sdy368.d <- contrast(em23_sdy368.d,alpha=0.05,interaction = "trt.vs.ctrl1")
em24_sdy368.d <- emmeans(mod21_sdy368.d,~study_time_collected|responder0)
p <- plot(em24_sdy368.d)
p+coord_flip()