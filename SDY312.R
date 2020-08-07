#############################################
# Building a local sqlite database instance #
# Quering the sqlite database               #
#############################################

library(RImmPort)
library(DBI)
#library(sqldf)
library(plyr)
#library("DBI")
#library(RMySQL)
library(RSQLite)
library("openxlsx")
#setwd("/Volumes/Samsung_T5/Immport_UH2/SDY416/Study/")

sdy <- "SDY312"
##Set tab_dir to the folder where the study files are located
studies_dir <- file.path(paste0("/Volumes/Samsung_T5/Immport_UH2/Companion_studies_3"),sdy,"Study")
tab_dir <- file.path(studies_dir,"Tab")
list.files(tab_dir)

##Set db dir to the folder where the database file "ImmPort.sqlite" should be stored
db_dir <- file.path(studies_dir,"Db")

##Build a local SQLite ImmPort database instance
buildNewSqliteDb(tab_dir, db_dir)
list.files(db_dir)

######
sqlite_conn <- dbConnect(SQLite(), dbname=file.path(db_dir, "ImmPort.sqlite"))
setImmPortDataSource(sqlite_conn)
getListOfStudies()


fcs_meta_Data <- dbGetQuery(sqlite_conn,"SELECT a.name, a.file_info_id, b.expsample_accession, b.file_info_id,
                            c.biosample_accession, c.expsample_accession, 
                            d.biosample_accession, d.study_time_collected, d.subject_accession,d.type,
                            e.subject_accession, e.max_subject_age
                            FROM file_info a
                            INNER JOIN expsample_2_file_info b
                            ON a.file_info_id = b.file_info_id
                            INNER JOIN expsample_2_biosample c
                            ON c.expsample_accession = b.expsample_accession
                            INNER JOIN biosample d
                            ON d.biosample_accession = c.biosample_accession
                            INNER JOIN arm_2_subject e
                            ON e.subject_accession = d.subject_accession
                            AND a.name LIKE '%.fcs%'
                            AND a.name NOT LIKE '%Compensation%';")
header_info <- read.xlsx("/Volumes/Samsung_T5/Immport_UH2/Companion_studies_3/SDY312/Header_SDY312.xlsx")                            
colnames(header_info)[1] <- "name" 

mm <- merge(fcs_meta_Data, header_info, by="name")
final_tab <- mm[which(mm$Panel == 3),c(1,2,3,5,8:10,12)]
write.csv(final_tab,"/Volumes/Samsung_T5/Immport_UH2/Companion_studies_3/SDY312/Tregs/MetaData_Tregs_SDY312.csv",row.names=FALSE)
ab_titer <- read.delim("/Volumes/Samsung_T5/Immport_UH2/Companion_studies_3/SDY312/SDY312-DR31_Tab/Tab/hai_result.txt",stringsAsFactors = FALSE,sep="\t",header=TRUE)
#ab_titer <- ab_titer[which(ab_titer$STUDY_TIME_COLLECTED == '28'),]
#ab_titer_d <- dcast(ab_titer,SUBJECT_ACCESSION~VIRUS_STRAIN_PREFERRED,value.var='VALUE_REPORTED')
#write.table(ab_titer_d,"/Volumes/Samsung_T5/Immport_UH2/Companion_studies_3/SDY312/Tregs/ab_titer_SDY312_D28.txt",sep = "\t",quote=FALSE,row.names = FALSE)
ab_titer_d <- dcast(ab_titer,SUBJECT_ACCESSION+STUDY_TIME_COLLECTED~VIRUS_STRAIN_PREFERRED,value.var='VALUE_REPORTED')
ab_titer_d$Geometric_mean <- apply(ab_titer_d[,3:9],1,function(x) geometric.mean(x,na.rm=TRUE))
resp <- dcast(ab_titer_d,SUBJECT_ACCESSION~STUDY_TIME_COLLECTED,value.var='Geometric_mean')
resp$FoldChange <- resp$`28`/resp$`0`
#write.table(resp,"/Volumes/Samsung_T5/Immport_UH2/Companion_studies_3/SDY312/Tregs/FoldChange.txt",sep="\t",quote=FALSE)


#########################################
Tregs_SDY312 <- read.xlsx("/Volumes/Samsung_T5/Immport_UH2/Resutls_Analysis/Tregs/SDY312_Tregs_proportions.xlsx",sheet=4)
colnames(resp)[1] <- 'subject_accession'
final_table <- na.omit(dplyr::left_join(Tregs_SDY312,resp,by=c('subject_accession')))
final_table <- final_table %>% mutate(Response=cut(FoldChange, breaks=c(-Inf,4,Inf),labels = c("NonResponder","Responder")))
#write.xlsx(final_table,"/Volumes/Samsung_T5/Immport_UH2/Resutls_Analysis/Tregs/SDY312_Final_Table.xlsx")
cor_data <- as.data.frame(final_Table[,c(2,6,22,23,24)],na.rm=TRUE)
rownames(cor_data) <- final_Table$subject_accession
cor(cor_data,method="s")
quartz()
ggplot(final_table)+geom_point(aes(x=Tregs,y=max_subject_age,colour = Response))+ylab("Age")+
  scale_fill_discrete(name = "Response")+ggtitle("Tregs proportion with CD3 as parent")

model <- glm(Response~Tregs+max_subject_age, data=final_table, family="binomial")

