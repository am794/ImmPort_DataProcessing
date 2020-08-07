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
#setwd("/Volumes/Samsung_T5/SDY67/Study/")

sdy <- "SDY739"
##Set tab_dir to the folder where the study files are located
studies_dir <- file.path(paste0("/Volumes/Samsung_T5/Immport_UH2/"),sdy,"Study")
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
                            d.biosample_accession, d.study_time_collected, d.subject_accession,
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
header_info <- read.xlsx("/Volumes/Samsung_T5/Immport_UH2/SDY739/Header_SDY739.xlsx",sheet=1)                            
colnames(header_info)[1] <- "name" 

mm <- merge(fcs_meta_Data, header_info, by="name")
final_tab <- mm[which(mm$Panel == 2),]
ab_titer <- read.delim("/Volumes/Samsung_T5/Immport_UH2/SDY739/Study/Tab/SDY739-DR32_Tab/Tab/neut_ab_titer_result.txt",stringsAsFactors = FALSE,sep="\t",header=TRUE)
#ab_titer <- ab_titer[which(ab_titer$STUDY_TIME_COLLECTED == '0' ),]
ab_titer <- ab_titer[which(ab_titer$STUDY_TIME_COLLECTED == '0' | ab_titer$STUDY_TIME_COLLECTED == '28'),]
ab_titer_d <- dcast(ab_titer,SUBJECT_ACCESSION+STUDY_TIME_COLLECTED~VIRUS_STRAIN_REPORTED,value.var='VALUE_REPORTED')
#write.table(ab_titer_d,"/Volumes/Samsung_T5/Immport_UH2/SDY739/Tregs/ab_titer_SDY739.txt",sep = "\t",quote=FALSE,row.names = FALSE)
ab_titer_d$Geometric_mean <- apply(ab_titer_d[,3:4],1,function(x) geometric.mean(x,na.rm=TRUE))
resp <- dcast(ab_titer_d,SUBJECT_ACCESSION~STUDY_TIME_COLLECTED,value.var='Geometric_mean')
resp$FoldChange <- resp$`28`/resp$`0`
#write.table(ab_titer_d,"/Volumes/Samsung_T5/Immport_UH2/Companion_studies_3/SDY312/Tregs/ab_titer_SDY312_D28.txt",sep = "\t",quote=FALSE,row.names = FALSE)

