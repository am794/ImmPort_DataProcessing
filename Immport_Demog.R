###############################
# Demographics table from     #
# ImmPort data with TSV input #
###############################

#Load libraries
install.packages("scales")
library("ggplot2")
library("reshape")

setwd("/Volumes/Samsung_T5/Immport_UH2/SDY144/")
studies <- list.files("/Volumes/Samsung_T5/Immport_UH2/",pattern="SDY*")
#tab <- list.files(paste0("/Volumes/Samsung_T5/Immport_UH2/",studies[i],"/"),pattern="SDY*")[3]
tab <- list.files(paste0("/Volumes/Samsung_T5/Immport_UH2/",studies[i],"/study/Tab"),,pattern=".zip")
## Function for merging the demographics data
subject <- read.table(paste0("/Volumes/Samsung_T5/Immport_UH2/",studies[i],"/",tab,"/Tab/subject.txt"),
                      header=TRUE,sep="\t",stringsAsFactors = FALSE)

arm_2_subject <- read.table(paste0("/Volumes/Samsung_T5/Immport_UH2/",studies[i],"/",tab,"/Tab/arm_2_subject.txt"),
                            header=TRUE,sep = "\t",stringsAsFactors = FALSE)

arm_or_cohort <- read.table(paste0("/Volumes/Samsung_T5/Immport_UH2/",studies[i],"/",tab,"/Tab/arm_or_cohort.txt"),
                            header=TRUE,sep = "\t",stringsAsFactors = FALSE)

merged1 <- merge(arm_or_cohort, arm_2_subject, by=intersect(colnames(arm_or_cohort),colnames(arm_2_subject)),all=TRUE)

subjects.merged <- merge(merged1, subject, by=c("SUBJECT_ACCESSION","WORKSPACE_ID"),all=TRUE)
subjects.merged <- subjects.merged[,c('SUBJECT_ACCESSION', 'ARM_ACCESSION', 'DESCRIPTION.x', 'GENDER', 'RACE', 'MIN_SUBJECT_AGE')]
subjects.merged <- cbind(subjects.merged,ARM_NAME=NA)
subjects.merged[which(subjects.merged$ARM_ACCESSION == "ARM894"),"ARM_NAME"] <- "Young"
subjects.merged[which(subjects.merged$ARM_ACCESSION == "ARM895"),"ARM_NAME"] <- "Old"

hai_result <- read.table(paste0("/Volumes/Samsung_T5/Immport_UH2/",studies[i],"/",tab,"/Tab/hai_result.txt"),
                            header=TRUE,sep = "\t",stringsAsFactors = FALSE)
hai.result.subset <- hai_result[,c("ARM_ACCESSION", "SUBJECT_ACCESSION", "STUDY_TIME_COLLECTED", "VALUE_PREFERRED", "VIRUS_STRAIN_REPORTED")]
hai.result.merged <- merge(hai.result.subset, subjects.merged, by=c("SUBJECT_ACCESSION","ARM_ACCESSION"))
hai.results <- cast(hai.result.merged, SUBJECT_ACCESSION + MIN_SUBJECT_AGE + GENDER + RACE + VIRUS_STRAIN_REPORTED + ARM_ACCESSION + ARM_NAME ~ STUDY_TIME_COLLECTED, value='VALUE_PREFERRED')
colnames(hai.results)[which(colnames(hai.results)=="0")] <- "Day0"
colnames(hai.results)[which(colnames(hai.results)=="28")] <- "Day28"
#hai.results.SDY <- cbind(hai.results, fold.change= hai.results[,"Day28"] / hai.results[,"Day0"])


####Mapping the data for input FCS files
experiment <- read.table(paste0("/Volumes/Samsung_T5/Immport_UH2/",studies[i],"/",tab,"/Tab/experiment.txt"),
                         header=TRUE,sep="\t",stringsAsFactors = FALSE)
experiment <- experiment[which(experiment$MEASUREMENT_TECHNIQUE == "Flow Cytometry"),]

expsample <- read.table(paste0("/Volumes/Samsung_T5/Immport_UH2/",studies[i],"/",tab,"/Tab/expsample.txt"),
                        header=TRUE,sep="\t",stringsAsFactors = FALSE)
#file_info <- read.table(paste0("/Volumes/Samsung_T5/Immport_UH2/",studies[i],"/",tab,"/Tab/file_info.txt"),
                        #header=TRUE,sep="\t",stringsAsFactors = FALSE)
#file_info <- file_info[which(file_info$DETAIL == "Flow cytometry result"),]
#expsample_2_file_info <- read.table(paste0("/Volumes/Samsung_T5/Immport_UH2/",studies[i],"/",tab,"/Tab/expsample_2_file_info.txt"),
                                    #header=TRUE,sep="\t",stringsAsFactors = FALSE)
#expsample_2_biosample <- read.table(paste0("/Volumes/Samsung_T5/Immport_UH2/",studies[i],"/",tab,"/Tab/expsample_2_biosample.txt"),
                                    #header=TRUE,sep="\t",stringsAsFactors = FALSE)

biosample <- read.table(paste0("/Volumes/Samsung_T5/Immport_UH2/",studies[i],"/",tab,"/Tab/biosample.txt"),
                        header=TRUE,sep="\t",stringsAsFactors = FALSE)
subject <- read.table(paste0("/Volumes/Samsung_T5/Immport_UH2/",studies[i],"/",tab,"/Tab/subject.txt"),
                      header=TRUE,sep="\t",stringsAsFactors = FALSE)
reagent <- read.table(paste0("/Volumes/Samsung_T5/Immport_UH2/",studies[i],"/",tab,"/Tab/reagent.txt"),
                      header=TRUE,sep="\t",stringsAsFactors = FALSE)
reagent_set_2_reagent <- read.table(paste0("/Volumes/Samsung_T5/Immport_UH2/",studies[i],"/",tab,"/Tab/reagent_set_2_reagent.txt"),
                                    header=TRUE,sep="\t",stringsAsFactors = FALSE)
fcs_analyzed_result <- read.table(paste0("/Volumes/Samsung_T5/Immport_UH2/",studies[i],"/",tab,"/Tab/fcs_analyzed_result.txt"),
           header=TRUE,sep="\t",stringsAsFactors = FALSE,na.strings = "NA",fill = TRUE)
study <- read.table(paste0("/Volumes/Samsung_T5/Immport_UH2/",studies[i],"/",tab,"/Tab/study.txt"),
                    header=TRUE,sep="\t",stringsAsFactors = FALSE)
expsample_2_reagent <- read.table(paste0("/Volumes/Samsung_T5/Immport_UH2/",studies[i],"/",tab,"/Tab/expsample_2_reagent.txt"),
                                  header=TRUE,sep="\t",stringsAsFactors = FALSE)
fcs_reagent_expsample_merged <- merge()
  
  
exp_merged <- merge(experiment,expsample,by=("EXPERIMENT_ACCESSION"),all=F)
expsample_merged <- merge(exp_merged,expsample_2_file_info,by=("EXPSAMPLE_ACCESSION"),all=F)
expsample_merged <- expsample_merged[which(expsample_merged$EXPERIMENT_ACCESSION == intersect(experiment$EXPERIMENT_ACCESSION,expsample$EXPERIMENT_ACCESSION)),]
#exp_fcs_analyzed_merged <- merge(expsample_merged,fcs_analyzed_result,by=("EXPSAMPLE_ACCESSION"),all = F) 
#exp_fcs_biosample <- merge()


exp_file_merged <- merge(expsample_merged,file_info,by=("FILE_INFO_ID"),all=T)
exp_file_merged <- exp_file_merged[]
expsample_2_bio_merged <- merge(exp_file_merged,expsample_2_biosample,by=("EXPSAMPLE_ACCESSION"),no.dups=T)
biosample_merged <- merge(expsample_2_bio_merged,biosample,by=("BIOSAMPLE_ACCESSION"),no.dups = T)
subject_biosample_merged <- merge(subject,biosample_merged,by=("SUBJECT_ACCESSION"),no.dups = T)
expsample_merged <- merge(subject_biosample_merged,expsample_2_reagent,by=("EXPSAMPLE_ACCESSION"),no.dups = T)
reagent_2_expsample_merged <- merge(expsample_merged ,expsample_2_reagent,by=("REAGENT_ACCESSION"),no.dups = T)
fcs_reagent_exp_merged <- merge(reagent_2_expsample_merged,fcs_analyzed_result,by=("SUBJECT_ACCESSION"),no.dups = T)


####RImmport : MySQL
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("RImmPort", version = "3.8")
library(RImmPort)
library("DBI")
#library(RMySQL)
library(RSQLite)

#buildNewSqliteDb(paste0("/Volumes/Samsung_T5/Immport_UH2/",studies[i],"/",tab,"/load"), paste0("/Volumes/Samsung_T5/Immport_UH2/",studies[i],"/",tab,"/SDY144_SQLite/"))
sqlite    <- dbDriver("SQLite")
exampledb <- dbConnect(sqlite,"/Volumes/Samsung_T5/Immport_UH2/SDY144/SDY144-DR28_Tab/SDY144_SQLite/ImmPort.sql")
dbListTables(exampledb)
dbListTables(exampledb,'study_file')
dbListFields(exampledb,'file_info')
dbSendQuery(exampledb,"select SUBJECT_ACCESSION from subject")
exampledb <- dbConnect(SQLite(),"/Volumes/Samsung_T5/Immport_UH2/SDY144/SDY144-DR28_Tab/SDY144_SQLite/ImmPort.sql")


##### RImmport: MySQL
mysql_conn <- dbConnect(MySQL(), user="root", password="G0davar!",host="localhost",port=3306)
setImmPortDataSource(mysql_conn)

# studies_dir <- system.file("extdata", "ImmPortStudies", package = "RImmPort") ##An example to load the Tab folder from RImmPort function
studies_dir <- system.file("/Volumes/Samsung_T5/Immport_UH2/SDY144/Study/Tab",mustWork = T)
studies_dir <- file.path("/Volumes/Samsung_T5/Immport_UH2/SDY144/Study")
tab_dir <- file.path("/Volumes/Samsung_T5/Immport_UH2/SDY144/Study/Tab/")
list.files(tab_dir)
db_dir <- file.path(studies_dir,"Db")
buildNewSqliteDb(tab_dir, db_dir)
list.files(db_dir)

sqlite_conn <- dbConnect(SQLite(), dbname=file.path(db_dir, "ImmPort.sqlite"))
setImmPortDataSource(sqlite_conn)
study_id <- 'SDY144'
sdy144 <- getStudy(study_id)

#demographics of the data
dm_df <- sdy144$special_purpose$dm_l$dm_df
head(dm_df)

#concomitant medications
cm_df <- sdy144$interventions$cm_l$cm_df
head(cm_df)

#Get a list of domains
getListOfDomains()
?"Demographics Domain"
# get domain code of Demographics domain
domain_name <- "Demographics"
getDomainCode(domain_name)
dm_l <- getDomainDataOfStudies(domain_name, "SDY144")
if (length(dm_l) > 0)
  names(dm_l)
head(dm_l$dm_df)

## Search and integrate functions 
domain_name <- "Cellular Quantification"
study_ids_l <- getStudiesWithSpecificDomainData(domain_name)
study_ids_l
domain_name <- "Cellular Quantification"
getDomainCode(domain_name)
domain_name <- "Cellular Quantification"
study_ids <- "SDY144"
zb_l <- getDomainDataOfStudies(domain_name, study_ids)
if (length(zb_l) > 0)
  names(zb_l)
head(zb_l$zb_df)

##FCM data
assay_type <- "Flow"
study_id = "SDY144"
flow_l <- getAssayDataOfStudies(study_id, assay_type)
if (length(flow_l) > 0)
  names(flow_l)
head(flow_l$zb_df)

dbGetQuery(sqlite_conn ,"SELECT sub.subject_accession sub.ethnicity sub.gender sub.race bio.biosample_accession bio.study_accession  
           FROM subject sub, biosample bio;
           AND bio.subject_accession = sub.subject_accession ;")

hai.data.SDY144 <- dbGetQuery(sqlite_conn,"SELECT hr.subject_accession, hr.arm_accession, hr.study_time_collected, hr.study_time_collected_unit, hr.value_preferred, hr.virus_strain_reported, s.race, s.gender, a2s.min_subject_age, a2s.age_unit
    FROM hai_result hr, subject s, arm_2_subject a2s
                WHERE hr.study_accession = 'SDY144'
                AND hr.subject_accession = s.subject_accession
    AND a2s.subject_accession = s.subject_accession;")

hai.data.SDY144 <- dbGetQuery(sqlite_conn,"SELECT fcs.biosample_accession, fcs.subject_accession, fcs.population_name_reported, fcs.experiment_accession,  b.biosample_accession,b.subject_accession,b.type
    FROM fcs_analyzed_result fcs, biosample b
                WHERE fcs.study_accession = 'SDY144'
                AND fcs.subject_accession = b.subject_accession
                AND fcs.biosample_accession = b.biosample_accession;")

##
fcs_subject_biosample <- dbGetQuery(sqlite_conn,"SELECT fcs.result_id, fcs.biosample_accession, fcs.arm_accession, fcs.EXPERIMENT_ACCESSION, fcs.EXPSAMPLE_ACCESSION, 
                                    fcs.POPULATION_NAME_REPORTED, fcs.POPULATION_NAME_REPORTED, fcs.STUDY_ACCESSION, fcs.SUBJECT_ACCESSION, s.SUBJECT_ACCESSION,
                                    s.ETHNICITY, s.GENDER, s.RACE, b.BIOSAMPLE_ACCESSION, b.SUBJECT_ACCESSION, b.STUDY_ACCESSION, h.fcs_file_name, h.fcs_header_id,
                                    h.fcs_header_text, h.fcs_file_name, h.file_info_id
                                    FROM fcs_analyzed_result fcs, subject s, biosample b, fcs_header h
                                    WHERE fcs.study_accession = 'SDY144'
                                    AND fcs.study_accession = b.study_accession
                                    AND fcs.biosample_accession = b.biosample_accession
                                    AND fcs.subject_accession = s.subject_accession
                                    AND fcs.expsample_accession = h.expsample_accession;")

## File mapping and arranging based on the reagent panel
fcs_file_mapping <- dbGetQuery(sqlite_conn,"SELECT fin.file_info_id, fin.detail, fin.name, 
                              fhd.file_info_id, fhd.fcs_file_name, fhd.expsample_accession, exs.expsample_accession,
                              exs.experiment_accession, exp.experiment_accession
                              FROM file_info fin, fcs_header fhd, expsample exs, experiment exp
                              WHERE fin.detail = 'Flow cytometry result'
                              AND fin.file_info_id = fhd.file_info_id
                              AND exs.expsample_accession = fhd.expsample_accession
                              AND exp.experiment_accession = exs.experiment_accession;")

## File mapping and arranging based on the reagent panel
fcs_file_mapping <- dbGetQuery(sqlite_conn,"SELECT fin.file_info_id, fin.detail, fin.name, 
                               fhd.file_info_id, fhd.fcs_file_name, exs.expsample_accession,
                               exs.experiment_accession, exp.description
                               FROM file_info fin, fcs_header fhd, expsample exs, experiment exp
                               WHERE fin.detail = 'Flow cytometry result'
                               AND fin.file_info_id = fhd.file_info_id
                               AND exs.expsample_accession = fhd.expsample_accession
                               AND exp.experiment_accession = exs.experiment_accession;")

fcs_file_mapping <- dbGetQuery(sqlite_conn,"SELECT a.expsample_accession, a.file_info_id, a.fcs_file_name,
                               b.experiment_accession, b.description,c.expsample_accession, c.experiment_accession
                               FROM fcs_header a
                               INNER JOIN expsample c
                               ON a.expsample_accession=c.expsample_accession
                               INNER JOIN experiment b
                               ON b.experiment_accession=c.experiment_accession;")
                               
          
###Merging
fcs_file_mapping <- dbGetQuery(sqlite_conn,"SELECT a.expsample_accession, a.file_info_id, a.fcs_file_name,
                               b.experiment_accession, b.measurement_technique,b.name,c.expsample_accession, 
                               c.experiment_accession, d.file_info_id, d.name
                               FROM expsample c
                               INNER JOIN fcs_header a
                               ON c.expsample_accession=a.expsample_accession
                               INNER JOIN experiment b
                               ON c.experiment_accession=b.experiment_accession
                               INNER JOIN file_info d
                               ON d.file_info_id=a.file_info_id;")

fcs_file_mapping <- dbGetQuery(sqlite_conn,"SELECT a.expsample_accession, a.file_info_id, a.fcs_file_name,
                               d.file_info_id, d.name
                               FROM fcs_header a
                               INNER JOIN file_info d
                               ON d.file_info_id=a.file_info_id
                               WHERE a.fcs_version = 3;")

fcs_file_mapping <- dbGetQuery(sqlite_conn,"SELECT a.expsample_accession, a.file_info_id, a.fcs_file_name,a.fcs_header_id,
                               d.file_info_id, d.name,c.expsample_accession,b.experiment_accession
                               FROM fcs_header a
                               INNER JOIN file_info d
                               ON d.file_info_id=a.file_info_id
                               INNER JOIN expsample c
                               ON c.expsample_accession=a.expsample_accession
                               INNER JOIN EXPERIment b
                               ON B.EXPERIMENT_ACCESSION=c.experiment_accession
                               WHERE c.result_schema='FCM'
                               AND a.fcs_version = 3
                               AND b.measurement_technique = 'Flow Cytometry';")
####Merging the files using sqlite
fcs_file_mapping <- dbGetQuery(sqlite_conn,"SELECT DISTINCT a.fcs_file_name,a.fcs_version, d.name,b.name,
                               a.file_info_id
                               FROM fcs_header a
                               INNER JOIN file_info d
                               ON d.file_info_id=a.file_info_id
                               INNER JOIN expsample c
                               ON c.expsample_accession=a.expsample_accession
                               INNER JOIN EXPERIment b
                               ON B.EXPERIMENT_ACCESSION=c.experiment_accession
                               WHERE c.result_schema='FCM'
                               AND b.measurement_technique = 'Flow Cytometry'
                               AND a.fcs_version=3
                               order BY a.fcs_file_name;")

fcs_file_mapping <- dbGetQuery(sqlite_conn,"SELECT DISTINCT a.fcs_file_name,a.fcs_version, d.name,b.name,
                               a.file_info_id,a.panel_reported,e.fcs_header_id,e.pnn_reported,e.pns_reported
                               FROM fcs_header a
                               INNER JOIN file_info d
                               ON d.file_info_id=a.file_info_id
                               INNER JOIN expsample c
                               ON c.expsample_accession=a.expsample_accession
                               INNER JOIN EXPERIment b
                               ON B.EXPERIMENT_ACCESSION=c.experiment_accession
                               INNER JOIN fcs_header_marker e
                               ON e.fcs_header_id = a.fcs_header_id
                               WHERE c.result_schema='FCM'
                               AND b.measurement_technique = 'Flow Cytometry'
                               order BY a.fcs_file_name;")

write.xlsx(fcs_file_mapping,"/Volumes/Samsung_T5/Immport_UH2/SDY144/Study/Tab/SDY144_panel_mapping_header_text.xlsx")
                              

fcs_file_mapping <- dbGetQuery(sqlite_conn,"SELECT DISTINCT a.fcs_file_name,a.fcs_version, d.name,b.name
                               FROM fcs_header a
                               INNER JOIN file_info d
                               ON d.file_info_id=a.file_info_id
                               INNER JOIN expsample c
                               ON c.expsample_accession=a.expsample_accession
                               INNER JOIN EXPERIment b
                               ON B.EXPERIMENT_ACCESSION=c.experiment_accession
                               WHERE c.result_schema='FCM'
                               AND b.measurement_technique = 'Flow Cytometry';")



INNER JOIN experiment_accession c
ON c.experiment_accession=a.experiment_accession



fcs_file_mapping <- dbGetQuery(sqlite_conn,"SELECT exp.experiment_accession, exp.description, fhd.expsample_accession,
                               exs.experiment_accession, exs.expsample_accession
                               FROM experiment exp
                               LEFT OUTER JOIN expsample as exs ON exs.experiment_accession = exp.experiment_accession
                               LEFT OUTER JOIN fcs_header as fhd ON fhd.expsample_accession = exs.expsample_accession;")


                               LEFT OUTER JOIN file_info as fin ON fin.file_info_id = fhd.file_info_id;")


                               FROM fcs_file_mapping as fhd 
                               LEFT OUTER JOIN experiment as exp ON
                               fhd.description = exp.description and fhd.expsample_accession = exp.expsample_accession
                               WHERE exs.expsample_accession = fhd.expsample_accession;")



