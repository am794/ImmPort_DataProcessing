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
#setwd("/Volumes/Samsung_T5/Immport_UH2/SDY80/Study/")

sdy <- "SDY80"
##Set tab_dir to the folder where the study files are located
studies_dir <- file.path(paste0("/Volumes/Samsung_T5/Immport_UH2/"),sdy,"/Study")
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
                               fhd.file_info_id, fhd.fcs_file_name, fhd.panel_reported, 
                               exs.expsample_accession,
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

### Visit information
fcs_file_mapping <- dbGetQuery(sqlite_conn,"SELECT e.biosample_accession, e.study_time_collected, e.type, 
                               f.biosample_accession, g.expsample_accession, g.file_info_id, g.result_schema,
                               a.expsample_accession,a.file_info_id,a.fcs_file_name,  b.experiment_accession
                               FROM expsample_2_file_info g
                               INNER JOIN expsample_2_biosample f
                               ON f.expsample_accession=g.expsample_accession
                               INNER JOIN biosample e
                               ON e.biosample_accession=f.biosample_accession
                               INNER JOIN fcs_header a
                               ON a.file_info_id=g.file_info_id
                               INNER JOIN expsample c
                               ON c.expsample_Accession = a.expsample_Accession
                               INNER JOIN experiment b
                               ON b.experiment_accession=c.experiment_accession
                               WHERE e.type='Other'
                               AND e.study_time_collected >='0'
                               AND g.result_schema='FCM';")

### Visit information SDY80
fcs_file_mapping <- dbGetQuery(sqlite_conn,"SELECT h.name,a.fcs_file_name,a.expsample_accession, 
                               e.study_time_collected, e.type, 
                               f.biosample_accession,b.name,
                               g.file_info_id,b.measurement_technique
                               FROM expsample_2_file_info g
                               INNER JOIN expsample_2_biosample f
                               ON f.expsample_accession=g.expsample_accession
                               INNER JOIN biosample e
                               ON e.biosample_accession=f.biosample_accession
                               INNER JOIN fcs_header a
                               ON a.file_info_id=g.file_info_id
                               INNER JOIN expsample c
                               ON c.expsample_Accession = a.expsample_Accession
                               INNER JOIN file_info h
                               ON h.file_info_id=a.file_info_id
                               INNER JOIN experiment b
                               ON b.experiment_accession=c.experiment_accession;")

fcs_file_mapping <- dbGetQuery(sqlite_conn,"SELECT DISTINCT a.comments,a.parent_population_reported,a.population_defnition_reported,
                               a.population_name_reported,a.study_time_collected
                               FROM fcs_analyzed_result a
                               WHERE a.experiment_accession = 'EXP14117';")

write.csv(fcs_file_mapping,"/Volumes/Samsung_T5/Immport_UH2/SDY80/Bcell/Cell_definition_SDY80.csv",row.names = FALSE)

##SDY 224: Bcell      
fcs_file_mapping <- dbGetQuery(sqlite_conn,"SELECT DISTINCT a.expsample_accession, a.name, a.experiment_accession,
                               b.biosample_accession, b.expsample_accession, c.expsample_accession, c.file_info_id,
                               d.file_info_id, d.name, d.original_file_name, e.biosample_accession, e.study_time_collected,
                               e.subject_accession, e.subtype
                               FROM expsample a
                               INNER JOIN expsample_2_biosample b
                               ON a.expsample_accession = b.expsample_accession
                               INNER JOIN expsample_2_file_info c
                               ON c.expsample_accession = a.expsample_accession
                               INNER JOIN file_info d
                               ON d.file_info_id = c.file_info_id
                               INNER JOIN biosample e
                               ON e.biosample_accession = b.biosample_accession
                               WHERE a.experiment_accession = 'EXP13737'
                               AND d.name LIKE '%.fcs%'
                               AND d.name NOT LIKE '%Control%'
                               AND d.name NOT LIKE 'Day%'
                               AND d.name LIKE 'L%';")

##SDY 224: B_Tcell      
fcs_file_mapping <- dbGetQuery(sqlite_conn,"SELECT DISTINCT a.expsample_accession, a.name, a.experiment_accession,
                               b.biosample_accession, b.expsample_accession, c.expsample_accession, c.file_info_id,
                               d.file_info_id, d.name, d.original_file_name, e.biosample_accession, e.study_time_collected,
                               e.subject_accession, e.subtype
                               FROM expsample a
                               INNER JOIN expsample_2_biosample b
                               ON a.expsample_accession = b.expsample_accession
                               INNER JOIN expsample_2_file_info c
                               ON c.expsample_accession = a.expsample_accession
                               INNER JOIN file_info d
                               ON d.file_info_id = c.file_info_id
                               INNER JOIN biosample e
                               ON e.biosample_accession = b.biosample_accession
                               WHERE a.experiment_accession = 'EXP10665'
                               AND d.name LIKE '%.fcs%'
                               AND d.name NOT LIKE '%Control%'
                               AND d.name NOT LIKE 'Day%'
                               AND d.name NOT LIKE 'XL%';")

write.xlsx(fcs_file_mapping,paste0("/Volumes/Samsung_T5/Immport_UH2/SDY224/B_Tcell/SDY224_BTcell_vaccine_info.xlsx"))
write.xlsx(fcs_file_mapping,paste0("/Volumes/Samsung_T5/Immport_UH2/",sdy,"/Study/Tab/SDY180_panel_mapping_header_text.xlsx"))
write.table(fcs_file_mapping,paste0("/Volumes/Samsung_T5/Immport_UH2/",
                                    sdy,"/Study/Tab/SDY180_panel_mapping_header_text.txt"),
                                    sep="\t",quote=F)
write.xlsx(fcs_file_mapping,"/Volumes/Samsung_T5/Immport_UH2/SDY144/Study/Tab/SDY144_panel_mapping_header_text.xlsx")

