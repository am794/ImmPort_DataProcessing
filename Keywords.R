###############################
# Read the FCS files and      #
# return .lst keyword files   #
###############################

library(flowCore)
#library(FCSTrans)
library(xfun)
library("openxlsx")
##enter FCS files path and the list of FCS files
raw_fcs_path <- "/Volumes/Samsung_T5/Immport_UH2/SDY180/FCS/"
files <- list.files(raw_fcs_path,pattern=".fcs")
#rm(key)

# Merge the header channel names and marker names
key <- c()
for(i in 1:length(files)){
if(file_ext(files[i]) == "fcs"){
  fcs_header <- read.FCS(paste0("/Volumes/Samsung_T5/Immport_UH2/SDY180/FCS/",files[i]))
  if(i == 1){
  key <- as.data.frame(cbind(paste(as.vector(fcs_header@parameters@data$name),collapse = "|"),
                paste(as.vector(fcs_header@parameters@data$desc),collapse = "|")))
  }else{
    key <- rbind(key,as.data.frame(cbind(paste(as.vector(fcs_header@parameters@data$name),collapse = "|"),
                                         paste(as.vector(fcs_header@parameters@data$desc),collapse = "|"))))
  }
}
}

# Rename the row names and column names
rownames(key) <- files
colnames(key) <- c("Channel","Marker")
write.xlsx(key,"/Volumes/Samsung_T5/Immport_UH2/SDY180/Headers_SDY180.xlsx",row.names=TRUE)
