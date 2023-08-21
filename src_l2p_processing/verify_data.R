source("./Code_Transformation/download_tools.R")
key<-Sys.getenv("key")
report<-list()
currentdir <- getwd()
rds_output <- paste0(currentdir,"/rds_output")
rid="ri.foundry.main.dataset.677fa906-0a50-4a98-9517-e670eb06bd36"
report["var_Clean_Raw_Counts_1"]<-'no comparison'
try({
branch="R3_Archive_Draft_20230821_Updated"
var_Clean_Raw_Counts_1files<-pullnidap_raw(key=key,rid=rid,branch=branch)
var_Clean_Raw_Counts_1_target<-figure_out_nidap_files(var_Clean_Raw_Counts_1files)
var_Clean_Raw_Counts_1_new<-readRDS(paste0(rds_output,"/var_Clean_Raw_Counts_1.rds"))
report["var_Clean_Raw_Counts_1"]<-report_differences(var_Clean_Raw_Counts_1_target,var_Clean_Raw_Counts_1_new)
},silent=TRUE)
print(report["var_Clean_Raw_Counts_1"])
###################################
rid="ri.foundry.main.dataset.fd9859ff-8dbd-4337-8082-4b0c4f0ef337"
report["var_DEG_Analysis_1"]<-'no comparison'
try({
branch="R3_Archive_Draft_20230821_Updated"
var_DEG_Analysis_1files<-pullnidap_raw(key=key,rid=rid,branch=branch)
var_DEG_Analysis_1_target<-figure_out_nidap_files(var_DEG_Analysis_1files)
var_DEG_Analysis_1_new<-readRDS(paste0(rds_output,"/var_DEG_Analysis_1.rds"))
report["var_DEG_Analysis_1"]<-report_differences(var_DEG_Analysis_1_target,var_DEG_Analysis_1_new)
},silent=TRUE)
print(report["var_DEG_Analysis_1"])
###################################
rid="ri.foundry.main.dataset.919bad60-3d70-4f8e-994e-f9937c0d156d"
report["var_Filtered_Counts_1"]<-'no comparison'
try({
branch="R3_Archive_Draft_20230821_Updated"
var_Filtered_Counts_1files<-pullnidap_raw(key=key,rid=rid,branch=branch)
var_Filtered_Counts_1_target<-figure_out_nidap_files(var_Filtered_Counts_1files)
var_Filtered_Counts_1_new<-readRDS(paste0(rds_output,"/var_Filtered_Counts_1.rds"))
report["var_Filtered_Counts_1"]<-report_differences(var_Filtered_Counts_1_target,var_Filtered_Counts_1_new)
},silent=TRUE)
print(report["var_Filtered_Counts_1"])
###################################
rid="ri.foundry.main.dataset.32433dfd-3ac6-4230-b864-01ddba459e4b"
report["var_Volcano_Enhanced_1"]<-'no comparison'
try({
branch="R3_Archive_Draft_20230821_Updated"
var_Volcano_Enhanced_1files<-pullnidap_raw(key=key,rid=rid,branch=branch)
var_Volcano_Enhanced_1_target<-figure_out_nidap_files(var_Volcano_Enhanced_1files)
var_Volcano_Enhanced_1_new<-readRDS(paste0(rds_output,"/var_Volcano_Enhanced_1.rds"))
report["var_Volcano_Enhanced_1"]<-report_differences(var_Volcano_Enhanced_1_target,var_Volcano_Enhanced_1_new)
},silent=TRUE)
print(report["var_Volcano_Enhanced_1"])
###################################
