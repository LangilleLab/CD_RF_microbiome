# Loop through RF summary objects (.rds files) and make summary table for disease and treatment response.

setwd("/home/gavin/github_repos/CD_RF_microbiome/")

source("BISCUIT_utility_code.R")

rds_files_16S <- list.files("RF_RDS_output/16S", full.names = TRUE)
rds_files_16S_CLR <- list.files("RF_RDS_output/16S_unrarified_CLR", full.names = TRUE)
rds_files_mgs <- list.files("RF_RDS_output/mgs", full.names = TRUE)

rds_files_16S_summary <- rf_summary_from_rds_files(rds_files_16S)
rds_files_16S_CLR_summary <- rf_summary_from_rds_files(rds_files_16S_CLR)
rds_files_mgs_summary <- rf_summary_from_rds_files(rds_files_mgs)

rds_files_16S_mgs_summary <- rbind(rds_files_16S_summary, rds_files_mgs_summary)

write.table( x=rds_files_16S_mgs_summary , file="raw_summary_out/RF_16S_mgs_summary.txt", 
             row.names=FALSE , quote = FALSE , sep="\t" )

write.table( x=rds_files_16S_CLR_summary , file="raw_summary_out/RF_16S_CLR_summary.txt", 
             row.names=FALSE , quote = FALSE , sep="\t" )
