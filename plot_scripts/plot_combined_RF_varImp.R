setwd("/home/gavin/github_repos/CD_RF_microbiome/")

source("BISCUIT_utility_code.R")

combined_disease_RF_out <- readRDS("RF_RDS_output/combined/combined_disease.rds")

combined_disease_RF_out_varImp <- read_in_varImp("RF_RDS_output/combined/combined_disease.rds")

combined_response_RF_out <- readRDS("RF_RDS_output/combined/combined_response.rds")

combined_response_RF_out_varImp <- read_in_varImp("RF_RDS_output/combined/combined_response.rds")

par(mar=c(5.1,6.1,4.1,2.1) , mfrow=c(1,2))
plot_varImp(imp = combined_disease_RF_out_varImp, input = combined_disease_RF_out$loocv_mod$trainingData, classes = c("CD", "CN"), x_val=c(0,2.5))
plot_varImp(imp = combined_response_RF_out_varImp, input = combined_response_RF_out$loocv_mod$trainingData, classes = c("RS", "NR"), x_val=c(0,3.5))