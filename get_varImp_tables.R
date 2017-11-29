# Output tables of varImp for significant models.

setwd("/home/gavin/github_repos/CD_RF_microbiome/")

source("BISCUIT_utility_code.R")

read_write_varImp("RF_RDS_output/16S/biscuit_16S_species_disease.rds", "raw_summary_out/varImp/disease/species_16S_disease_varImp.txt")
read_write_varImp("RF_RDS_output/16S/biscuit_16S_genus_disease.rds", "raw_summary_out/varImp/disease/genus_16S_disease_varImp.txt")
read_write_varImp("RF_RDS_output/16S/biscuit_16S_family_disease.rds", "raw_summary_out/varImp/disease/family_16S_disease_varImp.txt")
read_write_varImp("RF_RDS_output/16S/biscuit_16S_order_disease.rds", "raw_summary_out/varImp/disease/order_16S_disease_varImp.txt")
read_write_varImp("RF_RDS_output/16S/biscuit_16S_class_disease.rds", "raw_summary_out/varImp/disease/class_16S_disease_varImp.txt")
read_write_varImp("RF_RDS_output/16S/biscuit_16S_phylum_disease.rds", "raw_summary_out/varImp/disease/phylum_16S_disease_varImp.txt")
read_write_varImp("RF_RDS_output/16S/biscuit_16S_ko_disease.rds", "raw_summary_out/varImp/disease/KO_16S_disease_varImp.txt")

read_write_varImp("RF_RDS_output/mgs/biscuit_mgs_mgs_strain_disease_prep.rds", "raw_summary_out/varImp/disease/strain_mgs_disease_varImp.txt")
read_write_varImp("RF_RDS_output/mgs/biscuit_mgs_genus_disease_prep.rds", "raw_summary_out/varImp/disease/genus_mgs_disease_varImp.txt")
read_write_varImp("RF_RDS_output/mgs/biscuit_mgs_family_disease_prep.rds", "raw_summary_out/varImp/disease/family_mgs_disease_varImp.txt")
read_write_varImp("RF_RDS_output/mgs/biscuit_mgs_phylum_disease_prep.rds", "raw_summary_out/varImp/disease/phylum_mgs_disease_varImp.txt")
read_write_varImp("RF_RDS_output/mgs/biscuit_mgs_moduless_disease_prep.rds", "raw_summary_out/varImp/disease/module_mgs_disease_varImp.txt")


read_write_varImp("RF_RDS_output/16S/biscuit_16S_genus_response.rds", "raw_summary_out/varImp/response/genus_16S_response_varImp.txt")
read_write_varImp("RF_RDS_output/16S/biscuit_16S_family_response.rds", "raw_summary_out/varImp/response/family_16S_response_varImp.txt")
read_write_varImp("RF_RDS_output/16S/biscuit_16S_order_response.rds", "raw_summary_out/varImp/response/order_16S_response_varImp.txt")
read_write_varImp("RF_RDS_output/16S/biscuit_16S_class_response.rds", "raw_summary_out/varImp/response/class_16S_response_varImp.txt")

read_write_varImp("RF_RDS_output/mgs/biscuit_mgs_mgs_strain_response_prep.rds", "raw_summary_out/varImp/response/strain_mgs_response_varImp.txt")
read_write_varImp("RF_RDS_output/mgs/biscuit_mgs_genus_response_prep.rds", "raw_summary_out/varImp/response/genus_mgs_response_varImp.txt")
read_write_varImp("RF_RDS_output/mgs/biscuit_mgs_class_response_prep.rds", "raw_summary_out/varImp/response/class_mgs_response_varImp.txt")
read_write_varImp("RF_RDS_output/mgs/biscuit_mgs_ko_response_prep.rds", "raw_summary_out/varImp/response/KO_mgs_response_varImp.txt")
read_write_varImp("RF_RDS_output/mgs/biscuit_mgs_pathways_response_prep.rds", "raw_summary_out/varImp/response/pathway_mgs_response_varImp.txt")
