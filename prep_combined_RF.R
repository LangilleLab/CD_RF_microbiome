# Prep input files for combined RF models (based on top 3 features from all significant models).

setwd("/home/gavin/github_repos/CD_RF_microbiome/")

source("BISCUIT_utility_code.R")

# Read in metadata file.
biscuit_meta <- read.table("biscuit_metadata.txt",
                           header=T,
                           row.names=1,
                           stringsAsFactors=FALSE,
                           sep="\t")

# Remove samples S34/S38 from metadata due to low DNA quality.
biscuit_meta <- biscuit_meta[-which(rownames(biscuit_meta) %in% c("S34", "S38")),]

# Read in alpha diversity and GRS:
input_alpha <- read.table("raw_data/BISCUIT_16S/alpha_div_stats_BISCUIT.txt", header=T, sep="\t", row.names=1)
input_grs <- read.table("raw_data/BISCUIT_MGS/mgs_genetic_risk_scores.txt", header=T, sep="\t", row.names=1)

sample_order <- rownames(input_alpha)
input_grs <- input_grs[sample_order,, drop=FALSE]

# Read in top 3 features from all significant disease datasets.
species_16S_disease_top3 <- get_top3_sig_varImp_vals("RF_RDS_output/16S/biscuit_16S_species_disease.rds", 
                                                   sample_order, 
                                                   fix_ids=TRUE, 
                                                   "~/projects/crohns_RF/biscuit2randomIDmapping.txt")

genus_16S_disease_top3 <- get_top3_sig_varImp_vals("RF_RDS_output/16S/biscuit_16S_genus_disease.rds", 
                                                   sample_order, 
                                                   fix_ids=TRUE, 
                                                   "~/projects/crohns_RF/biscuit2randomIDmapping.txt")

family_16S_disease_top3 <- get_top3_sig_varImp_vals("RF_RDS_output/16S/biscuit_16S_family_disease.rds", 
                                                   sample_order, 
                                                   fix_ids=TRUE, 
                                                   "~/projects/crohns_RF/biscuit2randomIDmapping.txt")

order_16S_disease_top3 <- get_top3_sig_varImp_vals("RF_RDS_output/16S/biscuit_16S_order_disease.rds", 
                                                    sample_order, 
                                                    fix_ids=TRUE, 
                                                    "~/projects/crohns_RF/biscuit2randomIDmapping.txt")

class_16S_disease_top3 <- get_top3_sig_varImp_vals("RF_RDS_output/16S/biscuit_16S_class_disease.rds", 
                                                   sample_order, 
                                                   fix_ids=TRUE, 
                                                   "~/projects/crohns_RF/biscuit2randomIDmapping.txt")

phylum_16S_disease_top3 <- get_top3_sig_varImp_vals("RF_RDS_output/16S/biscuit_16S_phylum_disease.rds", 
                                                    sample_order, 
                                                   fix_ids=TRUE, 
                                                   "~/projects/crohns_RF/biscuit2randomIDmapping.txt")

ko_16S_disease_top3 <- get_top3_sig_varImp_vals("RF_RDS_output/16S/biscuit_16S_ko_disease.rds", 
                                                    sample_order, 
                                                    fix_ids=TRUE, 
                                                    "~/projects/crohns_RF/biscuit2randomIDmapping.txt")

strain_mgs_disease_top3 <- get_top3_sig_varImp_vals("RF_RDS_output/mgs/biscuit_mgs_mgs_strain_disease_prep.rds", 
                                                sample_order)

genus_mgs_disease_top3 <- get_top3_sig_varImp_vals("RF_RDS_output/mgs/biscuit_mgs_genus_disease_prep.rds", 
                                                    sample_order)

family_mgs_disease_top3 <- get_top3_sig_varImp_vals("RF_RDS_output/mgs/biscuit_mgs_family_disease_prep.rds", 
                                                   sample_order)

phylum_mgs_disease_top3 <- get_top3_sig_varImp_vals("RF_RDS_output/mgs/biscuit_mgs_phylum_disease_prep.rds", 
                                                    sample_order)

module_mgs_disease_top3 <- get_top3_sig_varImp_vals("RF_RDS_output/mgs/biscuit_mgs_moduless_disease_prep.rds", 
                                                    sample_order)

combined_disease_prep <- cbind(input_alpha,
                               input_grs,
                               species_16S_disease_top3,
                               genus_16S_disease_top3,
                               family_16S_disease_top3,
                               order_16S_disease_top3,
                               class_16S_disease_top3,
                               phylum_16S_disease_top3,
                               ko_16S_disease_top3,
                               strain_mgs_disease_top3,
                               genus_mgs_disease_top3,
                               family_mgs_disease_top3,
                               phylum_mgs_disease_top3,
                               module_mgs_disease_top3)
old_tmp <- combined_disease_prep

colnames(combined_disease_prep) <- c("num_OTUs",                                                                                                                                                               
                                    "GRS",
                                    "g__Akkermansia.s__muciniphila",                                               
                                    "g__Desulfovibrio.Unclassified",                                             
                                    "g__Prevotella.s__copri",                                                                       
                                    "g__Desulfovibrio",                                                          
                                    "g__Akkermansia",                                                              
                                    "g__Butyricimonas",                                                                         
                                    "f__Verrucomicrobiaceae",                                                                             
                                    "f__Prevotellaceae",                                                                                              
                                    "f__Odoribacteraceae",                                                                                          
                                    "o__Verrucomicrobiales",                                                                                                    
                                    "o__Actinomycetales",                                                                                                          
                                    "o__Rhizobiales",                                                                                                         
                                    "c__Verrucomicrobiae",                                                                                                                          
                                    "c__Deltaproteobacteria",                                                                                                                        
                                    "c__Erysipelotrichi",                                                                                                                                
                                    "p__Verrucomicrobia",                                                                                                                                              
                                    "k__Bacteria.p__Actinobacteria",                                                                                                                                               
                                    "k__Bacteria.p__Cyanobacteria",                                                                                                                                                
                                    "K03785",                                                                                                                                                                      
                                    "K09013",                                                                                                                                                                      
                                    "K03809",                                                                                                                                                                      
                                    "s__Alistipes_putredinis.t__GCF_000154465_MGS",                                         
                                    "t__Clostridium_symbiosum_unclassified_MGS",                    
                                    "t__Faecalibacterium_prausnitzii_unclassified_MGS",
                                    "g__Alistipes_MGS",                                                                                  
                                    "g__Oscillibacter_MGS",                                                                               
                                    "g__Dorea_MGS",                                                                                        
                                    "f__Rikenellaceae_MGS",                                                                                               
                                    "f__Oscillospiraceae_MGS",                                                                                                
                                    "f__Porphyromonadaceae_MGS",                                                                                          
                                    "p__Verrucomicrobia_MGS",                                                                                                                                            
                                    "p__Firmicutes_MGS",                                                                                                                                                   
                                    "p__Bacteroidetes_MGS",                                                                                                                                                
                                    "M00144_MGS",                                                                                                                                                                      
                                    "M00362_MGS",                                                                                                                                                                      
                                    "M00239_MGS")

write_RF_prep(as.data.frame(t(combined_disease_prep)), "prep_data/combined_prep/combined_disease_prep.txt", "disease", biscuit_meta)


### Get combined response RF prep.
response_samples <- rownames(biscuit_meta)[which(biscuit_meta$disease == "CD")]

genus_16S_response_top3 <- get_top3_sig_varImp_vals("RF_RDS_output/16S/biscuit_16S_genus_response.rds", 
                                                   response_samples, 
                                                   fix_ids=TRUE, 
                                                   "~/projects/crohns_RF/biscuit2randomIDmapping.txt")

family_16S_response_top3 <- get_top3_sig_varImp_vals("RF_RDS_output/16S/biscuit_16S_family_response.rds", 
                                                    response_samples, 
                                                    fix_ids=TRUE, 
                                                    "~/projects/crohns_RF/biscuit2randomIDmapping.txt")

order_16S_response_top3 <- get_top3_sig_varImp_vals("RF_RDS_output/16S/biscuit_16S_order_response.rds", 
                                                   response_samples, 
                                                   fix_ids=TRUE, 
                                                   "~/projects/crohns_RF/biscuit2randomIDmapping.txt")

class_16S_response_top3 <- get_top3_sig_varImp_vals("RF_RDS_output/16S/biscuit_16S_class_response.rds", 
                                                   response_samples, 
                                                   fix_ids=TRUE, 
                                                   "~/projects/crohns_RF/biscuit2randomIDmapping.txt")


strain_mgs_response_top3 <- get_top3_sig_varImp_vals("RF_RDS_output/mgs/biscuit_mgs_mgs_strain_response_prep.rds", 
                                                     response_samples)

genus_mgs_response_top3 <- get_top3_sig_varImp_vals("RF_RDS_output/mgs/biscuit_mgs_genus_response_prep.rds", 
                                                    response_samples)

class_mgs_response_top3 <- get_top3_sig_varImp_vals("RF_RDS_output/mgs/biscuit_mgs_class_response_prep.rds", 
                                                    response_samples)

ko_mgs_response_top3 <- get_top3_sig_varImp_vals("RF_RDS_output/mgs/biscuit_mgs_ko_response_prep.rds", 
                                                 response_samples)

pathway_mgs_response_top3 <- get_top3_sig_varImp_vals("RF_RDS_output/mgs/biscuit_mgs_pathways_response_prep.rds", 
                                                 response_samples)
combined_response_prep <- cbind(genus_16S_response_top3,
                                family_16S_response_top3,
                                order_16S_response_top3,
                                class_16S_response_top3,
                                strain_mgs_response_top3,
                                genus_mgs_response_top3,
                                class_mgs_response_top3,
                                ko_mgs_response_top3,
                                pathway_mgs_response_top3)

colnames(combined_response_prep) <- c("g__Dialister",                                                                                   
                                     "g__Bilophila",                                                               
                                     "g__Aggregatibacter",                                                                 
                                     "f__Desulfovibrionaceae",                                                                            
                                     "f__Porphyromonadaceae",                                                                                           
                                     "o__Clostridiales.f__UN",                                                                                                        
                                     "o__Erysipelotrichales",                                                                                                           
                                     "o__Desulfovibrionales",                                                                                                   
                                     "o__Pseudomonadales",                                                                                                      
                                     "c__Erysipelotrichi",                                                                                                                                 
                                     "c__Deltaproteobacteria",                                                                                                                         
                                     "c__Flavobacteriia",                                                                                                                               
                                     "t__Parabacteroides_merdae_UN_MGS", 
                                     "t__Sutterella_wadsworthensis_UN_MGS",
                                     "t__UN_Lachnospiraceae_GCF_000209405_MGS",                
                                     "g__Parabacteroides_MGS",                                                                        
                                     "g__Bacteroides_MGS",                                                                                
                                     "g__Lachnospiraceae_noname_MGS",                                                                        
                                     "c__Actinobacteria_MGS",                                                                                                                              
                                     "c__Viruses_noname_MGS",                      
                                     "c__Bacteroidia_MGS",                                                                                                                                  
                                     "K02954_MGS",                                                                                                                                                                       
                                     "K07259_MGS",                                                                                                                                                                       
                                     "K07793_MGS",                                                                                                                                                                       
                                     "ko00633_MGS",                                                                                                                                                                      
                                     "ko00250_MGS",                                                                                                                                                                      
                                     "ko00230_MGS")  

write_RF_prep(as.data.frame(t(combined_response_prep)), "prep_data/combined_prep/combined_response_prep.txt", "response", biscuit_meta)
