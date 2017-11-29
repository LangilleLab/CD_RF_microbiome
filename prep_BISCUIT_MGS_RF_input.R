### These commands were used to prepare the BISCUIT MGS tables to run with random forest.

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

response_samples <- rownames(biscuit_meta)[which(biscuit_meta$disease == "CD")]

# Read in HUMANN2 KO table.
mgs_ko <- read.table("raw_data/BISCUIT_MGS/biscuit_mgs_KOs.tsv", 
                         header=T, 
                         sep="\t", 
                         stringsAsFactors = FALSE, 
                         quote="", 
                         comment.char = "", 
                         row.names=1)

# Read in HUMANN2 modules table.
mgs_modules <- read.table("raw_data/BISCUIT_MGS/biscuit_mgs_KEGG_modules.tsv", 
                     header=T, 
                     sep="\t", 
                     stringsAsFactors = FALSE, 
                     quote="", 
                     comment.char = "", 
                     row.names=1)

# Read in HUMANN2 pathways table.
mgs_pathways <- read.table("raw_data/BISCUIT_MGS/biscuit_mgs_KEGG_pathways.tsv", 
                          header=T, 
                          sep="\t", 
                          stringsAsFactors = FALSE, 
                          quote="", 
                          comment.char = "", 
                          row.names=1)

# Remove "UNMAPPED" and "UNINTEGRATED"
mgs_modules <- mgs_modules[-which(rownames(mgs_modules) %in% c("UNMAPPED", "UNINTEGRATED")),]
mgs_pathways <- mgs_pathways[-which(rownames(mgs_pathways) %in% c("UNMAPPED", "UNINTEGRATED")),]

# Read in all collapsed taxonomy tables:
biscuit_mgs_strain <- read.table("raw_data/BISCUIT_MGS/metaphlan2_by_level/metaphlan2_out_merged_strain.tsv",
                                 sep="\t", header=T, row.names=1)
biscuit_mgs_species <- read.table("raw_data/BISCUIT_MGS/metaphlan2_by_level/metaphlan2_out_merged_species.tsv",
                                  sep="\t", header=T, row.names=1)
biscuit_mgs_genus <- read.table("raw_data/BISCUIT_MGS/metaphlan2_by_level/metaphlan2_out_merged_genus.tsv",
                                sep="\t", header=T, row.names=1)
biscuit_mgs_family <- read.table("raw_data/BISCUIT_MGS/metaphlan2_by_level/metaphlan2_out_merged_family.tsv",
                                 sep="\t", header=T, row.names=1) 
biscuit_mgs_order <- read.table("raw_data/BISCUIT_MGS/metaphlan2_by_level/metaphlan2_out_merged_order.tsv",
                                sep="\t", header=T, row.names=1)
biscuit_mgs_class <- read.table("raw_data/BISCUIT_MGS/metaphlan2_by_level/metaphlan2_out_merged_class.tsv",
                                sep="\t", header=T, row.names=1)
biscuit_mgs_phylum <- read.table("raw_data/BISCUIT_MGS/metaphlan2_by_level/metaphlan2_out_merged_phylum.tsv",
                                 sep="\t", header=T, row.names=1)

# Remove S34 and S38 from all taxa tables.
biscuit_mgs_strain <- biscuit_mgs_strain[,-which(colnames(biscuit_mgs_strain) %in% c("S34", "S38"))]
biscuit_mgs_species <- biscuit_mgs_species[,-which(colnames(biscuit_mgs_species) %in% c("S34", "S38"))]
biscuit_mgs_genus <- biscuit_mgs_genus[,-which(colnames(biscuit_mgs_genus) %in% c("S34", "S38"))]
biscuit_mgs_family <- biscuit_mgs_family[,-which(colnames(biscuit_mgs_family) %in% c("S34", "S38"))]
biscuit_mgs_order <- biscuit_mgs_order[,-which(colnames(biscuit_mgs_order) %in% c("S34", "S38"))]
biscuit_mgs_class <- biscuit_mgs_class[,-which(colnames(biscuit_mgs_class) %in% c("S34", "S38"))]
biscuit_mgs_phylum <- biscuit_mgs_phylum[,-which(colnames(biscuit_mgs_phylum) %in% c("S34", "S38"))]


# Remove features found in only <= 10% of samples and separate into disease and response tables.
# Also scale all tables.
mgs_ko_disease_filt <- as.data.frame(scale(remove_rare_rows(mgs_ko, 0.1)))
mgs_modules_disease_filt <- as.data.frame(scale(remove_rare_rows(mgs_modules, 0.1)))
mgs_pathways_disease_filt <- as.data.frame(scale(remove_rare_rows(mgs_pathways, 0.1)))

biscuit_mgs_strain_disease_filt <- as.data.frame(scale(remove_rare_rows(biscuit_mgs_strain, 0.1)))
biscuit_mgs_species_disease_filt <- as.data.frame(scale(remove_rare_rows(biscuit_mgs_species, 0.1)))
biscuit_mgs_genus_disease_filt <- as.data.frame(scale(remove_rare_rows(biscuit_mgs_genus, 0.1)))
biscuit_mgs_family_disease_filt <- as.data.frame(scale(remove_rare_rows(biscuit_mgs_family, 0.1)))
biscuit_mgs_order_disease_filt <- as.data.frame(scale(remove_rare_rows(biscuit_mgs_order, 0.1)))
biscuit_mgs_class_disease_filt <- as.data.frame(scale(remove_rare_rows(biscuit_mgs_class, 0.1)))
biscuit_mgs_phylum_disease_filt <- as.data.frame(scale(remove_rare_rows(biscuit_mgs_phylum, 0.1)))


mgs_ko_response_filt <- as.data.frame(scale(remove_rare_rows(mgs_ko[,response_samples], 0.1)))
mgs_modules_response_filt <- as.data.frame(scale(remove_rare_rows(mgs_modules[,response_samples], 0.1)))
mgs_pathways_response_filt <- as.data.frame(scale(remove_rare_rows(mgs_pathways[,response_samples], 0.1)))

biscuit_mgs_strain_response_filt <- as.data.frame(scale(remove_rare_rows(biscuit_mgs_strain[,response_samples], 0.1)))
biscuit_mgs_species_response_filt <- as.data.frame(scale(remove_rare_rows(biscuit_mgs_species[,response_samples], 0.1)))
biscuit_mgs_genus_response_filt <- as.data.frame(scale(remove_rare_rows(biscuit_mgs_genus[,response_samples], 0.1)))
biscuit_mgs_family_response_filt <- as.data.frame(scale(remove_rare_rows(biscuit_mgs_family[,response_samples], 0.1)))
biscuit_mgs_order_response_filt <- as.data.frame(scale(remove_rare_rows(biscuit_mgs_order[,response_samples], 0.1)))
biscuit_mgs_class_response_filt <- as.data.frame(scale(remove_rare_rows(biscuit_mgs_class[,response_samples], 0.1)))
biscuit_mgs_phylum_response_filt <- as.data.frame(scale(remove_rare_rows(biscuit_mgs_phylum[,response_samples], 0.1)))


# Write output file:
write_RF_prep(mgs_ko_disease_filt, "prep_data/BISCUIT_MGS/biscuit_mgs_ko_disease_prep.tsv", "disease", biscuit_meta)
write_RF_prep(mgs_modules_disease_filt, "prep_data/BISCUIT_MGS/biscuit_mgs_moduless_disease_prep.tsv", "disease", biscuit_meta)
write_RF_prep(mgs_pathways_disease_filt, "prep_data/BISCUIT_MGS/biscuit_mgs_pathways_disease_prep.tsv", "disease", biscuit_meta)

write_RF_prep(biscuit_mgs_strain_disease_filt, "prep_data/BISCUIT_MGS/biscuit_mgs_mgs_strain_disease_prep.tsv", "disease", biscuit_meta)
write_RF_prep(biscuit_mgs_species_disease_filt, "prep_data/BISCUIT_MGS/biscuit_mgs_species_disease_prep.tsv", "disease", biscuit_meta)
write_RF_prep(biscuit_mgs_genus_disease_filt, "prep_data/BISCUIT_MGS/biscuit_mgs_genus_disease_prep.tsv", "disease", biscuit_meta)
write_RF_prep(biscuit_mgs_family_disease_filt, "prep_data/BISCUIT_MGS/biscuit_mgs_family_disease_prep.tsv", "disease", biscuit_meta)
write_RF_prep(biscuit_mgs_order_disease_filt, "prep_data/BISCUIT_MGS/biscuit_mgs_order_disease_prep.tsv", "disease", biscuit_meta)
write_RF_prep(biscuit_mgs_class_disease_filt, "prep_data/BISCUIT_MGS/biscuit_mgs_class_disease_prep.tsv", "disease", biscuit_meta)
write_RF_prep(biscuit_mgs_phylum_disease_filt, "prep_data/BISCUIT_MGS/biscuit_mgs_phylum_disease_prep.tsv", "disease", biscuit_meta)

# Write output file:
write_RF_prep(mgs_ko_response_filt, "prep_data/BISCUIT_MGS/biscuit_mgs_ko_response_prep.tsv", "response", biscuit_meta)
write_RF_prep(mgs_modules_response_filt, "prep_data/BISCUIT_MGS/biscuit_mgs_moduless_response_prep.tsv", "response", biscuit_meta)
write_RF_prep(mgs_pathways_response_filt, "prep_data/BISCUIT_MGS/biscuit_mgs_pathways_response_prep.tsv", "response", biscuit_meta)

write_RF_prep(biscuit_mgs_strain_response_filt, "prep_data/BISCUIT_MGS/biscuit_mgs_mgs_strain_response_prep.tsv", "response", biscuit_meta)
write_RF_prep(biscuit_mgs_species_response_filt, "prep_data/BISCUIT_MGS/biscuit_mgs_species_response_prep.tsv", "response", biscuit_meta)
write_RF_prep(biscuit_mgs_genus_response_filt, "prep_data/BISCUIT_MGS/biscuit_mgs_genus_response_prep.tsv", "response", biscuit_meta)
write_RF_prep(biscuit_mgs_family_response_filt, "prep_data/BISCUIT_MGS/biscuit_mgs_family_response_prep.tsv", "response", biscuit_meta)
write_RF_prep(biscuit_mgs_order_response_filt, "prep_data/BISCUIT_MGS/biscuit_mgs_order_response_prep.tsv", "response", biscuit_meta)
write_RF_prep(biscuit_mgs_class_response_filt, "prep_data/BISCUIT_MGS/biscuit_mgs_class_response_prep.tsv", "response", biscuit_meta)
write_RF_prep(biscuit_mgs_phylum_response_filt, "prep_data/BISCUIT_MGS/biscuit_mgs_phylum_response_prep.tsv", "response", biscuit_meta)

