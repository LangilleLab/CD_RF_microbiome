### These commands were used to prepare the BISCUIT 16S tables to run with random forest.

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

# Read in PICRUSt KO table.
picrust_ko <- read.table("raw_data/BISCUIT_16S/ko.tsv", 
                         header=T, 
                         sep="\t", 
                         stringsAsFactors = FALSE, 
                         quote="", 
                         comment.char = "", 
                         row.names=1)

# Keep only KO id in rownames.
rownames(picrust_ko) <- gsub(":.+$", "", rownames(picrust_ko))

# Read in PICRUSt pathways tables.
picrust_pathways <- read.table("raw_data/BISCUIT_16S/ko_l3.spf", 
                               header=T, 
                               sep="\t", 
                               stringsAsFactors = FALSE, 
                               quote="", 
                               comment.char = "")

# Rename as pathway df rownames as concatenation of first 3 columns.
rownames(picrust_pathways) <- paste(picrust_pathways$Level_1, picrust_pathways$Level_2, picrust_pathways$Level_3, sep="|")

# Replace all spaces with underscores.
rownames(picrust_pathways) <- gsub(" ", "_", rownames(picrust_pathways))

# Remove the first 3 columns.
picrust_pathways <- picrust_pathways[ , -c(1:3)]

# Read in OTU table.
biscuit_otu <- read.table( "raw_data/BISCUIT_16S/otu_table_w_tax_BISCUIT.txt", 
                           header=T, 
                           sep="\t", 
                           stringsAsFactors = FALSE, 
                           quote="", 
                           comment.char = "",
                           row.names=1)

# Remove "taxonomy" column:
biscuit_otu <- biscuit_otu[,-which(colnames(biscuit_otu) %in% "taxonomy")]

# Read in all collapsed taxonomy tables:
biscuit_species <- read_spf_by_level("raw_data/BISCUIT_16S/otu_table_w_tax_BISCUIT.spf", "Species")
biscuit_genus <- read_spf_by_level("raw_data/BISCUIT_16S/otu_table_w_tax_BISCUIT.spf", "Genus")
biscuit_family <- read_spf_by_level("raw_data/BISCUIT_16S/otu_table_w_tax_BISCUIT.spf", "Family")
biscuit_order <- read_spf_by_level("raw_data/BISCUIT_16S/otu_table_w_tax_BISCUIT.spf", "Order")
biscuit_class <- read_spf_by_level("raw_data/BISCUIT_16S/otu_table_w_tax_BISCUIT.spf", "Class")
biscuit_phylum <- read_spf_by_level("raw_data/BISCUIT_16S/otu_table_w_tax_BISCUIT.spf", "Phylum")

biscuit_otu_unrarified <- read.table("raw_data/BISCUIT_16S/unrarified/otu_table_high_conf_BISCUIT.tsv", 
                                      header=T, 
                                      sep="\t", 
                                      stringsAsFactors = FALSE, 
                                      quote="", 
                                      comment.char = "",
                                      row.names=1)

# Remove "taxonomy" column:
biscuit_otu_unrarified <- biscuit_otu_unrarified[,-which(colnames(biscuit_otu_unrarified) %in% "taxonomy")]

biscuit_species_unrarified <- read_spf_by_level("raw_data/BISCUIT_16S/unrarified/otu_table_high_conf_BISCUIT.spf", "Species")
biscuit_genus_unrarified <- read_spf_by_level("raw_data/BISCUIT_16S/unrarified/otu_table_high_conf_BISCUIT.spf", "Genus")
biscuit_family_unrarified <- read_spf_by_level("raw_data/BISCUIT_16S/unrarified/otu_table_high_conf_BISCUIT.spf", "Family")
biscuit_order_unrarified <- read_spf_by_level("raw_data/BISCUIT_16S/unrarified/otu_table_high_conf_BISCUIT.spf", "Order")
biscuit_class_unrarified <- read_spf_by_level("raw_data/BISCUIT_16S/unrarified/otu_table_high_conf_BISCUIT.spf", "Class")
biscuit_phylum_unrarified <- read_spf_by_level("raw_data/BISCUIT_16S/unrarified/otu_table_high_conf_BISCUIT.spf", "Phylum")

# Remove features found in only <= 10% of samples and separate into disease and response tables.
# Scale non-taxonomic tables.
picrust_ko_disease_filt <- as.data.frame(scale(remove_rare_rows(picrust_ko, 0.1)))
picrust_pathways_disease_filt <- as.data.frame(scale(remove_rare_rows(picrust_pathways, 0.1)))

biscuit_otu_disease_filt <- remove_rare_rows(biscuit_otu, 0.1)
biscuit_species_disease_filt <- remove_rare_rows(biscuit_species, 0.1)
biscuit_genus_disease_filt <- remove_rare_rows(biscuit_genus, 0.1)
biscuit_family_disease_filt <- remove_rare_rows(biscuit_family, 0.1)
biscuit_order_disease_filt <- remove_rare_rows(biscuit_order, 0.1)
biscuit_class_disease_filt <- remove_rare_rows(biscuit_class, 0.1)
biscuit_phylum_disease_filt <- remove_rare_rows(biscuit_phylum, 0.1)

picrust_ko_response_filt <- as.data.frame(scale(remove_rare_rows(picrust_ko[,response_samples], 0.1)))
picrust_pathways_response_filt <- as.data.frame(scale(remove_rare_rows(picrust_pathways[,response_samples], 0.1)))

biscuit_otu_response_filt <- remove_rare_rows(biscuit_otu[,response_samples], 0.1)
biscuit_species_response_filt <- remove_rare_rows(biscuit_species[,response_samples], 0.1)
biscuit_genus_response_filt <- remove_rare_rows(biscuit_genus[,response_samples], 0.1)
biscuit_family_response_filt <- remove_rare_rows(biscuit_family[,response_samples], 0.1)
biscuit_order_response_filt <- remove_rare_rows(biscuit_order[,response_samples], 0.1)
biscuit_class_response_filt <- remove_rare_rows(biscuit_class[,response_samples], 0.1)
biscuit_phylum_response_filt <- remove_rare_rows(biscuit_phylum[,response_samples], 0.1)

# Remove rows from unrarified as well:
biscuit_otu_unrarified_disease_filt <- remove_rare_rows(biscuit_otu_unrarified, 0.1)
biscuit_species_unrarified_disease_filt <- remove_rare_rows(biscuit_species_unrarified, 0.1)
biscuit_genus_unrarified_disease_filt <- remove_rare_rows(biscuit_genus_unrarified, 0.1)
biscuit_family_unrarified_disease_filt <- remove_rare_rows(biscuit_family_unrarified, 0.1)
biscuit_order_unrarified_disease_filt <- remove_rare_rows(biscuit_order_unrarified, 0.1)
biscuit_class_unrarified_disease_filt <- remove_rare_rows(biscuit_class_unrarified, 0.1)
biscuit_phylum_unrarified_disease_filt <- remove_rare_rows(biscuit_phylum_unrarified, 0.1)

biscuit_otu_unrarified_response_filt <- remove_rare_rows(biscuit_otu_unrarified[,response_samples], 0.1)
biscuit_species_unrarified_response_filt <- remove_rare_rows(biscuit_species_unrarified[,response_samples], 0.1)
biscuit_genus_unrarified_response_filt <- remove_rare_rows(biscuit_genus_unrarified[,response_samples], 0.1)
biscuit_family_unrarified_response_filt <- remove_rare_rows(biscuit_family_unrarified[,response_samples], 0.1)
biscuit_order_unrarified_response_filt <- remove_rare_rows(biscuit_order_unrarified[,response_samples], 0.1)
biscuit_class_unrarified_response_filt <- remove_rare_rows(biscuit_class_unrarified[,response_samples], 0.1)
biscuit_phylum_unrarified_response_filt <- remove_rare_rows(biscuit_phylum_unrarified[,response_samples], 0.1)

# Impute 0 values and apply centred log-ratio transform for taxonomic tables. Scale output.
biscuit_otu_unrarified_disease_filt_clr <- as.data.frame(scale(CZM_and_clr(biscuit_otu_unrarified_disease_filt)))
biscuit_otu_unrarified_response_filt_clr <- as.data.frame(scale(CZM_and_clr(biscuit_otu_unrarified_response_filt)))

biscuit_species_unrarified_disease_filt_clr <- as.data.frame(scale(CZM_and_clr(biscuit_species_unrarified_disease_filt)))
biscuit_species_unrarified_response_filt_clr <- as.data.frame(scale(CZM_and_clr(biscuit_species_unrarified_response_filt)))

biscuit_genus_unrarified_disease_filt_clr <- as.data.frame(scale(CZM_and_clr(biscuit_genus_unrarified_disease_filt)))
biscuit_genus_unrarified_response_filt_clr <- as.data.frame(scale(CZM_and_clr(biscuit_genus_unrarified_response_filt)))

biscuit_family_unrarified_disease_filt_clr <- as.data.frame(scale(CZM_and_clr(biscuit_family_unrarified_disease_filt)))
biscuit_family_unrarified_response_filt_clr <- as.data.frame(scale(CZM_and_clr(biscuit_family_unrarified_response_filt)))

biscuit_order_unrarified_disease_filt_clr <- as.data.frame(scale(CZM_and_clr(biscuit_order_unrarified_disease_filt)))
biscuit_order_unrarified_response_filt_clr <- as.data.frame(scale(CZM_and_clr(biscuit_order_unrarified_response_filt)))

biscuit_class_unrarified_disease_filt_clr <- as.data.frame(scale(CZM_and_clr(biscuit_class_unrarified_disease_filt)))
biscuit_class_unrarified_response_filt_clr <- as.data.frame(scale(CZM_and_clr(biscuit_class_unrarified_response_filt)))

biscuit_phylum_unrarified_disease_filt_clr <- as.data.frame(scale(CZM_and_clr(biscuit_phylum_unrarified_disease_filt)))
biscuit_phylum_unrarified_response_filt_clr <- as.data.frame(scale(CZM_and_clr(biscuit_phylum_unrarified_response_filt)))

# Also scale data when it's non-CLR transformed:
biscuit_otu_disease_filt_NOclr <- as.data.frame(scale(biscuit_otu_disease_filt))
biscuit_otu_response_filt_NOclr <- as.data.frame(scale(biscuit_otu_response_filt))

biscuit_species_disease_filt_NOclr <- as.data.frame(scale(biscuit_species_disease_filt))
biscuit_species_response_filt_NOclr <- as.data.frame(scale(biscuit_species_response_filt))

biscuit_genus_disease_filt_NOclr <- as.data.frame(scale(biscuit_genus_disease_filt))
biscuit_genus_response_filt_NOclr <- as.data.frame(scale(biscuit_genus_response_filt))

biscuit_family_disease_filt_NOclr <- as.data.frame(scale(biscuit_family_disease_filt))
biscuit_family_response_filt_NOclr <- as.data.frame(scale(biscuit_family_response_filt))

biscuit_order_disease_filt_NOclr <- as.data.frame(scale(biscuit_order_disease_filt))
biscuit_order_response_filt_NOclr <- as.data.frame(scale(biscuit_order_response_filt))

biscuit_class_disease_filt_NOclr <- as.data.frame(scale(biscuit_class_disease_filt))
biscuit_class_response_filt_NOclr <- as.data.frame(scale(biscuit_class_response_filt))

biscuit_phylum_disease_filt_NOclr <- as.data.frame(scale(biscuit_phylum_disease_filt))
biscuit_phylum_response_filt_NOclr <- as.data.frame(scale(biscuit_phylum_response_filt))


# Write output file:
write_RF_prep(picrust_ko_disease_filt, "prep_data/BISCUIT_16S/biscuit_ko_disease_prep.tsv", "disease", biscuit_meta)
write_RF_prep(picrust_pathways_disease_filt, "prep_data/BISCUIT_16S/biscuit_pathways_disease_prep.tsv", "disease", biscuit_meta)
write_RF_prep(biscuit_otu_disease_filt_NOclr, "prep_data/BISCUIT_16S/biscuit_otu_disease_prep.tsv", "disease", biscuit_meta)
write_RF_prep(biscuit_species_disease_filt_NOclr, "prep_data/BISCUIT_16S/biscuit_species_disease_prep.tsv", "disease", biscuit_meta)
write_RF_prep(biscuit_genus_disease_filt_NOclr, "prep_data/BISCUIT_16S/biscuit_genus_disease_prep.tsv", "disease", biscuit_meta)
write_RF_prep(biscuit_family_disease_filt_NOclr, "prep_data/BISCUIT_16S/biscuit_family_disease_prep.tsv", "disease", biscuit_meta)
write_RF_prep(biscuit_order_disease_filt_NOclr, "prep_data/BISCUIT_16S/biscuit_order_disease_prep.tsv", "disease", biscuit_meta)
write_RF_prep(biscuit_class_disease_filt_NOclr, "prep_data/BISCUIT_16S/biscuit_class_disease_prep.tsv", "disease", biscuit_meta)
write_RF_prep(biscuit_phylum_disease_filt_NOclr, "prep_data/BISCUIT_16S/biscuit_phylum_disease_prep.tsv", "disease", biscuit_meta)

write_RF_prep(picrust_ko_response_filt, "prep_data/BISCUIT_16S/biscuit_ko_response_prep.tsv", "response", biscuit_meta)
write_RF_prep(picrust_pathways_response_filt, "prep_data/BISCUIT_16S/biscuit_pathways_response_prep.tsv", "response", biscuit_meta)
write_RF_prep(biscuit_otu_response_filt_NOclr, "prep_data/BISCUIT_16S/biscuit_otu_response_prep.tsv", "response", biscuit_meta)
write_RF_prep(biscuit_species_response_filt_NOclr, "prep_data/BISCUIT_16S/biscuit_species_response_prep.tsv", "response", biscuit_meta)
write_RF_prep(biscuit_genus_response_filt_NOclr, "prep_data/BISCUIT_16S/biscuit_genus_response_prep.tsv", "response", biscuit_meta)
write_RF_prep(biscuit_family_response_filt_NOclr, "prep_data/BISCUIT_16S/biscuit_family_response_prep.tsv", "response", biscuit_meta)
write_RF_prep(biscuit_order_response_filt_NOclr, "prep_data/BISCUIT_16S/biscuit_order_response_prep.tsv", "response", biscuit_meta)
write_RF_prep(biscuit_class_response_filt_NOclr, "prep_data/BISCUIT_16S/biscuit_class_response_prep.tsv", "response", biscuit_meta)
write_RF_prep(biscuit_phylum_response_filt_NOclr, "prep_data/BISCUIT_16S/biscuit_phylum_response_prep.tsv", "response", biscuit_meta)


# Write CLR input:
write_RF_prep(biscuit_otu_unrarified_disease_filt_clr, "prep_data/BISCUIT_16S/unrarified_CLR/biscuit_otu_unrarified_disease_prep.tsv", "disease", biscuit_meta)
write_RF_prep(biscuit_species_unrarified_disease_filt_clr, "prep_data/BISCUIT_16S/unrarified_CLR/biscuit_species_unrarified_disease_prep.tsv", "disease", biscuit_meta)
write_RF_prep(biscuit_genus_unrarified_disease_filt_clr, "prep_data/BISCUIT_16S/unrarified_CLR/biscuit_genus_unrarified_disease_prep.tsv", "disease", biscuit_meta)
write_RF_prep(biscuit_family_unrarified_disease_filt_clr, "prep_data/BISCUIT_16S/unrarified_CLR/biscuit_family_unrarified_disease_prep.tsv", "disease", biscuit_meta)
write_RF_prep(biscuit_order_unrarified_disease_filt_clr, "prep_data/BISCUIT_16S/unrarified_CLR/biscuit_order_unrarified_disease_prep.tsv", "disease", biscuit_meta)
write_RF_prep(biscuit_class_unrarified_disease_filt_clr, "prep_data/BISCUIT_16S/unrarified_CLR/biscuit_class_unrarified_disease_prep.tsv", "disease", biscuit_meta)
write_RF_prep(biscuit_phylum_unrarified_disease_filt_clr, "prep_data/BISCUIT_16S/unrarified_CLR/biscuit_phylum_unrarified_disease_prep.tsv", "disease", biscuit_meta)

write_RF_prep(biscuit_otu_unrarified_response_filt_clr, "prep_data/BISCUIT_16S/unrarified_CLR/biscuit_otu_unrarified_response_prep.tsv", "response", biscuit_meta)
write_RF_prep(biscuit_species_unrarified_response_filt_clr, "prep_data/BISCUIT_16S/unrarified_CLR/biscuit_species_unrarified_response_prep.tsv", "response", biscuit_meta)
write_RF_prep(biscuit_genus_unrarified_response_filt_clr, "prep_data/BISCUIT_16S/unrarified_CLR/biscuit_genus_unrarified_response_prep.tsv", "response", biscuit_meta)
write_RF_prep(biscuit_family_unrarified_response_filt_clr, "prep_data/BISCUIT_16S/unrarified_CLR/biscuit_family_unrarified_response_prep.tsv", "response", biscuit_meta)
write_RF_prep(biscuit_order_unrarified_response_filt_clr, "prep_data/BISCUIT_16S/unrarified_CLR/biscuit_order_unrarified_response_prep.tsv", "response", biscuit_meta)
write_RF_prep(biscuit_class_unrarified_response_filt_clr, "prep_data/BISCUIT_16S/unrarified_CLR/biscuit_class_unrarified_response_prep.tsv", "response", biscuit_meta)
write_RF_prep(biscuit_phylum_unrarified_response_filt_clr, "prep_data/BISCUIT_16S/unrarified_CLR/biscuit_phylum_unrarified_response_prep.tsv", "response", biscuit_meta)
