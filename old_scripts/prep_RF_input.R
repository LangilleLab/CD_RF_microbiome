setwd("/Users/Gavin/Dropbox (Langille Lab)/work/Langille/crohns_disease/random_forest/random_forest_analyses/")

cutoff_pro = 0.1

### read in metadata
metadata <- read.table( "../Biscuit_metadata.txt" , header=T , sep="\t" , stringsAsFactors = FALSE , quote = "" , comment.char = "" )

### This filter takes in a table and removes all rows that 
### have less than a specified proportion of non-zero values
remove_rare <- function( table , cutoff_pro ) {
  row2keep <- c()
  for ( i in 1:nrow(table) ) {
    row_nonzero <- length( which( table[ i , ]  > 0 ) ) 
    cutoff <- ceiling( min_pro * ncol(table) )
    if ( row_nonzero > cutoff ) {
      row2keep <- c( row2keep , i)
    }
  }
  return( table [ row2keep , , drop=F ])
}

### define function that takes in table, metadata table and output prefix
### outputs filtered and standardized table for disease AND treatment response
### filter out features that are non-zero in specified # of samples
### Note that this is done AFTER the samples are subsetted (i.e. only for CD patients in the case of treatment response)
### Standardization (mean centering and scaling) is done AFTER this above filtering 

output_prepped_tables <- function( input , map , out_prefix , cutoff_pro_nonzero ) {
   
   map_subset <- map[ which(map$X.SampleID %in% colnames( input ) ) , ]
   rownames(map_subset) <- map_subset$X.SampleID
   map_subset[ which(map_subset$IBDDiagnosis == "Crohn\'s disease") , "IBDDiagnosis"] <- "CD"
   map_subset[ which(map_subset$IBDDiagnosis == "control" ) , "IBDDiagnosis"] <- "CN"
   
   disease_samples <- map_subset[ , "X.SampleID" ]
   response_samples <- map_subset[ which( map_subset$response != "control" ) , "X.SampleID" ]
   
   input_disease <- input[ , disease_samples ]
   input_response <- input[ , response_samples ]
   
   input_disease_filtered <- data.frame( remove_rare( input_disease , cutoff_pro_nonzero ) ) 
   input_response_filtered <- data.frame( remove_rare( input_response , cutoff_pro_nonzero )  ) 
   
   input_disease_filtered_scaled <- scale( input_disease_filtered , center = TRUE , scale = TRUE )
   input_response_filtered_scaled <- scale( input_response_filtered , center = TRUE , scale = TRUE )
 
   input_disease_prep <- data.frame( t ( input_disease_filtered_scaled ) ) 
   input_response_prep <- data.frame( t ( input_response_filtered_scaled ) ) 
   
   input_disease_prep$disease <- map_subset[ rownames(input_disease_prep) , "IBDDiagnosis" ]
   input_response_prep$response <- map_subset[ rownames(input_response_prep) , "response" ]
   
   outfile_disease <- paste( out_prefix , "disease_input.txt" , sep="_")
   outfile_response <- paste( out_prefix , "response_input.txt" , sep="_")
     
   write.table( x=input_disease_prep , file=outfile_disease , row.names=T , col.names=NA, quote = FALSE , sep="\t" )
   write.table( x=input_response_prep , file=outfile_response , row.names=T , col.names=NA, quote = FALSE , sep="\t" )
   
}

collapse_taxa <- function( table , level ) {
 
  ranks <- c("Phylum" , "Class" , "Order" , "Family" , "Genus" , "Species" )
  
  rank_subset <- ranks[1:which(ranks %in% level)]
  
  table$taxa <- table$Kingdom
  
  for ( r in rank_subset) { 
    table$taxa <- paste( table$taxa , table[ , r ] , sep="." )
  }
  
  table <- table[ , -c(1:7) ]
  
  table_summed <- aggregate( . ~ taxa , data=table, FUN=sum )

  rownames(table_summed) <- table_summed$taxa
  
  table_summed <- table_summed[ , -which(colnames(table_summed) == "taxa") ]
  
  return( table_summed )
  
}

### read in 16S data
otus <- read.table( "../clean_16S/16S_otu_table_clean.txt"  , header=T , sep="\t" , stringsAsFactors = FALSE , quote = "" , comment.char = "" , row.names = 1 )
taxa_16S <- read.table( "../clean_16S/16S_taxa_table_clean.txt"  , header=T , sep="\t" , stringsAsFactors = FALSE , quote = "" , comment.char = "" )
picrust_ko <- read.table( "../clean_16S/16S_picrust_KO_clean.txt"  , header=T , sep="\t" , stringsAsFactors = FALSE , quote = "" , comment.char = "", row.names = 1 )
picrust_pathways <- read.table( "../clean_16S/16S_picrust_pathways_clean.txt"  , header=T , sep="\t" , stringsAsFactors = FALSE , quote = "" , comment.char = "", row.names = 1 )

### read in MGS data
strains <- read.table( "../clean_MGS/MGS_metaphlan_strains_clean.txt" ,  header=T , sep="\t" , stringsAsFactors = FALSE , quote = "" , comment.char = "" , row.names = 1 )
taxa_MGS <- read.table( "../clean_MGS/MGS_metaphlan_taxa_clean.txt" ,  header=T , sep="\t" , stringsAsFactors = FALSE , quote = "" , comment.char = ""  )
humann_ko <- read.table( "../clean_MGS/MGS_humann_KO_clean.txt" ,  header=T , sep="\t" , stringsAsFactors = FALSE , quote = "" , comment.char = "" , row.names = 1 )
humann_pathways <- read.table( "../clean_MGS/MGS_humann_pathways_clean.txt" ,  header=T , sep="\t" , stringsAsFactors = FALSE , quote = "" , comment.char = "" , row.names = 1 )
humann_modules <- read.table( "../clean_MGS/MGS_humann_modules_clean.txt" ,  header=T , sep="\t" , stringsAsFactors = FALSE , quote = "" , comment.char = "" , row.names = 1 )

### removing 2 samples from humann modules since they are all 0s:
humann_modules <- humann_modules[ , -which(colnames(humann_modules) %in% c("B21_CD" , "B79_CD" )) ]


### collapse taxa_16S and taxa_MGS to species, genus, family, class, order and phylum

### start out by creating new dataframes for each case
### then add column called "taxa" to each, which is just the level that should be collapsed to
### remove other taxonomic columns
### sum all rows with the same value in the taxa column
### rename rows to be taxa column
### remove taxa column

# species
taxa_16S_species_summed <- collapse_taxa( table = taxa_16S , level= "Species" )
taxa_MGS_species_summed <- collapse_taxa( table = taxa_MGS , level= "Species" )

taxa_16S_genus_summed <- collapse_taxa( table = taxa_16S , level= "Genus" )
taxa_MGS_genus_summed <- collapse_taxa( table = taxa_MGS , level= "Genus" )

taxa_16S_family_summed <- collapse_taxa( table = taxa_16S , level= "Family" )
taxa_MGS_family_summed <- collapse_taxa( table = taxa_MGS , level= "Family" )

taxa_16S_order_summed <- collapse_taxa( table = taxa_16S , level= "Order" )
taxa_MGS_order_summed <- collapse_taxa( table = taxa_MGS , level= "Order" )

taxa_16S_class_summed <- collapse_taxa( table = taxa_16S , level= "Class" )
taxa_MGS_class_summed <- collapse_taxa( table = taxa_MGS , level= "Class" )

taxa_16S_phylum_summed <- collapse_taxa( table = taxa_16S , level= "Phylum" )
taxa_MGS_phylum_summed <- collapse_taxa( table = taxa_MGS , level= "Phylum" )

### Output prepped tables for 16S OTUs, all 6 taxa levels or 16S and MGS, MGS strains, picrust KOs and pathways, HUMAnN KOs, pathways and modules 

output_prepped_tables( otus , metadata , "16S_otus" , cutoff_pro )
output_prepped_tables( strains , metadata , "MGS_strains" , cutoff_pro )

output_prepped_tables( picrust_ko , metadata , "16S_KOs" , cutoff_pro )
output_prepped_tables( picrust_pathways , metadata , "16S_pathways" , cutoff_pro )

output_prepped_tables( humann_ko , metadata , "MGS_KOs" , cutoff_pro  )
output_prepped_tables( humann_pathways , metadata , "MGS_pathways" , cutoff_pro )
output_prepped_tables( humann_modules , metadata , "MGS_modules" , cutoff_pro )

output_prepped_tables( taxa_16S_species_summed , metadata , "16S_species" , cutoff_pro )
output_prepped_tables( taxa_MGS_species_summed , metadata , "MGS_species" , cutoff_pro )

output_prepped_tables( taxa_16S_genus_summed , metadata , "16S_genus" , cutoff_pro )
output_prepped_tables( taxa_MGS_genus_summed , metadata , "MGS_genus" , cutoff_pro )

output_prepped_tables( taxa_16S_family_summed , metadata , "16S_family" , cutoff_pro )
output_prepped_tables( taxa_MGS_family_summed , metadata , "MGS_family" , cutoff_pro )

output_prepped_tables( taxa_16S_order_summed , metadata , "16S_order" , cutoff_pro )
output_prepped_tables( taxa_MGS_order_summed , metadata , "MGS_order" , cutoff_pro )

output_prepped_tables( taxa_16S_class_summed , metadata , "16S_class" , cutoff_pro )
output_prepped_tables( taxa_MGS_class_summed , metadata , "MGS_class" , cutoff_pro )

output_prepped_tables( taxa_16S_phylum_summed , metadata , "16S_phylum" , cutoff_pro )
output_prepped_tables( taxa_MGS_phylum_summed , metadata , "MGS_phylum" , cutoff_pro )

