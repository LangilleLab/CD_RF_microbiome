setwd( "/Users/Gavin/Dropbox/work/Langille/crohns_disease/random_forest/raw_MGS/")

humann_ko <- read.table( "hummann_kos.spf" , header=T , sep="\t", stringsAsFactors = FALSE, quote="", comment.char = "", row.names=1 )
rownames( humann_ko ) <- sub( ":.+$" , "", rownames(humann_ko) )
  
humann_modules <- read.table( "hummann_modules.spf" , header=T , sep="\t", stringsAsFactors = FALSE, quote="", comment.char = "", row.names=1 )
rownames( humann_modules ) <- sub( ":.+$" , "", rownames(humann_modules) )

humann_pathways <- read.table( "hummann_pathways.spf" , header=T , sep="\t", stringsAsFactors = FALSE, quote="", comment.char = "", row.names=1 )
rownames( humann_pathways ) <- sub( ":.+$" , "", rownames(humann_pathways) )

mgs_taxa <- read.table( "metaphlan_taxonomy_PhiX-screened.spf" , header=T , sep="\t" , stringsAsFactors=FALSE , quote="", comment.char = "" )
mgs_taxa_table <- mgs_taxa[ , -which(colnames(mgs_taxa) == "Strain" ) ]

rownames(mgs_taxa) <- paste( mgs_taxa$Kingdom , mgs_taxa$Phylum , mgs_taxa$Class , mgs_taxa$Order , mgs_taxa$Family , mgs_taxa$Genus , mgs_taxa$Species , mgs_taxa$Strain , sep="." )
mgs_taxa <- mgs_taxa[ , -c(1:8)]


write.table( x=humann_ko , file="../clean_MGS/MGS_humann_KO_clean.txt" , row.names=T , col.names=NA, quote = FALSE , sep="\t" )
write.table( x=humann_modules , file="../clean_MGS/MGS_humann_modules_clean.txt" , row.names=T , col.names=NA, quote = FALSE , sep="\t" )
write.table( x=humann_pathways , file="../clean_MGS/MGS_humann_pathways_clean.txt" , row.names=T , col.names=NA, quote = FALSE , sep="\t" )

write.table( x=mgs_taxa , file="../clean_MGS/MGS_metaphlan_strains_clean.txt" , row.names=T , col.names=NA, quote = FALSE , sep="\t" )
write.table( x=mgs_taxa_table , file="../clean_MGS/MGS_metaphlan_taxa_clean.txt" , row.names=F , col.names=T, quote = FALSE , sep="\t" )
