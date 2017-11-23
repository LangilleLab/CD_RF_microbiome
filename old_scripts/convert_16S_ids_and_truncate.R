setwd( "/Users/Gavin/Dropbox/work/Langille/crohns_disease/random_forest/raw_16S/")

ids <- read.table( "id_conversions.txt" , header=T , sep="\t" )
ids$X16S_ids <- as.character( ids$X16S_ids )
ids$genome_ids <- as.character( ids$genome_ids )
rownames(ids) <- ids$X16S_ids


picrust_ko <- read.table( "ko.spf" , header=T , sep="\t", stringsAsFactors = FALSE, quote="", comment.char = "", row.names=1 )

picrust_pathways <- read.table( "ko_L3.spf" , header=T , sep="\t", stringsAsFactors = FALSE, quote="", comment.char = "" )
rownames(picrust_pathways) <- paste( picrust_pathways$Level_1 , picrust_pathways$Level_2 , picrust_pathways$Level_3 , sep="_" )
picrust_pathways <- picrust_pathways[ , -c(1:3) ]

stamp_otu <- read.table( "otu_table.spf" , header=T , sep="\t", stringsAsFactors = FALSE, quote="", comment.char = "" )
otu <- read.table( "otu_table_w_tax.txt" , header=T , sep="\t", stringsAsFactors = FALSE, quote="", comment.char = "" )
rownames( otu ) <- paste( otu$OTU.ID , otu$taxonomy , sep="." )
rownames(otu) <- gsub( "; " , "." , rownames( otu ) )
otu <- otu[ , -c( which(colnames(otu) == "OTU.ID") , which(colnames(otu) == "taxonomy") ) ]

colnames(picrust_ko) <- sub( "^X" , "" , colnames(picrust_ko) )
colnames(picrust_pathways) <- sub( "^X" , "" , colnames(picrust_pathways) )
colnames(stamp_otu) <- sub( "^X" , "" , colnames(stamp_otu) )
colnames(otu) <- sub( "^X" , "" , colnames(otu) )


picrust_ko_sorted <- picrust_ko[ , rownames(ids) ]
picrust_pathways_sorted <- picrust_pathways[ , rownames(ids) ]
stamp_otu_sorted <- stamp_otu[ , c("Level_1" , "Level_2" ,"Level_3" ,"Level_4" ,"Level_5" ,"Level_6" ,"Level_7" , rownames(ids) ) ]
otu_sorted <- otu[ , rownames(ids) ]

colnames(picrust_ko_sorted) <- ids$genome_ids
colnames(picrust_pathways_sorted) <- ids$genome_ids
colnames(stamp_otu_sorted) <- c("Kingdom" , "Phylum" , "Class" , "Order" , "Family", "Genus", "Species" , ids$genome_ids )
colnames(otu_sorted) <- ids$genome_ids

# just get level 3 of pathways
picrust_pathways_sorted_rownames <- c()
picrust_pathways_rowsplit <- strsplit( rownames(picrust_pathways_sorted) , "_" )
for ( i in 1:nrow(picrust_pathways_sorted)){
  picrust_pathways_sorted_rownames <- c(picrust_pathways_sorted_rownames , picrust_pathways_rowsplit[[i]][3])
}
rownames( picrust_pathways_sorted ) <- picrust_pathways_sorted_rownames

# just get KEGG ID
rownames( picrust_ko_sorted ) <- sub( ":.+$" , "", rownames(picrust_ko_sorted))


write.table( x=picrust_ko_sorted , file="../clean_16S/16S_picrust_KO_clean.txt" , row.names=T , col.names=NA, quote = FALSE , sep="\t" )
write.table( x=picrust_pathways_sorted , file="../clean_16S/16S_picrust_pathways_clean.txt" , row.names=T , col.names=NA, quote = FALSE , sep="\t" )
write.table( x=stamp_otu_sorted , file="../clean_16S/16S_taxa_table_clean.txt" , row.names=F , col.names=T, quote = FALSE , sep="\t" )
write.table( x=otu_sorted , file="../clean_16S/16S_otu_table_clean.txt"  , row.names=T , col.names=NA, quote = FALSE , sep="\t" )


