setwd("/home/gavin/github_repos/CD_RF_microbiome/raw_data/")

source("../BISCUIT_utility_code.R")

raw_16S <- read.table("BISCUIT_16S/otu_table_w_tax_BISCUIT.spf", header=T, stringsAsFactors = FALSE, sep="\t")
colnames(raw_16S) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", colnames(raw_16S)[8:ncol(raw_16S)])

raw_mgs_genus <- read.table("BISCUIT_MGS/metaphlan2_by_level/metaphlan2_out_merged_genus.tsv", 
                            header=T, stringsAsFactors = FALSE, sep="\t", row.names=1)

mapfile <- read.table("../biscuit_metadata.txt", header=T , sep="\t", comment.char = "", quote = "", stringsAsFactors = FALSE)

mapfile <- mapfile[-which(mapfile$sample_id %in% c("S34", "S38")),]

rownames(mapfile) <- mapfile$sample_id


raw_16S_genus <- collapse_taxa(raw_16S, "Genus")

mapfile_sorted <- mapfile[colnames(raw_16S_genus) , ]

raw_16S_genus_summed_sweep <- sweep(raw_16S_genus, 2, colSums(raw_16S_genus) , '/') *100 
raw_mgs_genus_summed <- raw_mgs_genus[ , colnames(raw_16S_genus) ]

disease_state_col <- mapfile_sorted$disease
disease_state_col[ which(mapfile_sorted$disease == "CD")] <- "black"
disease_state_col[ which(mapfile_sorted$disease == "CN")] <- "white"

# Missing in MGS data so set all to 0
MGS_Desulfovibrio <- rep.int( 0 , times=38 )

par(mfrow=c(1,2))

boxplot( log(as.numeric( raw_16S_genus_summed_sweep[ grep( "g__Desulfovibrio" , rownames(raw_16S_genus_summed_sweep) ) , ]) + 1 ) , log( MGS_Desulfovibrio + 1 ) , outline=FALSE, ylim=c(0,1.2) , ylab="log(g__Desulfovibrio rel. abundance +1 )", names=c("16S" , "MGS")  )

stripchart( log(as.numeric( raw_16S_genus_summed_sweep[ grep( "g__Desulfovibrio" , rownames(raw_16S_genus_summed_sweep) ) , ]) + 1 )  , vertical = TRUE, jitter=0.2, method = "jitter", add = TRUE, pch = 21, bg = disease_state_col )

stripchart( log(MGS_Desulfovibrio + 1 )  , vertical = TRUE, jitter=0.2, method = "jitter", add = TRUE, pch = 21, bg = disease_state_col, at=2 )

legend( "topright" , legend=c("CD" , "CN") , fill=c("black" , "white") , cex=1)


boxplot( log(as.numeric( raw_16S_genus_summed_sweep[ grep( "Akkermansia" , rownames(raw_16S_genus_summed_sweep) ) , ]) + 1 ) , log( as.numeric( raw_mgs_genus_summed[ grep( "Akkermansia" , rownames(raw_mgs_genus_summed) ) ,  ] ) + 1 ) , outline=FALSE, ylim=c(0,4) , ylab="log(g__Akkermansia rel. abundance +1 )", names=c("16S" , "MGS") )

stripchart( log(as.numeric( raw_16S_genus_summed_sweep[ grep( "Akkermansia" , rownames(raw_16S_genus_summed_sweep) ) , ]) + 1 )  , vertical = TRUE, jitter=0.2, method = "jitter", add = TRUE, pch = 21, bg = disease_state_col )
 
stripchart( log(as.numeric( raw_mgs_genus_summed[ grep( "Akkermansia" , rownames(raw_mgs_genus_summed) ) , ]) + 1 )  , vertical = TRUE, jitter=0.2, method = "jitter", add = TRUE, pch = 21, bg = disease_state_col, at=2 )

legend( "topleft" , legend=c("CD" , "CN") , fill=c("black" , "white") , cex=1)
