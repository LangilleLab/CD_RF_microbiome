library(reshape2)
library(ggplot2)

setwd("/home/gavin/github_repos/CD_RF_microbiome/raw_data/")

source("../BISCUIT_utility_code.R")

raw_16S <- read.table("BISCUIT_16S/otu_table_w_tax_BISCUIT.spf", header=T, stringsAsFactors = FALSE, sep="\t")
colnames(raw_16S) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", colnames(raw_16S)[8:ncol(raw_16S)])

raw_mgs_class <- read.table("BISCUIT_MGS/metaphlan2_by_level/metaphlan2_out_merged_class.tsv", 
                            header=T, stringsAsFactors = FALSE, sep="\t", row.names=1)

mapfile <- read.table("../biscuit_metadata.txt", header=T , sep="\t", comment.char = "", quote = "", stringsAsFactors = FALSE)

samples_ordered <- c(1:40)
samples_ordered <- paste("S", samples_ordered, sep="")
samples_ordered <- samples_ordered[-which(samples_ordered %in% c("S34", "S38"))]

mapfile <- mapfile[-which(mapfile$sample_id %in% c("S34", "S38")),]

mapfile$sample_id <- factor(mapfile$sample_id, levels=c(samples_ordered))
mapfile$disease <- factor(mapfile$disease, levels=c("CN", "CD"))
mapfile$response <- as.factor(mapfile$response)

mapfile_sorted <- mapfile[with(mapfile, order(disease, response, sample_id)), ]

raw_16S_class <- collapse_taxa(raw_16S, "Class")

raw_16S_class_sorted <- raw_16S_class[, as.character(mapfile_sorted$sample_id)]
raw_mgs_class_sorted <- raw_mgs_class[, as.character(mapfile_sorted$sample_id)]


prep_16S_class <- raw_16S_class_sorted
prep_16S_class$taxa <- as.factor(rownames(raw_16S_class_sorted))
rownames(prep_16S_class) <- NULL
prep_16S_class_long <- melt(prep_16S_class, id.vars=c("taxa"))

prep_mgs_class <- raw_mgs_class_sorted
prep_mgs_class$taxa <- as.factor(rownames(raw_mgs_class_sorted))
rownames(prep_mgs_class) <- NULL
prep_mgs_class_long <- melt(prep_mgs_class, id.vars=c("taxa"))


qual_col <- c( "#a6cee3" ,  "#1f78b4" , "#b2df8a" , "#33a02c" , "#fb9a99" , "#e31a1c" ,
               "#fdbf6f" , "#ff7f00" , "#cab2d6" , "#6a3d9a" , "#ffff99" , "#b15928" )

qual_col_16S <- c(qual_col, "#a6cee3" ,  "#1f78b4" , "#b2df8a" , "#33a02c" , "#fb9a99" , "#e31a1c")

qual_col_mgs <- c("black", qual_col)

ggplot(prep_16S_class_long, aes(x=variable, y=value, fill=taxa)) +
         geom_bar(stat="identity") +
         scale_fill_manual(values=qual_col_16S) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  scale_y_continuous(expand = c(0.0, 0.0))

ggplot(prep_mgs_class_long, aes(x=variable, y=value, fill=taxa)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=qual_col_mgs) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  scale_y_continuous(expand = c(0.0, 0.0))