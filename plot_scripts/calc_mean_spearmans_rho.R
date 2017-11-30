setwd("/home/gavin/github_repos/CD_RF_microbiome/")

# First compare concordance between 16S and MGS samples.
mgs_to_16S_genera <- read.table("mgs_16S_genus_names.txt", header=T, sep="\t", stringsAsFactors = FALSE)

biscuit_genus <- read_spf_by_level("raw_data/BISCUIT_16S/otu_table_w_tax_BISCUIT.spf", "Genus")
rownames(biscuit_genus) <- gsub(".+g__", "g__" ,rownames(biscuit_genus)) 
biscuit_genus <- biscuit_genus[mgs_to_16S_genera$X16S,]

#biscuit_genus_CLR <- read.table("test_genus_CLR.txt", header=T, sep="\t", row.names=1)
#biscuit_genus_CLR <- data.frame(t(biscuit_genus_CLR[,-ncol(biscuit_genus_CLR)]))
#rownames(biscuit_genus_CLR) <- gsub(".+g__", "g__" ,rownames(biscuit_genus_CLR)) 
#biscuit_genus_CLR <- biscuit_genus_CLR[mgs_to_16S_genera$X16S,colnames(biscuit_genus)]
#genus_MGS_vs_16S_CLR <- c()

genus_MGS <- read.table("raw_data//BISCUIT_MGS/metaphlan2_by_level/metaphlan2_out_merged_genus.tsv",
                        header=T, sep="\t", stringsAsFactors = FALSE, row.names=1)
rownames(genus_MGS) <- gsub(".+g__", "g__" ,rownames(genus_MGS)) 
genus_MGS <- genus_MGS[mgs_to_16S_genera$mgs,colnames(biscuit_genus)]

genus_MGS_vs_16S <- c()

sampl_names <- c()
for(i in 1:ncol(genus_MGS)) {
  genus_MGS_vs_16S <- c(genus_MGS_vs_16S, cor.test(genus_MGS[,i], biscuit_genus[,i], method="spearman")$estimate)
  sampl_names <- c(sampl_names, colnames(genus_MGS)[i])
}

names(genus_MGS_vs_16S) <- sampl_names

print(genus_MGS_vs_16S)

mean(genus_MGS_vs_16S)
sd(genus_MGS_vs_16S)