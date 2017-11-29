library(stringi)

setwd("/home/gavin/github_repos/CD_RF_microbiome/")

rf_summary_tab <- read.table("raw_summary_out/RF_16S_mgs_summary.txt",
                             header=T,
                             stringsAsFactors = FALSE,
                             sep="\t")

# Capitalize all levels
rf_summary_tab$dataset <- stri_trans_totitle(rf_summary_tab$dataset)

# Fix up KO and OTU title and mgs -> MGS
rf_summary_tab[which(rf_summary_tab$dataset == "Ko") , "dataset"] <- "KO"
rf_summary_tab[which(rf_summary_tab$dataset == "Otu") , "dataset"] <- "OTU"
rf_summary_tab[which(rf_summary_tab$sequencing == "mgs") , "sequencing"] <- "MGS"

# Add in richness and GRS values.
num_otus_disease <- c("num otus manual",
                      1,
                      NA,
                      NA,
                      0.71,
                      NA,
                      NA,
                      0.009,
                      "16S",
                      "num_OTUs",
                      "disease")

num_otus_response <- c("num otus manual",
                      1,
                      NA,
                      NA,
                      0.33,
                      NA,
                      NA,
                      0.282,
                      "16S",
                      "num_OTUs",
                      "response")

grs_disease <- c("grs manual",
                      1,
                      NA,
                      NA,
                      0.625,
                      NA,
                      NA,
                      0.008,
                      "MGS",
                      "GRS",
                      "disease")

grs_response <- c("grs manual",
                 1,
                 NA,
                 NA,
                 0.55,
                 NA,
                 NA,
                 0.736,
                 "MGS",
                 "GRS",
                 "response")

rf_summary_tab_all <- rbind(rf_summary_tab,
                            num_otus_disease,
                            num_otus_response,
                            grs_disease,
                            grs_response)
                            

rf_summary_tab_all$sequencing <- factor(rf_summary_tab_all$sequencing, levels=c("MGS", "16S"))

rf_summary_tab_all$dataset <- factor(rf_summary_tab_all$dataset, 
                                     levels=rev(c("num_OTUs", "GRS", "OTU", "Strain", "Species", "Genus", "Family", "Order", "Class", "Phylum", "KO", "Pathway", "Module")))

rf_summary_disease <- rf_summary_tab_all[which(rf_summary_tab_all$trait == "disease") , ]
rf_summary_response <- rf_summary_tab_all[which(rf_summary_tab_all$trait == "response") , ]

rf_summary_disease_sorted <- arrange(rf_summary_disease , sequencing , dataset)
rf_summary_response_sorted <- arrange(rf_summary_response , sequencing , dataset)

# Increase left margin so name will fit and make 2 panels per plot
par(mar=c(5.1,6.1,4.1,2.1) , mfrow=c(1,2))
palette(c("green3" , "black"))

# Make disease plot
barplot(height = as.numeric(rf_summary_disease_sorted$LOOCV.Accuracy)*100, names.arg = rf_summary_disease_sorted$dataset , horiz = T , las=1, xlim=c(0,100), 
            col=rf_summary_disease_sorted$sequencing, xlab="Disease State Classification Accuracy" )

# Make response plot
barplot( height = as.numeric(rf_summary_response_sorted$LOOCV.Accuracy)*100, names.arg = rf_summary_response_sorted$dataset , horiz = T , las=1, xlim=c(0,100), 
            col=rf_summary_response_sorted$sequencing , xlab="Treatment Response Classification Accuracy" )