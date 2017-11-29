library(cowplot)

setwd("/home/gavin/github_repos/CD_RF_microbiome/")

otus_16S_orig_disease <- readRDS("RF_RDS_output/16S/biscuit_16S_otus_disease.rds")
species_16S_orig_disease <- readRDS("RF_RDS_output/16S/biscuit_16S_species_disease.rds")
genus_16S_orig_disease <- readRDS("RF_RDS_output/16S/biscuit_16S_genus_disease.rds")
family_16S_orig_disease <- readRDS("RF_RDS_output/16S/biscuit_16S_family_disease.rds")
order_16S_orig_disease <- readRDS("RF_RDS_output/16S/biscuit_16S_order_disease.rds")
class_16S_orig_disease <- readRDS("RF_RDS_output/16S/biscuit_16S_class_disease.rds")
phylum_16S_orig_disease <- readRDS("RF_RDS_output/16S/biscuit_16S_phylum_disease.rds")

otus_16S_unrarified_disease <- readRDS("RF_RDS_output/16S_unrarified_CLR/biscuit_otu_unrarified_disease_prep.rds")
species_16S_unrarified_disease <- readRDS("RF_RDS_output/16S_unrarified_CLR/biscuit_species_unrarified_disease_prep.rds")
genus_16S_unrarified_disease <- readRDS("RF_RDS_output/16S_unrarified_CLR/biscuit_genus_unrarified_disease_prep.rds")
family_16S_unrarified_disease <- readRDS("RF_RDS_output/16S_unrarified_CLR/biscuit_family_unrarified_disease_prep.rds")
order_16S_unrarified_disease <- readRDS("RF_RDS_output/16S_unrarified_CLR/biscuit_order_unrarified_disease_prep.rds")
class_16S_unrarified_disease <- readRDS("RF_RDS_output/16S_unrarified_CLR/biscuit_class_unrarified_disease_prep.rds")
phylum_16S_unrarified_disease <- readRDS("RF_RDS_output/16S_unrarified_CLR/biscuit_phylum_unrarified_disease_prep.rds")

otus_16S_orig_response <- readRDS("RF_RDS_output/16S/biscuit_16S_otus_response.rds")
species_16S_orig_response <- readRDS("RF_RDS_output/16S/biscuit_16S_species_response.rds")
genus_16S_orig_response <- readRDS("RF_RDS_output/16S/biscuit_16S_genus_response.rds")
family_16S_orig_response <- readRDS("RF_RDS_output/16S/biscuit_16S_family_response.rds")
order_16S_orig_response <- readRDS("RF_RDS_output/16S/biscuit_16S_order_response.rds")
class_16S_orig_response <- readRDS("RF_RDS_output/16S/biscuit_16S_class_response.rds")
phylum_16S_orig_response <- readRDS("RF_RDS_output/16S/biscuit_16S_phylum_response.rds")

otus_16S_unrarified_response <- readRDS("RF_RDS_output/16S_unrarified_CLR/biscuit_otu_unrarified_response_prep.rds")
species_16S_unrarified_response <- readRDS("RF_RDS_output/16S_unrarified_CLR/biscuit_species_unrarified_response_prep.rds")
genus_16S_unrarified_response <- readRDS("RF_RDS_output/16S_unrarified_CLR/biscuit_genus_unrarified_response_prep.rds")
family_16S_unrarified_response <- readRDS("RF_RDS_output/16S_unrarified_CLR/biscuit_family_unrarified_response_prep.rds")
order_16S_unrarified_response <- readRDS("RF_RDS_output/16S_unrarified_CLR/biscuit_order_unrarified_response_prep.rds")
class_16S_unrarified_response <- readRDS("RF_RDS_output/16S_unrarified_CLR/biscuit_class_unrarified_response_prep.rds")
phylum_16S_unrarified_response <- readRDS("RF_RDS_output/16S_unrarified_CLR/biscuit_phylum_unrarified_response_prep.rds")


disease_compare_df <- data.frame(matrix(NA, nrow=14, ncol=3))

colnames(disease_compare_df) <- c("Accuracy", "Workflow", "Level")

disease_compare_df$Accuracy <- 
c(as.numeric(as.character(otus_16S_orig_disease$summary$LOOCV.Accuracy)),
  as.numeric(as.character(species_16S_orig_disease$summary$LOOCV.Accuracy)),
  as.numeric(as.character(genus_16S_orig_disease$summary$LOOCV.Accuracy)),
  as.numeric(as.character(family_16S_orig_disease$summary$LOOCV.Accuracy)),
  as.numeric(as.character(order_16S_orig_disease$summary$LOOCV.Accuracy)),
  as.numeric(as.character(class_16S_orig_disease$summary$LOOCV.Accuracy)),
  as.numeric(as.character(phylum_16S_orig_disease$summary$LOOCV.Accuracy)),
  
  as.numeric(as.character(otus_16S_unrarified_disease$summary$LOOCV.Accuracy)),
  as.numeric(as.character(species_16S_unrarified_disease$summary$LOOCV.Accuracy)),
  as.numeric(as.character(genus_16S_unrarified_disease$summary$LOOCV.Accuracy)),
  as.numeric(as.character(family_16S_unrarified_disease$summary$LOOCV.Accuracy)),
  as.numeric(as.character(order_16S_unrarified_disease$summary$LOOCV.Accuracy)),
  as.numeric(as.character(class_16S_unrarified_disease$summary$LOOCV.Accuracy)),
  as.numeric(as.character(phylum_16S_unrarified_disease$summary$LOOCV.Accuracy)))

disease_compare_df$Workflow <- c(rep(x = "Rarified", 7), rep(x = "CLR", 7))

disease_compare_df$Level <- c("OTUs", "Species", "Genus", "Family", "Order", "Class", "Phylum",
                              "OTUs", "Species", "Genus", "Family", "Order", "Class", "Phylum")

disease_compare_df$Workflow <- as.factor(disease_compare_df$clr)
disease_compare_df$Level <- factor(disease_compare_df$Level, 
                                   levels=c("OTUs", "Species", "Genus", "Family", "Order", "Class", "Phylum"))

disease_barplot <- ggplot(disease_compare_df, aes(Level, Accuracy)) +   
  geom_bar(aes(fill = Workflow), position = "dodge", stat="identity") +
  scale_fill_manual(values=c("black", " dark grey")) +
  coord_cartesian(ylim = c(0.0, 1.0))

response_compare_df <- data.frame(matrix(NA, nrow=14, ncol=3))

colnames(response_compare_df) <- c("Accuracy", "Workflow", "Level")

response_compare_df$Accuracy <- 
  c(as.numeric(as.character(otus_16S_orig_response$summary$LOOCV.Accuracy)),
    as.numeric(as.character(species_16S_orig_response$summary$LOOCV.Accuracy)),
    as.numeric(as.character(genus_16S_orig_response$summary$LOOCV.Accuracy)),
    as.numeric(as.character(family_16S_orig_response$summary$LOOCV.Accuracy)),
    as.numeric(as.character(order_16S_orig_response$summary$LOOCV.Accuracy)),
    as.numeric(as.character(class_16S_orig_response$summary$LOOCV.Accuracy)),
    as.numeric(as.character(phylum_16S_orig_response$summary$LOOCV.Accuracy)),
    
    as.numeric(as.character(otus_16S_unrarified_response$summary$LOOCV.Accuracy)),
    as.numeric(as.character(species_16S_unrarified_response$summary$LOOCV.Accuracy)),
    as.numeric(as.character(genus_16S_unrarified_response$summary$LOOCV.Accuracy)),
    as.numeric(as.character(family_16S_unrarified_response$summary$LOOCV.Accuracy)),
    as.numeric(as.character(order_16S_unrarified_response$summary$LOOCV.Accuracy)),
    as.numeric(as.character(class_16S_unrarified_response$summary$LOOCV.Accuracy)),
    as.numeric(as.character(phylum_16S_unrarified_response$summary$LOOCV.Accuracy)))

response_compare_df$Workflow <- c(rep(x = "Rarified", 7), rep(x = "CLR", 7))

response_compare_df$Level <- c("OTUs", "Species", "Genus", "Family", "Order", "Class", "Phylum",
                              "OTUs", "Species", "Genus", "Family", "Order", "Class", "Phylum")

response_compare_df$Workflow <- as.factor(response_compare_df$clr)
response_compare_df$Level <- factor(response_compare_df$Level, 
                                   levels=c("OTUs", "Species", "Genus", "Family", "Order", "Class", "Phylum"))

response_barplot <- ggplot(response_compare_df, aes(Level, Accuracy)) +   
  geom_bar(aes(fill = Workflow), position = "dodge", stat="identity") +
  scale_fill_manual(values=c("black", " dark grey")) +
  coord_cartesian(ylim = c(0.0, 1.0))
  
plot_grid(disease_barplot, response_barplot, labels = c("A", "B"))
