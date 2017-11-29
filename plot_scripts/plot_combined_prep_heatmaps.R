library("Hmisc")
library("reshape2")
setwd("/home/gavin/github_repos/CD_RF_microbiome/")


plot_combined_corr_coefficent <- function( combined_table ) {
  
  combined_table_corr <- rcorr( as.matrix(combined_table) , type="spearman" )
  combined_table_corr_melt = melt(combined_table_corr$r)
  
  sig_corr <- c()
  
  for (cname in colnames(combined_table_corr$P)) {
    for (rname in rownames(combined_table_corr$P)) {
      feature_compare <- paste( sort( c(cname,rname) ) , collapse = ";" )
      if ( ! is.na(combined_table_corr$P[rname , cname]) && (combined_table_corr$P[rname , cname] < 0.05 ) ){
        sig_corr <- c( sig_corr , feature_compare )
      }
    }
  }
  sig_corr <- unique(sig_corr)
  
  past_comparisons <- c()
  combined_table_corr_melt_culled <- combined_table_corr_melt
  
  for(i in 1:nrow(combined_table_corr_melt)) {
    
    feature1 <- levels(combined_table_corr_melt$Var1)[combined_table_corr_melt[ i , "Var1" ]]
    feature2 <- levels(combined_table_corr_melt$Var2)[combined_table_corr_melt[ i , "Var2" ]]
    
    feature_compare <- paste( sort( c(feature1,feature2) ) , collapse = ";" )
    
    if ( ( feature1 == feature2 ) || ( feature_compare %in% past_comparisons ) || ( ! ( feature_compare %in% sig_corr) ) ) {
      combined_table_corr_melt_culled[ i, "value" ] <- NA
    } else {
      past_comparisons <- c(past_comparisons, feature_compare)
    }
  }
  
  g <- ggplot(combined_table_corr_melt_culled, aes(x=Var1, y=Var2, fill=value)) +
    geom_tile() + xlab("") + ylab("")+ theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1), panel.background = element_blank())+scale_fill_continuous(na.value="white", low="blue", high="red")
  
  return(g)
}

disease_combined_prep <- read.table("prep_data/combined_prep/combined_disease_prep.txt",
                                    header=T,
                                    row.names=1,
                                    stringsAsFactors = FALSE,
                                    sep="\t")

disease_combined_raw <- disease_combined_prep[, -which(colnames(disease_combined_prep) %in% "disease")]



response_combined_prep <- read.table("prep_data/combined_prep/combined_response_prep.txt",
                                    header=T,
                                    row.names=1,
                                    stringsAsFactors = FALSE,
                                    sep="\t")

response_combined_raw <- response_combined_prep[, -which(colnames(response_combined_prep) %in% "response")]

disease_combined_raw_heatmap <- plot_combined_corr_coefficent(disease_combined_raw)
response_combined_raw_heatmap <- plot_combined_corr_coefficent(response_combined_raw)

plot(response_combined_raw_heatmap)
plot(disease_combined_raw_heatmap)