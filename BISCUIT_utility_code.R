library("zCompositions")
library("plyr")

read_spf_by_level <- function(infile, level) {

  input_spf <- read.table(infile,
                          header=T,
                          sep="\t",
                          stringsAsFactors = TRUE)
  
  colnames(input_spf) <- c("Kingdom", "Phylum" , "Class" , "Order" , "Family" , "Genus" , "Species", 
                           colnames(input_spf)[8:ncol(input_spf)])
  
  
  return(collapse_taxa(input_spf, level))
  
}

collapse_taxa <- function(table , level) {
  
  ranks <- c("Phylum" , "Class" , "Order" , "Family" , "Genus" , "Species")
  
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


read_legacy_biom <- function(in_file) {

	return(read.table(in_file, 
                    header=T, 
                    sep="\t", 
                    stringsAsFactors = FALSE, 
                    quote="", 
                    comment.char = "",
                    row.names=1,
                    skip=1))
}

remove_rare_rows <- function(in_table, cutoff_pro) {
  return(in_table[which(rowSums(in_table > 0) > ceiling(cutoff_pro*ncol(in_table))), ])
}

CZM_and_clr <- function(in_table) {
  in_table_impute <- cmultRepl(t(in_table), label=0, method="CZM")
  return(as.data.frame(apply(in_table_impute, 1, function(x){log(x) - mean(log(x))})))
}

write_RF_prep <- function(in_table, out_file, column2add, meta_df) {
  
  in_table_t <- as.data.frame(t(in_table))
  
  if(column2add == "disease") {
    in_table_t$disease <- meta_df[rownames(in_table_t), "disease"]
  } else if(column2add == "response") {
    in_table_t$response <- meta_df[rownames(in_table_t), "response"]
  } else {
   stop("column2add needs to be disease or response") 
  }
  
  write.table(in_table_t,
              file=out_file, 
              row.names=TRUE,
              col.names=NA, 
              quote = FALSE , 
              sep="\t")
}

# Function to find which of one (and only one) substrings match a string.
multiple_match_single <- function(substrs, str2match) {
  matched <- NULL
  
  for(s in substrs) {
    if(length(grep(s, str2match) > 0)) {
      if(is.null(matched)) {
        matched <- s 
      } else {
        # If already made a match then return an error.
        matched_strings <- paste(c(matched, s), collapse = " ")
        stop(paste("Both of these substrings matched:", matched_strings, sep=" ")) 
      }
    }
  }
  if(! is.null(matched)) {
    return(matched)
  } else {
    # If none matched then return error.
    stop(paste("No substrings matched in", str2match, sep=" ")) 
  }
}

rf_summary_from_rds_files <- function(rds_files) {
  
  summary_table <- data.frame(matrix(c(NA), nrow = 0, ncol = 14))
  
  for(rds in rds_files) {
    rds_in <-  readRDS(rds)$summary
    rds_in$sequencing <- multiple_match_single(c("16S", "mgs"), rds_in$file)
    rds_in$dataset <- multiple_match_single(
      c("otu", "strain", "species", "genus", "family", "order", "class", "phylum", "ko", "pathway", "module"), 
      tolower(rds_in$file))
    rds_in$trait <- multiple_match_single(
      c("disease", "response"), 
      tolower(rds_in$file))
    
    summary_table <- rbind(summary_table, rds_in)
  }
  
  return(summary_table)
}

read_in_varImp <- function(inRDS) {
  
  imp_df <- as.data.frame(readRDS(inRDS)$mod$importance)
  
  imp_df$features <- rownames(imp_df)
  
  imp_df_sorted <- arrange(imp_df, desc(MeanDecreaseAccuracy))
  
  rownames(imp_df_sorted) <- imp_df_sorted$features
  
  imp_df_sorted <- imp_df_sorted[, -which(colnames(imp_df_sorted) %in% "features")]
  
  return(imp_df_sorted)
}


read_write_varImp <- function(in_rds, out_table) {
  
 input_varImp <- read_in_varImp(in_rds)
 
 
  write.table(input_varImp,
            file=out_table, 
            row.names=TRUE,
            col.names=NA, 
            quote = FALSE , 
            sep="\t")
  
}


get_top3_sig_varImp_vals <- function(in_rds, sample_order, fix_ids=FALSE, id_mapping) {
  
  in_results <- readRDS(in_rds)
  
  imp <- as.data.frame(in_results$mod$importance)
  
  imp$features <- rownames(imp)
  
  imp_sorted <- arrange(imp, desc(MeanDecreaseAccuracy))
  
  top3_features <- imp_sorted[1:3, "features"]
  
  tab <- in_results$loocv_mod$trainingData
  
  if(fix_ids) {
    old_ids <- rownames(tab)
    old_ids <- gsub("_.+$", "", old_ids)
    id_map <- read.table(id_mapping, header=T, sep="\t", stringsAsFactors = FALSE)
    rownames(id_map) <- id_map$BISCUIT
    rownames(tab) <- id_map[old_ids, "sample_id"]
  } 
  
  return(tab[sample_order, top3_features])

}


get_sig_and_direction <- function(input_table, two_classes, varImp_order) {
  
  col2return <- c()
  
  # get last column, which is what is being classified
  class_col <- input_table[ , ncol(input_table) , drop = F ]
  
  # remove from original table
  input_table_features <- input_table[ , varImp_order  ]
  
  # figure out rows that correspond to class1 and class2
  class1_rows <- which(class_col == two_classes[1] )
  class2_rows <- which(class_col == two_classes[2] )
  
  for ( i in 1:ncol(input_table_features) ) {
    
    feature_class1 <- input_table_features[ class1_rows , i ]
    feature_class2 <- input_table_features[ class2_rows , i ]
    
    sig <- wilcox.test( feature_class1 , feature_class2 )$p.value
    mean_diff <- mean( feature_class1 ) - mean( feature_class2 )
    
    if ( sig >= 0.05 ) {
      col2return <- c( col2return , "grey" )
    } else {
      
      if ( mean_diff > 0 ) {
        col2return <- c( col2return , "red" )
      } else if ( mean_diff < 0 ) {
        col2return <- c( col2return , "blue" )
      } else if ( mean_diff == 0 ) {
        col2return <- c( col2return , "black" ) ### unlikely, but just in case
      }
    }
  }
  
  return( col2return )
}

plot_varImp <- function( imp , input , classes, x_val , pdf=FALSE, pdf_name, output_table=FALSE , table_name ) {
  
  imp <- as.data.frame( imp )
  
  imp$features <- rownames( imp )
  
  imp_sorted <- arrange( imp  , MeanDecreaseAccuracy  )
  
  bar_col <- get_sig_and_direction( input_table=input , two_classes=classes , varImp_order = imp_sorted$features )
  
  if ( pdf == TRUE ) {
    
    pdf( pdf_name, width=6 , height=9 )
    
  }
  
  barplot( imp_sorted$MeanDecreaseAccuracy * 100 , xlab="Variable Importance\n(Mean Decrease in Accuracy)" , horiz=TRUE  , names.arg = imp_sorted$features , las=2 , xlim=x_val, cex.names=0.75, col=bar_col)
  
  abline(v=0 , lwd=1)
  
  if ( pdf == TRUE ) {
    dev.off()
  }
  
  if (output_table==TRUE){
    write.table(x=imp_sorted[, -which(colnames(imp_sorted)=="features")], file = table_name, row.names = imp_sorted$features , col.names=NA, sep="\t", quote=FALSE )
  }
  
}