library("zCompositions")

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
  return(in_table[which(rowSums(in_table > 0) > cutoff_pro*ncol(in_table)), ])
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