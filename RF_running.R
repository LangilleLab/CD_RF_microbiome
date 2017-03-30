### The below commands were used to run random forest on each input table. 
### Significantly low out-of-bag errors were identified based on a permutation test
### Accuracies were determined by re-running models with leave-one-out cross-validation

# Load required packages
library("randomForest") # for running RF and getting out of bag (OOB) error
library("rfUtilities") # for running OOB error significance test (requires re-running RF with permuted disease/response column)
library("caret") # for running leave-one-out cross-validation to estimate accuracy (wraps randomForest and a number of other packages)

library("doMC") # for multithreading when possible
registerDoMC( cores = 30 )

# This is where the input tables were located locally
setwd("./RF_input_tables/")

run_all_RF <- function( filename , num_tree , num_perm , set_seed = 712 , outfile ) {
  
  # Read in table. Column to classify is last column
  t <- read.table( filename , header=T , sep="\t" , row.names=1 )
  
  # Set random seed
  set.seed( set_seed )
  
  # Below is the default mtry parameter choice for randomForest
  mtry_count <- floor(sqrt(ncol(t[,1:(ncol(t)-1)])))
  
  # Run random forest, calculate importance and proximities as well
  rf_mod <- randomForest( x=t[,1:(ncol(t)-1)] , y=t[ , ncol(t)] , ntree=num_tree , 
                          proximity = TRUE , importance = TRUE , mtry=mtry_count )
  
  # Run leave-one-out cross-validation to estimate model accuracy
  fit_control <- trainControl( method = "LOOCV" )
  rf_mod_loocv <- train( x=t[,1:(ncol(t)-1)] , y=t[ , ncol(t)] , method="rf", ntree=num_tree , 
                  tuneGrid=data.frame( mtry=mtry_count ) , trControl=fit_control )
  
  # Test whether observed OOB error is significant
  rf_sig <- rf.significance( x=rf_mod ,  xdata=t[,1:(ncol(t)-1)] , nperm=num_perm , ntree=num_tree ) 

  # Make dataframe with summary data for RF
  results <- data.frame( matrix( c( filename , ncol( rf_mod_loocv$trainingData ) - 1 , 
                                      rf_mod$mtry , rf_mod$ntree , rf_mod_loocv$results$Accuracy , 
                                      rf_sig$test.OOB , median( rf_sig$RandOOB ) , rf_sig$pValue ) , nrow=1, ncol=8 ) )
  
  colnames(results) <- c( "file" , "num.features" , "mtry" , "num.trees" , "LOOCV.Accuracy" , "Model.OOB.error" , "Median.Rand.OOB.error" , "P.value" )
  
  
  rf_out <- list( "mod" = rf_mod , "sig" = rf_sig , "loocv_mod" = rf_mod_loocv , "summary" = results )
  
  save( rf_out , file = outfile )
  
  # Return list containing RF model, model result when running LOOCV and significance test results and summary of main metrics
  # return( rf_out )
  
}

### Looped through 3 different types of labels to loop through all input files automatically (except for KOs)

sequencing <- c( "16S" , "MGS" )

# Note that KOs are not included since they are run with a lower number of trees below
categories <- c( "phylum" , "class" , "order" , "family" , "genus" , "species", "otus" , "strains" , "modules" , "pathways" )

traits <- c( "disease" , "response" )

for ( s in sequencing ) {

  for ( c in categories ) {

    for ( t in traits ) {

     # specified input and output files based on these labels
     f <- paste( s , c , t , "input.txt" , sep="_" )
     out <- paste( "../RF_obj_out/", s , sep="" )
     out <- paste( out , c , t , "RF_out.rda" , sep="_"  )

     # if input file exists then ran run_all_RF function, which outputs an rf_out object for each input file    
     if ( file.exists( f ) ) {
        run_all_RF( filename = f , num_tree = 10001 , num_perm = 1000 , outfile=out  ) 
      }

    }

  }

}

### Run RFs for KO data with fewer trees since there are thousands of features
KOs_16S_response_rf_out <- run_all_RF( filename = "KOs/16S_KOs_response_input.txt" , num_tree = 501 , num_perm = 1000 , outfile= "../RF_obj_out/KOs/16S_KOs_response_RF_out.rda" )
KOs_MGS_response_rf_out <- run_all_RF( filename = "KOs/MGS_KOs_response_input.txt" , num_tree = 501 , num_perm = 1000 , outfile= "../RF_obj_out/KOs/MGS_KOs_response_RF_out.rda" )

KOs_16S_disease_rf_out <- run_all_RF( filename = "KOs/16S_KOs_disease_input.txt" , num_tree = 501 , num_perm = 1000 , outfile= "../RF_obj_out/KOs/16S_KOs_disease_RF_out.rda" )
KOs_MGS_disease_rf_out <- run_all_RF( filename = "KOs/MGS_KOs_disease_input.txt" , num_tree = 501 , num_perm = 1000 , outfile= "../RF_obj_out/KOs/MGS_KOs_disease_RF_out.rda" )

