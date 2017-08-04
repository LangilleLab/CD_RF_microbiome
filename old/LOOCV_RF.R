#!/usr/bin/env Rscript

### Written by Gavin Douglas (gavin.douglas@dal.ca)

version <- 1.0

# load optpase library
library(optparse)

option_list <- list(

		  make_option( c( "-i" , "--input"), type="character", default=NA, 
              help = "Input table (required). Tab-delimited table with rows as samples and columns as features. 
		One column should be what you want to classify using all other features. All pre-processing should be done in advance." , metavar="path"),

		  make_option( c( "-o" , "--output_prefix"), type="character", default=NA, 
              help = "Prefix for 3 output files (required). A table of accuracy values, a variable importance table and a file containing
        the observed median accuracy, median kappa and raw P-value based on the null distribution will be outputted" , 
              metavar="path"),

		  make_option( c(  "-c" , "--col" ) , type="character" , default=NA ,
			  help = "Column in table to classify (required)." ,
			  metavar = "string" ) ,

		  make_option( c(  "--num_rep" ) , type="numeric" , default=11 ,
			  help = "Number of RF replicates to take median of for individual accuracy estimates [default = %default].
		These RF replicates are run on the same input." ,
			  metavar = "number" ) ,

		  make_option( c(  "--num_rand" ) , type="numeric" , default=1000 ,
			  help = "Number of RF replicates to run on the randomized data to build null distribution [default = %default]" ,
			  metavar = "number" ) ,

		  make_option( c(  "--ntree" ) , type="numeric" , default=501 ,
			  help = "Number of bootstrapped trees to use in each RF model except when calculating variable importance [default = %default]" ,
			  metavar = "number" ) ,

		  make_option( c(  "--ntree_varImp" ) , type="numeric" , default=10001 ,
			  help = "Number of bootstrapped trees to use when calculating variable importance [default = %default]" ,
			  metavar = "number" ) ,

		  make_option( c(  "--mtry_pro" ) , type="numeric" , default=0.2 ,
			  help = "Proportion of variables to randomly sample as candidates at each split to a minimum of 1 [default = %default]" ,
			  metavar = "number" ) ,

		  make_option( c(  "--threads" ) , type="numeric" , default=1 ,
			  help = "Number of threads to use [default = %default]" ,
			  metavar = "number" ) ,

		  make_option( c(  "--seed" ) , type="numeric" , default=812 ,
			  help = "Starting random seed (important for reproducibility) [default = %default]" ,
			  metavar = "number" ) ,

		  make_option( c( "-v" , "--version" ) , action = "store_true" , type="logical" , default=FALSE ,
			  help = "Print out version number and exit", 
			  metavar = "string" ) 

		 )

opt_parser <- OptionParser( option_list=option_list , 
usage = "%prog [options] --input PATH --output_prefix PATH --col STRING

This Rscript is for running random forest (RF) on small datasets. Leave-one-out cross-validation is used on all samples to calculate accuracy. 
Three output files will be generated: a table of accuracies (per RF replicate), a table of variable importance (per input feature) and a file 
containing the observed median accuracy, median kappa and raw P-value based on the null distribution. Outputted accuracies are the median of a 
user-specified number of RF replicates. The actual accuracy will be outputted followed by accuracies inferred when the column of interest is 
permuted (which is used as the null distribution).", )

opt <- parse_args( opt_parser )

check_and_read_input <- function( opt_list , table ) {

	### function to check basic input parameters and options
	### also will read in input table and return it

	# check for version flag
	if( opt_list$version ) {
        cat( "Version" , version , "\n")
        options_tmp <- options(show.error.messages=FALSE)
    	on.exit(options(options_tmp))
        stop()
	}

	# check for required arguments
	if ( is.na( opt_list$input ) ) {
        stop( "Path to input table needs to be specified with --input\nType \"LOOCV_RF.R --help\" for help." )
	} else if ( is.na( opt_list$output_prefix )) {
        stop( "Prefix of output files needs to be specified with --output_prefix\nType \"LOOCV_RF.R --help\" for help." )
	} else if ( is.na( opt_list$col )) {
        stop( "Column to classify needs to be specified with --col\nType \"LOOCV_RF.R --help\" for help." )
    }

    table_in <-  read.table( table, header=T, sep="\t" , row.names=1 )

    # check that specified column is in table
    if ( ! ( opt_list$col %in% colnames(table_in) ) ) {
    	col_err <- paste( "Specificed column" , opt_list$col , "is not present in input file")
    	stop( col_err )
    }

    # check that mtry_pro is > 0 and <= 1
    if ( opt_list$mtry_pro <= 0 | opt_list$mtry_pro > 1 ) {
    	mtry_err <- paste( "Option --mtry_pro needs to be greater than 0 and less or equal to 1, not" , opt_list$mtry_pro , sep=" " )
    	stop( mtry_err )
    }

    return( table_in )
}

run_rep_RF <- function( table , n , set_seed , type )	{

	### This function runs a RF a set number of times and uses leave-one-out cross validation (LOOCV) to estimate accuracy. 
	### The mtry, median accuracy, median kappa and type ("actual" or "random") will be returned 

	rep_results <- data.frame( matrix( c(NA) , nrow=0 , ncol=3 ) )
	colnames(rep_results ) <- c( "mtry" , "Accuracy" , "Kappa" )

	for ( i in 1:n ) {

		# each replicate is based on a different random seed
		set.seed( set_seed + i ) 

		rf_fit <- train( classify ~ . , data=table, method="rf", trControl=fit_control , 
			ntree=opt$ntree , tuneGrid=data.frame( mtry=mtry_count  ) , importance=FALSE ) 

		# note if you did mtry tuning instead you would have to select the row that corresponded to the best tuning
		# in this case there is only 1 row of results anyway per replicate
		results_tmp <- rf_fit$results[ 1 , , drop=F ]
		rep_results <- rbind( rep_results , results_tmp )

	}

	median_result <- data.frame( matrix( c(NA) , nrow=1 , ncol=4 ) )
	colnames( median_result ) <- c( "mtry" , "Accuracy" , "Kappa" , "type" )
	median_result$mtry <- unique( rep_results$mtry )
	median_result$Accuracy <- median( rep_results$Accuracy )
	median_result$Kappa <- median( rep_results$Kappa )
	median_result$type <- type

	return( median_result )

}

# read in table and check parameters/options
input <- check_and_read_input( opt_list = opt , table = opt$input )

### load other required libraries (which take longer to load)
library(caret) 
library(doMC)
registerDoMC(cores = opt$threads )

# RF settings to be used for all models (except when calculating variable importance)
fit_control <- trainControl( method = "LOOCV" ) 

# make new column called "classify", which is column of interest
input$classify <- input[ , opt$col ] 

# remove original column of interest
input <- input[ , -c(which(colnames(input)== opt$col))]

# mtry parameter is taken to be mtry_pro*(# of features) OR 1, whichever is higher
mtry_count <- max( floor( ( ncol(input) - 1 ) * opt$mtry_pro ) , 1 )

#### First major step is to get table of variable importance, which can be produced by a large # of trees
# set random seed
set.seed( opt$seed ) 

rf_fit_varImp <- train( classify ~ . , data=input, method="rf", ntree=opt$ntree_varImp , 
 				tuneGrid=data.frame( mtry=mtry_count ) , importance=TRUE  )

varImp_scaled <- varImp(rf_fit_varImp , scale=TRUE )$importance[, 1, drop=F ]
varImp_unscaled <- varImp(rf_fit_varImp , scale=FALSE )$importance[, 1, drop=F ]

# just to make sure that the features are in the same order, although they should be anyway:
varImp_unscaled <- varImp_unscaled[ rownames( varImp_scaled ) , , drop=F ]

# only output 1 column of varImp (which are identical in caret output)
varImp_scaled <- varImp_scaled[ , 1 , drop=F ]
varImp_unscaled <- varImp_unscaled[ , 1 , drop=F ]

# rename column to reflect the type of varImp 
colnames( varImp_scaled ) <- "scaled"
colnames( varImp_unscaled ) <- "unscaled"

varImp <- cbind( varImp_scaled , varImp_unscaled )
varImp_output <- paste( opt$output_prefix , "varImp.txt" , sep="_" )
write.table( file = varImp_output , x=varImp, row.names = T, sep="\t", quote=F , col.names=NA )


#### Second major step is to get distribution of accuracy results
results <- data.frame( matrix( c(NA) , nrow=0 , ncol=4 ) )
colnames(results ) <- c( "mtry" , "Accuracy" , "Kappa" , "type" )

### Get actual results:
results <- rbind( results , run_rep_RF( table=input , n=opt$num_rep , set_seed=opt$seed , type="actual" ) )


### Get null distribution:
for (i in 1:opt$num_rand) {
  
  	# different random seed for each re-sampling
	opt$seed <- opt$seed + i*1000
	set.seed( opt$seed )

	# randomize "classify" column
	input_ran <- input
	input_ran$classify <- sample( input_ran$classify )

	results <- rbind( results , run_rep_RF( table=input_ran , n=opt$num_rep , set_seed=opt$seed , type="random" ) )

}

### Output full table of accuracies and kappa values
output_table <- paste( opt$output_prefix , "acc.txt" , sep="_" ) 
write.table(x = results, file = output_table , quote=FALSE , row.names=F  , sep=" " )

### Figure out median accuracy, median kappa and raw P-value
actual_result <- results[ which(results$type == "actual") , , drop=FALSE ]
random_results <- results[ which(results$type == "random") , , drop=FALSE ]

actual_result$raw_p <- length( which( random_results$Accuracy >= actual_result$Accuracy ) ) / nrow( random_results )

### Output just "actual" line with raw P-value as well:
actual_out <- paste( opt$output_prefix , "sig_test.txt" , sep="_" )
write.table(x = actual_result, file = actual_out , quote=FALSE , row.names=F  , sep=" " )
