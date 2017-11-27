# Read in package to read in command-line options.
library("optparse")

option_list <- list(
  
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Input table (required)." , metavar="path"),
  
  make_option(c("--seed"), type="integer", default=NULL,
              help="Random seed to make command reproducible.", metavar = "integer"),
  
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Output .rds file (required)." , metavar="path"),
  
  make_option(c("-t", "--threads"), type="integer", default=1,
              help="Number of threads (default: 1)",
              metavar="integer"),
  
  make_option(c("--num_tree"), type="integer", default=10001,
              help="Number of trees to use for random forest (default: 10001)",
              metavar="integer"),
  
  make_option(c("--num_perm"), type="integer", default=1000,
              help="Number of permutations to base null OOB error distribution on (default: 1000)",
              metavar="integer")
)

opt_parser <- OptionParser(
  option_list=option_list, 
  usage = "%prog [options] -i input --seed 213 --num_tree 10001 --num_perm 1000 -o outfile.rds",
  description = paste(
    "\nThis script runs random forest on an input table assuming",
    "that the last column is what you are classifying.",
    "Leave-one-out cross-validation will also be run to give an estimate of accuracy.",
    "Finally, model significance is also outputted based on the out-of-bag error test",
    "available in the rfUtilities package"))

opt <- parse_args(opt_parser)

# Check that required arguments are set:
if(is.null(opt$input)) {
  stop("path to input file needs to be set.")
}

if(is.null(opt$seed)) {
  stop("random seed needs to be set.")
}

if(is.null(opt$output)) {
  stop("path to output file needs to be set.")
}

# Load required packages
library("randomForest") # for running RF and getting out of bag (OOB) error
library("rfUtilities") # for running OOB error significance test (requires re-running RF with permuted classes)
library("caret") # for running leave-one-out cross-validation to estimate accuracy (wraps randomForest and a number of other packages)

if(opt$threads > 1) {
  library("doMC") # for multithreading when possible
  registerDoMC(cores = opt$threads)
}

run_all_RF <- function(filename, num_tree, num_perm, set_seed, outfile) {
  
  # Read in table. Column to classify is last column
  intable <- read.table( filename , header=T , sep="\t" , row.names=1 )
  
  # Set random seed
  set.seed(set_seed)
  
  # Below is the default mtry parameter choice for randomForest for classification
  mtry_count <- floor(sqrt(ncol(intable[,1:(ncol(intable)-1)])))
  
  # Run random forest, calculate importance and proximities as well
  rf_mod <- randomForest( x=intable[,1:(ncol(intable)-1)] , y=intable[ , ncol(intable)] , ntree = num_tree, 
                          proximity = TRUE, importance = TRUE, mtry = mtry_count)
  
  # Run leave-one-out cross-validation to estimate model accuracy
  fit_control <- trainControl(method = "LOOCV")
  rf_mod_loocv <- train(x = intable[,1:(ncol(intable)-1)], y = intable[ , ncol(intable)], method = "rf", ntree = num_tree, 
                  tuneGrid = data.frame(mtry = mtry_count), trControl = fit_control)
  
  # Test whether observed OOB error is significant
  rf_sig <- rf.significance(x=rf_mod,  xdata=intable[,1:(ncol(intable)-1)], nperm=num_perm, ntree=num_tree) 

  # Make dataframe with summary data for RF
  results <- data.frame( matrix( c( filename, ncol(rf_mod_loocv$trainingData) - 1,
                                      rf_mod$mtry, rf_mod$ntree, rf_mod_loocv$results$Accuracy, 
                                      rf_sig$test.OOB, median(rf_sig$RandOOB), rf_sig$pValue), nrow=1, ncol=8))
  
  colnames(results) <- c( "file" , "num.features" , "mtry" , "num.trees" , "LOOCV.Accuracy" , "Model.OOB.error" , "Median.Rand.OOB.error" , "P.value" )
  
  
  rf_out <- list("mod" = rf_mod , "sig" = rf_sig , "loocv_mod" = rf_mod_loocv , "summary" = results)
  
  saveRDS(rf_out, file = outfile)
  
}

# Run RF on input file.
run_all_RF(filename = opt$input, 
           num_tree = opt$num_tree, 
           num_perm = opt$num_perm, 
           set_seed = opt$seed, 
           outfile = opt$output)
