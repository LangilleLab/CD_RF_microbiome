# CD\_RF\_microbiome
Scripts used for the processing and analysis of human genetic,
metagenomic and 16S rRNA sequences from the biopsies of pediatric Crohn's
patients and healthy colon controls.

### Human Genotype Parsing
Most of the human genotype analyses were performed with publicly available
software. However, two scripts needed to be written to perform the imputation
pipeline.
* *splitVCFbyChr.pl* - script to split a VCF into different chromosomes.
* *impute2vcf.pl* - script to parse the IMPUTE2 output and to re-build a VCF
that includes the imputed variants.

### Random Forest Pre-processing

The three scripts in the *pre-processing* folder were used to prepare the input
tables for random forest.

* *convert_16S_ids_and_truncate.R* - commands required to convert ids to
consistent names used throughout the study. Also, cleaned up raw 16S tables.
* *clean_MGS.R* - commands to clean up the metagenomics raw files and in
particular to keep only function ids, not full descriptions.
* *prep_RF_input.R* - commands used to read in the cleaned tables generated
by the above two files, to collapse each table to different hierarchical levels,
to remove features that were present in less or equal to 10% of samples, and to
output a prepped input tables for both the disease state and treatment response
random forest models.


### Running Random Forest Models

* *RF_running.R* - the main script used to run the random forest models in a
loop. Each model was outputted as an R object.
* *get_summary_tables.R* - commands used to loop through the model R objects
to get summary information for each one.
* *old/LOOCV_RF.R* contains the original script that was written before the
rfUtilities package was determined to be a better method of determining model
significance.
