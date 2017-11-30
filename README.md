# CD\_RF\_microbiome
Scripts used for the processing and analysis of human genetic,
metagenomic and 16S rRNA sequences from the biopsies of pediatric Crohn's
patients and healthy colon controls.

### Pre-processing

```BISCUIT_utility_code.R``` contains useful functions and is sourced at the start
of several of scripts in this repository.

```raw_data/``` contains the raw abundance tables used as input for the random
forest analyses. ```bicuit_metadata.txt``` contains the relevant metadata for each
sample.

These tables were prepped for the RF analyses with these scripts:

* ```prep_BISCUIT_16S_RF_input.R```

* ```prep_BISCUIT_MGS_RF_input.R```

* ```prep_combined_RF.R```

### Running RF

The RF models were run using this command-line R-script:
```run_RF.R``` - you can type this script with the ```-h``` option to see the
input arguments.

The above script outputs RDS files for each RF model, which were parsed using
commands in these two scripts:

```get_summary_tables.R``` - Get summary metrics for each model.

```get_varImp_tables.R``` - Get features sorted by variable importance for each
model.

Several scripts contain the plotting commands used when analyzing this data and
are found in *plot_scripts/*.

### Human Genotype Parsing
Most of the human genotype analyses were performed with publicly available
software. However, two scripts needed to be written to perform the imputation
pipeline.
* ```splitVCFbyChr.pl``` - Split a VCF into different chromosomes.
* ```impute2vcf.pl``` - Parse the IMPUTE2 output and to re-build a VCF
that includes the imputed variants.
