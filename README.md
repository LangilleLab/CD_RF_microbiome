# IBD_microbiome
Scripts and workflows used for the prcessing and analysis of human genetic, metagenomic and 16S rRNA sequences from the biopsies of pediatric Crohn's patients and healthy colon controls.

LOOCV_RF.R is an Rscript (designed to be run in parallel on the command-line) written by Gavin Douglas to determine whether random forest accuracies estimated by leave-one-out cross-validation are statistically significant based on a permutation test. This test is much slower than testing for model significance based on out-of-bag error, which is why we instead opted to use the test implemented in the "rfUtilities" R package.
