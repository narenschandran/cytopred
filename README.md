# cytopred
The current version of the top scoring pair decision tree based signature to classify AML samples into "good" vs "intermediate/poor" cytogenetic risk groups using expression data. Running the master script should download most of the data required for the analysis. LAML_TCGA data and SRP050272 data should already be present in the data/input folder before the analysis is started.

The files required for LAML_TCGA are:
1. data_clinical.txt
2. LAML.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt

Both are from Broad Firehose.

The files required for SRP050272 are:
tx2gene.csv
salmon_output (folder)

The script to generate tx2gene.csv is in the rnaseq_scripts folder.

The current solution for generating the script for processing SRP050272 is by using generate_script.R in the rnaseq_scripts folder
