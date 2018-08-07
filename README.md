# cytopred
The current version of the top scoring pair decision tree based signature to classify AML samples into "good" vs "intermediate/poor" cytogenetic risk groups using expression data. Running the master script should download most of the data required for the analysis. LAML_TCGA data and SRP050272 data should already be present in the data/input folder before the analysis is started.

The files required for LAML_TCGA are:
1. data_clinical.txt
2. LAML.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt

Both are from Broad Firehose.

The files required for SRP050272 are:
tx2gene.csv
salmon_output (folder)

