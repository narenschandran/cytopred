### Setting up commonly used paths ###
data_path        <- "../data"
input_path       <- file.path(data_path, "input")

proc_path        <- file.path(data_path, "processed")
proc_sample_path <- file.path(proc_path, "sample_description")                   
proc_expr_path   <- file.path(proc_path, "expression_data")

analysis_root    <- "../analysis"

array_info_file          <- file.path("./misc_data/array_pkgname_info.csv")
GSE14468_surv_file       <- "./misc_data/OS_EFS_GSE6891.xls"
GSE14468_fixed_karyotype <- "./misc_data/GSE14468_fixed_karyotypes.xls"

rnaseq_datasets          <- c("SRP050272",
                              "LAML_TCGA")

GEO_microarray_datasets  <- c("GSE10358",
                              "GSE12417",
                              "GSE13159",
                              "GSE14468",
                              "GSE15434",
                              "GSE16015",
                              "GSE22845",
                              "GSE30285",
                              "GSE39730",
                              "GSE61804")


# I'm not using GPL97 because there is no overlap between the genes in that
# dataset and every other one, which means that pretty much most of the
# genes would not be available for common analysis
allowed_p_ids <- c("GPL96", 
                   "GPL570", 
                   "GPL5175", 
                   "rnaseq")

