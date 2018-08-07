# Ensuring that nothing else is loaded
rm(list=ls())
gc()

# Loading basic data
source("./setup_scripts/base_paths.R")
source("./setup_scripts/prognosis.R")                                         

clear_temp_data <- function(bdata) {
                    env_dat   <- ls(envir=globalenv())
                    clear_dat <- env_dat[!(env_dat %in% bdata)]
                    rm(list=clear_dat, envir=globalenv())
                    gc()
                   }

basic_data <- c(ls(), "basic_data")
# Loading the main scripts
print("Processing microarray .soft files")
source("./processing_scripts/01_process_microarray_sample_description.R")
clear_temp_data(basic_data)

print("Processing microarray expression data")
source("./processing_scripts/02_process_microarray_expression_data.R")
clear_temp_data(basic_data)

print("Processing rnaseq data")
source("./processing_scripts/03_process_rnaseq_data.R")
clear_temp_data(basic_data)

print("All done!")
