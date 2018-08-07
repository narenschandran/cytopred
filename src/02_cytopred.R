# Ensuring that nothing else is loaded                                           
rm(list=ls())                                                                    
gc()                                                                             
                                                                                 
# Loading basic data 

source("./setup_scripts/base_paths.R")
source("./gene_sig_scripts/00_cytopred_init.R")


clear_temp_data <- function(bdata) {
                    env_dat   <- ls(envir=globalenv())
                    clear_dat <- env_dat[!(env_dat %in% bdata)]
                    rm(list=clear_dat, envir=globalenv())
                    gc()
                   }

basic_data <- c(ls(), "basic_data")
# Loading the main scripts
source("./processing_scripts/04_working_data_processing.R")
clear_temp_data(basic_data)

source("./gene_sig_scripts/01_generate_trees.R")
#clear_temp_data(basic_data)

source("./gene_sig_scripts/02_select_trees.R")
#clear_temp_data(basic_data)



