analysis_name       <- "cytopred"

nthreads            <- 4 
label_type          <- "crisk_code"                                                   
ref_labels          <- c("good")                                                      
test_labels         <- c("intr", "intr_poor", "poor") # This is the phenotype which we are characterizing  

ref_group           <- "good"
test_group          <- "intr_poor"

analysis_path       <- file.path(analysis_root, analysis_name)
working_data_path   <- file.path(analysis_path, "working_data")
working_sample_path <- file.path(working_data_path, "sample_description")
working_expr_path   <- file.path(working_data_path, "expression_data")

results_path        <- file.path(analysis_path, "results")
working_desc_file   <- file.path(working_data_path, "sample_description.csv")
working_expr_file   <- file.path(working_data_path, "expression_data.csv")
analysis_info_file  <- file.path(working_data_path, "analysis_info.csv")
analysis_data_file  <- file.path(working_data_path, "analysis_data")
gpi_matrix_file     <- file.path(working_data_path, "gpi_matrix")
random_sample_path  <- file.path(working_data_path, "random_sample")
sample_table_file   <- file.path(working_data_path, "sample_table")
train_results_file  <- file.path(results_path, "train_results.csv")
test_results_file   <- file.path(results_path, "test_results.csv")
val_results_file    <- file.path(results_path, "val_results.csv")
aggregated_results_file <- file.path(results_path, "aggregated_results.csv")


analysis_dsets      <- c("GSE10358_GPL570",
                         "GSE12417_GPL96",
                         "GSE12417_GPL570",
                         "GSE13159_GPL570",
                         "GSE14468_GPL570",
                         "GSE15434_GPL570",
                         "GSE16015_GPL570",
                         "GSE22845_GPL570",
                         "GSE30285_GPL5175",
                         "GSE39730_GPL570",
                         "GSE61804_GPL570",
                         "LAML_TCGA_rnaseq",
                         "SRP050272_rnaseq")

# These are the datasets that will be chosen for training + validation split
select_dsets <- c("GSE10358_GPL570",
                  "GSE12417_GPL96",
                  "GSE12417_GPL570",
                  "GSE13159_GPL570",
                  "GSE14468_GPL570",
                  "GSE15434_GPL570",
                  "GSE22845_GPL570")

random_guesser_trials <- 100000
random_guesser_cutoff <- 0.9999 # Gene pairs tested must be better than 99.99% of random gene pairs

bucket_size <- 3000
nbuckets    <- 10000

# Common functions

gp_accuracy_check    <- function(true_labels, test_labels) {
                        	ptable <- rbind(c(sum(true_labels == F & test_labels == F),
                        	                  sum(true_labels == T & test_labels == F)),
                        	                c(sum(true_labels == F & test_labels == T),
                        	                  sum(true_labels == T & test_labels == T)))

                            main_diag <- (sum(diag(ptable))/sum(ptable)) >= 0.5
                            if (main_diag) {
                                retval <- diag(ptable)
                            } else {
                                retval <- c(ptable[2, 1], ptable[1, 2])
                            }
                            names(retval) <- c("ref_correct", "test_correct")
                            retval
                        }


random_class_guesser <- function(ll) {
                            pll <- sample(ll) # Permulted logic label
                            gp_accuracy_check(ll, sample(ll))
                        }

generate_gpi_matrix <- function(expr, gene_pairs) {
                        gmat <- do.call("cbind", lapply(gene_pairs, function(gp) {
                                    genes <- strsplit(gp, "_")[[1]]
                                    x <- expr[genes[1],] > expr[genes[2],]
                                    x[x==T] <- paste(genes[1], "_greater", sep="")
                                    x[x==F] <- paste(genes[2], "_greater", sep="")
                                    x
                                }))
                        colnames(gmat) <- gene_pairs
                        gmat
                       }

get_perf_stats  <- function(pred, l) {                                           
                    ptable <- table(pred, l)                                     
                    sens   <- diag(ptable)/colSums(ptable)                       
                    names(sens) <- paste(names(sens), "_sensitivity", sep="")    
                    pv     <- diag(ptable)/rowSums(ptable)                       
                    names(pv) <- paste(names(pv), "_predictive_value", sep="")    
                    round(c(sens, pv) * 100, 2)                                  
                   }                                                             
                                                                                 
get_prediction <- function(tr, gmat=NULL, cval, training_data=F) {                 
                    ptree <- prune(tr, cp=cval)                                  
                    if (training_data) {                                         
                        pred <- predict(ptree, type="class")                     
                    } else {                                                     
                        pred <- predict(ptree, type="class", newdata=gmat)         
                    }                                                            
                    list(ptree=ptree, pred=pred)                                 
                  }                                                              
                                                                                 
generate_results <- function(tdat, gmat=NULL, gene_pairs=NULL, l, training_data=F) {
                        if ((is.null(gmat)) & (training_data==F )) {             
                            print("Check input. Either gmat or gene pairs is missing, or training data should be set to TRUE")
                        }                                                        
                        cp_vals <- round(tdat$cp_vals, 8)                        
                        result  <- do.call("rbind", lapply(cp_vals, function(cp_val) {
                                    pred_dat <- get_prediction(tr=tdat$refit_tree, gmat=gmat, cval=cp_val, training_data=training_data)
                                    pred     <- pred_dat$pred                    
                                    ptree    <- pred_dat$ptree                   
                                    npairs <- ptree$frame$var                    
                                    npairs <- length(npairs[npairs!="<leaf>"])   
                                    c(bucket=tdat$bucket, cp_val=cp_val,         
                                      get_perf_stats(pred=pred, l=l), npairs=npairs)
                                   }))                                           
                    }                                                            
                                                                                 
get_tree_dat  <- function(sgpi_matrix, tlabels, bucket_name) {                                
                    train_mat  <- data.frame(sgpi_matrix, tree_labels=tlabels)
                    tree       <- rpart(formula=tree_labels ~ ., data=train_mat, method="class")
                    tree_frame <- tree$frame$var                                 
                    fit_vars   <- as.character(tree_frame[tree_frame != "<leaf>"])
                    fit_df     <- data.frame(sgpi_matrix[,fit_vars],       
                                             tree_labels=tree_labels)            
                    refit_tree <- rpart(tree_labels ~ ., data=fit_df, method="class")
                    tree_dat   <- list(refit_tree=refit_tree, fit_vars=fit_vars, 
                                       bucket=sub(pattern="bucket", x=bucket_name, replacement=""),
                                       cp_vals=refit_tree$cptable[-1,"CP"])      
                 }     

get_result_by_split <- function(rdf, sfactor_name) {                             
                        sfactor <- rdf[,sfactor_name]                            
                        sdf <- split(rdf, sfactor)                               
                        res <- lapply(sdf, function(df) {                        
                                ptable <- table(df$pred, df$sample_group)        
                                group_order <- colnames(ptable)                  
                                correct <- diag(ptable)                          
                                predicted <- table(df$pred)                      
                                tot     <- table(df$sample_group)                
                                sens    <- round(100 * diag(ptable)/colSums(ptable), 2)
                                pv      <- round(100 * diag(ptable)/rowSums(ptable), 2)
                                ret_val <- unlist(lapply(1:2, function(j) {      
                                            c(predicted[j], correct[j], tot[j], sens[j], pv[j])
                                           }))                                   
                                name_addon <- c("predicted", "correct", "total", "sensitivity", "predictor_val")
                                name_grid  <- expand.grid(name_addon, group_order)
                                retnames   <- apply(name_grid, 1, function(x) {  
                                                paste(x, collapse="_")           
                                              })                                 
                                names(ret_val) <- retnames                       
                                ret_val                                          
                        })                                                       
                        resnames <- names(res)                                   
                        res <- do.call("rbind", res)                             
                        res[is.na(res)] <- "not_applicable"                      
                        res <- data.frame(resnames, res)                         
                        colnames(res)[1] <- sfactor_name                         
                        res                                                      
                       }                                                         

ready_surv_data <- function(sdsets, rdf) {
                    sdesc_files <- list.files(working_sample_path, pattern=paste(sdsets, collapse="|"))
                    surv_dats   <- lapply(sdesc_files, function(sdesc_file) {
                                    sdesc_path <- file.path(working_sample_path,
                                                            sdesc_file)
                                    sdat <- data.frame(fread(sdesc_path),
                                                       row.names=1)

                                    if (grepl("GSE10358_GPL570", sdesc_file)) {
                                        selected_cols <- c("sample_id",
                                                           "os.months..3.31.10.",
                                                           "efs.months..3.31.10.",
                                                           "vital.status")
                                        sdat <- data.frame(sdat[,selected_cols],
                                                           time="months", efs_status="no_data")
                                        sdat[,"vital.status"] <- tolower(sdat[,"vital.status"])
                                        colnames(sdat) <- c("sample_id",
                                                            "os_time",
                                                            "efs_time",
                                                            "os_status",
                                                            "time_units",
                                                            "efs_status")
                                        unkn_i <- apply(sdat, 1, function(x) {
                                                    "extract_not_available" %in% x
                                                  })
                                        sdat <- sdat[!unkn_i,]

                                    } else if (grepl("GSE14468_GPL570", sdesc_file)) {
                                        selected_cols <- c("sample_id",
                                                           "os", "efs",
                                                           "efsi", "osi")
                                        sdat <- data.frame(sdat[,selected_cols],
                                                           time="months")
                                        colnames(sdat) <- c("sample_id",
                                                            "os_time",
                                                            "efs_time",
                                                            "efs_status",
                                                            "os_status",
                                                            "time_units")
                                    } else if (grepl("LAML_TCGA_rnaseq", sdesc_file)) {
                                        selected_cols <- c("sample_id",
                                                           "os_months",
                                                           "dfs_months",
                                                           "dfs_status",
                                                           "os_status")
                                        sdat <- data.frame(sdat[,selected_cols],
                                                           time="months")
                                        sdat$os_status[sdat$os_status=="DECEASED"] <- "dead"
                                        sdat$os_status[sdat$os_status=="LIVING"]   <- "alive"
                                        colnames(sdat) <- c("sample_id",
                                                            "os_time",
                                                            "efs_time",
                                                            "efs_status",
                                                            "os_status",
                                                            "time_units")
                                    } else {
                                        print("Surv data processing error")
                                    }
                                    sdat_i <- sdat[,"os_status"] %in% c("alive",
                                                                        "dead")
                                    sdat <- sdat[sdat_i,]
                                    sstatus <- rep(1, nrow(sdat))
                                    sstatus[sdat[,"os_status"]=="alive"] <- 0
                                    sdat <- data.frame(sdat,
                                                       simple_os_status=sstatus)
                                    i <- match(sdat[,"sample_id"],
                                               rdf$sample_id)
                                    retdat <- data.frame(sdat, rdf[i,])
                                    retdat[complete.cases(retdat),]
                   })}

