library(data.table)
library(rpart)
library(rpart.plot)
library(survival)
library(survminer)
load(analysis_data_file)
load(gpi_matrix_file)
num_conv <-  function(fval) { as.numeric(as.character(fval)) }
analysis_info <- read.csv(analysis_info_file, row.names=1, stringsAsFactors=F)
dset_df <- table(analysis_info$dset, analysis_info$class_labels)
write.csv(dset_df, file.path(results_path, "dset_description.csv"))


selected_buckets <- c(922, 3302)
selected_cp_vals <- c(0.01, 0.05325444)
do_surv_analysis <- TRUE
survival_dsets   <- c("GSE10358_GPL570", "GSE14468_GPL570", "LAML_TCGA_rnaseq")

tree_labels      <- factor(as.character(analysis_data$train$sample_group))
val_tree_labels  <- factor(as.character(analysis_data$validate$sample_group))
test_tree_labels <- factor(as.character(analysis_data$test$sample_group))

results <- lapply(1:length(selected_buckets), function(b_i) {            
           	i <- selected_buckets[b_i]
            bucket_name <- paste("bucket", i, sep="")
            cp_val      <- selected_cp_vals[b_i]
           	result_name <- paste(bucket_name, "__cp_val", cp_val, sep="")
           	outpath     <- file.path(results_path, result_name)
           	dir.create(outpath)
            load(file=file.path(random_sample_path, bucket_name)) # Provides the prestored random_sample
            tree_dat    <- get_tree_dat(gpi_matrix[,random_sample],
                                        tlabels=tree_labels,
                                        bucket_name=bucket_name)
            train_pred <- get_prediction(tr=tree_dat$refit_tree,             
                                         gmat=NULL,                        
                                         cval=cp_val,                    
                                         training_data=T) 
            val_pred    <- get_prediction(tr=tree_dat$refit_tree,
                                          gmat=data.frame(generate_gpi_matrix(analysis_data$validate$expr,
                                                                              tree_dat$fit_vars)),
                                                                              cval=cp_val,                    
                                                                              training_data=F)                
            test_pred   <- get_prediction(tr=tree_dat$refit_tree,
                                          gmat=data.frame(generate_gpi_matrix(analysis_data$test$expr,
                                                                              tree_dat$fit_vars)),
                                                                              cval=cp_val,                    
                                                                              training_data=F)                

            train_df=data.frame(sample_id=rownames(gpi_matrix), pred=train_pred$pred)
            val_df=data.frame(sample_id=names(val_pred$pred), pred=val_pred$pred)
            test_df=data.frame(sample_id=names(test_pred$pred), pred=test_pred$pred)
            all_predictions <- data.frame(do.call("rbind", list(train_df, val_df, test_df)))

            pred_match   <- match(as.character(all_predictions$sample_id), analysis_info$sample_id)
            subset_info  <- analysis_info[pred_match,]
            result_df    <- data.frame(all_predictions, 
                                       sample_group=subset_info$sample_group,
                                       dset=subset_info$dset,
                                       label_type=label_type,
                                       class_labels=subset_info$class_labels,
                                       analysis_set=subset_info$analysis_set)

            if (do_surv_analysis) {
            	surv_dats <- ready_surv_data(sdsets=survival_dsets, rdf=result_df)
            	names(surv_dats) <- survival_dsets
            	surv_main_outpath <- file.path(outpath, "survival_plots")
            	dir.create(surv_main_outpath)
            	splots <- lapply(names(surv_dats), function(sname) {
            	          	surv_dat   <- surv_dats[[sname]]
            	          	surv_outpath  <- file.path(surv_main_outpath, sname)
            	          	time_units <- unique(surv_dat[,"time_units"])
            	          	dir.create(surv_outpath)
            	          	png(filename=file.path(surv_outpath, paste(sname, "_prediction_survival_plot.png", sep="")),
            	          	    height=720, width=720)
            	          	pred_plot <- ggsurvplot(survfit(Surv(num_conv(os_time), 
            	          	                                     num_conv(simple_os_status))
            	          	                                ~ pred, data=surv_dat),
            	          	                        pval=T, title=paste(sname, "OS (gene signature prediction)"),
            	          	                        xlab=paste("Time in", time_units),
            	          	                        risk.table=T, data=surv_dat)
            	          	print(pred_plot)
            	          	dev.off()
            	          	png(filename=file.path(surv_outpath, paste(sname, "_actual_survival_plot.png", sep="")),
            	          	    height=720, width=720)
            	          	lab_plot  <- ggsurvplot(survfit(Surv(num_conv(os_time), 
            	          	                                     num_conv(simple_os_status))
            	          	                                ~ sample_group, data=surv_dat),
                                                    pval=T, title=paste(sname, "OS (actual sample grouping)"),
                                                    xlab=paste("Time in", time_units),
                                                    risk.table=T, data=surv_dat)
            	          	print(lab_plot)
            	          	dev.off()
            	          	png(filename=file.path(surv_outpath, paste(sname, "_both_survival_plot.png", sep="")),
            	          	    height=720, width=1280)
            	          	both_plot <- ggsurvplot(survfit(Surv(num_conv(os_time), 
            	          	                                     num_conv(simple_os_status))
            	          	                                ~ pred + sample_group, data=surv_dat),
                                                    pval=T, title=paste(sname, "OS (both)"),
                                                    xlab=paste("Time in", time_units),
                                                    risk.table=T, data=surv_dat)

            	          	print(both_plot)
            	          	dev.off()
            	          	list(pred_plot=pred_plot, label_plot=lab_plot,  
            	          	     both_plot=both_plot)
            	          })
            }
            comb_result <- get_result_by_split(result_df, "label_type")
            dset_result <- get_result_by_split(result_df, "dset")
            set_result <- get_result_by_split(result_df, "analysis_set")
            write.csv(result_df, file.path(outpath, paste("results_", result_name, ".csv", sep="")))
            write.csv(comb_result, file.path(outpath, paste("combined_result_", result_name, ".csv", sep="")))
            write.csv(dset_result, file.path(outpath, paste("dataset_results_", result_name, ".csv", sep="")))
            write.csv(set_result, file.path(outpath, paste("analysis_set_results_", result_name, ".csv", sep="")))
            png(filename=file.path(outpath, paste("training_", result_name, ".png", sep="")),
                width=1280,
                height=720)
            rpart.plot(train_pred$ptree, type=2)
            dev.off()
            png(filename=file.path(outpath, paste("validation_", result_name, ".png", sep="")),
                width=1280,
                height=720)
            rpart.plot(val_pred$ptree, type=2)
            dev.off()
            png(filename=file.path(outpath, paste("test_", result_name, ".png", sep="")),
                width=1280,
                height=720)
            rpart.plot(test_pred$ptree, type=2)
            dev.off()
            # c(rbind(x1, x2)) seems to interleave x1 and x2
            bplot_labels <- paste(c("correct", "total"), c(ref_group, ref_group, test_group, test_group), sep="_")
            bplot_df <- rbind(c(rbind(num_conv(dset_result[,bplot_labels[1]]),
                                      num_conv(dset_result[,bplot_labels[2]]))),
                              c(rbind(num_conv(dset_result[,bplot_labels[3]]), 
                                      num_conv(dset_result[,bplot_labels[4]]))))
            colnames(bplot_df) <- paste(c(rbind(as.character(dset_result$dset), as.character(dset_result$dset))), c("correct", "total"), sep="_")
            rownames(bplot_df) <- c(ref_group, test_group)
            png(filename=file.path(outpath, paste("dset_barplot_", result_name, ".png", sep="")),
                width=1280, height=980)
            par(mar=c(15, 5, 5, 5))
            barplot(bplot_df, ylim=c(0, max(colSums(bplot_df))),
                    las=2, ylab="Sample number", legend=rownames(bplot_df))
            dev.off()
            set_bplot_df <- rbind(c(rbind(num_conv(set_result[,bplot_labels[1]]),
                                          num_conv(set_result[,bplot_labels[2]]))),
                                  c(rbind(num_conv(set_result[,bplot_labels[3]]), 
                                          num_conv(set_result[,bplot_labels[4]]))))
            colnames(set_bplot_df) <- paste(c(rbind(as.character(set_result$analysis_set), as.character(set_result$analysis_set))), c("correct", "total"), sep="_")
            rownames(set_bplot_df) <- c(ref_group, test_group)
            png(filename=file.path(outpath, paste("analysis_set_barplot_", result_name, ".png", sep="")),
                width=800, height=600)
            par(mar=c(15, 5, 5, 5))
            barplot(set_bplot_df, ylim=c(0, max(colSums(set_bplot_df))),
                    las=2, ylab="Sample number", legend=rownames(set_bplot_df))
            dev.off()
            mosaic_df <- table(pred=result_df$pred, actual_group=result_df$sample_group)
            png(filename=file.path(outpath, "result_mosaic.png"),
                height=800, width=800)
            plot(mosaic_df, col=c("indianred3", "skyblue4"),
                 main="Mosaic plot: total results", 
                 xlab="Prediction", ylab="Actual labels", cex=1)
            dev.off()
            list(result_df=result_df, comb_result=comb_result, 
                 dset_result=dset_result, set_result=set_result,
                train_tree=train_pred$ptree, val_tree=val_pred$ptree, 
                test_tree=test_pred$ptree) 
           })



