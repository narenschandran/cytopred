library(parallel)                                                                            
library(rpart)
library(data.table)                                                              

working_desc <- data.frame(fread(working_desc_file, colClasses=c("character")), row.names=1)
working_group <- rep("skip", nrow(working_desc))
working_group[working_desc[,label_type]  %in% ref_labels] <- ref_group
working_group[working_desc[,label_type]  %in% test_labels] <- test_group
working_desc <- data.frame(working_desc, sample_group=working_group, stringsAsFactors=F)


working_expr <- data.frame(fread(working_expr_file, stringsAsFactors=F), row.names=1)
snames       <- colnames(working_expr)                                           
working_expr <- t(apply(working_expr, 1, function(x) {                           
                    as.numeric(as.character(x))                                  
                }))                                                              
colnames(working_expr) <- snames

# Sanity check. Should return true
# identical(colnames(working_expr), working_desc$sample_id)

test_select <- which((!(working_desc$dset %in% select_dsets)) & 
                     (working_desc[,"sample_group"] %in% c(ref_group, test_group)))


train_validation_select <- which((working_desc$dset %in% select_dsets) & 
                                (working_desc[,"sample_group"] %in% c(ref_group, 
                                                                    test_group)))

print("Splitting data into training, test and validation datasets")
set.seed(1)
train_i <- sample(1:length(train_validation_select), size=2*length(train_validation_select)/3)
train_select <- train_validation_select[train_i]
validation_select  <- train_validation_select[-train_i]

analysis_set <- rep("skip", nrow(working_desc))
analysis_set[validation_select] <- "validate"
analysis_set[train_select]      <- "train"
analysis_set[test_select]       <- "test"

analysis_info <- data.frame(sample_id=working_desc$sample_id, class_labels=working_desc[,label_type], analysis_set=analysis_set, dset=working_desc$dset, sample_group=working_desc$sample_group, stringsAsFactors=F)

write.csv(analysis_info, analysis_info_file)

analysis_data               <- list()

analysis_data$train         <- list()
analysis_data$validate      <- list()
analysis_data$test          <- list()

analysis_data$train$expr            <- working_expr[,train_select]
analysis_data$validate$expr         <- working_expr[,validation_select]
analysis_data$test$expr             <- working_expr[,test_select]

analysis_data$train$desc            <- working_desc[train_select,]
analysis_data$validate$desc         <- working_desc[validation_select,]
analysis_data$test$desc             <- working_desc[test_select,]

analysis_data$train$class_labels    <- analysis_data$train$desc[,label_type]
analysis_data$validate$class_labels <- analysis_data$validate$desc[,label_type]
analysis_data$test$class_labels     <- analysis_data$test$desc[,label_type]

analysis_data$train$sample_group    <- analysis_data$train$desc[,"sample_group"]
analysis_data$validate$sample_group <- analysis_data$validate$desc[,"sample_group"]
analysis_data$test$sample_group     <- analysis_data$test$desc[,"sample_group"]

analysis_data$train$logic_labels    <- analysis_data$train$class_labels %in% test_labels


print("Generating data for accuracy of random guesses")
set.seed(1)
random_guesses <- replicate(random_guesser_trials, random_class_guesser(analysis_data$train$logic_labels))
rownames(random_guesses) <- c("ref_correct", "test_correct")
analysis_data$train$random_guesses <- ceiling(apply(random_guesses, 1, function(x) {
                             	quantile(x, prob=random_guesser_cutoff)[[1]]
                            }))

print("Saving data to file")
save(analysis_data, file=analysis_data_file, compress=F)

expr     <- analysis_data$train$expr
texpr    <- t(expr)
rnames   <- rownames(expr)
cnames   <- colnames(expr)
rguesses <- analysis_data$train$random_guesses
tlabels  <- analysis_data$train$logic_labels

rm(working_expr)
rm(working_desc)
rm(analysis_data)
gc()


print("Generating gene pair indicator matrix")
cl <- makeCluster(nthreads, type="FORK")
lgpi_matrix <- do.call("cbind", parLapply(cl, 1:nrow(expr), function(i) {
               	g_index <- rnames[i] > rnames
               	if (sum(g_index) > 0) {
               		gmat <- matrix(expr[i,] > texpr[,g_index], 
               		               nrow=length(expr[i,]), 
               		               dimnames=list(cnames, rnames[g_index]))
               		colnames(gmat) <- paste(rnames[i], colnames(gmat), sep="_")
               		gp_select <- apply(gmat, 2, function(x) {
               		             	gp_accuracy_check(tlabels, x)
               		             })
              		select_index <- colSums(gp_select > rguesses)==2
              		if (sum(select_index) > 0) {
               			gmat[,select_index]
               		} else {
               			NULL
               		}
               	} else {
               		NULL
               	}
               }))
stopCluster(cl)

rm(expr)
rm(texpr)
gc()

print("Making gene pair indicator ready for making descriptive trees")
gpi_matrix <- do.call("cbind", lapply(1:ncol(lgpi_matrix), function(i) {
              	genes <- strsplit(colnames(lgpi_matrix)[i], "_")[[1]]
              	x <- lgpi_matrix[,i]
              	x[x==T] <- paste(genes[1], "_greater", sep="")
              	x[x==F] <- paste(genes[2], "_greater", sep="")
              	x
              }))
colnames(gpi_matrix) <- colnames(lgpi_matrix)

rm(lgpi_matrix)
gc()

save(gpi_matrix, file=gpi_matrix_file, compress=F)

print(paste("Number of gene pairs left after filtering:", ncol(gpi_matrix)))
bucket_index <- rep(F, ncol(gpi_matrix))
bucket_index[1:bucket_size] <- T

dir.create(random_sample_path)
set.seed(1)
sapply(1:nbuckets, function(i) {
	bucket_name <- paste("bucket", i, sep="")
	random_sample <- sample(bucket_index)
	save(random_sample, file=file.path(random_sample_path, bucket_name))
})


load(analysis_data_file)

# Regardless of the previous assignment to string, the data gets saved as factor
# This includes a factor column into it, which will pop up in tables
# To prevent that, I'm converting it into character, and then factor again.
# This removes "skip" from the factor levels. Must be mindful of the factor level
# manipulations order though for future
tree_labels      <- factor(as.character(analysis_data$train$sample_group))
val_tree_labels  <- factor(as.character(analysis_data$validate$sample_group))
test_tree_labels <- factor(as.character(analysis_data$test$sample_group))


print("Generating and testing trees")
cl <- makeCluster(4, type="FORK")
results <- parLapply(cl, 1:nbuckets, function(i) {
           	if(i %% 4 == 0) {
           		gc()
           	}
           	bucket_name <- paste("bucket", i, sep="")
           	load(file=file.path(random_sample_path, bucket_name)) # Provides the prestored random_sample
           	tree_dat    <- get_tree_dat(gpi_matrix[,random_sample], 
           	                            tlabels=tree_labels,
           	                            bucket_name=bucket_name)
           	train_results <- generate_results(tdat=tree_dat, 
           	                                  gmat=NULL, 
           	                                  l=tree_labels, 
           	                                  training_data=T)
           	val_results   <- generate_results(tdat=tree_dat,
           	                                  gmat=data.frame(generate_gpi_matrix(analysis_data$validate$expr,
           	                                                                      tree_dat$fit_vars)), 
           	                                                                      l=val_tree_labels)
           test_results  <- generate_results(tdat=tree_dat,
                                             gmat=data.frame(generate_gpi_matrix(analysis_data$test$expr, 
                                                                                 tree_dat$fit_vars)), 
                                                                                 l=test_tree_labels)
           list(train_results=train_results,
                val_results=val_results,
                test_results=test_results)
})
stopCluster(cl)

result_to_df <- function(res) {                                                  
                    cnames <- colnames(res)                                      
                    rownames(res) <- NULL                                        
                    res <- data.frame(res)                                       
                    res[,1] <- as.numeric(as.character(res[,1]))
                    res[,2] <- as.numeric(as.character(res[,2]))                 
                    res[,3] <- as.numeric(as.character(res[,3]))                 
                    res[,4] <- as.numeric(as.character(res[,4]))                 
                    res[,5] <- as.numeric(as.character(res[,5]))                 
                    res[,6] <- as.numeric(as.character(res[,6]))                 
                    res[,7] <- as.numeric(as.character(res[,7]))                 
                    res                                                          
                }                                                                

train_results <- result_to_df(do.call("rbind", 
                              lapply(results, function(r) { r$train_results })))
test_results  <- result_to_df(do.call("rbind", 
                              lapply(results, function(r) { r$test_results })))
val_results   <- result_to_df(do.call("rbind", 
                              lapply(results, function(r) { r$val_results })))

write.csv(train_results, train_results_file)
write.csv(test_results, test_results_file)
write.csv(val_results, val_results_file)

all_results <- list(train_results, test_results, val_results)

categories <- c("good_sensitivity", "intr_poor_sensitivity", "good_predictor_value", "intr_poor_predictor_value")
min_result_table <- do.call("cbind", lapply(categories, function(cat) {
                        cat_table <- do.call("cbind", lapply(all_results, function(res) {
                            res[,cat]
                        }))
                        vals <- apply(cat_table, 1, function(x) { min(x) })
                    }))
colnames(min_result_table) <- paste("min_", categories, sep="")
aggregated_results <- data.frame(bucket_id=train_results[,"bucket"], cp_val=train_results[,"cp_val"], min_result_table, npairs=train_results[,"npairs"])

write.csv(aggregated_results, aggregated_results_file)
