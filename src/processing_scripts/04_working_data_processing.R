library(data.table)

dir.create(analysis_root)
dir.create(analysis_path)
dir.create(working_data_path)
dir.create(working_sample_path)
dir.create(working_expr_path)
dir.create(results_path)

all_expr_files <- list.files(proc_expr_path)
all_expr_dsets <- sub("_expression.csv", "", all_expr_files)

all_sdesc_files <- list.files(proc_sample_path)
all_sdesc_dsets <- sub("_sample_desc.csv", "", all_sdesc_files)

# I'm doing this just in case either expression data or sample desc has more 
# datasets (like in the case of one being processed and the other not)
common_dsets <- Reduce("intersect", list(all_expr_dsets, all_sdesc_dsets, analysis_dsets))

expr_i  <- grepl(paste(common_dsets, collapse="|"), all_expr_files) &
           grepl(paste(allowed_p_ids, collapse="|"), all_expr_files)

sdesc_i <- grepl(paste(common_dsets, collapse="|"), all_sdesc_files) &
           grepl(paste(allowed_p_ids, collapse="|"), all_sdesc_files)

expr_files  <- all_expr_files[expr_i]
sdesc_files <- all_sdesc_files[sdesc_i]

sdesc <- lapply(sdesc_files, function(sdesc_file) {
         	print(paste("Reading in sample description for", sdesc_file))
         	data.frame(fread(file.path(proc_sample_path, sdesc_file), 
         	           stringsAsFactors=F), row.names=1)
         })
names(sdesc) <- sdesc_files

ssamples <- unlist(lapply(sdesc, function(desc) {
         		i <- desc$isAML==T & desc$karyotype_available==T
         		desc$sample_id[i]
         	}))

expr <- lapply(expr_files, function(expr_file) {
         	print(paste("Reading in expression data for", expr_file))
        	data.frame(fread(file.path(proc_expr_path, expr_file), 
        	           stringsAsFactors=F), row.names=1)
        })
names(expr) <- expr_files

common_genes <- Reduce(intersect, lapply(expr, function(edat) {
                	rownames(edat)
                }))

esamples <- unlist(lapply(expr, function(edat) {
            	colnames(edat)
         	}))

common_samples <- intersect(ssamples, esamples)

final_desc <- lapply(sdesc, function(desc) {
              	i <- desc$sample_id %in% common_samples
              	fdesc <- desc[i,]
              	fdesc <- fdesc[order(fdesc$sample_id),]
              })
 
final_expr <- lapply(expr, function(edat) {
              	i <- rownames(edat) %in% common_genes
              	j <- colnames(edat) %in% common_samples
              	fedat <- edat[i,j]
              	fedat <- fedat[,order(colnames(fedat))]
              	fedat <- fedat[order(rownames(fedat)),]
              })

# Sanity check: Should come up TRUE for all entries
#sapply(1:length(final_expr), function(i) {
#    print(i)
#	check <- colnames(final_expr[[i]]) == final_desc[[i]]$sample_id
#	!(FALSE %in% check)
#})

lapply(names(final_desc), function(fname) {
	desc <- final_desc[[fname]]
	write.csv(desc, file.path(working_sample_path, fname))
})

lapply(names(final_expr), function(fname) {
	edat <- final_expr[[fname]]
	write.csv(edat, file.path(working_expr_path, fname))
})

working_desc <- do.call("rbind", lapply(1:length(final_desc), function(i) {
                	dset <- sub("_sample_desc.csv", "", names(final_desc)[i])
                	
                	desc <- final_desc[[i]]
                	data.frame(sample_id=desc$sample_id,
                	           dset=rep(dset, nrow(desc)),
                	           isAML=desc$isAML, isCN=desc$isCN,
                	           karyotype_available=desc$karyotype_available,
                	           karyotype=desc$karyotype,
                	           crisk_code=desc$crisk_code)
                }))

names(final_expr) <- NULL
working_expr <- do.call("cbind", final_expr)

#write.csv(working_desc, file.path(analysis_path, "sample_description.csv"))
#write.csv(working_expr, file.path(analysis_path, "expression_data.csv"))
write.csv(working_desc, working_desc_file)
write.csv(working_expr, working_expr_file)

