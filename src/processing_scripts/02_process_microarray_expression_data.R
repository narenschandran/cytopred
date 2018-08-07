library(GEOquery)                                                                
library(oligo)
library(hgu133plus2hsentrezg.db)                                                 
library(huex10sthsentrezg.db)                                                    
library(hgu133bhsentrezg.db)                                                     
library(hgu133ahsentrezg.db)  

dir.create(input_path)
dir.create(proc_path)
dir.create(proc_expr_path)

probe_annotations <- list()
probe_annotations[["GPL96"]] <- toTable(hgu133ahsentrezgSYMBOL)
probe_annotations[["GPL97"]] <- toTable(hgu133bhsentrezgSYMBOL)
probe_annotations[["GPL570"]] <- toTable(hgu133plus2hsentrezgSYMBOL)
probe_annotations[["GPL5175"]]  <- toTable(huex10sthsentrezgSYMBOL)

array_info <- read.csv(file=array_info_file, colClasses=c("character"))

exp_data <- lapply(GEO_microarray_datasets, function(dset) {
            	dset_path <- file.path(input_path, dset)
            	exp_path  <- file.path(dset_path, "RAW_files")
            	dir.create(dset_path)
            	dir.create(exp_path)
            	desc_files <- list.files(path=proc_sample_path, pattern=dset)
            	desc       <- do.call("rbind", lapply(desc_files, function(desc_file) {
   	        	              	dfile_path <- file.path(proc_sample_path, desc_file)
   	        	              	desc      <- read.csv(dfile_path, colClasses="character")
   	        	              }))

            	fname  <- paste(file.path(exp_path, paste(dset, "_RAW.tar", sep="")))
            	if (!file.exists(fname)) {
            		getGEOSuppFiles(GEO=dset, makeDirectory=F, baseDir=exp_path)
            	} else {
            		print(paste("Using already downloaded file for", dset))
            	}
            	untar(tarfile=fname, exdir=exp_path)

            	all_celfiles <- file.path(exp_path, list.files(exp_path, 
            	                                               pattern="*.cel.gz", 
            	                                               ignore.case=T))
            	p_ids <- unique(desc$platform_id)
            	edat <- lapply(p_ids, function(p_id) {
            	        	pkgname  <- array_info$pkgname[array_info$platform_id==p_id]
            	        	pr_annot <- probe_annotations[[p_id]]
            	        	s_id_ind <- desc$platform_id==p_id & as.logical(desc$isAML)==TRUE
            	        	s_ids    <- desc$sample_id[s_id_ind]
            	        	s_index <- sapply(all_celfiles, function(cfile) {
            	        	           	cf_yes <- sapply(s_ids, function(s_id) {
            	        	           	          	grepl(s_id, cfile)
            	        	           	          })
            	        	           	TRUE %in% cf_yes
            	        	           })
            	        	cfiles <- all_celfiles[s_index]
            	        	dat    <- exprs(rma(read.celfiles(cfiles,
            	        	                                  pkgname=pkgname)))
            	        	annot <- sapply(rownames(dat), function(pr_id){
            	        	       	 	if (pr_id %in% pr_annot$probe_id) {
            	        	       	 		pr_annot$symbol[pr_annot$probe_id==pr_id]
            	        	       	 	} else {
            	        	       	 		pr_id
            	        	       	 	}
            	        	       	 })
            	        	rownames(dat) <- annot
            	        	dup_genes <- rownames(dat)[duplicated(rownames(dat))]     
            	        	dat <- dat[!(rownames(dat) %in% dup_genes),]   
            	        })
            	names(edat) <- p_ids
            	edat
            })
names(exp_data) <- GEO_microarray_datasets
#save(exp_data, file="exp_data")

subset_edat <- lapply(names(exp_data), function(dset) {
               	edat   <- exp_data[[dset]]
               	p_ids  <- names(edat)
               	s_edat <- lapply(p_ids, function(p_id) {
               	          	e      <- edat[[p_id]]
               	          	cnames <- colnames(e)
               	          	if (dset=="GSE10358" & p_id=="GPL570") {
               	          		cnames1 <- sapply(cnames[1:25], function(cn) {
               	          		           	strsplit(cn, "_")[[1]][1]
               	          		           })
               	          		cnames2 <- sapply(cnames[26:209], function(cn) {
               	          		           	strsplit(cn, "[.]")[[1]][1]
               	          		           })
               	          		cnames3 <- sapply(cnames[210:300], function(cn) {
               	          		           	strsplit(cn, "_")[[1]][1]
               	          		           })
               	          		final_cn <- c(cnames1, cnames2, cnames3)
               	          	} else if ((dset=="GSE12417" & p_id=="GPL96")    |
               	          	           (dset=="GSE12417" & p_id=="GPL97")    |
               	          	           (dset=="GSE12417" & p_id=="GPL570")   |
               	          	           (dset=="GSE13159" & p_id=="GPL570")   |
               	          	           (dset=="GSE14468" & p_id=="GPL570")   |
               	          	           (dset=="GSE15434" & p_id=="GPL570")   |
               	          	           (dset=="GSE16015" & p_id=="GPL570")   |
               	          	           (dset=="GSE22845" & p_id=="GPL570")   |
             	          	           (dset=="GSE30285" & p_id=="GPL5175")) {
               	          		final_cn <- sapply(cnames, function(cn) {
               	          		            	strsplit(cn, "[.]")[[1]][1]
               	          		            })
               	          	} else if ((dset=="GSE39730" & p_id=="GPL570")  |
               	          	           (dset=="GSE61804" & p_id=="GPL570")) {
               	          		final_cn <- sapply(cnames, function(cn) {
               	          		            	strsplit(cn, "_")[[1]][1]
               	          		            })
               	          	} else {
               	          		print(paste("No special formatting done for",
               	          		            dset, p_id))
               	          		final_cn <- cnames
               	          	}
               	          	colnames(e) <- final_cn
               	          	e
               	          })
               	print(dset)
               	names(s_edat) <- paste(dset, p_ids, sep="_")
               	s_edat
               })
names(subset_edat) <- names(exp_data)

lapply(names(subset_edat), function(dset) {
	s_edat <- subset_edat[[dset]]
	lapply(names(s_edat), function(sdset) {
		fname <- file.path(proc_expr_path, paste(sdset, 
		                                         "_expression.csv", sep=""))
		write.csv(s_edat[[sdset]], fname)
	})
})
		









               	          		
               	          		
               	          		
               	          

