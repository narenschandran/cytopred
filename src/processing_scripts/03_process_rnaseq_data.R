library("tximport")

all_dsets <- list.files(input_path)
dsets <- all_dsets[all_dsets %in% rnaseq_datasets]

rnaseq_dat <- lapply(dsets, function(dset) {
              	dset_path <- file.path(input_path, dset)
              	if (dset=="LAML_TCGA") {
              		sdesc_file <- file.path(dset_path, "data_clinical.txt")
              		edat_file  <- file.path(dset_path, "LAML.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt")
              		sdesc <- read.csv(sdesc_file, sep="\t", skip=5,
              		                  stringsAsFactors=F)
              		colnames(sdesc) <- tolower(colnames(sdesc))
              		sdesc$sample_id <- make.names(sdesc$sample_id)
              		isAML <- rep(T, nrow(sdesc))
              		isCN  <- sdesc$cytogenetic_code_other=="Normal Karyotype"
              		karyotype_available <- !(sdesc$cytogenetic_code_other=="N.D.")
              		platform_id <- rep("sequencing", nrow(sdesc))
              		crisk <- rep("unknown", nrow(sdesc))
              		crisk[sdesc$risk_cyto=="Good"]         <- "good"
              		crisk[sdesc$risk_cyto=="Intermediate"] <- "intr"
              		crisk[sdesc$risk_cyto=="Poor"]         <- "poor"
              		sdesc <- data.frame(sdesc, platform_id=platform_id,
              		                    karyotype=sdesc$cytogenetic_code_other,
              		                    crisk_code=crisk,
              		                    isAML=isAML, isCN=isCN,
              		                    karyotype_available=karyotype_available)
              		# I'm removing the first row because it just has entries with
              		# the string "gene_id" in the first column and "normalized count"
              		# in all the rest
              		edat  <- read.csv(edat_file, sep="\t", row.names=1, header=T,
              		                  stringsAsFactors=F)[-1,]
              		genes <- sapply(rownames(edat), function(rn) {
              		         	strsplit(rn, "[|]")[[1]][1]
              		         })
              		d_genes <- genes[duplicated(genes)]
              		removal_index <- genes %in% d_genes
              		edat  <- edat[!removal_index,]
              		rownames(edat) <- genes[!removal_index]
              		colnames(edat) <- sapply(colnames(edat), function(cn) {
              		                  	paste(strsplit(cn, "[.]03[A/B]")[[1]][1],
              		                  	      ".03", sep="")
              		                  })
              	} else if (dset=="SRP050272") {
              		tx2gene <- read.csv(file.path(dset_path, "tx2gene.csv"))
              		epath   <- file.path(dset_path, "salmon_output")
              		efiles  <- file.path(epath, 
              		                     list.files(path=epath, pattern="SRX"),
              		                     "quant.sf")
              		txi     <- tximport(efiles, type="salmon", tx2gene=tx2gene)
              		edat    <- txi$abundance 
              		colnames(edat) <- list.files(path=epath, pattern="SRX")
              		rownames(edat) <- make.names(rownames(edat))

              		platform_id <- rep("sequencing", ncol(edat))
              		karyotype   <- rep("Normal karyotype", ncol(edat))
              		crisk_code  <- rep("intr", ncol(edat))
              		isAML       <- rep(TRUE, ncol(edat))
              		isCN        <- rep(TRUE, ncol(edat))
              		karyotype_available <- rep(TRUE, ncol(edat))
              		sdesc   <- data.frame(sample_id=colnames(edat),
              		                      platform_id=platform_id,
              		                      karyotype=karyotype,
              		                      crisk_code=crisk_code,
              		                      isAML=isAML, isCN=isCN,
              		                      karyotype_available=karyotype_available)

              	} else {
              		print(paste("No instructions on how to process", dset))
              	}
              	write.csv(edat, 
              	          file.path(proc_expr_path, 
              	                    paste(dset, "_rnaseq_expression.csv", 
              	                          sep="")))
              	write.csv(sdesc, 
              	          file.path(proc_sample_path, 
              	                    paste(dset, "_rnaseq_sample_desc.csv", 
              	                          sep="")))
 
              	list(expression_data=edat, sample_desc=sdesc)

              })

names(rnaseq_dat) <- dsets
