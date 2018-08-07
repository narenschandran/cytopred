desc_extract <- function(gse, extraction_feature) {
                	gsmlist <- GSMList(gse)
                	sapply(gsmlist, function(gsm) {
                		Meta(gsm)[[extraction_feature]]
                	})
                }

library(GEOquery)
library(readxl)

dir.create(input_path)
dir.create(proc_path)
dir.create(proc_sample_path)


# For isAML, only samples that can proven to be not MDS is taken. This means that
# samples with unknown FAB are sometimes not selected because some datasets don't 
# demarcate between MDS and AML, and unknown ones can fall in either. Better be safe
# than sorry and remove unknown FAB which might contain RAEB (MDS) samples

sample_description <- sapply(GEO_microarray_datasets, function(dset) {
                      	dset_path    <- file.path(input_path, dset)
                      	destdir_path <- file.path(dset_path, "soft_file")
                      	dir.create(dset_path)
                      	dir.create(destdir_path)
                      	print(paste("Getting soft file for:", dset))
                      	gse <- getGEO(dset, GSEMatrix=F, destdir=destdir_path)
                      	sample_id   <- desc_extract(gse, "geo_accession")
                      	platform_id <- desc_extract(gse, "platform_id")

                      	if (dset=="GSE10358") {
                      		extr   <- desc_extract(gse, "characteristics_ch1")
                      		cnames <- unique(sapply(unlist(extr), function(x) {
                      		          	strsplit(x, ":")[[1]][1]
                      		          }))
                      		cnames <- cnames[!grepl("Tumor", cnames)]
                      		dat    <- do.call("rbind", lapply(extr, function(ext) {
                      		          	datc <- sapply(ext, function(x) {
                      		          	        	strsplit(x, ": ")[[1]][1]
                      		          	        })
                      		          	sapply(cnames, function(cname) {
                      		          		i <- datc %in% cname
                      		          		s <- if (TRUE %in% i) {
                      		          		     	strsplit(ext[i], ": ")[[1]][2]
                      		          		     } else {
                      		          		     	"extract_not_available"
                      		          		     }
                      		          	})
                      		          }))
                      		# Karyotype data is in three columns. Trying to
                      		# merge them here
                      		k_ind <- colnames(dat) %in% c("cytogenetic change", 
                      		                              "cytogenetics", 
                      		                              "karyotype")
                      		k_dat <- apply(dat[,k_ind], 1, function(x) {
                      		         	gsub(pattern="extract_not_available", 
                      		         	     replacement="", 
                      		         	     paste(x, collapse=""))
                      		         })
                      		a_ind <- colnames(dat) %in% c("age on study",
                      		                              "age")
                      		a_dat <- apply(dat[,a_ind], 1, function(x) {
                      		         	gsub(pattern="extract_not_available", 
                      		         	     replacement="", 
                      		         	     paste(x, collapse=""))
                      		         })
                      		n_ind <- colnames(dat) %in% c("npm1 ins",
                      		                              "nmp1 ins")
                      		n_dat <- apply(dat[,n_ind], 1, function(x) {
                      		         	sub(pattern="extract_not_available",
                      		         	    replacement="",
                      		         	    paste(x, collapse=""))
                      		         })
                      		isAML   <- rep(TRUE, nrow(dat))
                      		CNkars  <- c("Normal", "normal")
                      		isCN <- k_dat %in% CNkars
                      		karyotype_available <- !(k_dat %in% prognosis[[dset]][["unknown"]])
                      		desc  <- data.frame(sample_id=sample_id,
                      		                    platform_id=platform_id,
                      		                    age=a_dat,
                      		                    dat[,!(k_ind | a_ind | n_ind)],
                      		                    npm1_ins=n_dat, 
                      		                    karyotype=k_dat,
                      		                    isAML=isAML,
                      		                    isCN=isCN,
                      		                    karyotype_available=karyotype_available)
                      		colnames(desc)[colnames(desc)=="X.bm.blast"] <- "bm_blast_percentage"
                      		colnames(desc)[colnames(desc)=="X.pb.blast"] <- "pb_blast_percentage"
                      		colnames(desc)[colnames(desc)=="pb wbc@ presentation"] <- "pb_wbc_at_presentation"
                      	} else if (dset=="GSE12417") {
                      		extr <- desc_extract(gse, "characteristics_ch1")
                      		dat  <- t(sapply(extr, function(ex) {
                      		        	x <- strsplit(ex, "; ")[[1]]
                      		        	c(strsplit(x[1], "[,]")[[1]][1],
                      		        	  strsplit(x[1], ") ")[[1]][2],
                      		        	  strsplit(x[2], " =")[[1]][2],
                      		        	  strsplit(x[3], "= ")[[1]][2],
                      		        	  strsplit(x[4], ": ")[[1]][2],
                      		        	  strsplit(x[1], " [(]")[[1]][1]
                      		        	  )
                      		        }))
                      		colnames(dat) <- c("disease", "fab", "age", "os_days", 
                      		                   "simple_os_status_1_is_dead", "karyotype")
                      		isAML <- sapply(dat[,"disease"], function(d) {
                      		         	grepl("AML", d)
                      		         })
                      		isCN <- sapply(dat[,"karyotype"], function(d) {
                      		           	grepl("normal", d)
                      		           })
                      		crisk   <- rep("unknown", nrow(dat))
                      		crisk[isCN] <- "intr"
                      		karyotype_available <- rep(TRUE, nrow(dat))
                      		
                      		desc <- data.frame(sample_id=sample_id,
                      		                   platform_id=platform_id,
                      		                   dat, isAML=isAML, 
                      		                   isCN=isCN,
                      		                   crisk_code=crisk,
                      		                   karyotype_available=karyotype_available)
                      		                   
                      	} else if (dset=="GSE13159") {
                      		extr <- desc_extract(gse, "characteristics_ch1")
                      		cnames <- unique(sapply(unlist(extr), function(x) {
                      		          	strsplit(x, ":")[[1]][1]
                      		          }))
                      		dat    <- apply(extr, 1, function(ext) {
                      		          	sapply(ext, function(ex) {
                      		          		strsplit(ex, ": ")[[1]][2]
                      		          	})
                      		          })
                        	colnames(dat) <- make.names(cnames)
                        	isAML <- grepl("AML", dat[,"leukemia.class"])
                        	isCN  <- rep(F, nrow(dat))
                        	karyotype_available <- rep(T, nrow(dat))
                        	desc   <- data.frame(sample_id=sample_id,
                        	                     platform_id=platform_id,
                        	                     dat, karyotype=dat[,"leukemia.class"],
                        	                     isAML=isAML, isCN=isCN,
                        	                     karyotype_available=karyotype_available)
                      	} else if (dset=="GSE14468") {
                      		extr   <- desc_extract(gse, "characteristics_ch1")
                      		cnames <- unique(unlist(lapply(extr[1:461], function(ext) {
                      		          	sapply(ext, function(ex) {
                      		          		strsplit(ex, ": ")[[1]][1]
                      		          	})
                      		          })))
                      		dat1   <- do.call("rbind", lapply(extr[1:461], function(ext) {
                      		          	sapply(cnames, function(cname) {
                      		          		i <- grepl(cname, ext)
                      		          		if (TRUE %in% i) {
                      		          			strsplit(ext[i], ": ")[[1]][2]
                      		          		} else {
                      		          			"extract_not_available"
                      		          		}
                      		          	})
                      		          }))
                      		dat2   <- t(sapply(extr[462:526], function(ext) {
                      		          	c("Acute myeloid leukemia",
                      		          	  rep("extract_not_available", 3),
                      		          	  strsplit(ext, ", ")[[1]],
                      		          	  rep("extract_not_available", 8))
                      		          }))
                      		title_extr <- desc_extract(gse, "title")
                      		vnummer    <- sapply(title_extr, function(ext) {
                      		              	if (grepl("AML", ext)) {
                      		              		strsplit(ext, "AML ")[[1]][2]
                      		              	} else {
                      		              		ext
                      		              	}
                      		              })
                      		dat <- data.frame(sample_id=sample_id,
                      		                  volgnummer=vnummer,
                      		                  platform_id=platform_id,
                      		                  rbind(dat1, dat2))
                      		dat[is.na(dat)] <- "extract_not_available"
                      		dat[,"karyotype"] <- as.character(dat[,"karyotype"])
                      		dat[,"risk"] <- as.character(dat[,"risk"])
                      		fixed_kdat      <- data.frame(read_excel(GSE14468_fixed_karyotype))
                      		fixed_i <- match(fixed_kdat[,"sample_id"], dat[,"sample_id"])
                      		dat[fixed_i, "karyotype"] <- fixed_kdat[,"karyotype"]
                      		dat[fixed_i, "risk"]      <- fixed_kdat[,"risk"]
                      		isAML <- !(
                      		          grepl("RAEB", as.character(dat[,"score"])) 
                      		          |
                      		          grepl("unknown", as.character(dat[,"score"])) 
                      		          |
                                      grepl("MDS", as.character(dat[,"karyotype"]))
                                      )
                      		isCN  <- grepl("NN", as.character(dat[,"karyotype"])) 
                      		crisk <- rep("unknown", nrow(dat))
                      		crisk[grepl("good", dat[,"risk"])] <- "good"
                      		crisk[grepl("intermediate", dat[,"risk"])] <- "intr"
                      		crisk[grepl("poor", dat[,"risk"])] <- "poor"
                      		karyotype_available <- !(dat[,"karyotype"]=="failure" |
                      		                         dat[,"karyotype"]=="extract_not_available")
                      		sdat <- data.frame(read_excel(GSE14468_surv_file))
                      		match_i <- match(dat$volgnummer, sdat$volgnummer)
                      		sdat    <- sdat[match_i,]                          
                      		sdat$volgnummer <- mapply(function(sval, cval) {
                      		                   	if (is.na(sval)) {
                      		                   		cval
                      		                   	} else {
                      		                   		if(sval==cval) {
                      		                   			sval
                      		                   		} else {
                      		                   			"mismatch"
                      		                   		}
                      		                   	}
                      		                   }, as.character(sdat$volgnummer), 
                      		                      as.character(dat$volgnummer))
                      		sdat[is.na(sdat)] <- "extract_not_available"                        
                      		colnames(sdat)[1] <- "matched_volgnummer"
                      		desc <- data.frame(dat, sdat, crisk_code=crisk,
                      		                   isAML=isAML, isCN=isCN,
                      		                   karyotype_available=karyotype_available)
                      		colnames(desc)[colnames(desc)=="score"] <- "fab"
                      	} else if (dset=="GSE15434") {
                      		extr <- t(desc_extract(gse, "characteristics_ch1"))
                      		cnames <- apply(extr, 2, function(ext) {
                      		          	unique(unlist(sapply(ext, function(ex) {
                      		          		strsplit(ex, ": ")[[1]][1]
                      		          	})))
                      		          })
                      		dat   <- t(apply(extr, 1, function(ext) {
                      		         	sapply(cnames, function(cname) {
                      		         		i <- grepl(cname, ext)
                      		         		strsplit(ext[i], ": ")[[1]][2]
                      		         	})
                      		         }))
                      		isAML <- !grepl("MDS", as.character(dat[,"disease state"]))
                      		isCN  <- grepl("normal", as.character(dat[,"diagnosis"]))
                      		crisk <- rep("unknown", nrow(dat))
                      		crisk[isCN] <- "intr"
                      		karyotype_available <- rep(TRUE, nrow(dat))
                      		desc <- data.frame(sample_id=sample_id,
                      		                   platform_id=platform_id,
                      		                   dat,
                      		                   karyotype=dat[,"diagnosis"],
                      		                   crisk_code=crisk,
                      		                   isAML=isAML,
                      		                   isCN=isCN,
                      		                   karyotype_available=karyotype_available)
                      	} else if (dset=="GSE16015") {
                      		extr <- t(desc_extract(gse, "characteristics_ch1"))
                      		cnames <- apply(extr, 2, function(ext) {
                      		          	unique(unlist(sapply(ext, function(ex) {
                      		          		strsplit(ex, ": ")[[1]][1]
                      		          	})))
                      		          })
                      		dat   <- t(apply(extr, 1, function(ext) {
                      		         	sapply(cnames, function(cname) {
                      		         		i <- grepl(cname, ext)
                      		         		strsplit(ext[i], ": ")[[1]][2]
                      		         	})
                      		         }))
                      		isAML <- rep(TRUE, nrow(dat))
                      		isCN  <- grepl("normal", dat[,"karyotype"])
                      		crisk <- rep("unknown", nrow(dat))
                      		crisk[isCN] <- "intr"
                      		karyotype_available <- rep(TRUE, nrow(dat))
                      		desc <- data.frame(sample_id=sample_id,
                      		                   platform_id=platform_id,
                      		                   dat, crisk_code=crisk,
                      		                   isAML=isAML, isCN=isCN,
                      		                   karyotype_available=karyotype_available)
                      	} else if (dset=="GSE22845") {
                       		extr   <- desc_extract(gse, "characteristics_ch1")
                       		cnames <- unique(sapply(extr, function(ext) {
                       		          	strsplit(ext, ":")[[1]][1]
                       		          }))
                       		dat    <- t(sapply(extr, function(ext) {
                       		          	if (grepl(cnames, ext, fixed=T)) {
                       		          		strsplit(strsplit(ext, ": ")[[1]][2],
                       		          		         ",")[[1]]
                       		            } else {
                       		            	rep("check_extraction_code", 2)
                       		            }
                       		          }))
                       		colnames(dat) <- c("cebpa_single_mutatnt",
                       		                   "cebpa_double_mutant")
                       		title_dat <- desc_extract(gse, "title")
                       		# check treatment_protocol_ch1 in Meta for proof of:
                       		isAML <- rep(TRUE, nrow(dat))
                       		isCN  <- grepl("normal karyotype", title_dat)
                       		crisk <- rep("unknown", nrow(dat))
                       		crisk[isCN] <- "intr"
                       		karyotype_available <- rep(TRUE, nrow(dat))
                       		desc    <- data.frame(sample_id=sample_id,
                       		                      platform_id=platform_id,
                       		                      dat, karyotype=title_dat,
                       		                      crisk_code=crisk,
                       		                      isAML=isAML,
                       		                      isCN=isCN,
                       		                      karyotype_available=karyotype_available)
                      	} else if (dset=="GSE30285") {
                      		extr <- t(desc_extract(gse, "characteristics_ch1"))
                      		cnames <- apply(extr, 2, function(ext) {
                      		          	unique(unlist(sapply(ext, function(ex) {
                      		          		strsplit(ex, ": ")[[1]][1]
                      		          	})))
                      		          })
                      		dat   <- t(apply(extr, 1, function(ext) {
                      		         	sapply(cnames, function(cname) {
                      		         		i <- grepl(cname, ext, fixed=T)
                      		         		if (TRUE %in% i) {
                      		         			strsplit(ext[i], ": ")[[1]][2]
                      		         		} else {
                      		         			"extract_not_available"
                      		         		}
                      		         	})
                      		         }))
                      		isAML <- grepl("AML", as.character(dat[,"subtype"]))
                      		# Check the dat[,"abnormality"] for proof of:
                      		isCN  <- rep(FALSE, nrow(dat))
                      		karyotype_available <- rep(TRUE, nrow(dat))
                      		desc   <- data.frame(sample_id=sample_id,
                      		                     platform_id=platform_id,
                      		                     dat,
                      		                     karyotype=dat[,"abnormality"],
                      		                     isAML=isAML, isCN=isCN,
                      		                     karyotype_available=karyotype_available)
                      	} else if (dset=="GSE39730") {
                      		extr <- desc_extract(gse, "characteristics_ch1")
                      		dat  <- apply(extr, 1, function(ext_row) {
                      		        	sapply(ext_row, function(x) {
                      		        		strsplit(x, ": ")[[1]][2]
                      		        	})
                      		        })
                      		colnames(dat) <- c("disease", "karyotype",
                      		                    "TP53_mutation_status",
                      		                    "mir34a_expression")
                      		isAML <- grepl("AML", dat[,"disease"])
                      		isCN  <- rep(FALSE, nrow(dat))
                      		crisk <- rep("poor", nrow(dat))
                      		karyotype_available <- rep(TRUE, nrow(dat))
                      		desc <- data.frame(sample_id, platform_id, dat,
                      		                   isAML=isAML, isCN=isCN,
                      		                   crisk_code=crisk,
                      		                   karyotype_available=karyotype_available)
                      	} else if (dset=="GSE61804") {
                      		extr <- t(desc_extract(gse, "characteristics_ch1"))
                      		cnames <- apply(extr, 2, function(ext) {
                      		          	unique(unlist(sapply(ext, function(ex) {
                      		          		strsplit(ex, ": ")[[1]][1]
                      		          	})))
                      		          })
                      		dat   <- t(apply(extr, 1, function(ext) {
                      		         	sapply(cnames, function(cname) {
                      		         		i <- grepl(cname, ext, fixed=T)
                      		         		if (TRUE %in% i) {
                      		         			strsplit(ext[i], ": ")[[1]][2]
                      		         		} else {
                      		         			"extract_not_available"
                      		         		}
                      		         	})
                      		         }))
                      		isAML <- grepl("AML", as.character(dat[,"condition"]))
                      		isCN  <- grepl("normal karyotype", as.character(dat[,"condition"]))
                      		karyotype_available <- rep(TRUE, nrow(dat))
                      		desc   <- data.frame(sample_id=sample_id,
                      		                     platform_id=platform_id,
                      		                     dat,
                      		                     karyotype=dat[,"condition"],
                      		                     isAML=isAML, isCN=isCN,
                      		                     karyotype_available=karyotype_available)
                      	} else {
                      		"No information given on how to process dataset"
                      	}
                      	if (dset %in% names(prognosis)) {
                      		dset_prog <- prognosis[[dset]]
                      		crisk     <- rep("unknown", nrow(desc))
                      		for (dsp in names(dset_prog)) {
                      			p_i        <- desc$karyotype %in% dset_prog[[dsp]]
                      			crisk[p_i] <- dsp
                      		}
                      	desc <- data.frame(desc, crisk_code=crisk)
                      	}
                      	desc$platform_id <- as.character(desc$platform_id)
                      	lapply(unique(desc$platform_id), function(p_id) {
                      		i <- desc$platform_id==p_id
                      		outdat <- desc[i,]
                      		outfile <- file.path(proc_sample_path, paste(dset, "_", p_id, "_sample_desc.csv", sep=""))
                      		write.csv(outdat, file=outfile)
                      	})
                      })
