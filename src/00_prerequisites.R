# For reading in files
library(readxl)                                                                                                                                                                                 
library(data.table)

# For downloading, processing and reading in microarray data
library(GEOquery)                                                                
library(oligo)
library(hgu133plus2hsentrezg.db)  # File is from Brain Array Custrom CDFv21
library(huex10sthsentrezg.db)     # File is from Brain Array Custrom CDFv21
library(hgu133bhsentrezg.db)      # File is from Brain Array Custrom CDFv21
library(hgu133ahsentrezg.db)      # File is from Brain Array Custrom CDFv21

# For reading in salmon RNAseq data
library("tximport")

# For survival analysis
library(survival)
library(survminer)

# For generating decision trees
library(rpart)                                                                   
library(rpart.plot)

# For parallelizing work
library(parallel)                                                                            


