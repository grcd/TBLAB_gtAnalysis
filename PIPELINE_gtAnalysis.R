
# PIPELINE gtAnalysis

rm(list=ls())

################################################################
#                                                              #
#                         PARAMETERS                           #
#                                                              #
################################################################

# WORKING_DIR: root directory
WORKING_DIR = "/Users/danielegreco/Desktop/CNR_TCGA"

# EXPRESSION_DIR: output directory resulting from buildTCGASets the given expression dataset 
EXPRESSION_DIR = "datasets/TCGA_BRCA_Batch93/buildTCGASets/T"  # (T = Tumor set)

# PUTATIVE_FILE: Rdata file containing predicted putative target mRNAs for miRNAs listed in the expression dataset
PUTATIVE_FILE = "datasets/TBLAB_PI_BRCA_Batch93/BRCA_Batch93_PI_T_Rdata"

# HGNC_FILE: EnsembleId-EntrezId mapping file for protein coding genes
HGNC_FILE = "datasets/NCBI_HGNC/hgnc_protein_coding_reduced.txt"

# FDR_ALPHA: confidence level for pvalues FDR adjustament/filtering 
FDR_ALPHA = 0.01

# RESULTS_DIR: output directory for this analysis
RESULTS_DIR = ""

# DEBUG
DEBUG=NULL


################################################################
#                                                              #
#                               CODE                           #
#                                                              #
################################################################

setwd(WORKING_DIR)

# Perform global test direct (original) analysis: Y=miRNAs, X=predicted target mRNAs
source("TBLAB_gtAnalysis_original.R")

# Perform global test reverse (original) analysis: Y=mRNAs, X=predicted potentially regulating miRNAs
source("TBLAB_gtAnalysis_reverse.R")

# Get miRNA:mRNA significant pairs for both
load("results_CNR_TCGA/rd_FDR_significant.Rdata") #rd_adj
load("results_CNR_TCGA/rr_FDR_significant.Rdata") #rr_adj
total_pairs <- merge(rd_adj, rr_adj, by=c("miRNA", "mRNA", "Association"), suffixes=c(".rd", ".rr"))
write.table(total_pairs,
            file=file.path(results.PATH, "FINAL_significant_pairs.tsv"),
            quote=F, sep='\t', row.names=F)

