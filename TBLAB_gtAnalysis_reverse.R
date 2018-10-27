
#
#   TBLAB Globaltest REVERSE model Analysis (y=mRNA, X= putative miRNAs regulators)
#

#rm(list=ls())
#WORKINGDIR = "/Users/danielegreco/Desktop/CNR_TCGA"
#setwd(WORKINGDIR)
library(globaltest)
                                                       # utils
source("lib/utils/readExpressionData.R")
source("lib/utils/entrezEnsembl.R")                    # not used anymore
source("lib/utils/retrievePutativeTargets_SQLite.R")   # not used anynore

source("lib/miRNAmRNA/R/mirnamrna.R")                  # globaltest related
source("lib/miRNAmRNA/R/db_functions.R")
source("lib/miRNAmRNA/R/parsers.R")


################################################################
#                                                              #
#                         SETUP:                               #
#                                                              #
################################################################

# 1.1. Matrici di espressione ottenute da BuildTCGASets.R:
#
#   X: mRNA   [mRNAs]  x [samples]
#   Z: miRNA  [miRNAs] x [samples]
#
#   Le matrici sono accoppiate: le colonne di X e Z provengono dagli gli stessi individui.

X.file 	= file.path(EXPRESSION_DIR, "X")                  # mRNA expression matrix  
Z.file 	= file.path(EXPRESSION_DIR, "Z")                  # miRNA expression matrix 
X.names = file.path(EXPRESSION_DIR, "mrnas.txt")          # rownames(X), mRNA labels
Z.names	= file.path(EXPRESSION_DIR, "mirnas.txt")         # rownames(Z), miRNA labels
samples.names = file.path(EXPRESSION_DIR, "samples.txt")  # colnames(X), sample unique ID

# 1.2. Interazioni putative
#
# Le interazioni putative fornite da TBLAB sono state riversate in un db 
# SQLite3 (sl3) locale e successivamente spostate in un data-frame (PI) di R.
#
# Per rendere il processing ancora più veloce, è stato conveniente pre-estrarre
# le interazioni putative relative ad un particolare Batch (93) e ad una
# particolare condizione clinica (T = Tumor) e costruire la matrice di incidenza
# PI.
#
# PI[i, j] = TRUE se il target mRNA j (colonne) è un target predetto per il
#            miRNA regolatore i (righe); FALSE altrimenti.
#
# PI.path è il percorso al file R contenente la matrice PI.
#PI.path = file.path(getwd(), "datasets/TBLAB_PI_BRCA_Batch93/BRCA_Batch93_PI_T_Rdata")   # Tumor   samples
#PI.path = file.path(getwd(), "TBLAB_PI_BRCA_Batch93/BRCA_Batch93_PI_N_Rdata")  # Normal  samples
PI.path = PUTATIVE_FILE   # Tumor samples


# 1.3  Mapping EntrezID <--> EnsemblID for protein-coding mRNAs
#
# I mRNA di TCGA (matrice X) sono dati con nomenclatura EntrezID.
# I mRNA di TBLAB (PI) sono dati con nomenclatura EnsemblID.
# Un mapping tra i due mondi è dato da HGNC (http://www.genenames.org/).
# HGNC.proteincoding è una versione "ridotta all'osso" del file originario
# scaricabile da HGNC, generato dallo script ausiliario lib/utils/getProteinCodingGenes.R .
#HGNC.proteincoding 	= "datasets/NCBI_HGNC/hgnc_protein_coding_reduced.txt"  
HGNC.proteincoding 	= HGNC_FILE  


# 1.4  False Discovery Rate for pvalues adjustments
#FDR.alpha = 0.01
FDR.alpha = FDR_ALPHA  # FDR threshold

# 1.5  Risultati in output
if (RESULTS_DIR == "") {
  results.PATH = sprintf("results_%s", basename(getwd()))  # output directory
} else {
  results.PATH = RESULTS_DIR
}

################################################################
#                                                              #
#                               CODE                           #
#                                                              #
################################################################

cat(sprintf("------------------------------------------------------\n"))
cat(sprintf("globaltest ANALYSIS: REVERSE model (y=mRNAs, X=miRNAs)\n"))
cat(sprintf("    response:   single mRNA\n"))
cat(sprintf("    covariates: predicted regulating miRNAs\n"))
cat(sprintf("------------------------------------------------------\n\n"))

# Se non esiste, crea la cartella dove scrivere i risultati
if (!dir.exists(results.PATH)) {
    dir.create(results.PATH)
    if (!dir.exists(results.PATH))
        stop("** ERROR: unable to create results path. Permission denied.")
}

# Inizializzo un log file nel quale annoto, passo passo, le modifiche ai dati
# operate dallo script (se ve ne sono)
logFile = "log_file_REVERSE.txt"
ll = sprintf("# This is a log file for %s globaltest analysis (reverse model).\n", basename(getwd()))
cat(ll, file=file.path(results.PATH, logFile), append=FALSE, sep = "\n")

#	1. Read expression data (readExpressionData.R)
cat(sprintf("-- Reading expression data...\n"))
data = readExpressionData(X.file, Z.file, X.names, Z.names, samples.names)
ll = sprintf("Reading expression data...\n+ Number of samples: %d\n+ Number of mRNAs: %d\n+ Number of miRNAs: %d\n", ncol(data$X), nrow(data$X), nrow(data$Z))
cat(ll, file=file.path(results.PATH, logFile), append=TRUE, sep = "\n")

#	2.	Extract *RNA and miRNAs identifiers*
mrnas.entrez = rownames(data$X)
mirnas = rownames(data$Z) 

#	3.	Filter only	mRNAs that are: 
#       (a) currently marked as 'alive' by HGNC 
#       (b) protein-coding 
#
#  NOTA:  
#  Il sequenziatore estrae mRNA coding e non coding. Dato che vogliamo correlare
#  in un secondo momento le quantità di mRNA espresso con le quantità di proteine
#  corrispodenti, ho selezionato e filtrato solo gli mRNA alive e protein-coding.

# Rimuovi righe contenenti NAs o entry vuote
hgnc = read.delim2(HGNC.proteincoding, stringsAsFactors=FALSE, 	colClasses=c("character", "character"))
remove_nas      = union( which( is.na(hgnc$entrez_id) ), which( is.na(hgnc$ensembl_gene_id) ))
remove_blanks   = union( which( hgnc$entrez_id == "" ), which( hgnc$ensembl_gene_id == "" ))
remove_all      = union(remove_nas, remove_blanks)                    
hgnc = hgnc[ -remove_all, ]

# Prendo solo i coding dai dati di espressione
coding = intersect(mrnas.entrez, hgnc$entrez_id)

# Quelli non-coding scrivili in RNA_noncoding.csv per successiva consultazione
write.csv(file=file.path(results.PATH, "RNA_noncoding.csv"), 
          data.frame(Entrez=setdiff(mrnas.entrez, coding)), quote=FALSE, row.names=FALSE)

# Riduci la matrice di espressione dei mRNA per avere solo i mRNA coding
data$X = data$X[ coding, ]
mrnas.entrez = coding
rm(coding)

# 4a.  Conversione EntrezID --> EnsemblID
rownames(hgnc) = hgnc$entrez_id
mrnas.ensembl =  hgnc[ mrnas.entrez, "ensembl_gene_id"]		

# 4b.  Retrieve putative interactions (PI)
# Il file PI.path è un data.frame R contenente le interazioni putative pre-estratte.
cat(sprintf("-- Retrieving PI (putative interactions) boolean matrix from TBLAB...\n"))
load(file=PI.path)

#	5.  GLOBALTEST (gt) ANALYSIS 
#
# Il gt viene calcolato dalla runrgt.TBLAB di miRNAmRNA/R/mirnamrna.R.
# Il codice originario (dell'articolo) è stato modificato per utilizzare le
# interazioni putative fornite in PI. Il funzionamento è descritto nel file.
#
# Nota il nome della funzione: run r gt.TBLAB (r=reverse model).
#
# I risultati raw dell'analisi sono scritti in rr_raw che viene salvato in
# un file R, prima di ulteriori modifiche.

rownames(data$X) = mrnas.ensembl  
t1=Sys.time()
cat("-- Running globaltest (reversed)...") # 19 mins!
rr_raw  <- runrgt.TBLAB(trgs=rownames(data$X), 
                        X=data$X, 
                        Y=data$Z, 
                        path="", 
                        dbName="", 
                        credentials=NULL, 
                        putatives=PI,
                        debug=DEBUG)
save(rr_raw, file=file.path(results.PATH, "rr_raw.Rdata") )
cat("OK!\n")
t2=Sys.time(); difftime(t2, t1)

# 6.  Adjust p-values using an FDR (BH) approach (use FDR.alpha parameter)
#
# Per aggiustare i p-value usando la FDR threshold viene chiamata toTable.TBLAB.
# Ancora una volta, si tratta di una versione modificata di toTable fornita con
# l'articolo (miRNAmRNA/R/mirnamrna.R).
#
cat(sprintf("-- Adjusting p-values by BH (FDR < %3.2f)\n", FDR.alpha))
rr_adj = toTabler.TBLAB(rr_raw, alpha=FDR.alpha)
save(rr_adj, file=file.path(results.PATH, "rr_FDR_all.Rdata") )

#   7.  Seleziona interazioni significative e ordina per p-value non-decrescenti
rr_adj = subset(rr_adj, Pvalue_global < FDR.alpha)  # select significant
rr_adj = subset(rr_adj, Pvalue_pair < FDR.alpha)

ii = order(rr_adj$Pvalue_global, rr_adj$Pvalue_pair) # order by
rr_adj = rr_adj[ii, ]

save(rr_adj, file=file.path(results.PATH, "rr_FDR_significant.Rdata") ) # save

# 8.  Scriviamo i dati formattati
#
# Scriviamo due file tsv (tab-separated value) separati: 
#    rr_FDR_FINAL_negative.Rdata:  coppie mRNA-miRNA anticorrelate 
#    rr_FDR_FINAL_positive.Rdata:  coppie mRNA-miRNA correlate
cat("-- Writing out results: reverse model...\n")

hgnc_r = hgnc
rownames(hgnc_r) = hgnc_r$ensembl_gene_id
rr_adj = cbind(hgnc_r[ rr_adj$mRNA, ]$entrez_id, rr_adj)

colnames(rr_adj)[1] = "mRNA.entrez"
colnames(rr_adj)[3] = "mRNA.ensembl"
rr_adj = rr_adj[, c(2,3,1,4:6)]

if ( length(which(rr_adj$Association == 0)) > 0)
  rr_adj[ which(rr_adj$Association == 0), ]$Association = "negative"

if ( length(which(rr_adj$Association == 1)) > 0)
  rr_adj[ which(rr_adj$Association == 1), ]$Association = "positive"

rr_adj_negative = subset(rr_adj, Association == "negative") # anti-correlati
rr_adj_positive = subset(rr_adj, Association == "positive") # correlati

save(rr_adj_negative, file=file.path(results.PATH, "rr_FDR_FINAL_negative.Rdata"))
save(rr_adj_positive, file=file.path(results.PATH, "rr_FDR_FINAL_positive.Rdata"))

write.table(rr_adj_negative,
            file=file.path(results.PATH, "rr_FDR_FINAL_negative.tsv"),
            quote=F, sep='\t', row.names=F)

write.table(rr_adj_positive,
            file=file.path(results.PATH, "rr_FDR_FINAL_positive.tsv"),
            quote=F, sep='\t', row.names=F)


# 9. Statistiche: numero di coppie putative vs numero di coppie significative
#
# Scriviamo qualche statistica sul file di log.
# Il numero di coppie filtrate utilizzando i dati di espressione è calcolato
# sommando correlate ed anti-correlate.

putative_pairs = length(which(PI == TRUE))
narrowed_pairs = nrow(rr_adj_positive)+nrow(rr_adj_negative)

ll = sprintf("\n* REVERSEMODEL: PUTATIVE predicted pairs: %d", putative_pairs)
cat(ll, file=file.path(results.PATH, logFile), append=TRUE, sep = "\n")

ll = sprintf("* REVERSEMODEL: FILTERED significant pairs: %d", narrowed_pairs)
cat(ll, file=file.path(results.PATH, logFile), append=TRUE, sep = "\n")
