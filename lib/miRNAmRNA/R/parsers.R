##' Mapping between various gene identifiers using .db-package.
##'
##' Details follow.
##' 
##' @title Map identifiers
##' @param x vector of ids type 'from' 
##' @param from name of the current identifier
##' @param to name of identifier mapped to
##' @param organism organism
##' @return vector of ids type 'to' 
##' @author Maarten van Iterson
##' @export
##' @importClassesFrom AnnotationDbi
##' @importMethodsFrom AnnotationDbi mget
mapIds <- function(x, from, to="ENTREZ", organism=c("Hs", "Mm"))
{
  
  x <- as.character(x)
  
  require(paste("org.", organism, ".eg.db", sep=""), character.only=TRUE) #Library needed for mapping  
  if(from != "ENTREZ")
    {
      xx <- mget(x, get(paste("org.", organism, ".eg", from, "2EG", sep="")), ifnotfound=NA)
      x <- sapply(xx, function(x) unlist(x)[1], USE.NAMES=FALSE) #in case of multiple just use one    
    }
  
  if(to == "ENTREZ")
    return(x)

  xx <- mget(x[!is.na(x)], get(paste("org.", organism, ".eg", to, sep="")), ifnotfound=NA)
  y <- sapply(xx, function(x) unlist(x)[1], USE.NAMES=FALSE) #in case of multiple just use one    
  x[!is.na(x)] <- y 
  x
}

pita <- function(file, Org, full=TRUE)
{
  datafull <- read.table(file, sep="\t", header = TRUE, comment.char="", as.is=TRUE)
  
  columns <- c("microRNA", "RefSeq" ,"Score")
  cid <- colnames(datafull) %in% columns
  data <- datafull[, columns]
  
  data$Target <- sub(";.*", "", data$RefSeq) #required regex to reduce multiple ID's to one ID    

  data$mRNA <- mapIds(data$Target, from="REFSEQ", to="ENTREZ", organism=Org)   
  data$Symbol <- mapIds(data$Target, from="REFSEQ", to="SYMBOL", organism=Org) 
  
  columns <- c("microRNA", "mRNA", "Symbol", "Target", "Score")
  data <- data[, columns]
  
  colnames(data) <- c("miRNA", "mRNA", "Symbol", "Target", "Score")
  data$Score <- -1*as.numeric(data$Score) ##binding energy the more negative the better
  
  ##add all information
  if(full)
    data <- cbind(data, datafull[, !cid])
  
  data <- data[order(data$Score),]  
  data
}

targetscan <- function(file, Org, full=TRUE)
{
    datafull <- read.table(file, sep="\t", header = TRUE, comment.char="", as.is=TRUE, na.string="NULL")

    columns <- c("mirbase_id", "Gene.ID") #version 4.1
    ##columns <- c("miRNA", "Gene.ID", "context..score") #version 6.1
    ##columns <- c("miRNA", "Gene.ID", "context_score") #version 6.0
   
    contextscore <- colnames(datafull)[grep("context.*score", colnames(datafull))]
    
    columns <- c("miRNA", "Gene.ID", contextscore)
    
    cid <- colnames(datafull) %in% columns
    
    Species <- colnames(datafull)[grep("Gene.Tax.ID|Species.ID|UTR_tax_id", colnames(datafull))]
    
    rows <- datafull[,Species] == ifelse(Org == "Mm", 10090, 9606)
    
    data <- datafull[rows,  columns]
   
    data$Symbol <- mapIds(as.character(data$Gene.ID), from="ENTREZ", to="SYMBOL", organism="Mm") 
    
    columns <- c("miRNA", "Gene.ID", "Symbol", contextscore)
    
    data <- data[, columns]
    
    colnames(data) <- c("miRNA", "mRNA", "Symbol", "Score")
        
    data$Score <- -1*as.numeric(data$Score) #reverse score orientation because targetscan has a negative score scaling

    ##add all information
    if(full)
      data <- cbind(data, datafull[rows, !cid])
    
    data <- data[order(data$Score), 1:4]    
    data
}


gffRead <- function(gffFile, nrows = -1)
{

  gff <- read.table(gffFile, sep="\t", as.is=TRUE, quote="", header=FALSE, comment.char="#", nrows = nrows, 
                   colClasses=c("character", #SEQ
                                "character", #METHOD
                                "character", #FEATURE
                                "integer",   #START 
                                "integer",   #END
                                "numeric",   #SCORE
                                "numeric",   #PVALUE_OG
                                "character", #STRAND 
                                "character", #PHASE 
                                "character")) #ATTRIBUTES
  
  colnames(gff) <- c("seqname", "method", "feature", "start", "end",  "score", "pvalue_OG", "strand", "phase", "attributes")

  stopifnot(!any(is.na(gff$start)), !any(is.na(gff$end)))
  return(gff)
}

microcosm <- function(file, Org, full=TRUE)
  {
    #   Legge correttamente il file gff e riassegna i nomi delle colonne
    datafull <- gffRead(file)    
    
    #   Se l'organismo (Org) è Mm (MusMusculus) i miRNA hanno il prefisso mmu
    #   e preleva soltant questi
    rows <- grep(ifelse(Org=="Mm", "mmu", "hsa"), datafull$seqname)       
    data <- datafull[rows, ]
    colnames(data)[10] <- "Ensembl"
    
    extractID <- function(x) gsub('\\"', "", gsub(".*gene=", "", x))
    data$Target <- sapply(data$Ensembl, extractID, USE.NAMES=FALSE)

    data$mRNA <- mapIds(data$Target, from="ENSEMBLTRANS", to="ENTREZ", organism=Org) 
    data$Symbol <- mapIds(data$Target, from="ENSEMBLTRANS", to="SYMBOL", organism=Org) 
    
    columns <- c("seqname", "mRNA", "Symbol", "Target", "score")
    cid <- colnames(datafull) %in% columns
    
    data <- data[, columns]

     ##add all information
    if(full)
      data <- cbind(data, datafull[rows, !cid])
    
    colnames(data)[1:5] <- c("miRNA", "mRNA", "Symbol", "Target", "Score")
    data$Score <- as.numeric(data$Score)
    data <- data[order(data$Score),]    
    data
  }

