
library(RSQLite) 

retrievePutativeInteractions <- function(mirnas, mrnas, targetDB) {
    # Queries targetDB and returns a |mirnas| x |mrnas| incidence matrix
    #
    # Args:
    #	mirnas: vector of mirnas identifiers suited for targetDB
    #   targetDB: (filename) an edge-list of (mirna, putative_target) pairs
    #
    # Returns:
    #	C: 	a binary incidence matrix having mirnas in rows, mrnas in columns
    #		having C(i,j) = 1 if and only if mrnas[j] is a predicted target for
    #		mirnas[i]
    #
    # Notes:
    #	Assumes targetDB has exactly three columns: mrna, mirna, interaction_score
    #
    
    if (! file.exists(targetDB))
        stop("targetDB does not exists or permission denied.")
    
    
    #   Establish db connection
    con = dbConnect(RSQLite::SQLite(), dbname=targetDB)
    
    
    #	Initialize incidence matrix		
    C = matrix(FALSE, nrow=length(mirnas), ncol=length(mrnas))
    rownames(C) = mirnas
    colnames(C) = mrnas
    
    
    #	Query putative interactions
    query = sprintf("SELECT DISTINCT mRNA, LOWER(miRNA) as miRNA FROM predictions WHERE mRNA IN ( %s ) INTERSECT SELECT DISTINCT mRNA, LOWER(miRNA) as miRNA FROM predictions WHERE LOWER(miRNA) IN ( %s )", paste(shQuote(mrnas), collapse=","), paste(shQuote(mirnas), collapse=",") )
    putatives = dbGetQuery( con, query )

    
    #	*Magic*: in a single line we turn ON putative interactions
    C[ cbind(putatives$miRNA, putatives$mRNA) ] = TRUE
    #colnames(C) = ensembl2Entrez(colnames(C))			 # columns (mRNAs) are translated back to the EntrezId format
    
    invisible(dbDisconnect(con))
    return(C)
}


retrievePutativeTargets.gt <- function(mirnas, targetDB,
                                    tau=FALSE) {
    # Queries targetDB and returns a |mirnas| x |mrnas| incidence matrix
    #
    # Args:
    #	mirnas: vector of mirnas identifiers suited for targetDB
    #   targetDB: (filename) an edge-list of (mirna, putative_target) pairs
    #	tau:		a threshold for filtering putative interactions ($score, if any)
    #
    #
    # Returns:
    #	C: 	a binary incidence matrix having mirnas in rows, mrnas in columns
    #		having C(i,j) = 1 if and only if mrnas[j] is a predicted target for
    #		mirnas[i]
    #
    # Notes:
    #	Assumes targetDB has exactly three columns: mrna, mirna, interaction_score
    #
    
    if (! file.exists(targetDB))
        stop("targetDB does not exists or permission denied.")
    
    
    #   Establish db connection
    con = dbConnect(RSQLite::SQLite(), dbname=targetDB)
    
    
    #	Initialize incidence matrix		
    C = matrix(FALSE, nrow=length(mirnas), ncol=length(mrnas))
    rownames(C) = mirnas
    colnames(C) = mrnas
    
    
    #	Query putative interactions
    #   query1 = sprintf("SELECT DISTINCT lower(miRNA), mRNA FROM predictions WHERE (lower(miRNA) IN ( %s )) ",  paste(shQuote(mirnas), collapse=","))
    queryPutative = sprintf("SELECT DISTINCT mRNA, miRNA FROM predictions WHERE (lower(miRNA) IN ( %s ))", paste(shQuote(mrnas), collapse=","))
    
    #print(queryPutative)	# DEBUG
    putatives = dbGetQuery( con, queryPutative )
    putatives$miRNA = tolower(putatives$miRNA)
    
    
    #	*Magic*: in a single line we turn ON putative interactions
    C[ cbind(putatives$miRNA, putatives$Target) ] = TRUE
    #colnames(C) = ensembl2Entrez(colnames(C))			 # columns (mRNAs) are translated back to the EntrezId format
    
    
    return(as.numeric(C))
}

retrievePutativeTargets <- function(mirnas, mrnas, targetDB,
									tau=FALSE) {
	# Queries targetDB and returns a |mirnas| x |mrnas| incidence matrix
	#
	# Args:
	#	mirnas: vector of mirnas identifiers suited for targetDB
	#   mrnas:  vector of mrnas identifiers  suited for targetDB
	#   targetDB: (filename) an edge-list of (mirna, putative_target) pairs
	#	tau:		a threshold for filtering putative interactions ($score, if any)
	#
	#
	# Returns:
	#	C: 	a binary incidence matrix having mirnas in rows, mrnas in columns
	#		having C(i,j) = 1 if and only if mrnas[j] is a predicted target for
	#		mirnas[i]
	#
	# Notes:
	#	Assumes targetDB has exactly three columns: mrna, mirna, interaction_score
	#
	
	if (! file.exists(targetDB))
		stop("targetDB does not exists or permission denied.")
	
	
    #   Establish db connection
    con = dbConnect(RSQLite::SQLite(), dbname=targetDB)
    
    
	#	Initialize incidence matrix		
	C = matrix(FALSE, nrow=length(mirnas), ncol=length(mrnas.ensembl))
	rownames(C) = mirnas
	colnames(C) = mrnas.ensembl
	
	
	#	Query putative interactions
	queryPutative = sprintf("SELECT DISTINCT mRNA, miRNA FROM predictions WHERE (lower(miRNA) IN ( %s ))", paste(shQuote(mirnas), collapse=","))
	
	#print(queryPutative)	# DEBUG
	putatives = dbGetQuery( con, queryPutative )

		
	#	*Magic*: in a single line we turn ON putative interactions
	C[ putatives$miRNA, putatives$Target ] = TRUE
	#colnames(C) = ensembl2Entrez(colnames(C))			 # columns (mRNAs) are translated back to the EntrezId format
	
	
	return(as.numeric(C))
}