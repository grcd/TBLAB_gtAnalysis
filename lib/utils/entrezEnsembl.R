
mapping = "NCBI_HGNC/hgnc_reduced_set.txt"			# required, complete_set downloaded frm: http://www.genenames.org/


entrez2Ensembl <- function(entrezs) {
	# Returns ensembl_gene_ids and missing of the given entresz
	
	if (! file.exists(mapping))
		stop("No mapping file available.")
		
	mp = read.delim(mapping, stringsAsFactor=FALSE)[, c( "ensembl_gene_id", "entrez_id" )]
	
	ensembls = c()
	missing  = c();
	
	for (i in 1:length(entrezs)) {
		element = which(mp$entrez_id == entrezs[i], arr.ind=TRUE)
		if (length(element) == 0) {
			missing = c(missing, entrezs[i])
		} else {
			ensembls = c(ensembls, mp$ensembl_gene_id[ element ] )			
		}
	}
	
	return( list(ensembles=ensembls, missing=missing) )
}


ensembl2Entrez <- function(ensembls) {
	# Returns ensembls of the given entresz and missing
	
	if (! file.exists(mapping))
		stop("No mapping file available.")

	mp = read.delim(mapping, stringsAsFactor=FALSE)[, c( "ensembl_gene_id", "entrez_id" )]
	
	entrezs = c()
	missing = c()
	
	for (i in 1:length(ensembls)) {
		element = which(mp$ensembl_gene_id == ensembls[i], arr.ind=TRUE)
		if (length(element) == 0) {
			missing = c(missing, ensembls[i])
		} else {
			entrezs = c(entrezs, mp$entrez_id[ element ] )			
		}
	}
	
	return( list(entrezs=entrezs, missing=missing) )
}
