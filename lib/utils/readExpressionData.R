
readExpressionData <- function( X.mat, Z.mat, 
								X.names="", Z.names="", samples.names="") {
	# Reads sample-matched expression data from X.mat (mRNAs), Z.mat (miRNAs)
	#
	# Args:
	#   X.mat: mRNA  expression matrix [mRNAs x samples] filename
	#   Z.mat: miRNA expression matrix [miRNAs x samples] filename
	# Optional:
	#	X.names: rownames for X.mat (mRNA ids/names) filename
	#	Z.names: rownames for Z.mat (miRNA ids/names) filename
	#	samples.names: colnames for X and Z (sample matched)
	#
	#
	# Returns:
	#	A list containing the two sample-matched expression matrices X and Z
	
	if (! file.exists(X.mat))
		stop("No X.mat file")
	
	if (! file.exists(Z.mat))
		stop("No Z.mat file")
	
	X = read.delim(X.mat, sep='\t', stringsAsFactors=F)
	if (X.names != "") {
		rows = read.table(X.names, header=TRUE, stringsAsFactors=F)
		rownames(X) = rows[, 1]
	}

	Z = read.delim(Z.mat, sep='\t', stringsAsFactors=F)
	if (Z.names != "") {
		rows = read.table(Z.names, header=TRUE, stringsAsFactors=F)
		rownames(Z) = rows[, 1]
	}

	if (samples.names != "") {
		samples = read.table(samples.names, header=TRUE, stringsAsFactors=F)
		colnames(X) = samples[, 1]
		colnames(Z) = samples[, 1]
	}

	return(list(X=X, Z=Z))
		
}