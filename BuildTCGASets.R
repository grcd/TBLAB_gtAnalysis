
#	BUILDTCGASETS.R
#
#	Una volta scaricato il Batch da TCGA, decomprimerlo e far partire lo script nella stessa cartella. 
#	Verrà creata una directory buildTCGASets/ al cui interno troviamo:
#		
#		buildTCGASets/	file_manifest_T.txt				file_manifest.txt filtered sample	Tumor
#						file_manifest_N.txt												    Normal
#						file_manifest.TN.txt												Tumor-Normal matched			
#						samples_missing.mirbase20.txt		campioni esclusi per mancanza di miRNASeq data
#						samples_missing.rsem.txt			campioni esclusi per mancanza di RNASeqV2 data
#
#		buildTCGASets/T/	X_T				matrice [mRNAs]  x [samples]
#						Z_T				matrice [miRNAs] x [samples]	
#						T_mirnas.txt		lista dei mirnas
#						T_mrnas.txt		lista degli mrnas
#						T_samples.txt	lista dei campioni/SampleID
#						
#					 N/	{come sopra}
#					 TN/	{come sopra}
#	

rm(list=ls())
setwd("/Users/danielegreco/Desktop/CNR_TCGA/TCGA_BRCA_Batch93")
library(dplyr)


#	PARAMETERS
MANIFEST 		= "file_manifest.txt"	
outputDIR 		= "buildTCGASets"


#	CODE
cat("Building TCGA datasets...\n")

if (! dir.exists(outputDIR))
	dir.create(outputDIR, recursive=FALSE)

if (dir.exists(outputDIR)) {
	dir.create( paste(outputDIR, "T",  sep=.Platform$file.sep), recursive=FALSE )
	dir.create( paste(outputDIR, "N",  sep=.Platform$file.sep), recursive=FALSE )
	dir.create( paste(outputDIR, "TN", sep=.Platform$file.sep), recursive=FALSE )
} else {
	stop("** ERROR: unable to create output directory, permission denied.")
}


#	------------------------------------------------------------------------------
#
#	Phase0:	PREPROCESSING		setup the stage...
#
#	------------------------------------------------------------------------------

#	1. Open file_manifest, filter out useless rows
manifest = read.delim(MANIFEST, stringsAsFactor=FALSE)
manifest = subset( manifest, Platform.Type=="miRNASeq" | Platform.Type=="RNASeqV2")


#	2. Filter out $Sample(s) with missing mirbase20 expression data
#	a. ottengo la lista dei campioni che hanno mirbase20
mirbase20_samples = unique( (manifest[ grepl('mirbase20.mirna.quantification', manifest$File.Name) ,])$Sample  )

#	b. li sottraggo dalla lista di tutti i campioni e quindi ottengo i campioni che **NON** hanno mirbase20; e li scrivo su un file esterno.
missing.mirbase20 = setdiff( unique(manifest$Sample), mirbase20_samples )
write.table(missing.mirbase20, 
			sprintf("%s%s%s", outputDIR, .Platform$file.sep, "samples_missing.mirbase20.txt"), 
			sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

#	c. tolgo tutte le righe che si riferiscono ai campioni missing.mirbase20 
if (length(missing.mirbase20) > 0)		# ATTENZIONE: se la tua espressione regolare è vuota, lui cancella tutta la matrice
	manifest = manifest[ !grepl( pattern=paste(missing.mirbase20, collapse='|'), manifest$Sample), ]
rm(mirbase20_samples)
rm(missing.mirbase20)


#	3. Repeat (2) with 'rsem.genes.results' (RNASeqV2 expression data)

rsem_samples = unique( (manifest[ grepl('rsem.genes.results', manifest$File.Name) ,])$Sample  )
missing.rsem = setdiff( unique(manifest$Sample), rsem_samples )
write.table(missing.rsem, 
			sprintf("%s%s%s", outputDIR, .Platform$file.sep, "samples_missing.rsem.txt"), 
			sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

if (length(missing.rsem) > 0)
	manifest = manifest[ !grepl( pattern=paste(missing.rsem, collapse='|'), manifest$Sample), ]
rm(rsem_samples)
rm(missing.rsem)


#	4. Filter out other expression data (junction, exon, isoforms)
manifest = manifest[ grepl('mirbase20.mirna.quantification|rsem.genes.results', manifest$File.Name) ,]


#	------------------------------------------------------------------------------
#
#	Phase1:	BUILDING SETS		manifest.TN, manifest.T, manifest.N
#
#	------------------------------------------------------------------------------

manifest$normal 		= as.logical(0)
manifest$Participant 	= ""

manifest.TN = data.frame();
manifest.T 	= data.frame();
manifest.N 	= data.frame();

#	mark all the normal tissues, and extract up to Partipant 'Aliquot-slice'
for (i in 1:nrow(manifest)) {	# TODO: 	I GUESS it can be done in single code-line with some lapply stuff
	tmp = unlist(strsplit(manifest$Sample[i], split="-"))
	manifest$Participant[i] = sprintf("%s-%s-%s", tmp[1], tmp[2], tmp[3])
}

manifest$normal = FALSE
manifest[ grep('TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-1[0-9]{1}' , manifest$Sample) ,]$normal = TRUE		# selects -0*
#manifest[ grep('TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-1[0-9]1' , manifest$Sample) ,]$normal = TRUE		# selects -01 only 


#	BUILDING manifest.T		tumor  (only) samples
manifest.T = manifest[ manifest$normal == FALSE, ]
for (sample in unique(manifest.T$Sample)) { # double-check
	if ( nrow( manifest.T[ manifest.T$Sample == sample, ] ) != 2) {
		cat(sprintf("** [manifest.T] WARNING: SAMPLE=%s could not contain all the data needed, please check it out. \n", sample  ))		
	}
}

#	BUILDING manifest.N		normal (only) samples
manifest.N = manifest[ manifest$normal == TRUE, ]
for (sample in unique(manifest.N$Sample)) { # double-check
	if ( nrow( manifest.N[ manifest.N$Sample == sample, ] ) != 2) {
		cat(sprintf("** [manifest.N] WARNING: SAMPLE=%s could not contain all the data needed, please check it out. \n", sample  ))		
	}
}

#	BUILDING manifest.TN	tumor-normal matched samples
matchedParticipants = unique(manifest.N$Participant)
manifest.TN = rbind( manifest.N, 
                     manifest.T[ manifest.T$Participant %in% matchedParticipants, ] )

for (participant in unique(manifest.TN$Participant)) { # double-check
	if ( nrow( manifest.TN[ manifest.TN$Participant == participant, ] ) != 4) { # 4=2 (tumor, rsem+mirbase) + 2 (normal, same as before) files 
		cat(sprintf("** [manifest.TN] WARNING: PARTICIPANT=%s could not contain all the data needed, please check it out. \n", participant  ))		
	}
}

#	Write out manifest(s) and a report
write.table(manifest.T, 
			sprintf("%s%s%s", outputDIR, .Platform$file.sep, "file_manifest_T.txt"), 
			quote=FALSE, sep='\t', row.names=FALSE)

write.table(manifest.N, 
			sprintf("%s%s%s", outputDIR, .Platform$file.sep, "file_manifest_N.txt"), 
			quote=FALSE, sep='\t', row.names=FALSE)

write.table(manifest.TN, 
			sprintf("%s%s%s", outputDIR, .Platform$file.sep, "file_manifest_TN.txt"), 
			quote=FALSE, sep='\t', row.names=FALSE)

generated_set = c("file_manifest_T.txt","file_manifest_N.txt","file_manifest_TN.txt")
description = c("Tumor only samples", "Normal only samples", "Tumor-Normal matched samples")
number_of_samples = c(length(unique(manifest.T$Sample)), length(unique(manifest.N$Sample)), length(unique(manifest.TN$Sample)))
report = data.frame(generated_set, description, number_of_samples)
write.table(report, 
			sprintf("%s%s%s", outputDIR, .Platform$file.sep, "buildTCGASets_REPORT.txt"), 
			quote=FALSE, sep='\t', row.names=FALSE)

rm(generated_set); rm(description); rm(number_of_samples); rm(report)


#	------------------------------------------------------------------------------
#
#	Phase2:	BUILDING MATRICES		
#
#	------------------------------------------------------------------------------


#	------------------------------------------------------------------------------
#	Tumor set:	X_T		mRNA 		mrnas  x samples
#				Z_T		miRNA		mirnas x samples
#	------------------------------------------------------------------------------

cat("* Generating set T: ")
samples = sort(unique(manifest.T$Sample))
mat.directory = sprintf("%s%s%s", outputDIR, .Platform$file.sep, "T")

#
#	Z_T, miRNAs
#

#	Generate the intersection of mirnas across all the samples (maybe not needed)
mirnas = c()
for (sample in samples) {
	
	#	Open the miRNASeq expression file
	row = subset( manifest.T, Platform.Type == "miRNASeq" & Sample == sample )
	mirbase20.file = sprintf("%s/%s/%s__%s/Level_%s/%s", 
						getwd(),
						row$Platform.Type,
						row$Center,
						row$Platform,
						row$Level,
						row$File.Name)
	mr = (read.delim(mirbase20.file, stringsAsFactor=FALSE))$miRNA_ID

	#	Intersect miRNAs in this file with miRNAs so far collected	
	if (length(mirnas) == 0) {
		mirnas = mr
	} else {
		mirnas = intersect(mirnas, mr);	
	}
	
	# cat(sprintf("SAMPLE=%s, |mirnas|=%d\n", sample, length(mirnas)))
}

#	Now, we can safely construct the matrix
Z_T = matrix(0.0, nrow=length(mirnas), ncol=length(samples))
rownames(Z_T) = mirnas
colnames(Z_T) = samples

#	For each sample, we open the mirbase20 expression file again and take **only** expression values for miRNAs in mirnas vector
for (sample in samples) {
	
	#	Open the miRNASeq expression file
	row = subset( manifest.T, Platform.Type == "miRNASeq" & Sample == sample )
	mirbase20.file = sprintf("%s/%s/%s__%s/Level_%s/%s", 
						getwd(),
						row$Platform.Type,
						row$Center,
						row$Platform,
						row$Level,
						row$File.Name)

	#	Read the miRNA expression data
	mr = (read.delim(mirbase20.file, stringsAsFactor=FALSE))[, c("miRNA_ID", "reads_per_million_miRNA_mapped")]

	#	Label rows with miRNA_IDs, needed for the subsequent step
	rownames(mr) = mr$miRNA_ID
	
	#	In a single step:	select only counts relative to the common set of mirnas obtained above 
	Z_T[, sample] = mr[mirnas, ]$reads_per_million_miRNA_mapped

}

write.table(Z_T, 
			paste(mat.directory, "Z_T", sep=.Platform$file.sep),
			quote=FALSE, sep='\t')

write.table(data.frame(miRNA_ID=mirnas), 
			paste(mat.directory, "T_mirnas.txt", sep=.Platform$file.sep), 
			quote=FALSE, row.names=FALSE)	
					
rm(mirnas); rm(mr)


#
#	X_T, mRNAs
#

#	Generate the intersection of mrnas across all the samples
mrnas = c()
for (sample in samples) {
	
	#	Open the RNASeqV2 expression file
	row = subset( manifest.T, Platform.Type == "RNASeqV2" & Sample == sample )
	rsem.file = sprintf("%s/%s/%s__%s/Level_%s/%s", 
						getwd(),
						row$Platform.Type,
						row$Center,
						row$Platform,
						row$Level,
						row$File.Name)
	mr = (read.delim(rsem.file, stringsAsFactor=FALSE))$gene_id				
	mr = mr[ !grepl('\\?', mr) ]								# filter out unknown ? gene entries
	mr = sapply(strsplit(as.character(mr), "\\|"), "[[", 2)		# take only the gene Entrez id (unique)

	#	Intersect mRNAs in this file with mRNAs so far collected	
	if (length(mrnas) == 0) {
		mrnas = mr
	} else {
		mrnas = intersect(mrnas, mr);	
	}
	
	#cat(sprintf("SAMPLE=%s, |mrnas|=%d\n", sample, length(mrnas)))
	
}

#	Now, we can safely construct the matrix
X_T = matrix(0.0, nrow=length(mrnas), ncol=length(samples))
rownames(X_T) = mrnas
colnames(X_T) = samples

#	For each sample, we open the rsem.genes expression file again and take **only** expression values for miRNAs in mirnas
for (sample in samples) {
	
	#	Open the RNASeqV2 expression file
	row = subset( manifest.T, Platform.Type == "RNASeqV2" & Sample == sample )
	rsem.file = sprintf("%s/%s/%s__%s/Level_%s/%s", 
						getwd(),
						row$Platform.Type,
						row$Center,
						row$Platform,
						row$Level,
						row$File.Name)
	
	#	Read the mRNA expression data
	mr = (read.delim(rsem.file, stringsAsFactor=FALSE))[, c("gene_id", "scaled_estimate")]
	mr = mr[ !grepl('\\?', mr$gene_id), ]		
	mr$gene_id = sapply(strsplit(as.character(mr$gene_id), "\\|"), "[[", 2)

	#	Label rows with mRNA_IDs, needed for the subsequent step
	rownames(mr) = mr$gene_id
	
	#	In a single step:	select only counts relative to the common set of mrnas obtained above 
	X_T[, sample] = mr[mrnas, ]$scaled_estimate * 10^6      # TPM=scaled_estimate*10^6

}

write.table(X_T, 
			paste(mat.directory, "X_T", sep=.Platform$file.sep),
			quote=FALSE, sep='\t')

write.table(data.frame(mRNA_EntrezId=mrnas), 
			paste(mat.directory, "T_mrnas.txt", sep=.Platform$file.sep), 
			quote=FALSE, row.names=FALSE)			
			
write.table(data.frame(SampleID=samples), 
			paste(mat.directory, "T_samples.txt", sep=.Platform$file.sep), 
			quote=FALSE, row.names=FALSE)	

cat(sprintf("DONE. number_of_samples=%d, mRNAs=%d, miRNAs=%d\n", length(samples), nrow(X_T), nrow(Z_T)))
rm(mrnas); rm(mr); rm(X_T); rm(Z_T);



#	------------------------------------------------------------------------------
#	Normal   set:	X_N		mRNA 		mrnas  x samples
#				 	Z_N		miRNA		mirnas x samples
#	------------------------------------------------------------------------------

cat("* Generating set N: ")
samples = sort(unique(manifest.N$Sample))
mat.directory = sprintf("%s%s%s", outputDIR, .Platform$file.sep, "N")

#
#	Z_N, miRNAs
#

#	Generate the intersection of mirnas across all the samples
mirnas = c()
for (sample in samples) {
	
	#	Open the miRNASeq expression file
	row = subset( manifest.N, Platform.Type == "miRNASeq" & Sample == sample )
	mirbase20.file = sprintf("%s/%s/%s__%s/Level_%s/%s", 
						getwd(),
						row$Platform.Type,
						row$Center,
						row$Platform,
						row$Level,
						row$File.Name)
	mr = (read.delim(mirbase20.file, stringsAsFactor=FALSE))$miRNA_ID

	#	Intersect miRNAs in this file with miRNAs so far collected	
	if (length(mirnas) == 0) {
		mirnas = mr
	} else {
		mirnas = intersect(mirnas, mr);	
	}
	
	# cat(sprintf("SAMPLE=%s, |mirnas|=%d\n", sample, length(mirnas)))
	
}

#	Now, we can safely construct the matrix
Z_N = matrix(0.0, nrow=length(mirnas), ncol=length(samples))
rownames(Z_N) = mirnas
colnames(Z_N) = samples

#	For each sample, we open the mirbase20 expression file again and take **only** expression values for miRNAs in mirnas vector
for (sample in samples) {
	
	#	Open the miRNASeq expression file
	row = subset( manifest.N, Platform.Type == "miRNASeq" & Sample == sample )
	mirbase20.file = sprintf("%s/%s/%s__%s/Level_%s/%s", 
						getwd(),
						row$Platform.Type,
						row$Center,
						row$Platform,
						row$Level,
						row$File.Name)

	#	Read the entire miRNA expression file
	mr = (read.delim(mirbase20.file, stringsAsFactor=FALSE))[, c("miRNA_ID", "reads_per_million_miRNA_mapped")]

	#	Label rows with miRNA_IDs, needed for the subsequent step
	rownames(mr) = mr$miRNA_ID
	
	#	In a single step:	select only counts relative to the common set of mirnas obtained above 
	Z_N[, sample] = mr[mirnas, ]$reads_per_million_miRNA_mapped

}

write.table(Z_N, 
			paste(mat.directory, "Z_N", sep=.Platform$file.sep),
			quote=FALSE, sep='\t')

write.table(data.frame(miRNA_ID=mirnas), 
			paste(mat.directory, "N_mirnas.txt", sep=.Platform$file.sep), 
			quote=FALSE, row.names=FALSE)		

rm(mirnas); rm(mr)


#
#	X_N, mRNAs
#

#	Generate the intersection of mrnas across all the samples
mrnas = c()
for (sample in samples) {
	
	#	Open the RNASeqV2 expression file
	row = subset( manifest.N, Platform.Type == "RNASeqV2" & Sample == sample )
	rsem.file = sprintf("%s/%s/%s__%s/Level_%s/%s", 
						getwd(),
						row$Platform.Type,
						row$Center,
						row$Platform,
						row$Level,
						row$File.Name)
	mr = (read.delim(rsem.file, stringsAsFactor=FALSE))$gene_id				
	mr = mr[ !grepl('\\?', mr) ]									# filter out unknown ? gene entries
	mr = sapply(strsplit(as.character(mr), "\\|"), "[[", 2)		# take only the gene Entrez id (unique)

	#	Intersect mRNAs in this file with mRNAs so far collected	
	if (length(mrnas) == 0) {
		mrnas = mr
	} else {
		mrnas = intersect(mrnas, mr);	
	}
	
	#cat(sprintf("SAMPLE=%s, |mrnas|=%d\n", sample, length(mrnas)))
}

#	Now, we can safely construct the matrix
X_N = matrix(0.0, nrow=length(mrnas), ncol=length(samples))
rownames(X_N) = mrnas
colnames(X_N) = samples

#	For each sample, we open the rsem.genes expression file again and take **only** expression values for miRNAs in mirnas
for (sample in samples) {
	
	#	Open the miRNASeq expression file
	row = subset( manifest.N, Platform.Type == "RNASeqV2" & Sample == sample )
	rsem.file = sprintf("%s/%s/%s__%s/Level_%s/%s", 
						getwd(),
						row$Platform.Type,
						row$Center,
						row$Platform,
						row$Level,
						row$File.Name)
	
	#	Read the entire mRNA expression file
	mr = (read.delim(rsem.file, stringsAsFactor=FALSE))[, c("gene_id", "scaled_estimate")]
	mr = mr[ !grepl('\\?', mr$gene_id), ]		
	mr$gene_id = sapply(strsplit(as.character(mr$gene_id), "\\|"), "[[", 2)

	#	Label rows with mRNA_IDs, needed for the subsequent step
	rownames(mr) = mr$gene_id
	
	#	In a single step:	select only counts relative to the common set of mrnas obtained above 
	X_N[, sample] = mr[mrnas, ]$scaled_estimate * 10^6

}

write.table(X_N, 
			paste(mat.directory, "X_N", sep=.Platform$file.sep),
			quote=FALSE, sep='\t')

write.table(data.frame(mRNA_EntrezId=mrnas), 
			paste(mat.directory, "N_mrnas.txt", sep=.Platform$file.sep), 
			quote=FALSE, row.names=FALSE)		

write.table(data.frame(SampleID=samples), 
			paste(mat.directory, "N_samples.txt", sep=.Platform$file.sep),
			quote=FALSE, row.names=FALSE)	

cat(sprintf("DONE. number_of_samples=%d, mRNAs=%d, miRNAs=%d\n", length(samples), nrow(X_N), nrow(Z_N)))
rm(mrnas); rm(mr); rm(X_N); rm(Z_N);



#	------------------------------------------------------------------------------
#	Tumor-Normal   set:	X_TN		mRNA 		mrnas  x samples
#				     	Z_TN		miRNA		mirnas x samples
#	------------------------------------------------------------------------------

cat("* Generating set TN: ")
samples = sort(unique(manifest.TN$Sample))
samples_T = samples[ grepl('TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}-01', samples) ]
samples_N = setdiff( samples, samples_T )
mat.directory = sprintf("%s%s%s", outputDIR, .Platform$file.sep, "TN")

#
#	Z_TN, miRNAs
#

#	Generate the intersection of mirnas across all the samples
mirnas = c()
for (sample in samples) {
    
    #	Open the miRNASeq expression file
    row = subset( manifest.TN, Platform.Type == "miRNASeq" & Sample == sample )
    mirbase20.file = sprintf("%s/%s/%s__%s/Level_%s/%s", 
                             getwd(),
                             row$Platform.Type,
                             row$Center,
                             row$Platform,
                             row$Level,
                             row$File.Name)
    mr = (read.delim(mirbase20.file, stringsAsFactor=FALSE))$miRNA_ID
    
    #	Intersect miRNAs in this file with miRNAs so far collected	
    if (length(mirnas) == 0) {
        mirnas = mr
    } else {
        mirnas = intersect(mirnas, mr);	
    }
    
    # cat(sprintf("SAMPLE=%s, |mirnas|=%d\n", sample, length(mirnas)))
    
}

#	Now, we can safely construct the matrix
Z_TN = matrix(0.0, nrow=length(mirnas), ncol=length(samples))
rownames(Z_TN) = mirnas
colnames(Z_TN) = samples

#	For each sample, we open the mirbase20 expression file again and take **only** expression values for miRNAs in mirnas vector
for (sample in samples) {
    
    #	Open the miRNASeq expression file
    row = subset( manifest.TN, Platform.Type == "miRNASeq" & Sample == sample )
    mirbase20.file = sprintf("%s/%s/%s__%s/Level_%s/%s", 
                             getwd(),
                             row$Platform.Type,
                             row$Center,
                             row$Platform,
                             row$Level,
                             row$File.Name)
    
    #	Read the entire miRNA expression file
    mr = (read.delim(mirbase20.file, stringsAsFactor=FALSE))[, c("miRNA_ID", "reads_per_million_miRNA_mapped")]
    
    #	Label rows with miRNA_IDs, needed for the subsequent step
    rownames(mr) = mr$miRNA_ID
    
    #	In a single step:	select only counts relative to the common set of mirnas obtained above 
    Z_TN[, sample] = mr[mirnas, ]$reads_per_million_miRNA_mapped
    
}

#   Now splits Z_TN into Z_T, Z_N
write.table(Z_TN[, samples_T], 
             paste(mat.directory, "Z_T", sep=.Platform$file.sep),
             quote=FALSE, sep='\t')

write.table(Z_TN[, samples_N], 
            paste(mat.directory, "Z_N", sep=.Platform$file.sep),
            quote=FALSE, sep='\t')


#
#	X_TN, mRNAs
#

#	Generate the intersection of mrnas across all the samples
mrnas = c()
for (sample in samples) {
    
    #	Open the RNASeqV2 expression file
    row = subset( manifest.TN, Platform.Type == "RNASeqV2" & Sample == sample )
    rsem.file = sprintf("%s/%s/%s__%s/Level_%s/%s", 
                        getwd(),
                        row$Platform.Type,
                        row$Center,
                        row$Platform,
                        row$Level,
                        row$File.Name)
    mr = (read.delim(rsem.file, stringsAsFactor=FALSE))$gene_id				
    mr = mr[ !grepl('\\?', mr) ]									# filter out unknown ? gene entries
    mr = sapply(strsplit(as.character(mr), "\\|"), "[[", 2)		# take only the gene Entrez id (unique)
    
    #	Intersect mRNAs in this file with mRNAs so far collected	
    if (length(mrnas) == 0) {
        mrnas = mr
    } else {
        mrnas = intersect(mrnas, mr);	
    }
    
    #cat(sprintf("SAMPLE=%s, |mrnas|=%d\n", sample, length(mrnas)))
    
}

#	Now, we can safely construct the matrix
X_TN = matrix(0.0, nrow=length(mrnas), ncol=length(samples))
rownames(X_TN) = mrnas
colnames(X_TN) = samples

#	For each sample, we open the rsem.genes expression file again and take **only** expression values for miRNAs in mirnas
for (sample in samples) {
    
    #	Open the miRNASeq expression file
    row = subset( manifest.TN, Platform.Type == "RNASeqV2" & Sample == sample )
    rsem.file = sprintf("%s/%s/%s__%s/Level_%s/%s", 
                        getwd(),
                        row$Platform.Type,
                        row$Center,
                        row$Platform,
                        row$Level,
                        row$File.Name)
    
    #	Read the entire mRNA expression file
    mr = (read.delim(rsem.file, stringsAsFactor=FALSE))[, c("gene_id", "scaled_estimate")]
    mr = mr[ !grepl('\\?', mr$gene_id), ]		
    mr$gene_id = sapply(strsplit(as.character(mr$gene_id), "\\|"), "[[", 2)
    
    #	Label rows with mRNA_IDs, needed for the subsequent step
    rownames(mr) = mr$gene_id
    
    #	In a single step:	select only counts relative to the common set of mrnas obtained above 
    X_TN[, sample] = mr[mrnas, ]$scaled_estimate * 10^6
    
}

write.table(X_TN[, samples_T], 
             paste(mat.directory, "X_T", sep=.Platform$file.sep),
             quote=FALSE, sep='\t')

write.table(X_TN[, samples_N], 
            paste(mat.directory, "X_N", sep=.Platform$file.sep),
            quote=FALSE, sep='\t')

write.table(data.frame(mRNA_EntrezId=mrnas), 
             paste(mat.directory, "TN_mrnas.txt", sep=.Platform$file.sep), 
             quote=FALSE, row.names=FALSE)	

write.table(data.frame(miRNA_ID=mirnas), 
            paste(mat.directory, "TN_mirnas.txt", sep=.Platform$file.sep),
            quote=FALSE, row.names=FALSE)	

write.table(data.frame(SampleID=samples_T), 
             paste(mat.directory, "T_samples.txt", sep=.Platform$file.sep), 
             quote=FALSE, row.names=FALSE)	

write.table(data.frame(SampleID=samples_N), 
            paste(mat.directory, "N_samples.txt", sep=.Platform$file.sep), 
            quote=FALSE, row.names=FALSE)	


cat(sprintf("DONE. number_of_paired_samples=%d, mRNAs=%d, miRNAs=%d\n", length(samples_T), nrow(X_TN), nrow(Z_TN)))
#rm(mirnas); rm(mrnas); rm(mr); rm(X_TN); rm(Z_TN);

cat("DONE.")

# X_T = X_TN[, samples_T]
# X_N = X_TN[, samples_N]
# Z_T = Z_TN[, samples_T]
# Z_N = Z_TN[, samples_N]
# 
# data_T = readExpressionData("buildTCGASets/TN/X_T", "buildTCGASets/TN/Z_T", "buildTCGASets/TN/TN_mrnas.txt", "buildTCGASets/TN/TN_mirnas.txt", "buildTCGASets/TN/T_samples.txt")
# data_N = readExpressionData("buildTCGASets/TN/X_N", "buildTCGASets/TN/Z_N", "buildTCGASets/TN/TN_mrnas.txt", "buildTCGASets/TN/TN_mirnas.txt", "buildTCGASets/TN/N_samples.txt")
