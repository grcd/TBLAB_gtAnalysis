
#   getProteinCodingGenes.R
#       Scarica la lista completa HGNC che sono protein-coding
#       e crea un mapping entrez_id|ensembl_gene_id.

HGNC.filename = "NCBI_HGNC/HGNC_gene_with_protein_product.txt"
HGNC.reduced  = "NCBI_HGNC/hgnc_protein_coding_reduced_2.txt"

#   1.  Download 'gene with protein product' genes, http://www.genenames.org/cgi-bin/statistics .
cat("Downloading HGNC gene-with-protein-product file...")
download.file( "ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/locus_types/gene_with_protein_product.txt",
               HGNC.filename)

#   2.  Select only 'entrez_id', 'ensembl_gene_id' columns
hgnc = read.delim2(HGNC.filename, stringsAsFactors=FALSE)
hgnc = hgnc[, c("entrez_id", "ensembl_gene_id")]

#   3.  Write a reduced version to be used with TBLAB script(s)
write.table(hgnc, HGNC.reduced, quote=TRUE, row.names=FALSE, sep='\t')