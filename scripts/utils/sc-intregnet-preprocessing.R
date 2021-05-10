#!/usr/bin/Rscript

preprocess_exp_sample <- function(gene_names, gene_format=c('ensembl_gene_id', 'hgnc_symbol'), values){
  library(biomaRt)
  hsmart <- useMart(dataset = "hsapiens_gene_ensembl", biomart = "ensembl")

  if (gene_format == "ensembl_gene_id"){
    gene_names <- sapply(gene_names, function(x){gsub("\\.[0-9]*","",x)})
  }

  mapping <- getBM(
    attributes = c('ensembl_gene_id', 'hgnc_symbol'), 
    filters = match.arg(gene_format),
    values = gene_names,
    mart = hsmart
  )

  data = data.frame(gene = gene_names, TPM = values)
  colnames(data)[1] = match.arg(gene_format)
  
  mapping = mapping[!duplicated(mapping$hgnc_symbol),]

  data = merge(data, mapping, by=match.arg(gene_format))
  return(data[c('ensembl_gene_id', 'hgnc_symbol', 'TPM')])
}