#!/usr/bin/env Rscript

library(parallel)

WORKDIR <- Sys.getenv("WORKDIR")
MATLAB_path <- file.path(WORKDIR, "scripts", "main", "Matlab")
setwd(WORKDIR)

sample_name <- Sys.getenv("SAMPLE_NAME") # sample dir as well
sample_geo <- "GSM000000"

INPUT_DIR <- Sys.getenv("INPUT_DIR")
OUTPUT_DIR <- Sys.getenv("OUTPUT_DIR")

## INPUT

target_exp_path <- Sys.getenv("target_exp_path")

## LOAD DEFAULT FILES
EnsemblToTF_path <- Sys.getenv("EnsemblToTF_path")
EnsemblToTF <- read.table(EnsemblToTF_path, stringsAsFactors = FALSE, header = TRUE, sep = "\t")

# extract TFs sub marix fro JSD
Ens_TF_path = Sys.getenv("Ens_TF_path")
Ens_TF = read.table(Ens_TF_path, sep = "\t", stringsAsFactors = F) # V1 is symbol V2 is EnsID
homo_sapiens.TFs = as.character(Ens_TF[, 2])

description_file_TPM_clust_path = Sys.getenv("description_file_TPM_clust_path")
load(description_file_TPM_clust_path)
description_TPM <- description_file_TPM_clust

# load background matrix of TPMs, we need this matrix to compute correlation between new sample (we want to add) and the background.
Recount_TPM_GSM_Bg_path = Sys.getenv("Recount_TPM_GSM_Bg_path")
load(Recount_TPM_GSM_Bg_path) # 57992  8952 (dim of bg data).
dataset_Recount_TPM = Recount_TPM_GSM_Bg
rm(Recount_TPM_GSM_Bg)
raw_TPM_data = TRUE

# load raw TPMs for TFs only (a subset of original background matrix)
dataset_raw_TPM_TFs_path = Sys.getenv("dataset_raw_TPM_TFs_path")
load(dataset_raw_TPM_TFs_path)

cor_TPM_path = Sys.getenv("cor_TPM_path")
load(cor_TPM_path)

JSD_ranked_TFs_Cor_path = Sys.getenv("JSD_ranked_TFs_Cor_path")
jsd <- read.table(JSD_ranked_TFs_Cor_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

## OUTPUT

data_Bool_path = file.path(OUTPUT_DIR, paste0(sample_name, "_Bool.txt"))
data_JSD_path = file.path(OUTPUT_DIR, paste0("JSD_", sample_name, "_Query.tsv"))
data_JSD_TF_path = file.path(OUTPUT_DIR, paste0("JSD_", sample_name, "_Query_TF.tsv"))
data_core_TF_path = file.path(OUTPUT_DIR, paste0(sample_name, "_TF.tsv"))
cor_query_TMP_path = file.path(OUTPUT_DIR, paste0("cor_TPM_", sample_name, "_query.RData"))
inclus_list_path = file.path(OUTPUT_DIR, paste0("inclusionList_", sample_name, "_0.75_v2.RData"))
inclus_lists_path = file.path(OUTPUT_DIR, paste0("inclusionLists_", sample_name, "_0.75_v2.RData"))
finalTFList_path = file.path(OUTPUT_DIR, paste0(sample_name, "_mRNA_TFs.txt"))

# RUN

print("Starting booleanization")
command <- paste0('matlab -nodisplay -nosplash -r "cd ', MATLAB_path, '; BooleanizeFromFile(\'', target_exp_path, '\', \'', data_Bool_path, '\');quit"')
system(command)
print("Booleanization done, starting JSD")

# add sample name (e.g, "PC3_WT1") and geo accession (e.g. "GSM1563053") to the existing metadata file
description_TPM[nrow(description_TPM) + 1,] = list(NA, sample_name, sample_geo, NA)

# load query matrix
# To add new query samples, get their TPM or process relevant SRP from recount and merge with this df to have same rows of genes.
data <- read.table(target_exp_path, stringsAsFactors = FALSE, header = TRUE)
data <- data[, c("ensembl_gene_id", "TPM")]
colnames(data) <- c("ensembl_gene_id", sample_geo)

common_genes = intersect(data$ensembl_gene_id, rownames(dataset_Recount_TPM))
Recount_TPM_GSM_Query <- data[which(data$ensembl_gene_id %in% common_genes),]
dataset_Recount_TPM <- dataset_Recount_TPM[which(rownames(dataset_Recount_TPM) %in% common_genes),]

Recount_TPM_GSM_Query <- Recount_TPM_GSM_Query[order(Recount_TPM_GSM_Query$ensembl_gene_id),]
rownames(Recount_TPM_GSM_Query) <- Recount_TPM_GSM_Query$ensembl_gene_id
Recount_TPM_GSM_Query$ensembl_gene_id <- NULL

# NOTE: Query (Recount_TPM_GSM_Query) and background data (dataset_Recount_TPM) should have identical row.names

# extract TFs sub marix from the Query TPM matrix
common_tfs = intersect(homo_sapiens.TFs, rownames(Recount_TPM_GSM_Query))
Recount_TPM_GSM_Query$ensembl_gene_id = rownames(Recount_TPM_GSM_Query)
dataset_raw_TPM_TFs$ensembl_gene_id = rownames(dataset_raw_TPM_TFs)

Recount_TPM_GSM_Query_TF <- Recount_TPM_GSM_Query[Recount_TPM_GSM_Query$ensembl_gene_id %in% common_tfs,]
dataset_raw_TPM_TFs <- dataset_raw_TPM_TFs[dataset_raw_TPM_TFs$ensembl_gene_id %in% Recount_TPM_GSM_Query_TF$ensembl_gene_id,]

rownames(Recount_TPM_GSM_Query) = Recount_TPM_GSM_Query$ensembl_gene_id
rownames(Recount_TPM_GSM_Query_TF) = Recount_TPM_GSM_Query_TF$ensembl_gene_id
rownames(dataset_raw_TPM_TFs) = dataset_raw_TPM_TFs$ensembl_gene_id

Recount_TPM_GSM_Query$ensembl_gene_id <- NULL
Recount_TPM_GSM_Query_TF$ensembl_gene_id <- NULL
dataset_raw_TPM_TFs$ensembl_gene_id <- NULL


final_output <- c()
cor_TPM_query <- data.frame()
for (ind_ct in 1:ncol(Recount_TPM_GSM_Query)) {
  selected_sample <- cor(dataset_Recount_TPM, Recount_TPM_GSM_Query[, ind_ct]) # correlation of query against bg matrix
  colnames(selected_sample) <- colnames(Recount_TPM_GSM_Query[ind_ct])

  ifelse(ind_ct == 1, cor_TPM_query <- selected_sample, cor_TPM_query <- cbind(cor_TPM_query, selected_sample)) # save query's cor columns

  selected_sample[selected_sample < 0.75] <- NA # threshold for similar sample is cor >= .75
  print(colnames(selected_sample))
  selected_sample <- na.omit(selected_sample) # remove samples which are NA (cor < 0.75)
  excluded_clustersIds <- unique(rownames(selected_sample))

  excluded_clusters_col <- which(colnames(dataset_raw_TPM_TFs) %in% excluded_clustersIds) # get column number (indices) of samples to be excluded
  if (length(excluded_clusters_col) == 0) {
    TF.dataset <- dataset_raw_TPM_TFs
  } else {
    TF.dataset <- dataset_raw_TPM_TFs[, - excluded_clusters_col] # remove the samples having cor > 0.75 from background
  }
  TF.dataset <- cbind(TF.dataset, Recount_TPM_GSM_Query_TF[, ind_ct]) # bind the query sample with background df
  colnames(TF.dataset)[ncol(TF.dataset)] <- colnames(Recount_TPM_GSM_Query[ind_ct]) # rename the newly added column
  considered_TFs <- rownames(TF.dataset)
  cat("dimensions of background data. \n")
  print(dim(TF.dataset))

  selected_samples <- description_TPM[description_TPM[, "geo_accession"] == colnames(selected_sample), c("mix.cell", "geo_accession")]
  selected_samples.colnames <- colnames(selected_samples)
  selected_samples <- cbind(selected_samples, which(colnames(TF.dataset) %in% selected_samples[, 2])) # get column numbers of selected samples
  colnames(selected_samples) <- c(selected_samples.colnames, "Column_number")
  num_selected_samples <- dim(selected_samples)[1]

  # cat("JSD computation is started.\n") 
  ## make cluster for parallel computations
  cl <- makePSOCKcluster(detectCores())
  clusterExport(cl, "considered_TFs")
  clusterExport(cl, "TF.dataset")
  clusterExport(cl, "raw_TPM_data")

  normalized_expr_for_JSD <- parSapply(cl, considered_TFs, function(tf) {
    distr_vector <- TF.dataset[tf,]
    distr_vector <- distr_vector / sum(distr_vector)
  })

  clusterExport(cl, "selected_samples")
  clusterExport(cl, "normalized_expr_for_JSD")

  JSD_ranked_TFs <- parSapply(cl, 1:num_selected_samples, function(i) {
    target_mix.cell_col <- selected_samples[i, 3]
    exclude_from_background <- selected_samples[-i, 3]

    ideal_distr_vector <- rep(0, dim(TF.dataset)[2])
    ideal_distr_vector[target_mix.cell_col] <- 1

    JSD_based_expr_score <- sapply(considered_TFs, function(tf) {
      distr_vector <- as.numeric(normalized_expr_for_JSD[, tf])
      # Exclude target cluster from background.
      distr_vector[exclude_from_background] <- 0
      distr_vector <- distr_vector / sum(distr_vector)

      mean_distr_vector <- 0.5 * (ideal_distr_vector + distr_vector)
      KL1 <- log2(1 / mean_distr_vector[target_mix.cell_col])
      KL2 <- sum(sapply(1:length(distr_vector), function(i) {
        if (distr_vector[i] == 0)
          0
        else
          distr_vector[i] * log2(distr_vector[i] / mean_distr_vector[i])
      }))

      0.5 * (KL1 + KL2)
    })

    sJSD_expr_score <- sort(JSD_based_expr_score)

    labels(sJSD_expr_score)
  })
  stopCluster(cl)

  #colnames(JSD_ranked_TFs) <- selected_samples[,1] # To assign sample name as the column name
  colnames(JSD_ranked_TFs) <- selected_samples[, 2] # To assign GSM id as the column name

  ifelse(ind_ct == 1, final_output <- JSD_ranked_TFs, final_output <- cbind(final_output, JSD_ranked_TFs))
}

write.table(final_output, data_JSD_path, sep = "\t") # JSD for query samples with Ensembl TF names
cat("JSD computation is done and file saved.\n")

# save correlation columns of every individual query sample in a RData 
cor_TPM_query <- round(cor_TPM_query, 2)
save(cor_TPM_query, file = cor_query_TMP_path)

# convert the ensembl ids to HGNC symbols
JSD <- read.table(data_JSD_path, sep = "\t", header = T, stringsAsFactors = F)
lookUp1 <- setNames(as.character(Ens_TF$V1), Ens_TF$V2)
res <- data.frame(lapply(JSD, function(i) lookUp1[i]))

# save core TFs and JSD_TFs list
write.table(res, file = data_JSD_TF_path, sep = "\t", row.names = F, quote = F)
write.table(res[1:10, 1], file = data_core_TF_path, sep = "\t", row.names = F, quote = F, col.names = FALSE)

print("JSD Done, starting inclusionList")

tfToRemove <- jsd[which(!(jsd[, 1] %in% res[, 1])), 1]
jsd <- as.data.frame(apply(jsd, 2, function(x) { x[!(x %in% tfToRemove)] }))
jsd <- cbind.data.frame(jsd, res, stringsAsFactors = FALSE)
TFs <- jsd[, 1]
TFs <- TFs[order(TFs)]
TFs <- TFs[!duplicated(TFs)]
TFs <- TFs[!is.na(TFs)]

jsdRanks <- apply(jsd, 2, function(x) { match(TFs, x) })

ztest <- function(x, m, s) {
  (x - m) / (s)
}

incl_list <- vector("list", ncol(jsdRanks))
i = ncol(jsdRanks)
incl_samples <- rownames(cor_TPM_query)[which(cor_TPM_query <= 0.75)]
incl_samples <- which(colnames(jsd) %in% incl_samples)
z <- ztest(jsdRanks[, i], sapply(seq(1, length(TFs)), function(x) { mean(jsdRanks[x, incl_samples]) }), sapply(seq(1, length(TFs)), function(x) { sd(jsdRanks[x, incl_samples]) }))
incl <- TFs[which(z <= -1.5)]
incl_list[[colnames(jsd)[i]]] <- incl

save(incl_list, file = inclus_lists_path)
incl <- as.data.frame(incl)
save(incl, file = inclus_list_path)

boolExp <- read.table(data_Bool_path, header = TRUE, sep = ",", stringsAsFactors = FALSE)
boolExp <- boolExp[which(boolExp[, 2] == 1),]
finalTFList <- boolExp$output1[which(boolExp$output1 %in% incl$incl)]
write.table(finalTFList, finalTFList_path, sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
