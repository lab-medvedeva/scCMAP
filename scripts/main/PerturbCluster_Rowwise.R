#!/usr/bin/env Rscript

library(dplyr)
options(stringsAsFactors = FALSE)

sample_name <- Sys.getenv("SAMPLE_NAME") # sample dir as well
initial_celltype <- Sys.getenv("initial_celltype")
number_tfs_in_combo <- as.numeric(Sys.getenv("number_tfs_in_combo")) # number TFs in combo

WORKDIR <- Sys.getenv("WORKDIR")
setwd(WORKDIR)

INPUT_DIR <- Sys.getenv("INPUT_DIR")
OUTPUT_DIR <- Sys.getenv("OUTPUT_DIR")
TFS_PRED_DIR <- Sys.getenv("TFS_PRED_DIR")

## INPUT

networkFile_path <- file.path(OUTPUT_DIR, paste0(sample_name, "_Network_v1.txt"))
prism_output_file_path <- file.path(OUTPUT_DIR, paste0(sample_name, "_Prism_Output"))
Sample_Shell1_Pro_Acc_OnlyTF <- file.path(OUTPUT_DIR, paste0(sample_name, "_Shell1_Pro_Acc_OnlyTF.bed"))

H3K27ac_file <- Sys.getenv("target_H3K27ac")
# H3K4me3_file <- Sys.getenv("target_H3K4me3")

initial_celltype_expression <- Sys.getenv("initial_celltype_expression")
initial_celltype_H3K27ac_file <- Sys.getenv("initial_celltype_H3K27ac")
initial_celltype_H3K4me3_file <- Sys.getenv("initial_celltype_H3K4me3")

GeneHancer_GRCh38_path <- file.path(WORKDIR, "default_files", "Processed_GeneHancer_GRCh38.bed")

## OUTPUT

final_evals_path <- file.path(TFS_PRED_DIR, paste0(initial_celltype, "_to_", sample_name, "_FinalVals.tsv"))

## RUN

# Read Network
network <- read.table(networkFile_path, header = TRUE, sep = ",", stringsAsFactors = FALSE)

# Read prior (p-values, correspond to probability of being zero)
prior <- read.table(initial_celltype_expression, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Estimate ecdf -  empirical cumulative distribution function
ECD_func <- ecdf(prior$TPM)
prior$ecdf <- ECD_func(prior$TPM)

prior <- prior[prior$hgnc_symbol %in% network$Gene,]
prior <- prior[order(match(prior$hgnc_symbol, network$Gene)),]
tmp1 <- tempfile(tmpdir = "./")
tmp2 <- tempfile(tmpdir = "./")
tmp3 <- tempfile(tmpdir = "./")

# if initial_celltype_H3K4me3_file is provided
if (initial_celltype_H3K4me3_file != "") {
  system(paste0("cut -f6,7,8,9 ", Sample_Shell1_Pro_Acc_OnlyTF, " | sort | uniq | intersectBed -u -b ", initial_celltype_H3K4me3_file, " -a - > ", tmp1))
} else {
  system(paste0("cut -f6,7,8,9 ", Sample_Shell1_Pro_Acc_OnlyTF, " | sort | uniq > ", tmp1))
}

# if H3K27ac_file is provided
if (H3K27ac_file != "") {
  system(paste0("intersectBed -u -a ", GeneHancer_GRCh38_path, " -b ", H3K27ac_file, " > ", tmp2))
} else {
  system(paste0("cat ", GeneHancer_GRCh38_path, " > ", tmp2))
}

# if initial_celltype_H3K27ac_file is provided
if (initial_celltype_H3K27ac_file != "") {
  system(paste0("intersectBed -u -a ", GeneHancer_GRCh38_path, " -b ", initial_celltype_H3K27ac_file, " > ", tmp3))
} else {
  system(paste0("cat ", GeneHancer_GRCh38_path, " > ", tmp3))
}


Final_Prom_Active_In_Start <- read.table(tmp1, header = FALSE, stringsAsFactors = FALSE)
Final_Enhancers <- read.table(tmp2, header = FALSE, stringsAsFactors = FALSE)
Final_Enhancers <- Final_Enhancers[Final_Enhancers$V7 %in% network$Gene,]
Start_Enhancers <- read.table(tmp3, header = FALSE, stringsAsFactors = FALSE)
Start_Enhancers <- Start_Enhancers[Start_Enhancers$V7 %in% network$Gene,]
totalEnhancers <- sapply(network$Gene, function(x) { length(unique(union(Start_Enhancers$V5[Start_Enhancers$V7 == x], Final_Enhancers$V5[Final_Enhancers$V7 == x]))) })
toClose <- sapply(network$Gene, function(x) { length(unique(setdiff(Start_Enhancers$V5[Start_Enhancers$V7 == x], Final_Enhancers$V5[Final_Enhancers$V7 == x]))) })
toOpen <- sapply(network$Gene, function(x) { length(unique(setdiff(Final_Enhancers$V5[Final_Enhancers$V7 == x], Start_Enhancers$V5[Start_Enhancers$V7 == x]))) })
promToOpen <- sapply(seq(1, length(network$Gene)), function(x) { if (network$Gene[x] %in% Final_Prom_Active_In_Start$V4) { 0 } else { 1 }})
commonEnhancers <- sapply(network$Gene, function(x) { length(unique(intersect(Start_Enhancers$V5[Start_Enhancers$V7 == x], Final_Enhancers$V5[Final_Enhancers$V7 == x]))) })
system(paste0("rm -f ", tmp1))
system(paste0("rm -f ", tmp2))
system(paste0("rm -f ", tmp3))

print("Starting TFs combo perturbations...")

con = file(prism_output_file_path, "r")
print(prism_output_file_path)
perturbagens <- data.frame(V1 = network$Gene)
pp <- combn(perturbagens$V1, number_tfs_in_combo)
length(pp)
nominator <- rep(0, ncol(pp))
denominator <- rep(0, ncol(pp))
#counter <- 0
while (TRUE) {
  #counter <- counter + 1
  #print(counter)
  line = readLines(con, n = 1)
  if (length(line) == 0) {
    break
  }
  #Break it down to a vector
  line_vec <- sapply(strsplit(line, "\\s+")[[1]], as.numeric)
  line_vec <- unname(line_vec)
  #Check which perturbation fullfills the criteria
  pert_idx <- which(sapply(seq(1, ncol(pp)), function(x) { all(line_vec[which(perturbagens$V1 %in% pp[, x])] == 1) }))
  for (idx in pert_idx) {
    zeros <- which(line_vec == 0)
    ones <- setdiff(which(line_vec == 1), which(perturbagens$V1 %in% pp[, idx]))
    stateProb <- prod(1 - prior$ecdf[zeros]) * prod(prior$ecdf[ones])
    chromChangeProb <- prod(1 - commonEnhancers[zeros] / totalEnhancers[zeros]) * prod(commonEnhancers[ones] / totalEnhancers[ones])
    prob <- ifelse(is.na(tail(line_vec, 1)), 0, 1 / tail(line_vec, 1))
    nominator[idx] <- nominator[idx] + stateProb * chromChangeProb * prob
    denominator[idx] <- denominator[idx] + stateProb * chromChangeProb
  }
}

# save(nominator, denominator, pp, file = final_evals_path)
tfs_score <- nominator / denominator
pp <- t(pp)

tfs_combo_predictions <- data.frame()
for (i in 1:dim(pp)[1]) {
  tfs_combo_i <- c(tfs_score[i], nominator[i], denominator[i], as.character(pp[i,]))
  tfs_combo_predictions <- rbind(tfs_combo_predictions, tfs_combo_i)
}
colnames(tfs_combo_predictions) <- c("tfs_score", "nominator", "denominator", sprintf("TF%d", seq(1:ncol(pp))))
tfs_combo_predictions <- tfs_combo_predictions[order(tfs_combo_predictions$tfs_score, decreasing = TRUE),]
write.table(tfs_combo_predictions, file = final_evals_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

close(con)