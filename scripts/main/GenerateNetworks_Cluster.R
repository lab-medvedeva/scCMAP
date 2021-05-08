#!/usr/bin/env Rscript

library(plyr)
library(igraph)

WORKDIR <- Sys.getenv("WORKDIR")
setwd(WORKDIR)

INPUT_DIR <- Sys.getenv("INPUT_DIR")
OUTPUT_DIR <- Sys.getenv("OUTPUT_DIR")

sample_name <- Sys.getenv("SAMPLE_NAME") # sample dir as well

# INPUT FILES

ppiInteractionFile <- file.path(WORKDIR, "default_files", "ppi_unmod.txt")

promoterFile <- file.path(OUTPUT_DIR, paste0(sample_name, "_Shell1_Pro_Acc_OnlyTF.bed"))
enhancerFile <- file.path(OUTPUT_DIR, paste0(sample_name, "_Shell1_Enh_Acc_OnlyTF.bed"))

# OUTPUT FILES
promoter_network_path <- file.path(OUTPUT_DIR, paste0(sample_name, "_Promoter_Network_v1.txt"))
enhancer_network_path <- file.path(OUTPUT_DIR, paste0(sample_name, "_Enhancer_Network_v1.txt"))
promoter_regulators_path <- file.path(OUTPUT_DIR, paste0(sample_name, "_PromRegulators.txt"))
enhancer_regulators_path <- file.path(OUTPUT_DIR, paste0(sample_name, "_EnhRegulators.txt"))
networkFile_path <- file.path(OUTPUT_DIR, paste0(sample_name, "_Network_v1.txt"))

# TEMPORAL FILES
data_new_path <- file.path(OUTPUT_DIR, paste0(sample_name, "_data_new.bed"))
data_enh_new_path <- file.path(OUTPUT_DIR, paste0(sample_name, "_enh_data_new.bed"))

# RUN

#####
# See below, there are some hard-coded paths and filenames that are only used for temporary storage. Change them once
# and you can reuse them (unless you run multiple R session with this script at the same time)
####

writeLogicAsPolynomial <- function(df, ppi, directed, mode) {
  targetGenes <- unique(df$V9)
  logicRegulation_df <- data.frame(Gene = targetGenes, Logic = rep("", length(targetGenes)), stringsAsFactors = FALSE)
  for (i in 1:length(targetGenes)) {
    logic <- c()
    ints <- df[which(df$V9 == targetGenes[i]), c(4, 16)]
    graph <- graph_from_data_frame(ints, directed = directed)
    SCC <- clusters(graph, mode = mode)
    #Get all clusters and obtain gene level clusters
    clusts <- SCC$membership
    names(clusts) <- gsub("_[0-9]+", "", names(clusts))
    clusts <- cbind.data.frame(names(clusts), clusts, stringsAsFactors = FALSE)
    clusts <- clusts[!duplicated(clusts),]
    #Treat each cluster separately. Iterate over them
    for (c in 1:length(SCC$csize)) {
      cur_clust <- clusts[clusts[, 2] == c, 1]
      if (length(cur_clust) == 1) {
        #A single gene in the cluster does not make a difference in terms of logic. Include it as a singleton
        logic <- c(logic, cur_clust[1])
      } else {
        #There is more than one gene in the cluster. Subset ppi network
        ppi_subs <- ppi[ppi$symbol1 %in% cur_clust & ppi$symbol2 %in% cur_clust,]
        #Identify connected components
        ppi_graph <- graph_from_data_frame(ppi_subs[, 1:2], directed = FALSE)
        ppi_graph_comp <- clusters(ppi_graph)
        #If more than one gene participates in connected component connect with and otherwise include as singleton
        if (length(ppi_graph_comp$csize) > 0) {
          for (ppi_c in 1:length(ppi_graph_comp$csize)) {
            if (ppi_graph_comp$csize[ppi_c] == 1) {
              logic <- c(logic, names(ppi_graph_comp$membership)[which(ppi_graph_comp$membership == ppi_c)])
            } else {
              logic <- c(logic, paste0("( ", paste(names(ppi_graph_comp$membership)[which(ppi_graph_comp$membership == ppi_c)], collapse = " * "), " )"))
            }
          }
        }

        #Now remove all genes in componentes/with ppi's and add the remaining genes as singletons
        remaining <- setdiff(cur_clust, names(ppi_graph_comp$membership))
        if (length(remaining) > 0) {
          for (j in 1:length(remaining)) {
            logic <- c(logic, remaining[j])
          }
        }
      }
    }
    #Add everything that is remaining in 'data' for 'targetgene' as seperate entities. Should not happen but let's be sure
    regMotifs <- data[which(data$V9 == targetGenes[i]), 4]
    regMotifs <- setdiff(regMotifs, names(SCC$membership))
    regMotifs <- gsub("_[0-9]+", "", regMotifs)
    regMotifs <- regMotifs[!duplicated(regMotifs)]
    if (length(regMotifs) > 0) {
      logic <- c(logic, regMotifs)
    }

    #Collapse logic and add to logicRegulation_df
    logic <- unique(logic)
    if (length(logic) > 1) {
      logic <- sapply(logic, function(x) { paste0("( 1- ", x, " )") })
      logicRegulation_df[which(logicRegulation_df$Gene == targetGenes[i]), 2] <- paste0("(1-( ", paste(unique(logic), collapse = " * "), " ))")
    } else {
      logicRegulation_df[which(logicRegulation_df$Gene == targetGenes[i]), 2] <- logic[1]
    }

  }
  return(logicRegulation_df)
}

writeLogic <- function(df, ppi, directed, mode) {
  targetGenes <- unique(df$V9)
  logicRegulation_df <- data.frame(Gene = targetGenes, Logic = rep("", length(targetGenes)), stringsAsFactors = FALSE)
  for (i in 1:length(targetGenes)) {
    logic <- c()
    ints <- df[which(df$V9 == targetGenes[i]), c(4, 16)]
    graph <- graph_from_data_frame(ints, directed = directed)
    SCC <- clusters(graph, mode = mode)
    #Get all clusters and obtain gene level clusters
    clusts <- SCC$membership
    names(clusts) <- gsub("_[0-9]+", "", names(clusts))
    clusts <- cbind.data.frame(names(clusts), clusts, stringsAsFactors = FALSE)
    clusts <- clusts[!duplicated(clusts),]
    #Treat each cluster separately. Iterate over them
    for (c in 1:length(SCC$csize)) {
      cur_clust <- clusts[clusts[, 2] == c, 1]
      if (length(cur_clust) == 1) {
        #A single gene in the cluster does not make a difference in terms of logic. Include it as a singleton
        logic <- c(logic, cur_clust[1])
      } else {
        #There is more than one gene in the cluster. Subset ppi network
        ppi_subs <- ppi[ppi$symbol1 %in% cur_clust & ppi$symbol2 %in% cur_clust,]
        #Identify connected components
        ppi_graph <- graph_from_data_frame(ppi_subs[, 1:2], directed = FALSE)
        ppi_graph_comp <- clusters(ppi_graph)
        #If more than one gene participates in connected component connect with and otherwise include as singleton
        if (length(ppi_graph_comp$csize) > 0) {
          for (ppi_c in 1:length(ppi_graph_comp$csize)) {
            if (ppi_graph_comp$csize[ppi_c] == 1) {
              logic <- c(logic, names(ppi_graph_comp$membership)[which(ppi_graph_comp$membership == ppi_c)])
            } else {
              logic <- c(logic, paste0("( ", paste(names(ppi_graph_comp$membership)[which(ppi_graph_comp$membership == ppi_c)], collapse = " & "), " )"))
            }
          }
        }

        #Now remove all genes in componentes/with ppi's and add the remaining genes as singletons
        remaining <- setdiff(cur_clust, names(ppi_graph_comp$membership))
        if (length(remaining) > 0) {
          for (j in 1:length(remaining)) {
            logic <- c(logic, remaining[j])
          }
        }
      }
    }
    #Add everything that is remaining in 'data' for 'targetgene' as seperate entities. Should not happen but let's be sure
    regMotifs <- data[which(data$V9 == targetGenes[i]), 4]
    regMotifs <- setdiff(regMotifs, names(SCC$membership))
    regMotifs <- gsub("_[0-9]+", "", regMotifs)
    regMotifs <- regMotifs[!duplicated(regMotifs)]
    if (length(regMotifs) > 0) {
      logic <- c(logic, regMotifs)
    }

    #Collapse logic and add to logicRegulation_df
    logicRegulation_df[which(logicRegulation_df$Gene == targetGenes[i]), 2] <- paste(unique(logic), collapse = " | ")

  }
  return(logicRegulation_df)
}

# load promoter data
data <- read.table(promoterFile, header = FALSE, sep = '\t', stringsAsFactors = FALSE)
# data <- data[which(data$V4 %in% unique(data$V9)),]
# data <- data[order(data$V6,data$V7,data$V8),]

# load enhancer data
data_enh <- read.table(enhancerFile, header = FALSE, sep = '\t', stringsAsFactors = FALSE)
# data_enh <- data_enh[data_enh$V4 %in% unique(data_enh$V12),]
# data_enh <- data_enh[which(data_enh$V4 %in% unique(data$V9)),]
# data_enh <- data_enh[which(data_enh$V12 %in% unique(data$V9)),]
# data_enh <- data_enh[data_enh$V4 %in% unique(data_enh$V12),]
# data_enh <- data_enh[order(data_enh$V6,data_enh$V7,data_enh$V8),]

while (length(unique(data_enh$V12)) != length(unique(data$V9))) {
  data <- data[which(data$V4 %in% unique(data$V9)),]
  data_enh <- data_enh[which(data_enh$V12 %in% unique(data$V9)),]
  data_enh <- data_enh[which(data_enh$V4 %in% unique(data_enh$V12)),]
  data <- data[which(data$V9 %in% unique(data_enh$V12)),]
  data <- data[which(data$V4 %in% unique(data$V9)),]
}
dim(data)
dim(data_enh)

data <- data[order(data$V6, data$V7, data$V8),]
data_enh <- data_enh[order(data_enh$V6, data_enh$V7, data_enh$V8),]

if (!exists("ppi_unmod")) {
  ppi_unmod <- read.table(ppiInteractionFile, header = TRUE, stringsAsFactors = FALSE, sep = '\t')
}


uniqueRegulators <- data[, 1:4]
uniqueRegulators <- uniqueRegulators[!duplicated(uniqueRegulators),]
uniqueRegulators <- cbind.data.frame(uniqueRegulators, seq(1, nrow(uniqueRegulators)))
colnames(uniqueRegulators)[5] <- "V5"
uniqueRegulators$V5 <- paste0(uniqueRegulators$V4, "_", uniqueRegulators$V5)

data <- join(data, uniqueRegulators, by = c("V1", "V2", "V3", "V4"))
data$V4 <- data[, 13]
data <- data[, -13]
write.table(data, data_new_path, append = FALSE, quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)

system(paste0("intersectBed -loj -r -f 0.622 -a ", data_new_path, " -b ", data_new_path, " > relation.txt"))

relations <- read.table("relation.txt", sep = '\t', header = FALSE, stringsAsFactors = FALSE)
relations$V25 <- abs(relations$V5 - relations$V17)
#relations <- relations[which(relations$V4 != relations$V16 & relations$V9 == relations$V21),]
#relations <- relations[which(relations$V25 <= 50),]
targetGenes <- unique(relations$V9)
print("promoter file processed")
print(targetGenes)

#Do the same for enhancers
# data_enh <- read.table(enhancerFile, header = FALSE, sep = '\t', stringsAsFactors = FALSE)
# data_enh <- data_enh[which(data_enh$V4 %in% unique(data$V9)),]
# data_enh <- data_enh[which(data_enh$V12 %in% unique(data$V9)),]
# data_enh <- data_enh[order(data_enh$V6,data_enh$V7,data_enh$V8),]

uniqueRegulators_enh <- data_enh[, 1:4]
uniqueRegulators_enh <- uniqueRegulators_enh[!duplicated(uniqueRegulators_enh),]
uniqueRegulators_enh <- cbind.data.frame(uniqueRegulators_enh, seq(1, nrow(uniqueRegulators_enh)))
colnames(uniqueRegulators_enh)[5] <- "V5"
uniqueRegulators_enh$V5 <- paste0(uniqueRegulators_enh$V4, "_", uniqueRegulators_enh$V5)

data_enh <- join(data_enh, uniqueRegulators_enh, by = c("V1", "V2", "V3", "V4"))
data_enh$V4 <- data_enh[, 15]
data_enh <- data_enh[, -15]
write.table(data_enh, data_enh_new_path, append = FALSE, quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)

system(paste0("intersectBed -loj -r -f 0.622 -a ", data_enh_new_path, " -b ", data_enh_new_path, " > relation_enh.txt"))

relations_enh <- read.table("relation_enh.txt", sep = '\t', header = FALSE, stringsAsFactors = FALSE)
targetGenes <- unique(relations_enh$V12)
relations_enh <- relations_enh[, - c(9, 10, 11, 23, 25)]
relations_enh <- cbind.data.frame(relations_enh[, 1:9], data.frame(A = rep(1, nrow(relations_enh))), relations_enh[, 10:21], data.frame(B = rep(1, nrow(relations_enh))), relations_enh[, 22:23], stringsAsFactors = FALSE)
relations_enh$V25 <- abs(relations_enh$V5 - relations_enh$V19)
# relations_enh <- relations_enh[which(relations_enh$V4 != relations_enh$V18 & relations_enh$V12 == relations_enh$V26),]
# relations_enh <- relations_enh[which(relations_enh$V25 <= 50),]
colnames(relations_enh) <- colnames(relations)
print("enhancer file processed")

# Do Network evolvement.
mode <- "weak"
directed <- FALSE
allEnhancersReq <- FALSE
ppi <- ppi_unmod
logicRegulation_df <- writeLogicAsPolynomial(relations, ppi, directed = directed, mode = mode)
ppi <- ppi_unmod

if (allEnhancersReq) {
  print("Done")
  uniqEnhancers <- unique(relations_enh$V21)
  tmp_df <- data.frame(Gene = vector(mode = "character", length = length(uniqEnhancers)), Logic = vector(mode = "character", length = length(uniqEnhancers)), stringsAsFactors = FALSE)
  for (i in 1:length(uniqEnhancers)) {
    tmp_df[i,] <- writeLogicAsPolynomial(relations_enh[relations_enh$V21 == uniqEnhancers[i],], ppi, directed = directed, mode = mode)
  }
  uniqGenes <- unique(tmp_df$Gene)
  logicRegulation_enh_df <- data.frame(Gene = vector(mode = "character", length = length(uniqGenes)), Logic = vector(mode = "character", length = length(uniqGenes)), stringsAsFactors = FALSE)
  for (i in 1:length(uniqGenes)) {
    logicRegulation_enh_df[i, 1] <- uniqGenes[i]
    logicRegulation_enh_df[i, 2] <- paste(tmp_df[tmp_df$Gene == uniqGenes[i], 2], collapse = " * ")
  }
} else {
  logicRegulation_enh_df <- writeLogicAsPolynomial(relations_enh, ppi, directed = directed, mode = mode)
}

colnames(logicRegulation_enh_df) <- c("Gene", "EnhLogic")
commonGenes <- intersect(logicRegulation_df$Gene, logicRegulation_enh_df$Gene)
uniqueGenes <- union(setdiff(logicRegulation_df$Gene, commonGenes), setdiff(logicRegulation_enh_df$Gene, commonGenes))
logicRegulation_df <- logicRegulation_df[which(logicRegulation_df$Gene %in% commonGenes),]
logicRegulation_enh_df <- logicRegulation_enh_df[which(logicRegulation_enh_df$Gene %in% commonGenes),]

if (length(uniqueGenes) > 0) {
  for (i in 1:length(uniqueGenes)) {
    logicRegulation_df$Logic <- apply(logicRegulation_df, 1, function(x) { gsub(paste0("\\( 1- ", uniqueGenes[i], " \\) \\* "), "", x[2]) })
    logicRegulation_df$Logic <- apply(logicRegulation_df, 1, function(x) { gsub(paste0(" \\* \\( 1- ", uniqueGenes[i], " \\)"), "", x[2]) })
    logicRegulation_df$Logic <- apply(logicRegulation_df, 1, function(x) { gsub(paste0(uniqueGenes[i], " \\* "), "", x[2]) })
    logicRegulation_df$Logic <- apply(logicRegulation_df, 1, function(x) { gsub(paste0(" \\* ", uniqueGenes[i]), "", x[2]) })

    logicRegulation_enh_df$EnhLogic <- apply(logicRegulation_enh_df, 1, function(x) { gsub(paste0("\\( 1- ", uniqueGenes[i], " \\) \\* "), "", x[2]) })
    logicRegulation_enh_df$EnhLogic <- apply(logicRegulation_enh_df, 1, function(x) { gsub(paste0(" \\* \\( 1- ", uniqueGenes[i], " \\)"), "", x[2]) })
    logicRegulation_enh_df$EnhLogic <- apply(logicRegulation_enh_df, 1, function(x) { gsub(paste0(uniqueGenes[i], " \\* "), "", x[2]) })
    logicRegulation_enh_df$EnhLogic <- apply(logicRegulation_enh_df, 1, function(x) { gsub(paste0(" \\* ", uniqueGenes[i]), "", x[2]) })
  }
}

logics <- join(logicRegulation_df, logicRegulation_enh_df, by = c("Gene"))
logics$EnhLogic[which(is.na(logics$EnhLogic))] <- "FALSE"
logics$Logic[which(is.na(logics$Logic))] <- "FALSE"

write.table(logics[, c(1, 2)], promoter_network_path, sep = ',', quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(logics[, c(1, 3)], enhancer_network_path, sep = ',', quote = FALSE, row.names = FALSE, col.names = TRUE)
logics$Logic <- apply(logics, 1, function(x) { paste0("( ", x[2], " ) * ( ", x[3], " )") })
logics <- logics[, -3]
targetGenes <- logics$Gene
colnames(logics) <- c("Gene", "Logic")
write.table(logics, networkFile_path, sep = ',', quote = FALSE, row.names = FALSE, col.names = TRUE)

tmp <- read.table(promoter_network_path, header = TRUE, sep = ',', stringsAsFactors = FALSE)
tmp$Logic <- apply(tmp, 1, function(x) { gsub("1-", "", x[2]) })
tmp$Logic <- apply(tmp, 1, function(x) { gsub("\\*", "", x[2]) })
tmp$Logic <- apply(tmp, 1, function(x) { gsub("\\(", "", x[2]) })
tmp$Logic <- apply(tmp, 1, function(x) { gsub("\\)", "", x[2]) })
tmp$Logic <- apply(tmp, 1, function(x) { paste(unique(strsplit(x[2], "\\s+")[[1]]), collapse = ",") })
tmp$Logic <- apply(tmp, 1, function(x) { if (startsWith(x[2], ",")) { substring(x[2], 2) } else { x[2] }})
write.table(tmp, promoter_regulators_path, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")

tmp <- read.table(enhancer_network_path, header = TRUE, sep = ',', stringsAsFactors = FALSE)
tmp$EnhLogic <- apply(tmp, 1, function(x) { gsub("1-", "", x[2]) })
tmp$EnhLogic <- apply(tmp, 1, function(x) { gsub("\\*", "", x[2]) })
tmp$EnhLogic <- apply(tmp, 1, function(x) { gsub("\\(", "", x[2]) })
tmp$EnhLogic <- apply(tmp, 1, function(x) { gsub("\\)", "", x[2]) })
tmp$EnhLogic <- apply(tmp, 1, function(x) { paste(unique(strsplit(x[2], "\\s+")[[1]]), collapse = ",") })
tmp$EnhLogic <- apply(tmp, 1, function(x) { if (startsWith(x[2], ",")) { substring(x[2], 2) } else { x[2] }})
write.table(tmp, enhancer_regulators_path, quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")


if (file.exists("relation_enh.txt")) file.remove("relation_enh.txt")
if (file.exists("relation.txt")) file.remove("relation.txt")
if (file.exists(data_new_path)) file.remove(data_new_path)
if (file.exists(data_enh_new_path)) file.remove(data_enh_new_path)


