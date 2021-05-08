#!/usr/bin/env Rscript

WORKDIR <- Sys.getenv("WORKDIR")
setwd(WORKDIR)

sample_name <- Sys.getenv("SAMPLE_NAME") # sample dir as well

INPUT_DIR <- Sys.getenv("INPUT_DIR")
OUTPUT_DIR <- Sys.getenv("OUTPUT_DIR")

# INPUT
networkFile_path <- file.path(OUTPUT_DIR, paste0(sample_name, "_Network_v1.txt"))

# OUTPUT
prism_networkFile_path <- file.path(OUTPUT_DIR, paste0(sample_name, "_Network_Prism_v1.txt"))

# RUN
cat("network_file:", networkFile_path, "\n")
cat("prism network - output file:", prism_networkFile_path, "\n")

#args <- commandArgs(trainlingOnly=TRUE)
#prism_networkFile_path <- paste0("/Volumes/LCSB_EPI_CBG/Data/Evan_Harvard/",ct,"/",ct,"_Network_Prism_v1.txt")
#network <- read.table(paste0("/Volumes/LCSB_EPI_CBG/Data/Evan_Harvard/",ct,"/",ct,"_Network_v1.txt"), header = TRUE, sep=',', stringsAsFactors = FALSE)

network <- read.table(networkFile_path, header = TRUE, sep = ",", stringsAsFactors = FALSE)
#prism_networkFile_path <- args[2]

sink(prism_networkFile_path)
cat("dtmc", sep = "\n")

for (i in 1:nrow(network)) {
  cat(paste0("module M", network[i, 1]), sep = "\n")
  cat("\n", sep = "\n")
  cat(paste0(network[i, 1], " : [0..1];"), sep = "\n")
  cat(paste0("[] (", network[i, 2], "=1) -> 1: (", network[i, 1], "'=1);"), sep = "\n")
  cat(paste0("[] (", network[i, 2], "=0) -> 1: (", network[i, 1], "'=0);"), sep = "\n")
  cat("endmodule", sep = "\n")
  cat("\n", sep = "\n")
}

cat("init", sep = "\n")
cat("true", sep = "\n")
cat("endinit", sep = "\n")

cat("rewards", sep = "\n")
cat("[] true:1;", sep = "\n")
cat("endrewards", sep = "\n")

sink()



