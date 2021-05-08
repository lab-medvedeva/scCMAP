#!/usr/bin/env Rscript

WORKDIR <- Sys.getenv("WORKDIR")
setwd(WORKDIR)

sample_name <- Sys.getenv("SAMPLE_NAME") # sample dir as well

INPUT_DIR <- Sys.getenv("INPUT_DIR")
OUTPUT_DIR <- Sys.getenv("OUTPUT_DIR")

# INPUT

networkFile_path <- file.path(OUTPUT_DIR, paste0(sample_name, "_Network_v1.txt"))
cat("network_file:", networkFile_path, "\n")

# OUTPUT

reward_file_path <- file.path(OUTPUT_DIR, paste0(sample_name, "_Reward.txt"))
reachability_file_path <- file.path(OUTPUT_DIR, paste0(sample_name, "_Reachability.txt"))

# RUN

writeReachabilityProperty <- function(network, outfile) {
  sink(outfile)
  cat(paste0(
    "filter(avg, P=? [ (F ",
    paste(network$Gene, collapse = "+"),
    "=",
    nrow(network),
    ") ]",
    ")"
  ), sep = "\n")
  sink()
}

writeRewardProperty <- function(network, outfile) {
  sink(outfile)
  cat(paste0(
    "filter(avg, R=? [ (F ",
    paste(network$Gene, collapse = "+"),
    "=",
    nrow(network),
    ") ]",
    ")"
  ), sep = "\n")
  sink()
}

network <- read.table(networkFile_path, header = TRUE, sep = ",", stringsAsFactors = FALSE)
print("read network")
writeRewardProperty(network, reward_file_path)
print("reward written")
writeReachabilityProperty(network, reachability_file_path)
print("reachability written")
