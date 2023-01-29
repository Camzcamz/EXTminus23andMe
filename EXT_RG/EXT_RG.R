#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = T)
file <- as.character(arguments[1])
output_prefix <- as.character(arguments[2])
ncores <- as.numeric(arguments[3])

output_file = paste0(output_prefix, "_GSEM_EXT_RG.txt")

library(data.table)
source("./EXT_genetic_cor_func.R")

INPUT <- fread(file, header=FALSE)
cat(nrow(INPUT), " traits \n")

FILES <- INPUT[[1]]
NAMES <- INPUT[[2]]
if (ncol(INPUT)==4){
  cat("Sample and Population prevalence specified \n")
  SAMPLE_PREV <- INPUT[[3]]
  POP_PREV <- INPUT[[4]]  
} else {
  SAMPLE_PREV <- NULL
  POP_PREV <- NULL
}

OUT = RUN_GSEM_RG_MAIN(FILES, NAMES, SAMPLE_PREV, POP_PREV, ncores)

cat("Finished. The results are saved to: ", output_file, "\n")
fwrite(OUT, output_file, sep="\t")