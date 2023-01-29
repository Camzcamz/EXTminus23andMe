### Example that uses RUN_GSEM_RG_MAIN within own R session ###

# load the wrapper functions
source("./EXT_genetic_cor_func.R")

# prepare inputs
FILES <- list.files(path = "./example_data", pattern = ".sumstats.gz", full.names = TRUE)
NAMES <- gsub("./example_data/|.sumstats.gz", "", FILES)

# Run
RUN_GSEM_RG_MAIN(FILES, NAMES, cores = 3)

