
RUN_GSEM_RG <- function(trait, trait.name, sample.prev = NA, population.prev = NA){

    # Run LDSC
    LDSCoutput <- ldsc_EXT(trait = trait, trait.name = trait.name, sample.prev = sample.prev, population.prev = population.prev)

    # Set model 
    EXTmodel <- 'EXTERNALIZING =~ ADHD+ AGE_FIRST_SEX_REV_CODED+ ALCOHOL_PROB+ EVER_CANNABIS+ EVER_SMOKER+ NUM_SEX_PARTNER+ RISK_TOL
                    EXTERNALIZING ~~ EXTERNALIZING
                    ALCOHOL_PROB ~~ EVER_SMOKER
                    EVER_CANNABIS ~~ AGE_FIRST_SEX_REV_CODED
                    EXTERNALIZING ~~ TRAIT'

    # Run GSEM 
    capture.output(output <- usermodel(LDSCoutput, estimation = "DWLS", model = gsub("TRAIT", trait.name, EXTmodel), std.lv = TRUE))

    # Prepare results
    RES <- data.table(output$results)
    RES <- RES[lhs=="EXTERNALIZING" & rhs==trait.name, .(rhs, STD_Genotype, STD_Genotype_SE, p_value)]
    names(RES) <- c("TRAIT", "RG", "SE", "P")
    cat("\n", trait.name, " done \n")
    RES
}


RUN_GSEM_RG_MAIN <- function(file_vec, trait.name_vec, sample.prev_vec=NULL, population.prev_vec=NULL, cores=NULL){

    library(GenomicSEM)
    library(data.table)
    library(pbmcapply)
    source("EXT_LDSC_GSEM.R")

    cat("Reading inputs \n")
    # READ pre-computed inputs
    load("./data/LDSC_EXT_intermediate_res.RData", envir = .GlobalEnv) 
    # READ pre-merged EXT sumstats
    .GlobalEnv$all_y <- lapply(trait.names_EXT, function(x) fread(paste0("./data/", x, ".sumstats_merged.gz")))

    if (is.null(sample.prev_vec)) sample.prev_vec <- rep(NA, length(file_vec))
    if (is.null(population.prev_vec)) population.prev_vec <- rep(NA, length(file_vec))

    RUN_GSEM_RG_PAR <- function(t){
        RUN_GSEM_RG(file_vec[t], trait.name_vec[t], sample.prev_vec[t], population.prev_vec[t])
    }
    if (is.null(cores) | length(file_vec)==1) cores = 1
    
    cat("Running GSEM to compute genetic correlation with EXT factor, using ", cores, "threads \n")
    OUT = pbmclapply(1:length(file_vec), RUN_GSEM_RG_PAR, mc.cores = cores, ignore.interactive = TRUE)
    print("Result check (first 3 traits)")
    print(head(OUT,3))
    return( rbindlist(OUT) )
}
