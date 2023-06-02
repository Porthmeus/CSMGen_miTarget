# Porthmeus
# 10.05.21
dbg <- FALSE
if(!dbg){
    log <- file(snakemake@log[[1]], open = "wt")
    sink(log)
    sink(log, type = "message")
}
require(doMC)
require(foreach)
require(vegan)
require(data.table)


if(dbg){
    #distTar <- "results/data/FVA/MatjesAbsorption.colormore22_zeroBiomass_Rxnjacc.tar.gz"
    distTar <- "temp/temp_jacc.tar.gz"
    meta <- "resources/META_emed_future.csv"
    clinic <- "resources/ClinicalData.csv"
    sampleBlacklist <- "resources/sampleBlacklist.csv"
    temp <- tempdir()
    formula <- "~Remission:Time_seq"
    permutations <- "how(nperm=2, blocks = data$PatientID, minperm=1)"
    out <- file.path("temp","results.csv")
    threads <- parallel::detectCores()
} else {
    distTar <- snakemake@input[["distTar"]]
    meta <- snakemake@input[["meta"]]
    clinic <- snakemake@input[["clinic"]]
    sampleBlacklist <- snakemake@input[["sampleBlacklist"]]
    temp <- snakemake@params[["tempdir"]]
    formula <- snakemake@params[["formula"]]
    permutations <- snakemake@params[["permutations"]]
    threads <- snakemake@threads[[1]]
    if(temp == "NULL"){
        temp <- tempdir()
    }
    out <- snakemake@output[["results"]]
}

if (!dir.exists(temp)){
    dir.create(temp, recursive=TRUE)
}
print(temp)

# extract the files - this is somewhat confusing now as the cluster was not handling the files correctly and some major hacks had to be made. Here is what is done:
# first a temporary directory is generated - L41-44
# secondary the tar file is unzipped using a pipe to a file of the same name but localized in the temporary file - L47-53
# the files are listed afterwards - L55-56
# each file is than processed in parallel to read and perform the PERMANOVA
distTarOri <- distTar
distTar <- gsub(".gz","",distTar)
distTar <- file.path(temp,basename(distTar))
system(paste("gunzip -c",distTarOri, ">",distTar ))
if(!file.exists(distTar)){
    stop(paste0(distTar, " was not generated - cannot perform task"))
}
files <- untar(distTar, list=TRUE)
files <- grep(".csv$",files, value=TRUE)
if(dbg){
    files <- grep("\\(",files, value=TRUE)
}

# read the meta data
meta <- fread(meta)
meta[,PatientID_Time := paste(PatientID, Time_seq, sep ="_")]
meta <- meta[Time_seq != -1,]
clinic <- fread(clinic, drop=1)
meta <- merge(meta,clinic, by="PatientID_Time", all.x=TRUE, suffix = c(".meta",""))

# remove degraded samples
blcklst <- fread(sampleBlacklist)
blckSmpls <- blcklst[relevance < 2, Sample]
meta <- meta[!(SeqID %in% blckSmpls),]

#
setkey(meta, "SeqID")
dat <- as.data.frame(meta[Tissue=="Biopsy",])
row.names(dat) <- as.character(meta[Tissue == "Biopsy",SeqID])
print(permutations)

getVariablesFromFormula <- function(LHS_form){
    # the function basically splits the LHS of the formula for adonis into the
    # final coeffecients which are calculated - this is needed to construct the
    # summary table
    facs <- gsub("^~| ","", LHS_form)
    facs <- strsplit(facs, split="\\+")[[1]]
    facs <- sapply(facs, gsub, pattern ="(.*)\\*(.*)", replacement = "\\1\\+\\2\\+\\1:\\2")
    facs <- unlist(strsplit(facs, split="\\+"))
    return(facs)
}

testAndSumAdonis <- function(distMat, formula, data, permutations = 999){
    # basically a reformatting of the adonis output - the function performs the
    # an adonis test
    # @distMat - a csv file (sep = "," header = TRUE, row.names = 1, decimal
    # =".") containing the squared distance measure. The rownames/colnames
    # should be identifieable in the rownames of the table given in data.
    # @formula - the LHS of the formula to perform the adonis on
    # @data - the table which containst he variables to be tested - the factors
    # of formula should be available as columns in the data.frame. The rownames
    # should be contain colnames/rownames of the distMat. The table must be
    # indexable by rownames/colnames of distMat.
    # @permutations - the number of permutations to be executed to determine
    # the p-vale
    
    # match distance matrix and meta data
    #distMat <- read.csv(distMatFile,row.names=1)
    sel <- rownames(distMat) %in% rownames(data)
    distMat <- distMat[sel,sel]
    data <- data[rownames(distMat),]
    if(nrow(data) == 0 | nrow(distMat) == 0){
        stop("Rownames of the distance matrix is not matching the rownames of the meta data")
    }
    distMat <- as.dist(distMat)
    
    # get the permutations information if provided as character string through snakemake
    if(class(permutations) == "character"){
        permutations <- eval(parse(text=permutations))
    }
    
    # extract the reaction name from the file
    #rxn <- gsub(".csv$","",basename(distMatFile))

    # get the number of factors which have been tested
    facs <- getVariablesFromFormula(formula)
    i <- length(facs)
    
    if(var(distMat > 0)){
        # do the actual test
        fom <- as.formula(paste0("distMat",formula))
        adn <- adonis2(formula = fom, data=data, permutations=permutations, parallel = 1)
        
        # reformat the results
        df <- data.frame( #Rxn = rxn, 
                         Variable = facs,
                         SSq = adn[["SumOfSqs"]][1:i],
                         R2 = adn[["R2"]][1:i],
                         F = adn$F[1:i],
                         pval = adn[["Pr(>F)"]][1:i])
    } else {
        df <- data.frame(#Rxn = rxn, 
                         Variable = facs,
                         SSq = 0,
                         R2 = NaN,
                         F = NaN,
                         pval = NA)
    }

    return(df)

}

UntarAdonis <- function(file, distTar, formula, data, permutations){
    require(data.table)
    # extract the reaction ID
    rxn <- gsub(".csv$","",basename(file))

    file <- gsub("\\((.)\\)","\\\\(\\1\\\\)",file)
    distMat <- fread(cmd= paste("tar xOf", distTar, file))
    nms <- unlist(distMat[,1])
    distMat <- as.matrix(data.frame(distMat)[,-1])
    rownames(distMat) <- nms
    df <- testAndSumAdonis(distMat = distMat,
                     formula = formula,
                     data = data,
                     permutations = permutations)
    df <- cbind(Rxn = rxn, df)
    return(df)
}


registerDoMC(cores = threads)
resultsTab <- foreach(distFile = files, .combine = "rbind") %dopar%
    UntarAdonis(file = distFile,
                distTar = distTar,
                     formula = formula,
                     data = dat,
                     permutations = permutations)

resultsTab <- data.table(resultsTab)
resultsTab[,padj := p.adjust(pval, method = "BH"), by = .(Variable)]
                     
write.csv(resultsTab, file = out, row.names=FALSE) # safe the results

# remove files which might be still available
distTar <- distTarOri
if (dir.exists(gsub("tar.gz$","", distTar))){
    unlink(gsub("tar.gz$","", distTar), recursive = TRUE)
}
if (file.exists(gsub(".gz$","", distTar))){
    unlink(gsub(".gz$","", distTar), recursive = TRUE)
}
if (dir.exists(temp)){
    unlink(temp, recursive=TRUE)
}


