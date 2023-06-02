# Porthmeus
# 26.04.21

debug <- FALSE

if(!debug){
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")
}

require(data.table)

if(!debug){
files <- snakemake@input[["fva_results"]]
mod <- snakemake@input[["model"]]
condition <- snakemake@params[["condition"]]
model <- snakemake@params[["model"]]
diet <- snakemake@params[["diet"]]
samples <- snakemake@params[["samples"]]
threads <- snakemake@threads[[1]]
out <- snakemake@output[["matrix"]]
out2 <- snakemake@output[["rxnMatrix"]]
RDS <- snakemake@input[["fva_results"]]
}

### debug ###
if(debug){
    debug <- TRUE
    mod <- "results/data/consistentModels/MatjesAbsorption.colormore22_Consistent.csv"
    samples <- fread("resources/RC_emed_future.csv")
    meta <- fread("resources/META_emed_future.csv")
    samples <- colnames(samples)[colnames(samples) %in% meta[,SeqID]]
    files <- file.path("results/data/FVA", paste0("MatjesAbsorption.colormore22-", samples,"_FVAzeroBiomass.csv"))
    RDS <- "results/FVA/FVA_MatjesAbsorption.colormore22-zeroBiomass.RDS"
}

mod <- fread(mod)

#minMax <- array(0, dim = c(nrow(mod), length(samples), 2), 
#        dimnames = list(rxns = mod[["rxn"]], samples = samples, bounds = c("minimum","maximum")))
#
#for(smpl in samples){
#    f <- grep(smpl, files, value = TRUE)
#    tbl <- fread(f)
#    minMax[tbl[,V1],smpl,"minimum"] <- tbl[,minimum]
#    minMax[tbl[,V1],smpl,"maximum"] <- tbl[,maximum]
#}
minMax <- readRDS(RDS)
minMax <- minMax[,,1:2]

calcRxnDist <- function(rxn1Bounds, rxn2Bounds, normalize = TRUE){
    # the function calculates the distance between two reactions for which an fva was calculated
    # params:
    # @rxn1Bounds - a vector of length two giving the lower and upper boundaries from the FVA
    # @rxn2Bounds - same as rxn1Bounds for the second reaction
    # @normalize - whether the result should be normalized by the range covered by the two reactions
    
    dist <- max(    rxn1Bounds[1] - rxn2Bounds[2],
                    rxn2Bounds[1] - rxn1Bounds[2],
                    rxn1Bounds[1] - rxn1Bounds[2],
                    rxn2Bounds[1] - rxn2Bounds[2]
                    )
    # normalize distance by mean range of the two reactions
    if(normalize){
        fac <- mean(c(rxn1Bounds[2] - rxn1Bounds[1],
                    rxn2Bounds[2] - rxn2Bounds[1]))
        if(fac == 0){
            fac <- 1
            dist <- dist-1
        }
        # add +1 to transform the distance back to a positive realm
        dist <- (dist/fac) +1
    }
    return(dist)
}

minMaxNorm <- function(x){
    # a function to normalize the values between the highest and the lowest value of a vector - resulting in values distributed between 0 and 1
    # @x - the vector of values to be normalized

    val <- (x - min(x))/(max(x)-min(x))
    val[is.na(val)] <- 0
    return(val)
}


calcSampleDistance <- function(boundsArray, normalize=TRUE, minMaxNorm = TRUE, filename = NULL){
    # a function which calculates the fva-distance over all reactions of different samples/simulations 
    # @boundsArray - a 3-dimensional array with reactions, samples and min-max bounds in the first, second, and third dimension
    # @normalize - whether to normalize the reaction distance to be spaced between 0 and 1
    # value:
    # a 2 dimensional distance matrix for all sample pairs

    if(!is.null(filename)){
        if(!grepl("tar.gz$", filename)){
            dirnm <- filename
            filename <- paste0(filename, ".tar.gz")
        } else {
            dirnm <- gsub("tar.gz$","", filename)
        }
        dir.create(dirnm)
    }

    mat <- matrix(0,
                 ncol = dim(boundsArray)[2],
                 nrow = dim(boundsArray)[2],
                 dimnames = list(dimnames(boundsArray)[[2]],
                                 dimnames(boundsArray)[[2]])
                 )
    pb <- txtProgressBar(min = 0, max = dim(boundsArray)[1],style = 3)
    i <- 0
    for(rxn in dimnames(boundsArray)[[1]]){
        i <- i+1
        if(!is.null(filename)){
            fn <- file.path(dirnm, paste0(rxn,".csv"))
        } else {
            fn <- NULL
        }
        dist <- calcDistRxnSamples(boundsArray, rxn, normalize =normalize, minMaxNorm = minMaxNorm, filename = fn)
        mat <- mat+dist
        setTxtProgressBar(pb,i)
    }
    close(pb)
    if(!is.null(filename)){
        # compress all the files into a single tar.gz archive
        tar(tarfile = filename, files = dirnm, compression = "gzip")
        unlink(dirnm, recursive = TRUE)
    }
    mat <- mat/dim(boundsArray)[1]
    return(mat)
}          


calcSampleDistancePar <- function(boundsArray, normalize = TRUE, minMaxNorm=TRUE, filename = NULL, threads = parallel::detectCores()){
    # a function which calculates the fva-distance over all reactions of different samples/simulations 
    # @boundsArray - a 3-dimensional array with reactions, samples and min-max bounds in the first, second, and third dimension
    # @normalize - whether to normalize the reaction distance to be spaced between 0 and 1
    # @filename - a filepath to save the compressed (tar.gz) rxn distances between samples
    # value:
    # a 2 dimensional distance matrix for all sample pairs

    require(foreach)
    require(doMC)
    registerDoMC(threads)

    # sort the filenames
    if(!is.null(filename)){
        if(!grepl("tar.gz$", filename)){
            dirnm <- filename
            filename <- paste0(filename, ".tar.gz")
        } else {
            dirnm <- gsub("tar.gz$","", filename)
        }
        dir.create(dirnm)
        
        # run the calculation and save the results for each rxn in an extra file
        dist <- foreach(rxn = dimnames(boundsArray)[[1]], .combine = "+") %dopar%
            calcDistRxnSamples(boundsArray, rxn, normalize = normalize, minMaxNorm = minMaxNorm, filename = file.path(dirnm, paste0(rxn,".csv")))

        # compress all the files into a single tar.gz archive
        tar(tarfile = filename, files = dirnm, compression = "gzip")
        unlink(dirnm, recursive = TRUE)

    } else {

        # otherwise just calculate without any savings of temporary results
        dist <- foreach(rxn = dimnames(boundsArray)[[1]], .combine = "+") %dopar%
            calcDistRxnSamples(boundsArray, rxn, normalize = normalize, minMaxNorm = minMaxNorm)
    }
    
    # normalize the distance by the number of reactions
    dist <- dist/dim(boundsArray)[1]
    return(dist)
}

calcDistRxnSamples <- function(boundsArray, rxn, normalize=TRUE, minMaxNorm=TRUE, filename = NULL){
    mat_check <- matrix(-1,
             ncol = dim(boundsArray)[2],
             nrow = dim(boundsArray)[2],
             dimnames = list(dimnames(boundsArray)[[2]],
                             dimnames(boundsArray)[[2]])
             )
    dist <- matrix(0,
             ncol = dim(boundsArray)[2],
             nrow = dim(boundsArray)[2],
             dimnames = list(dimnames(boundsArray)[[2]],
                             dimnames(boundsArray)[[2]])
             )

    for(s1 in dimnames(boundsArray)[[2]]){
        for(s2 in dimnames(boundsArray)[[2]]){
            if(s1 != s2 & (mat_check[s1,s2]==-1 & mat_check[s2,s1] == -1)){
            mat_check[s1,s2] <- 1
            mat_check[s2,s1] <- 1
            d <- calcRxnDist(boundsArray[rxn,s1,],
                                boundsArray[rxn,s2,], 
                                normalize=normalize)
            dist[s1,s2] <- d
            dist[s2,s1] <- d
            }
        }
    }
    if(minMaxNorm){
        dist <- minMaxNorm(dist)
    }
    if(!is.null(filename)){
        write.csv(dist, file = filename)
    }
    return(dist)
}

if(!debug){
    if(out2 %in% c("NULL","","null","Null")){
        out2 <- NULL
    }
    distMat <- calcSampleDistancePar(minMax, normalize=TRUE, minMaxNorm = TRUE,filename =out2)
    write.csv(distMat, file = out)
}
