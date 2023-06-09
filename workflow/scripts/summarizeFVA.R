# Porthmeus
# 01.03.21

# read all the FVA analysis and combine the range of each reaction into one matrix - also write a matrix for min, max and center value
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

getSamplenName <- function(sample){
    # extract the name of the sample
    nm <- basename(sample)
    nm <- strsplit(nm, split = "-")[[1]]
    nm <- gsub("_FVA.*.csv","",paste(nm[-1], collapse = "-"))
    return(nm)
}

getRangeOfSample <- function(dat){
    # get the range of each reaction

    rng <- dat[,"maximum"] - dat[,"minimum"]
    names(rng) <- row.names(dat)
    return(rng)
}

samples <- snakemake@input[["samples"]]
rxns_file <- snakemake@input[["rxns"]]
out_range <- snakemake@output[["range"]]
out_min <- snakemake@output[["min"]]
out_max <- snakemake@output[["max"]]
out_center <- snakemake@output[["center"]]
out_RDS <- snakemake@output[["RDS"]]
sampleSel <- snakemake@params[["sampleSel"]]
sampleIDs <- snakemake@params[["sampleIDs"]]

# filter the samples to match the model/diet/FVA_condition
sel <- sapply(samples, function(x){strsplit(basename(x),"-")[[1]][1] == sampleSel})
samples <- samples[sel]

# define the matrix
rxns <- read.csv(rxns_file)


rangeMat <- matrix(0, ncol = length(samples),
              nrow = nrow(rxns),
              dimnames = list(rxn = rxns[["rxn"]], sample = sampleIDs))
minMat <- rangeMat
maxMat <- rangeMat
centerMat <- rangeMat

# populate the matrix
for(smplID in sampleIDs){
    smpl <- grep(smplID, samples, value = TRUE, fixed=TRUE)
    dat <- read.csv(smpl, row.names=1)
    vals<-getRangeOfSample(dat)
    rangeMat[names(vals),smplID] <- vals
    minMat[rownames(dat),smplID] <- dat[,"minimum"]
    maxMat[rownames(dat),smplID] <- dat[,"maximum"]
    centerMat[rownames(dat),smplID] <- apply(dat,1,mean)
}

write.csv(rangeMat, file = out_range)
write.csv(minMat, file = out_min)
write.csv(maxMat, file = out_max)
write.csv(centerMat, file = out_center)
saveRDS(simplify2array(list(min= minMat, max = maxMat, range = rangeMat, center = centerMat)), file = out_RDS)
