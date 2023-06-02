# Porthmeus
# 02.06.21

log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")


require(data.table)

dbg <- FALSE
if(dbg){
    dr <- "results/data/extractedModels"
    samples <- file.path(dr, list.files(dr))
    smpl <- samples[1]
    dat <- fread("temp/temp.csv")
    consistent <- "results/data/consistentModels/MatjesAbsorption.colormore22_Consistent.csv"
    out <- "temp/temp_fluxMat.csv"
}


samples <- snakemake@input[["samples"]]
consistent <- snakemake@input[["consistent"]]
out <- snakemake@output[["out"]]

# read the consistent model and create a matrix
cnsstnt <- fread(consistent)
mat <- matrix(0, ncol = length(samples), nrow= nrow(cnsstnt),
              dimnames = list(rxns = cnsstnt[,rxn], sample = samples))

# populate the flux matrix
for(smpl in samples){
    dat <- fread(smpl)
    mat[unlist(dat[,1]),smpl] <- unlist(dat[,2])
}

# correct the colnames of the matrix
samplenames <- sapply(basename(colnames(mat)), function(x) gsub(".csv$","",strsplit(x, split ="-")[[1]][2]))
colnames(mat) <- samplenames

write.csv(mat, file = out)
    


