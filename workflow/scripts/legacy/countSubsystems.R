# Porthmeus
# 19.01.21

# count the presence and absence of rxns in the extraced models and do the same for the subsystems

log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

require(sybil)
if(!("sybilSBML" %in% installed.packages()[,"Package"])){
    devtools::install_github("cran/sybilSBML")
}
require(sybilSBML)

models <- snakemake@input[["models"]]
cstModel <- snakemake@input[["cstModel"]]

# models <- file.path("results/data/extractedModels/",
#                     list.files("results/data/extractedModels/", pattern = "Matjes-"))
# cstModel <- file.path("results/data/consistentModels/",
#                       list.files("results/data/consistentModels/", pattern = "Matjes_"))


# hack snakemake here, otherwise R will load all models for each consistent model
models_origin <- gsub("extractedModels","consistentModels",paste0(sapply(models, function(x) strsplit(x, split="-")[[1]][1]), "_Consistent.xml"))
sel <- models_origin == cstModel
models <- models[sel]


# define a small helper to convert the names back correctly
convertAsciiCoding <- function(strings){
    
    ascii_conv <- c(40,41,42,58)
    names(ascii_conv) <- strsplit(rawToChar(as.raw(ascii_conv)), split = "")[[1]]
    conversionChars <- paste0("__",ascii_conv,"__")
    names(conversionChars) <- names(ascii_conv)

    rtnString <- c()
    for(x in strings){
        for(i in 1:length(conversionChars)){
            ascii <- conversionChars[i]
            x <- gsub(ascii,names(ascii),x)
        }
        rtnString <- c(rtnString,x)
    }
    return(rtnString)
}

# read the original model from which the reactions where extracted
cstModel <- readSBMLmod(cstModel)
react_id(cstModel) <- convertAsciiCoding(react_id(cstModel))

# create an ID-matrix for each reaction in the extracted models
samples <- sapply(models, function(x) gsub("_SBML.xml","",strsplit(basename(x), split = "-")[[1]][2]))
names(models) <- samples
rxnIdMat <- matrix(0, ncol = length(samples), nrow = length(react_id(cstModel)),
                   dimnames=list(rxns = react_id(cstModel), sample = samples))
for(i in 1:length(models)){
    mod <- models[i]
    exMod <- readSBMLmod(mod)
    react_id(exMod) <- convertAsciiCoding(react_id(exMod))
    rxnIdMat[,names(mod)] <- rownames(rxnIdMat) %in% react_id(exMod)
}


# create a matrix for each subsystem and the number of rxns in this subsystem
subSysMat <- matrix(0, nrow=ncol(subSys(cstModel)), ncol = length(samples),
                    dimnames= list(subsystem = colnames(subSys(cstModel)), sample = samples))

subSysCst <- as.matrix(subSys(cstModel))
for(sbss in colnames(subSysCst)){
    for(smpl in colnames(rxnIdMat)){
        subSysMat[sbss,smpl] <- sum(react_id(cstModel)[subSysCst[,sbss]] %in% rownames(rxnIdMat)[as.logical(rxnIdMat[,smpl])])
    }
}

# write the matrices to disk
write.csv(file = snakemake@output[["subSysMat"]], subSysMat)
write.csv(file = snakemake@output[["rxnIdMat"]], rxnIdMat)
