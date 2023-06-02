# Porthmeus
# 31.03.21

if(!exists("dbg")){
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")
dbg <- FALSE
}

if(dbg){

    setClass("snakemake", representation(input = "list", output = "list", threads="list", params = "list"))
    snakemake <- new("snakemake", 
                        input = list(
                            rangeMat = "results/data/FVA/MatjesAbsorption.colormore22-rangeFVA_zeroBiomass.csv",
                            meta = "resources/META_emed_future.csv",
                            clinic = "resources/ClinicalData.csv",
                            sampleBlacklist = "resources/sampleBlacklist.csv"
                        ),
                        output = list(
                            results ="results/data/FVA/MatjesAbsorption.colormore22-rangeFVA_zeroBiomass_LQMMmodels.csv"
                        ),
                        threads = list(1),
                        params = list(
                                      dependent = "HB_Mayo_impu",
                                      grp = "PatientID",
                                      random = "~1",
                                      covariable = "+Age"
                        )
                    )

}


print(snakemake@output[["results"]])

require(data.table)
if(!require(lqmm)){
    install.packages("lqmm")
    require(lqmm)
}
require(doMC)
require(parallel)
require(foreach)

# register threads
threads <- snakemake@threads[[1]]

# get the params for model generation
dependent <- snakemake@params[["dependent"]]
group <- snakemake@params[["grp"]]
random <- snakemake@params[["random"]]
covariable <- snakemake@params[["covariable"]]
if(length(covariable) > 1){
    covariable <- paste(covariable, collapse = "+")
}


# load the data
fva <- read.csv(snakemake@input[["rangeMat"]], row.names=1)

meta <- fread(snakemake@input[["meta"]])
meta[,PatientID_Time := paste(PatientID, Time_seq, sep ="_")]
meta <- meta[Time_seq != -1,]
clinic <- fread(snakemake@input[["clinic"]], drop=1)
meta <- merge(meta,clinic, by="PatientID_Time", all.x=TRUE, suffix = c(".meta",""))

# remove degraded samples
blcklst <- fread(snakemake@input[["sampleBlacklist"]])
blckSmpls <- blcklst[relevance < 2, Sample]

meta <- meta[!(SeqID %in% blckSmpls),]
meta <- meta[SeqID %in% colnames(fva),]
fva <- fva[,meta[,SeqID]]


# reformat the data

fva_rxn <- t(fva)

renaming <- data.frame(new= c("_hyphen_","_plus_","_dot_","_LPAREN_", "_RPAREN_"),
                       old = c("\\-","\\+","\\.","\\(","\\)"))
nms <- paste0("rxn_",colnames(fva_rxn))
for(i in 1:nrow(renaming)){
    nms <- gsub(renaming[i,"old"], renaming[i,"new"], nms)
}
colnames(fva_rxn) <-nms# gsub("-","_hyphen_",gsub("[.|+| |(|)]","_", paste0("rxn_",colnames(fva_rxn))))
mergeDat <- data.table(fva_rxn, keep.rownames=TRUE)
mergeDat <- merge(mergeDat, meta, by.x = "rn", by.y = "SeqID")

# do the stats only on the biopsy data
mergeDat <- mergeDat[Tissue == "Biopsy",]

rxn_sel <- (apply(fva_rxn[mergeDat[,rn],],2,var)>1e-6 &
            apply(fva_rxn[mergeDat[,rn],],2,
                  function(x) {
                      uniqv <- unique(x)
                      mod <- uniqv[which.max(tabulate(match(x, uniqv)))]
                      sum(mod == x)/length(x)
                  } < 0.9
                  ) 
            )
rxns <- colnames(fva_rxn)[rxn_sel]

getVariablesFromFormula <- function(RHS_form){
    # the function basically splits the LHS of the formula for adonis into the
    # final coeffecients which are calculated - this is needed to construct the
    # summary table
    facs <- gsub("^~| ","", RHS_form)
    facs <- strsplit(facs, split="\\+")[[1]]
    facs <- sapply(facs, gsub, pattern ="(.*)\\*(.*)", replacement = "\\1\\+\\2\\+\\1:\\2")
    facs <- unlist(strsplit(facs, split="\\+"))
    return(facs)
}

test <- function(a=3){
    Call <- match.call()
    print(Call[["a"]])
    return(a)
}


# short function to do the actual test
lqmmShort <- function(dat, formula_LHS, formula_RHS, random = as.formula("~1"), group="PatientID"){
    require(data.table)
   # 
    dat <- data.table(dat)
    # first check the formula and what to extract
    facs <- getVariablesFromFormula(formula_RHS)
    if(class(random) == "character"){
        random <- as.formula(random)
    }
    # hack to make the lqmm work in an automated fashion
    dat[,group := get(group)]
#
    fom <- paste0(formula_LHS,"~",formula_RHS)
    mod <- lqmm(data = dat,
            fixed = as.formula(fom),
            random = random,
            group = group)

    # another hack to get the correct formula into the summary function
    Call <- match.call()
    dd <- Call[["dat"]]
    mod$call <- str2lang(paste0("lqmm(data = ",dd, ", fixed = ", fom, ", random = ", paste(random, collapse = ""),", group = ",group,")"))
    sum_mod <- summary(mod)
    res <- data.frame(formula_LHS = formula_LHS,
                formula_RHS = formula_RHS,
                variable = facs,
                coefficient = sum_mod[["tTable"]][facs,1],
                stdErr = sum_mod[["tTable"]][facs,2],
                lw_bnd = sum_mod[["tTable"]][facs,3],
                up_bnd = sum_mod[["tTable"]][facs,4],
                pval = sum_mod[["tTable"]][facs,5],
                random = paste(as.character(random),collapse = ""),
                group = group
    )
    return(res)
}




# init parallel
print(threads)
registerDoMC(cores = threads)

# fit the models and get the stats
if(dbg){
    rxns <- tail(rxns)
}
results <- foreach(rxn = rxns, .combine = "rbind") %dopar% 
                lqmmShort(dat=mergeDat,
                          formula_LHS= dependent,
                          formula_RHS = paste0(rxn,covariable),
                          random = random,
                          group = group
                          )

# save the results
results <- data.table(results)
results[,padj := p.adjust(pval, method ="BH"),by = .(variable)]
vrbl <- gsub("^rxn_","",results[,variable])
for(i in 1:nrow(renaming)){
    vrbl <- gsub(renaming[i,"new"],renaming[i,"old"], vrbl)
}
results[,variable := vrbl]
if(!dbg){
    write.csv(as.matrix(results), file = snakemake@output[["results"]], row.names=FALSE)
}
### old unparallalized code
#i <- 0
#for(rxn in rxns){
#    print(rxn)
#    i <- i+1 
#    print(i)
#    mod <- lqmm(data = mergeDat,
#                fixed = as.formula(paste0("HB_Mayo_impu ~",rxn)),
#                random = ~1,
#                group = PatientID)
#    print(rxn)
#    sum_mod <- summary(mod)
#    pvals[rxn] <- sum_mod[["tTable"]][2,5]
#    coeffs[rxn] <- sum_mod[["tTable"]][2,1]
#}

#results <- data.table(Rxn = rxns,
#                      Rxn_ori = rownames(fva)[rxn_sel],
#                      coefficient = coeffs,
#                      pvalue = pvals,
#                      padj = p.adjust(pvals, method = "BH"))


