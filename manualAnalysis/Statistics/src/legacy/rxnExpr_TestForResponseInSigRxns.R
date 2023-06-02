# Porthmeus
# 17.03.22

require(data.table)
require(lme4)
require(lmerTest)
require(car)
require(parallel)
require(doMC)
require(foreach)


# load data
rxnExpr <- read.csv("../../../results/RxnExpression/rxnExpr.GL10-L50-GU90.colormore3D.csv", row.names = 1)

clinic <- fread("../../../resources/ClinicalData.csv")
meta <- fread("../../../resources/META_emed_future.csv")
blcklst <- fread("../../../resources/sampleBlacklist.csv")

# merge meta and clincal data
meta[,PatientID_Time := paste(PatientID, Time_seq, sep ="_")]
meta <- meta[Time_seq != -1,]
meta <- merge(meta,clinic, by="PatientID_Time", all.x=TRUE, suffix = c(".meta",""))

# remove degraded samples
blckSmpls <- blcklst[relevance < 2, Sample]
meta <- meta[!(SeqID %in% blckSmpls),]
rxnExpr <- rxnExpr[,!(colnames(rxnExpr) %in% blckSmpls)]
meta <- meta[SeqID %in% colnames(rxnExpr),]

# remove control samples
meta <- meta[Response != "C",]
rxnExpr <- rxnExpr[, meta[,SeqID]]

# load results from the HB/Mayo association
sigs.rxnExpr <- fread("../results/rxnExpr_coefficients_fullTable.csv")
sigs.rxnExpr <- sigs.rxnExpr[padj < 0.05,]

# reduce the rxnExpr matrix to those reactions with HB/Mayo association
rxns <- sigs.rxnExpr[,rxn]
rxnExpr <- rxnExpr[rxns.sig,]

# scale the values
rxnExpr <- t(scale(t(rxnExpr)))

# melt the data and merge it with the meta data
rxnExpr.melt <- melt(data.table(rxnExpr, keep.rownames = TRUE),
                     id.vars = "rn",
                     variable.name = "SeqID",
                     value.name = "Expression")

rxnExpr.melt <- merge(rxnExpr.melt, meta, by = "SeqID")


# fit the models
threads <- detectCores() - 1
registerDoMC(threads)

rxns <- unique(rxnExpr.melt[,rn])

mods <- foreach(rxn = rxns)%dopar%{
    dat <- rxnExpr.melt[rn == rxn,]
    glmer(data = dat,
         factor(Response) ~ HB_Mayo_impu + Expression + (1|PatientID),
         family = "binomial")
}
names(mods) <- rxns

stats <- list() 
for(rxn in names(mods)){
    mod <- mods[[rxn]]
    stats[[rxn]] <- cbind(
                          rxn = rxn, # add the reaction name to the table
                          as.data.frame(t(coef(summary(mod))["Expression",]))
                          )
}
stats <- do.call(rbind, stats)
stats <- data.table(stats)
stats[,padj := p.adjust(`Pr(>|z|)`, method = "BH")]
stats[,OR := round(exp(Estimate),3)]
stats[padj < 0.05,]

# save the model and the table
write.csv(stats, row.names=FALSE, file = "../results/rxnExpr_GLMMsResponse_coefs.csv")
saveRDS(mods, file = "../results/rxnExpr_GLMMsResponse_models.RDS")

