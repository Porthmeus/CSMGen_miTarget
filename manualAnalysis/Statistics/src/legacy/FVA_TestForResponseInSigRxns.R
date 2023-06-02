# Porthmeus
# 17.03.22


# required packages
require(data.table)
require(lme4)
require(lmerTest)
require(foreach)
require(parallel)
require(doMC)



# load FVA results
FVA.min <- read.csv("../../../results/FVA/minFVA.GL10-L50-GU90.MatjesAbsorption.colormore3D.csv", row.names=1)
FVA.max <- read.csv("../../../results/FVA/maxFVA.GL10-L50-GU90.MatjesAbsorption.colormore3D.csv", row.names=1)
# convert the values to ranges and centers
FVA.range <- FVA.max-FVA.min
FVA.center <- (FVA.max+FVA.min)/2

# load metadata
clinic <- fread("../../../resources/ClinicalData.csv")
meta <- fread("../../../resources/META_emed_future.csv")
blcklst <- fread("../../../resources/sampleBlacklist.csv")

# merge meta and clincal data
meta[,PatientID_Time := paste(PatientID, Time_seq, sep ="_")]
meta <- meta[Time_seq != -1,]
meta <- merge(meta,clinic, by="PatientID_Time", all.x=TRUE, suffix = c(".meta",""))

# remove controls
meta <- meta[Response != "C",]

# remove degraded samples
blckSmpls <- blcklst[relevance < 2, Sample]
meta <- meta[!(SeqID %in% blckSmpls),]
meta <- meta[SeqID %in% colnames(FVA.min),]
sel<-!(colnames(FVA.min) %in% blckSmpls)
FVA.range <- FVA.range[,sel]
FVA.center <- FVA.center[,sel]

# load results from the HB/Mayo association
sigs.FVA <- fread("../results/FVA_LMMsRandomPatient_fullCoefTable.csv")
sigs.FVA[,padj := p.adjust(`Pr(>|t|)`, method = "BH")]
sigs.FVA <- sigs.FVA[padj < 0.05,]

# remove non-significant reactions
rxns <- unique(sigs.FVA[,rxn])
FVA.range <- FVA.range[rxns,]
FVA.center <- FVA.center[rxns,]

# scale FVA
FVA.range <- t(scale(t(FVA.range)))
FVA.center <- t(scale(t(FVA.center)))


# melt the FVA data to one data.table to make the associations
FVA.range.melt <- melt(
                       data.table(
                                  FVA.range,
                                  keep.rownames = TRUE
                                  ),
                       id.vars = "rn",
                       variable.name = "SeqID",
                       value.name = "range"
                       )
FVA.center.melt <- melt(
                       data.table(
                                  FVA.center,
                                  keep.rownames = TRUE
                                  ),
                       id.vars = "rn",
                       variable.name = "SeqID",
                       value.name = "center"
                       )
FVA <- merge(FVA.center.melt, FVA.range.melt)

# merge the meta data
FVA <- merge(FVA, meta, by="SeqID")


# fit the model
threads <- detectCores()-1
registerDoMC(threads)

mods <- foreach(reac = rxns) %dopar% {
    dat <- FVA[rn == reac,]
    form <- as.formula(paste0("factor(Response) ~ ",paste(sigs.FVA[rxn == reac, coef], collapse = "+"), "+ (1|PatientID)"))
    glmer(data = dat,
          formula = form,
          family = "binomial")
}
names(mods) <- rxns

# get the stats
stats <- list() 
for(rxn in names(mods)){
    mod <- mods[[rxn]]
    # get all the interesting coefficients
    form_split <- strsplit(as.character(formula(mod))[3], split = " \\+ ")[[1]]
    form_split <- form_split[-(length(form_split))]
    if(length(form_split) >1){
        dat <- as.data.frame(coef(summary(mod))[form_split,])
    } else {
        dat <- as.data.frame(t(coef(summary(mod))[form_split,]))
    }
    stats[[rxn]] <- cbind(
                          rxn = rxn, # add the reaction name to the table
                          coef = form_split,
                          dat
                          )
}
stats <- do.call(rbind, stats)
stats <- data.table(stats)
stats[,padj := p.adjust(`Pr(>|z|)`, method = "BH")]
stats[,OR := round(exp(Estimate),3)]



# save the models and the table
write.csv(stats, row.names=FALSE, file = "../results/FVA_GLMMsResponse_coefs.csv")
saveRDS(mods, file = "../results/FVA_GLMMsResponse_models.RDS")

