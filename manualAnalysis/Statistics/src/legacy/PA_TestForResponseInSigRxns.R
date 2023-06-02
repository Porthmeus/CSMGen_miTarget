# Porthmeus
# 17.03.22

# Test for Response in those reactions which have a association to HB/Mayo

require(data.table)
require(lme4)
require(lmerTest)
require(car)
require(parallel)
require(foreach)
require(doMC)

# start with the PA analysis

rxnCount <- read.csv("../../../results/ModelAnalysis/rxnIdMat.GL10-L50-GU90.MatjesAbsorption.colormore3D.csv", row.names=1)

clinic <- fread("../../../resources/ClinicalData.csv")
meta <- fread("../../../resources/META_emed_future.csv")
blcklst <- fread("../../../resources/sampleBlacklist.csv")

# make rxnCounts a integer matrix
rxnCount2 <- rxnCount
rxnCount2[,] <- 0
rxnCount2[rxnCount == "True"] <- 1
rxnCount <- rxnCount2

# merge meta and clincal data
meta[,PatientID_Time := paste(PatientID, Time_seq, sep ="_")]
meta <- meta[Time_seq != -1,]
meta <- merge(meta,clinic, by="PatientID_Time", all.x=TRUE, suffix = c(".meta",""))


# remove degraded samples
blckSmpls <- blcklst[relevance < 2, Sample]
meta <- meta[!(SeqID %in% blckSmpls),]
rxnCount <- rxnCount[,!(colnames(rxnCount) %in% blckSmpls)]
meta <- meta[SeqID %in% colnames(rxnCount),]

# remove control samples
meta <- meta[Response != "C",]
rxnCount <- rxnCount[, meta[,SeqID]]

# load results from the HB/Mayo association
sigs.PA <- fread("../results/presenceAbsence_LMMsRandomPatient_fullCoefTable.csv")
sigs.PA <- sigs.PA[padj < 0.05,]
sigs.PA[,rxn := rn]


# reduce the rxnCount matrix to those reactions with HB/Mayo association
rxns.sig <- sigs.PA[,rxn]
rxnCount <- rxnCount[rxns.sig,]

# sort meta data to fit the samples
meta <- data.frame(meta)
rownames(meta) <- meta[,"SeqID"]
meta <- meta[colnames(rxnCount),]

# create models where Response is a factor in the model
threads <- detectCores()-1
registerDoMC(threads)

mods <- foreach(rxn = rownames(rxnCount)) %do% {
    dat <- data.frame(Response = as.factor(meta[["Response"]]),
                      HB_Mayo_impu = meta[["HB_Mayo_impu"]],
                      rxnCount = as.factor(unlist(rxnCount[rxn,])),
                      Time = meta[["Time_seq"]],
                      PatientID = meta[["PatientID"]])
    mod <- tryCatch({ 
        glmer(data = dat,
         Response ~ HB_Mayo_impu+rxnCount+(1|PatientID),
         family = "binomial")
    }, error = function(e) {
        print(e)
        print(paste0("Model expression for ", rxn, " failed. Will fall back on only rxnCounts association with response"))
        glmer(data = dat,
         Response ~ rxnCount+(1|PatientID),
         family = "binomial")
    })
}
names(mods) <- rownames(rxnCount)

# create the statistics
stats <- list() 
for(rxn in names(mods)){
    mod <- mods[[rxn]]
    stats[[rxn]] <- cbind(
                          rxn = rxn, # add the reaction name to the table
                          HBMayo_included = grepl("HB_Mayo_impu",formula(mod)[3]), # indicate if the model fell back to the more simple version because of complete seperation
                          as.data.frame(t(coef(summary(mod))["rxnCount1",]))
                          )
}
stats <- do.call(rbind, stats)
stats <- data.table(stats)
stats[,padj := p.adjust(`Pr(>|z|)`, method = "BH")]
stats[,OR := round(exp(Estimate),3)]

# save the models and the table
write.csv(stats, row.names=FALSE, file = "../results/presenceAbsence_GLMMsResponse_coefs.csv")
saveRDS(mods, file = "../results/presenceAbsence_GLMMsResponse_models.RDS")

#for(rxn in stats[padj < 0.05, rxn]){
#    dat <- data.frame(Response = as.factor(meta[["Response"]]),
#                      HB_Mayo_impu = meta[["HB_Mayo_impu"]],
#                      rxnCount = as.factor(unlist(rxnCount[rxn,])),
#                      Time = meta[["Time_seq"]],
#                      PatientID = meta[["PatientID"]])
#    print(rxn)
#    print(table(dat[,c("rxnCount","Time","Response")])) 
#}
