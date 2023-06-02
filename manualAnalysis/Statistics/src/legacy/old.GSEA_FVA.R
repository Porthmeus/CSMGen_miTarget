# Porthmeus
# 28.02.22

# make a GSEA for FVA results

# load libraries
require(data.table)
require(clusterProfiler)


# load data
subsystems <- fread("../dat/subsystems.csv")
cluster <- fread("../results/FVA_DBSCAN_Cluster.csv")
stats <- fread("../results/FVA_LMMsRandomPatient_fullCoefTable.csv")
stats.padj <- stats[coef != "(Intercept)",]
stats.padj[,padj := p.adjust(`Pr(>|t|)`, method = "BH")]



# first get all reactions
rxns.estimate <- cluster[,Rxn]

# get the estimate with the smallest p.value attached to it, if there are several for center and range
stats.padj[,.(Estimate = Estimate[which.min(padj)]), by=rxn]
