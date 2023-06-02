# Porthmeus 
# 04.11.20

log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

require(data.table)

# define parsed variables

TPMmat <- snakemake@input[["TPM"]] #"~/Work/miTarget/FUTURE/PipelineAnalysis/results/data/TPM_emed_future.csv"
outDir <- snakemake@params[["outDir"]] #"~/Work/miTarget/FUTURE/PipelineAnalysis/results/data/TPM_emed_future_SPLIT"

# load data and split it columnwise

if (!dir.exists(outDir)){
    dir.create(outDir, recursive = TRUE)
}


df <- fread(TPMmat)
colnames(df)[1] <- "ENSEMBL"
for(clm in 2:ncol(df)){
    write.csv(df[,c(1,get('clm'))],
              file= file.path(outDir,paste0("SPLIT_",colnames(df)[clm],".csv")), 
              row.names = FALSE)
}

