# Porthmeus 
# 04.11.20

log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

require(data.table)

# define parsed variables

Mat <- snakemake@input[["matrix"]] #
#Mat <- "~/Work/miTarget/eMed_FUTURE/Pipeline_MM/results/data/coreRxns/coreRxns.GL25|L50|GU75_colormore3D.csv"
outDir <- snakemake@params[["outDir"]]
#outDir <- "temp/SPLIT"
prefix <- snakemake@params[["prefix"]]
#prefix <- "SPLIT_GL25|L50|GU75"


# load data and split it columnwise

if (!dir.exists(outDir)){
    dir.create(outDir, recursive = TRUE)
}


df <- fread(Mat)
colnames(df)[1] <- "ID"
for(clm in 2:ncol(df)){
    write.csv(df[,c(1,get('clm'))],
              file= file.path(outDir,paste0(prefix,colnames(df)[clm],".csv")), 
              row.names = FALSE)
}

