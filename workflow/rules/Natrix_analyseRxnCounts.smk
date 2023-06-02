# Porthmeus
# 19.02.21

# run the analysis on the counts of the extracted reactions

rule Natrix_analyseRxnCounts:
    input:
        rxnCount = "results/data/ModelAnalysis/{diet}.{model}-rxnIdMat.csv",
        subSysCount = "results/data/ModelAnalysis/{diet}.{model}-SubSysMat.csv",
        meta = "resources/META_emed_future.csv",
        clinic = "resources/ClinicalData.csv",
        subsysTable = "resources/Subsystems.csv",
        sampleBlacklist = "resources/sampleBlacklist.csv"
    output:
        html = "results/reports/{diet}.{model}-NatrixAnalyseRxnCounts.html"
    log:
        "logs/NatrixAnalyseRxnCounts-{diet}.{model}.log"
    conda: "../envs/StatsRxnCounts.yaml"
    script: "../Rmarkdown/Natrix_analyseRxnCounts.Rmd"
