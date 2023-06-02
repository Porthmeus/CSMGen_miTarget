# Porthmeus
# 30.03.21

# create a report for the analysis of the FVA range data

rule analyseFVArange:
    input:
        rangeMat = "results/data/FVA/{diet}.{model}-rangeFVA_{condition}.csv",
        lqmmModels = "results/data/FVA/{diet}.{model}-LQMMmodelsRangeFVA_{condition}.csv",
        meta = "resources/META_emed_future.csv",
        clinic = "resources/ClinicalData.csv",
        subsysTable = "resources/Subsystems.csv",
        sampleBlacklist = "resources/sampleBlacklist.csv"
    output:
        html = "results/reports/{diet}.{model}-rangeFVA_{condition}_analysis.html"
    log: "logs/{diet}.{model}-rangeFVA_{condition}_analysis.log"
    conda:
        "../envs/StatsRxnCounts.yaml"
    script: "../Rmarkdown/Natrix_analyseFVArange.Rmd"
