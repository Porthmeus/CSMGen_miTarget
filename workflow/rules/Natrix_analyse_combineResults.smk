# Porthmeus
# 11.06.21

rule combineResults:
    input:
        rxnCount = "results/data/ModelAnalysis/{diet}.{model}-rxnIdMat.csv",
        subSysCount = "results/data/ModelAnalysis/{diet}.{model}-SubSysMat.csv",
        lqmmModels = "results/data/FVA/{diet}.{model}-LQMMmodelsRangeFVA_{condition}.csv",
        rxn_stats = "results/data/FVA/{diet}.{model}-distFVA_{condition}vs~HB_Mayo_impu-PERMANOVAmodels.csv",
        fva = "results/data/FVA/{diet}.{model}-FVA_{condition}.RDS",
        distFiles = "results/data/FVA/{diet}.{model}_{condition}_Rxn{distance}.tar.gz",
        meta = "resources/META_emed_future.csv",
        clinic = "resources/ClinicalData.csv",
        subsysTable = "resources/Subsystems.csv",
        sampleBlacklist = "resources/sampleBlacklist.csv",
        GPR = "resources/{model}_GPR.csv"
    output:
        html = "results/reports/{diet}.{model}-{condition}FVA{distance}-Natrix_analyse_combineResults.html"
    log:
        "logs/{diet}.{model}-{condition}FVA{distance}-Natrix_analyse_combineResults.log"
    conda: "../envs/StatsRxnCounts.yaml"
    script: "../Rmarkdown/Natrix_analyse_combineResults.Rmd"

