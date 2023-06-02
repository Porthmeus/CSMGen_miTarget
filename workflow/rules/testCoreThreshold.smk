# Porthmeus
# 17.12.21

# build different machine learning models to test the predictive capacity of the different core reaction thresholds

rule testCoreThresholds:
    input:
        meta = "resources/META_emed_future.csv",
        clinic = "resources/ClinicalData.csv",
        rxns = "results/ModelAnalysis/rxnIdMat.{thrsld}.{diet}.{model}.csv",
        fva_center = "results/FVA/centerFVA.{thrsld}.{diet}.{model}.csv",
        fva_range = "results/FVA/rangeFVA.{thrsld}.{diet}.{model}.csv",
    output:
        stats = "results/testCoreThresholds/MLstats.{thrsld}.{diet}.{model}.csv",
    conda: "../envs/R_ML.yaml"
    log: "logs/testCoresThresholds.{thrsld}.{diet}.{model}.log"
    threads: 32
    params: n_repeats = 5 # how often should the data set be split into train and test set to obtain the measures for the regression
    script: "../scripts/testCoreThreshold.R"
