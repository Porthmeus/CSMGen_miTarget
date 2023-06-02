# Porthmeus
# 10.05.21

# fit lqmm models to each of the reactions and save the results in an object
rule Analysis_PERMANOVAonRxnDist:
    input:
        distTar = "results/FVA/Rxn{distance}_{diet}.{model}_{condition}.tar.gz",
        meta = "resources/META_emed_future.csv",
        clinic = "resources/ClinicalData.csv",
        sampleBlacklist = "resources/sampleBlacklist.csv"
    params:
        formula = "{formula}",
        tempdir = "NULL",
        #tempdir = "temp_PERMANOVA{diet}.{model}_{condition}vs{formula}-{distance}",
        permutations = "how(nperm = 999, blocks = data$PatientID)"
    threads: 32
    output:
        results = "results/FVA/PERMANOVAmodels_{diet}.{model}-{distance}FVA_{condition}vs{formula}.csv"
    log: "logs/PERMANOVAmodels_{diet}.{model}-{distance}FVA_{condition}_{formula}.log"
    conda: "../envs/Rparallel.yaml"
    script: "../scripts/Analysis_PERMANOVAonRxnDist.R"
