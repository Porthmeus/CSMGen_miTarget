# Porthmeus
# 31.03.21

# fit lqmm models to each of the reactions and save the results in an object
rule createLQMMmodels:
    input:
        rangeMat = "results/FVA/rangeFVA_{diet}.{model}-{condition}.csv",
        meta = "resources/META_emed_future.csv",
        clinic = "resources/ClinicalData.csv",
        sampleBlacklist = "resources/sampleBlacklist.csv"
    output:
        results = "results/FVA/LQMMmodelsRangeFVA_{diet}.{model}-{condition}.csv"
    params:
        dependent = "HB_Mayo_impu", # the fixed effect
        grp = "PatientID", # the random intercept
        random = "~1", # the random slope as RHS formula
        covariable = "" # covariates which should be controlled for - in R formula style, only the RHS, leaving the first parameter out, as this will be the rxn to be tested
    threads: 32
    log: "logs/LQMMmodelsRangeFVA_{diet}.{model}-{condition}.log"
    conda: "../envs/Rparallel.yaml"
    script: "../scripts/Analysis_createLQMMmodels.R"
