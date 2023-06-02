# Porthmeus
# 31.03.21

# fit lmm models to each of the reactions and save the results in an object
rule createLMMmodels:
    input:
        rangeMat = "results/FVA/rangeFVA_{diet}.{model}-{condition}.csv",
        meta = "resources/META_emed_future.csv",
        clinic = "resources/ClinicalData.csv",
        sampleBlacklist = "resources/sampleBlacklist.csv"
    output:
        results = "results/FVA/LMMmodelsRangeFVA_{diet}.{model}-{covariable}_{condition}.csv"
    params:
        dependent = "HB_Mayo_impu", # the fixed effect{covariable}.
        random = "+(1|PatientID)", # the random slope|intercept as RHS formula
        covariable = "{covariable}" # covariates which should be controlled for - in R formula style, only the RHS, leaving the first parameter out, as this will be the rxn to be tested
    threads: 32
    log: "logs/LMMmodelsRangeFVA_{diet}.{model}-{covariable}.{condition}.log"
    conda: "../envs/Rparallel.yaml"
    script: "../scripts/Analysis_createLMMmodels.R"
