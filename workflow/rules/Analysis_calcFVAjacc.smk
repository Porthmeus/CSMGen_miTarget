# Porthmeus
# 26.04.21

# calculate the sample distances between FVAs
import os
#import numpy as np

#conditions = os.listdir("results/data/FVA/")
#conditions = np.unique([x.split("_")[-1].replace(".csv","") for x in conditions if x.endswith(".csv")])
#conditions = [x for x in conditions if x.startswith("FVA")]

smpls = ["{diet}.{model}-" + x +"_FVA{condition}.csv" for x in samples]

rule Analysis_calcFVAjacc:
    input:
        fva_results = "results/FVA/FVA_{diet}.{model}-{condition}.RDS",
        model = "results/data/consistentModels/Consistent_{diet}.{model}.csv"

    output:
        matrix = "results/FVA/Jacc_{diet}.{model}_{condition}.csv",
        rxnMatrix = "results/FVA/Rxnjacc_{diet}.{model}_{condition}.tar.gz"

    params:
        diet = "{diet}",
        model = "{model}",
        condition = "{condition}",
        samples = samples

    threads: 32
    conda: "../envs/Rparallel.yaml"
    log:
        "logs/calcFVAjacc_{diet}.{model}_{condition}.log"
    script:
        "../scripts/Analysis_calcFVAjacc.R"


