# Porthmeus
# 26.04.21

# calculate the sample distances between FVAs
import os
#import numpy as np

#conditions = os.listdir("results/data/FVA/")
#conditions = np.unique([x.split("_")[-1].replace(".csv","") for x in conditions if x.endswith(".csv")])
#conditions = [x for x in conditions if x.startswith("FVA")]

smpls = ["{diet}.{model}-" + x +"_FVA{condition}.csv" for x in samples]

rule Analysis_calcFVAdist:
    input:
        fva_results = "results/FVA/FVA_{diet}.{model}-{condition}.RDS",
        # fva_results = [os.path.join("results/data/FVA/", x) for x in smpls],
        #expand("{diet}.{model}-{sample}_{condition}.csv",
         #       diet = "{diet}", model = "{model}", sample = samples, condition = "{conditions}")
         model = "results/data/consistentModels/Consistent_{diet}.{model}.csv"

    output:
        matrix = "results/FVA/dist_{diet}.{model}_{condition}.csv",
        rxnMatrix = "results/FVA/Rxndist_{diet}.{model}_{condition}.tar.gz"

    params:
        diet = "{diet}",
        model = "{model}",
        condition = "{condition}",
        samples = samples
    threads: 32
    conda: "../envs/Rparallel.yaml"
    log:
        "logs/calcFVAdist_{diet}.{model}_{condition}.log"
    script:
        "../scripts/Analysis_calcFVAdist.R"


