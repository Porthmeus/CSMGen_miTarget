# Porthmeus
# 11.02.21

# calculate a the FVA for each extracted model
rule calcFVA:
    input:
        model = "resources/models/{model}.xml",
        diet = "resources/diets/{diet}.csv",
        extract = "results/data/extractedModels/{diet}.{model}-{sample}.csv"
    output:
        out = temp("results/data/FVA/FVAconstr_{diet}.{model}-{sample}.csv")
    log:
        "logs/FVAconstr_{diet}.{model}-{sample}.log"
    conda:
        "../envs/python.yaml"
    params:
        frac_opt = 0.9,
        pfba_fac = 1.1,
        setMaxBounds = False
    resources:
        mem_mb = "16000",
        time = "08:00:00",
        task_per_node = 1,
    threads: 1
    conda: "../envs/python.yaml"
    script: "../scripts/calcFVA.py"

        
