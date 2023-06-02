# Porthmeus
# 11.02.21

# calculate a the FVA for each extracted model
rule sampleModel:
    input:
        model = "resources/models/{model}.xml",
        diet = "resources/diets/{diet}.csv",
        extract = "results/data/extractedModels/{diet}.{model}-{sample}.csv"
    output:
        out = "results/data/sampling/{diet}.{model}-{sample}_optgp.csv"
    params:
        n = 500
    log:
        "logs/{diet}.{model}-{sample}_optgp.log"
    conda:
        "../envs/python.yaml"
    resources:
        mem_mb = "16000",
        time = "08:00:00",
        task_per_node = 1,
    threads: 16
    conda: "../envs/python.yaml"
    script: "../scripts/sampleModel.py"

        
