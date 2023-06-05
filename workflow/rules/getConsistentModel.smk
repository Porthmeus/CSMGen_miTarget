# Porthmeus
# 11.10.21

rule getConsistentModel:
    input:
        model = "resources/models/{model}.xml",
        diet = "resources/diets/{model}_{diet}.csv",
    output:
        cnst_rxns = "results/data/consistentModels/CnstMod.{diet}.{model}.csv"

    log: "logs/getConsistentModel_{diet}.{model}.log"
    conda: "../envs/python.yaml"
    threads: 12
    resources:
        mem_mb = 32000
    script: "../scripts/getConsistentModel.py"
