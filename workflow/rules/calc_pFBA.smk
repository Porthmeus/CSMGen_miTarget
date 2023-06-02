# Porthmeus
# 02.06.21

# calculate a pFBA for the extracted models

rule calc_pFBA:
    input:
        model = "resources/models/{model}.xml",
        diet = "resources/diets/{diet}.csv",
        extract = "results/data/extractedModels/{diet}.{model}-{sample}.csv"
    output:
        out = temp("results/data/FBA/fluxes_{diet}.{model}-{sample}.csv")
    params:
        frac_opt = 1,
        setMaxBounds = True
    conda: "../envs/python.yaml"
    log: "logs/calc_pFBA_{diet}.{model}-{sample}.log"
    script: "../scripts/calc_pFBA.py"
        
