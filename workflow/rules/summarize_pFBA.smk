# Porthmeus
# 02.06.21

# summarize the pFBA results into a concise matrix

rule summarize_pFBA:
    input:
        samples =["results/data/FBA/fluxes_{diet}.{model}-"+ x +".csv" for x in samples],
        consistent = "results/data/consistentModels/Consistent_{diet}.{model}.csv",
    output:
        out = "results/FBA/fluxes_{diet}.{model}.csv"
    conda: "../envs/R.yaml"
    log: "logs/summarize_FBA_{diet}.{model}.log"
    script: "../scripts/summarize_FBA.R"

