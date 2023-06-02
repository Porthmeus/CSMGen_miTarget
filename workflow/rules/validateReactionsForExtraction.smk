# Porthmeus
# 29.10.21

# this is another layer of attempts to get fastcore running by testing all reactions whether they cause the inconsistency error

rule validateReactionsForExtraction:
    input:
        model = "resources/models/{model}.xml",
        diet = "resources/diets/{model}_{diet}.csv",
        cnst_rxns = "results/data/consistentModels/CnstMod.{diet}.{model}.csv"
    params:
        rxn = "{reaction}"
    output:
        out = temp("results/data/invalidFastcoreRxns/invalRxn.{diet}.{model}.{reaction}.txt")
    threads: 1
    conda: "../envs/python.yaml"
    log: "logs/invalRxn.{diet}.{model}.{reaction}.log"
    script: "../scripts/validateReactionsForExtraction.py"
