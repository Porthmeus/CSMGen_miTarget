# Porthmeus
# 01.11.21

# this is soo stupid and I apologize to everyone reading this. I will simply combine the output for all validated reactions into one file which I can further work on. 

import pandas as pd
import os

def getInputFiles(model, diet):
    rxns = pd.read_csv("results/data/consistentModels/CnstMod.{diet}.{model}.csv".format(diet = diet, model = model)).index.to_list()
    files = ["results/data/invalidFastcoreRxns/invalRxn.{diet}.{model}.{reaction}.txt".format(diet = diet, model = model, reaction = reaction) for reaction in rxns]
    return(files)


rule validateReactionsForExtraction_concat:
    input:
        # load the consistent model and 
        inval_rxns = lambda wildcards: getInputFiles(wildcards.model, wildcards.diet)

    output:
        out = "results/data/invalidFastcoreRxns/invalRxns.{diet}.{model}.csv"
    threads: 1
    log: "logs/invalRxns.{diet}.{model}.log"
    shell: "echo 'rxn' > {output.out} 2> {log}; cat {input.inval_rxns} >> {output.out} 2> {log}"
