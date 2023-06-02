# Porthmeus
# 01.03.21

# combine the results of the FVA into a single matrix of ranges

rule rangeFVA_constr:
    input:
        rxns = "results/data/consistentModels/{diet}.{model}_Consistent.csv",
        samples = expand("results/data/FVA/{diet}.{model}-{sample}_FVAconstr.csv",
                diet = dietNames,
                model = modelNames,
                sample = samples)
    output:
        mat = "results/data/FVA/{diet}.{model}-rangeFVA_constr.csv"
    log : "logs/{diet}.{model}_rangeFVA_constr.csv"
    conda : "../envs/R.yaml"
    script : "../scripts/rangeFVA.R"
