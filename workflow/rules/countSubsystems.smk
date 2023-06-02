# Porthmeus
# 19.01.21

# load all the extraced models and count the reactions in the subsystem as well as their presence in general

rule countSubsystems:
    input:
        extModels = expand("results/data/extractedModels/{diet}.{model}-{sample}.csv",
                sample = samples, diet=dietNames, model = modelNames),
        cstModel = "results/data/consistentModels/Consistent_{diet}.{model}.csv",
        model = "resources/models/{model}.xml"
    output:
        rxnIdMat = "results/data/ModelAnalysis/rxnIdMat_{diet}.{model}.csv",
        subSysMat = "results/data/ModelAnalysis/SubSysMat_{diet}.{model}.csv"
    log: "logs/countSubsystems_{diet}.{model}.log"
    resources:
        mem_mb = 16000,
        time = "05:00:00",
        cpus_per_task = 8,
        task_per_node = 1
    conda: "../envs/python.yaml"
    script: "../scripts/countSubsystems.py"

