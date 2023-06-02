# Porthmeus
# 04.11.20

# take the TPM matrix and split it by columns, so that each columen can be used separately in the model extraction method

# some extra code to make snakemake find the right files and produce the right
# output - change the input file parameter here

rule splitCoreRxnsByColumn:
    input:
        matrix = "results/data/coreRxns/coreRxns.{thrld}.{model}.csv"
    output:
        Out = temp(
                expand("results/data/SPLIT_CoreRxns/SPLIT.{{thrld}}.{{model}}.{sample}.csv",
                    sample = samples) # samples are defined in the snakmake main file
                )
    params:
        outDir = "results/data/SPLIT_CoreRxns/",
        prefix = "SPLIT{thrld}.{model}."
    log: "logs/splitCoreRxns_{thrld}.{model}.log"
    resources:
        mem_mb = 16000,
        time = "01:00:00",
        cpus_per_task = 1,
        task_per_node = 1
    conda: "../envs/R.yaml"
    script: "../scripts/splitMatrixByColumn2.R"
