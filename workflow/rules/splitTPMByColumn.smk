# Porthmeus
# 04.11.20

# take the TPM matrix and split it by columns, so that each columen can be used separately in the model extraction method

# some extra code to make snakemake find the right files and produce the right
# output - change the input file parameter here

rule splitTPMByColumn:
    input:
        TPM = "results/data/TPM_emed_future.csv"
    output:
        Out = temp(
                expand("results/data/SPLIT_TPM_emed_future/SPLIT_{sample}.csv",
                    sample = samples) # samples are defined in the snakmake main file
                )


    params:
        outDir = "results/data/SPLIT_TPM_emed_future"
    log: "logs/splitTPMByColumn.log"
    resources:
        mem_mb = 16000,
        time = "01:00:00",
        cpus_per_task = 1,
        task_per_node = 1
    conda: "../envs/R.yaml"
    script: "../scripts/splitMatrixByColumn.R"
