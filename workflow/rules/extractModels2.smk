# Porthmeus
# 03.11.20

# extract tissue/sample specific model


rule extractModels:
    input:
        TPM = "results/data/SPLIT_CoreGene_emed_future/SPLIT_{thrld}_{model}-{sample}.csv",
        model = "results/data/CBT_models/{model}.mat",
        diet = "resources/diets/{diet}.csv",
        CT = "resources/ensemble2HGNC.csv",
        consRxn = "results/data/consistentModels/Consistent_{diet}.{model}.csv",
    output:
        "results/data/extractedModels/{diet}.{thrld}_{model}-{sample}.csv"
    log: "logs/extractModels_{diet}.{thrld}_{model}-{sample}.log"
    params:
        cobraPath = config["cobraToolbox_location"], # define the location of the cobra toolbox which should be initialized
        objective = "true" # this can be set to true if one is interested in running an FBA on the output
    resources:
        mem_mb = 16000,
        time = "02:00:00",
        cpus_per_task = 1,
        task_per_node = 1
    shell:
        '''echo "Starting matlab" &> {log}; {config[matlab_command]} -nodisplay -nodesktop -nosplash -nojvm -r " try; addpath(genpath(fullfile('workflow','scripts'))); extractModels('{input.model}', '{output}', '{input.diet}', '{input.consRxn}', '{input.TPM}', '{input.CT}', {params.thrsld}, '{params.cobraPath}', '{params.objective}'); catch ME; rethrow(ME); exit; end; exit;" &>> {log}'''

