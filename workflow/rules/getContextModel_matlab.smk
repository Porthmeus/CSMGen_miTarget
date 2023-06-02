# Porthmeus
# 03.11.20

# extract tissue/sample specific model


rule getContextModel:
    input:
        model = "results/data/CBT_models/{model}.mat",
        diet = "resources/diets/{model}_{diet}.csv",
        core = "results/data/coreRxns/coreRxns.{thrld}.{model}.csv",
        cnst_rxns = "results/data/consistentModels/CnstMod.{diet}.{model}.csv",
    output:
        "results/data/ContextSpecificModels/CntxtRxns.{thrld}.{sample}.{diet}.{model}.csv"
    log: "logs/extractModels_{diet}.{thrld}_{model}-{sample}.log"
    params:
        cobraPath = config["cobraToolbox_location"], # define the location of the cobra toolbox which should be initialized
        objective = "true", # this can be set to true if one is interested in running an FBA on the output
        sample = "{sample}", # current sample to extract
        zero_cutoff = "10e-4" # cutoff for which flux is considered zero
    resources:
        mem_mb = 16000,
        time = "02:00:00",
        cpus_per_task = 1,
        task_per_node = 1
    shell:
        '''echo "Starting matlab" &> {log};
        {config[matlab_command]} -nodisplay -nodesktop -nosplash -nojvm -r "try; addpath(genpath(fullfile('workflow','scripts'))); fastcore('{input.model}', '{output}', '{input.diet}', '{input.cnst_rxns}', '{input.core}', '{params.sample}', '{params.cobraPath}', '{params.zero_cutoff}'); catch ME; disp(ME); exit; end; exit;" &>> {log}'''

