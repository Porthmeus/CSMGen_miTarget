# Porthmeus
# 06.11.20

# extract a consistent model (no dead ends in the reactions) from the initial
# metabolic model

rule extractConsistentModel:
    input: 
        model = 'results/data/CBT_models/{model}.mat',
        diet = "resources/diets/{diet}.csv"
    output: 'results/data/consistentModels/Consistent_{diet}.{model}.csv'
    log: 'logs/extractConsistent_{diet}.{model}.log'
    params:
        cobraPath = config["cobraToolbox_location"]
    resources:
        mem_mb = 16000,
        time = "01:00:00",
        cpus_per_task = 1,
        task_per_node = 1
    shell:
        '''echo "Starting matlab" &> {log}; {config[matlab_command]} -nodisplay -nodesktop -nosplash -nojvm -r "try; addpath(genpath('workflow/scripts')); extractConsistentModel('{input.model}', '{input.diet}', '{output}','{params.cobraPath}'); exit; catch ME; disp(ME); exit; end;" &>>{log}'''
        
