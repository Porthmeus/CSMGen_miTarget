# Porthmeus
# 02.02.21

# takes an sbml model imports it to matlab's cobratoolbox and exports it as a .mat file, in that way import to matlab is much faster

rule generateMatModel:
    input: 
        model = "resources/models/{model}.xml"
    output:
        model = "results/data/CBT_models/{model}.mat"
    params:
        cobraPath = config["cobraToolbox_location"]
    threads: 1
    log: "logs/generateMatModel_{model}.log"
    shell:
        '''echo "Starting matlab" &> {log}; {config[matlab_command]} -nodisplay -nodesktop -nosplash -nojvm -r "try; addpath(genpath(fullfile('workflow','scripts'))); generateMatModel('{input.model}','{output.model}', '{params.cobraPath}'); catch ME; rethrow(ME); exit; end; exit" &>>{log}'''
