# Porthmeus
# 08.06.21

# do some basic associations on the fluxes of the models

rule Natrix_analyse_FBA:
    input:
        fba_results = "results/data/FBA/{diet}.{model}-fluxes.csv",
        subSysCount = "results/data/ModelAnalysis/{diet}.{model}-SubSysMat.csv",
        meta = "resources/META_emed_future.csv",
        clinic = "resources/ClinicalData.csv",
        subsysTable = "resources/Subsystems.csv",
        sampleBlacklist = "resources/sampleBlacklist.csv"
    output:
        html = "results/reports/{diet}.{model}-Natrix_analyse_FBA.html"
    wildcard_constraints:
        model = "[0-9a-zA-Z]+"
    conda:
        "../envs/R_analyseFBA.yaml"
    log:"logs/Analyse_FBA_{diet}.{model}.log"
    script:
        "../Rmarkdown/Natrix_analyse_FBA.Rmd"
        
