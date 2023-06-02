# Porthmeus
# 09.06.21

# automatically plot the results of the distance based stats

rule anlyseFVAdist:
    input:
        fva_dist = "results/data/FVA/{diet}.{model}_{condition}_{distance}.csv",
        rxn_stats = "results/data/FVA/{diet}.{model}-{distance}FVA_{condition}vs{formula}-PERMANOVAmodels.csv",
        meta = "resources/META_emed_future.csv",
        clinic = "resources/ClinicalData.csv",
        sampleBlacklist = "resources/sampleBlacklist.csv",
        subsysTable = "resources/Subsystems.csv"
    output:
        html = "results/reports/{diet}.{model}-Natrix_analyseFV{distance}-{condition}vs{formula}.html"
    conda:
        "../envs/R_analyseFBA.yaml"
    params:
        variable = "{formula}"
    log: "logs/{diet}.{model}-Natrix_analyseFV{distance}-{condition}vs{formula}.log"
    script:
        "../Rmarkdown/Natrix_analyseFVAdist.Rmd"
