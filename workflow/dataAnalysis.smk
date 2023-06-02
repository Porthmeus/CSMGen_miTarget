# Porthmeus
# 19.02.21

# seperate snakemake script for data analysis

models_path =  "resources/models"
diets_path = "resources/diets"

# read the header of the initial read count matrix to obtain the sample names
# and compare it to the samples which have meta data associated
import pandas as pd
import itertools as iter
import os
meta = pd.read_csv("resources/META_emed_future.csv")
samples = pd.read_csv("resources/RC_emed_future.csv")
samples = list(iter.compress(list(samples.columns),samples.columns.isin(meta.SeqID)))

models = [os.path.join(models_path,x) for x in os.listdir(models_path) if x.endswith(".xml")]
modelNames = [os.path.splitext(os.path.basename(model))[0] for model in models]

diets = [os.path.join(diets_path,x) for x in os.listdir(diets_path) if x.endswith(".csv")]
dietNames = [os.path.splitext(os.path.basename(diet))[0] for diet in diets]
conditions = ["zeroBiomass","constr"]
formulas = ["~HB_Mayo_impu","~Remission:Time_seq"]
distances = ["dist"]

wildcard_constraints:
    model = "[0-9a-zA-Z]+"

rule main:
    input:
        [expand("results/reports/{diet}.{model}-NatrixAnalyseRxnCounts.html",
                diet = dietNames, model = modelNames),
         expand("results/reports/{diet}.{model}-rangeFVA_{condition}_analysis.html",
             diet = dietNames, model = modelNames, condition = ["zeroBiomass","constr"]),
         expand("results/reports/{diet}.{model}-Natrix_analyse_FBA.html",
             diet = dietNames, model = modelNames),
         expand("results/reports/{diet}.{model}-Natrix_analyseFV{distance}-{condition}vs{formula}.html",
             diet = dietNames, model = modelNames, formula = formulas, condition = conditions, distance = distances),
         expand("results/reports/{diet}.{model}-{condition}FVA{distance}-Natrix_analyse_combineResults.html",
             diet = dietNames, model = modelNames, condition = conditions, distance = distances)

        ]

include:"rules/Natrix_analyseRxnCounts.smk"
include:"rules/Natrix_analyseFVArange.smk"
include:"rules/Natrix_analyse_FBA.smk"
include:"rules/Natrix_analyseFVAdist.smk"
include:"rules/Natrix_analyse_combineResults.smk"
