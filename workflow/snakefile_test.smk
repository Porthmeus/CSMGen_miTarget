# A snakemake file to pipe the analysis of different transcriptomes in one framework to analyse the metabolic state of the data.

models_path =  "resources/models"
diets_path = "resources/diets"

# read the header of the initial read count matrix to obtain the sample names
# and compare it to the samples which have meta data associated
import pandas as pd
import itertools as iter
import os
from pip._internal import main as pip
try:
    import corpse
except ImportError:
    pip(["install","git+"+"https://github.com/Porthmeus/CORPSE.git"])
    import corpse

meta = pd.read_csv("resources/META_emed_future.csv")
samples = pd.read_csv("resources/RC_emed_future.csv")
samples = list(iter.compress(list(samples.columns),samples.columns.isin(meta.SeqID)))

models = [os.path.join(models_path,x) for x in os.listdir(models_path) if x.endswith(".xml")]
modelNames = [os.path.splitext(os.path.basename(model))[0] for model in models]

diets = [os.path.join(diets_path,x) for x in os.listdir(diets_path) if x.endswith(".csv")]
dietNames = [os.path.splitext(os.path.basename(diet))[0] for diet in diets]
fva_conditions = ["constr","zeroBiomass"]
distance_measures = ["dist","jacc"]
formulas = ["~HB_Mayo_impu","~Remission:Time_seq","~HB_Mayo_impu*Diagnosis"]
covariables = [" ","*Diagnosis"]
thrshlds = ["GL50","GL25|L50", "GL25|L50|GU90"] 

#### test version ####
samples = samples[0:10]
models = ["colormore22"]
diets = ["MatjesAbsorption"]
fva_conditions = ["zeroBiomass"]
distance = ["dist"]
formulas = ["HB_Mayo_impu"]
covariables = [" "]
#######################

rule main:
    input:
        expand("results/data/coreGenes/coreGenes.{thrld}_{model}.csv",
                thrld = thrshlds, model = models),
    '''
    input:
        [
                expand("results/data/ModelAnalysis/SubSysMat_{diet}.{model}.csv",
                    model = modelNames, diet = dietNames),
                expand("results/FVA/PERMANOVAmodels_{diet}.{model}-{distance}FVA_{condition}vs{formula}.csv",
                     diet = dietNames,
                     model = modelNames,
                     condition = fva_conditions,
                     distance = distance_measures,
                     formula = formulas),

                 expand("results/FVA/LQMMmodelsRangeFVA_{diet}.{model}-{condition}.csv",
                     diet = dietNames,
                     model = modelNames,
                     condition = fva_conditions),
                 expand("results/FVA/LMMmodelsRangeFVA_{diet}.{model}-{covariable}_{condition}.csv",
                     diet = dietNames,
                     model = modelNames,
                     condition = fva_conditions,
                     covariable = covariables),
                expand("results/FBA/fluxes_{diet}.{model}.csv",
                    diet = dietNames,
                    model = modelNames)
        ]
        '''

include:"rules/generateMatModel.smk"
include:"rules/createTPMMatrix.smk"
include:"rules/extractModels.smk"
include:"rules/splitTPMByColumn.smk"
include:"rules/extractConsistentModel.smk"
include:"rules/countSubsystems.smk"
include:"rules/sampleModel.smk"
include:"rules/calcFVA.smk"
include:"rules/summarizeFVA.smk"
include:"rules/calcFVA_zeroBiomass.smk"
include:"rules/Analysis_calcFVAdist.smk"
include:"rules/Analysis_calcFVAjacc.smk"
include:"rules/Analysis_PERMANOVAonRxnDist.smk"
include:"rules/Analysis_createLQMMmodels.smk"
include:"rules/Analysis_createLMMmodels.smk"
include:"rules/calc_pFBA.smk"
include:"rules/summarize_pFBA.smk"
include:"rules/getCoreGenes.smk"
