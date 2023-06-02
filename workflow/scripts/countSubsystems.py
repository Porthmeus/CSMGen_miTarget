# Porthmeus
# 25.01.21

# count the reactions and the reactions in the subsystems for each extracted model
sys.stdout = sys.stderr = open(snakemake.log[0], 'w')

import pandas as pd
import cobra
import os
import numpy as np

extModels = snakemake.input["extModels"]
cstModel = snakemake.input["cstModel"]
model = snakemake.input["model"]


# hack snakemake here, otherwise all models in the extractedModels folder are loaded - select only those models, which where derived from the parsed consistent model
sel = [os.path.basename(x).split("-")[0] == os.path.basename(cstModel).split("_")[0] for x in extModels]
extModels = [x for x,y in zip(extModels,sel) if y]



# read the original model from which the reactions where extracted
cstModel = pd.read_csv(cstModel)
mod_ori= cobra.io.read_sbml_model(model)
mod_ori.remove_reactions([x for x in mod_ori.reactions if x.id not in cstModel.rxn.to_list()])


# create data frame with zeros for each reaction and sample to fill it with ids afterwards

samples = [os.path.basename(x).strip(".csv").split("-")[1] for x in extModels]
rxns = cstModel.rxn
rxnIDMat = np.zeros((len(rxns),len(samples)), dtype = bool)
rxnIDMat = pd.DataFrame(rxnIDMat, columns = samples, index = rxns)

# read the extracted models and check the reactions therein
for extModel in extModels:
    smpl = os.path.basename(extModel).strip(".csv").split("-")[1]
    extMod = pd.read_csv(extModel)
    rxns = extMod.rxn
    rxnIDMat.loc[rxns, smpl] = True

# remove the reactions with only zeros
rxnIDMat = rxnIDMat.loc[np.sum(rxnIDMat,1) != 0,:]
    

# count the reactions in the subsystems
groups = [x.id for x in mod_ori.groups]
subSysMat = np.zeros((len(groups), len(samples)),dtype = int)
subSysMat = pd.DataFrame(subSysMat, columns = samples, index = groups)

for group in groups:
    rxns = [x.id for x in mod_ori.groups.get_by_id(group).members]
    for extModel in extModels:
        smpl = os.path.basename(extModel).strip(".csv").split("-")[1]
        subSysMat.loc[group, smpl] = sum(rxnIDMat.loc[rxns,smpl])


# write the matrices to disk
rxnIDMat.to_csv(snakemake.output["rxnIdMat"])
subSysMat.to_csv(snakemake.output["subSysMat"])

