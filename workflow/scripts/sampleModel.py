# Porthmeus
# 09.02.21

# sample the possible flux space for each model from the extraction

sys.stdout = sys.stderr = open(snakemake.log[0], 'w')

import cobra
import pandas as pd
import numpy as np
import warnings

# load data
model = cobra.io.read_sbml_model(snakemake.input["model"])
extract = pd.read_csv(snakemake.input["extract"])
diet = pd.read_csv(snakemake.input["diet"])
threads = snakemake.threads
output = snakemake.output["out"]
n = snakemake.params["n"]

# reduce model
modRed = model.copy()
rxns = [x.id for x in modRed.reactions]
rm_rxns = [x for x in rxns if x not in extract.iloc[:,0].to_list()]

# security check before removing all or none of the reactions
if (extract.shape[0] + len(rm_rxns) != len(rxns)) or len(rm_rxns) == 0 or extract.shape[0] == 0:
    raise Exception("Something went wrong with the model extraction, please check the extracted reactions file and the model, compare namespaces!")

rm_rxns = [modRed.reactions.get_by_id(x) for x in rm_rxns]
modRed.remove_reactions(rm_rxns)


# adjust diet
diet = diet.set_index("ex_rxn",drop = False)
for ex_rxn in modRed.exchanges:
    ex_rxn.lower_bound = 0
    if ex_rxn.id in list(diet.index):
        if not pd.isna(diet.loc[ex_rxn.id][1]):
            ex_rxn.lower_bound = -1*diet.loc[ex_rxn.id][1]

# run the sampling
try:
    sampling = cobra.sampling.sample(modRed, n = n*threads, processes = threads, method = "optgp")
except Exception as optErr:
    err = optErr
    sampling = np.full(shape = [n*threads, len(modRed.reactions)], fill_value= np.nan)
    sampling = pd.DataFrame(sampling, columns = [x.id for x in modRed.reactions])
    warnings.warn(str(optErr) + "\n Will return a table with only NaN")

# write the output
sampling.to_csv(output)
