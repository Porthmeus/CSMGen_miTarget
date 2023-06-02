# Porthmeus
# 07.10.21

sys.stdout = sys.stderr = open(snakemake.log[0],"w")

import cobra as cb
import pandas as pd
import numpy as np
import warnings
import joblib
import multiprocessing
from corpse import simpleFastcore

dbg = False

if dbg:
    model = "resources/models/colormore22.xml"
    diet = "resources/diets/colormore22_MatjesAbsorption.csv"
    #core = "results/data/SPLIT_CoreRxns/SPLITGL25|L50|GU75_colormore22-F02243_L1_S10_L002.csv"
    out = "temp/inconsistentRxns.csv"
    cnst = "results/data/consistentModels/CnstMod.MatjesAbsorption.colormore22.csv"
    ncores = 3
    rxn_ix = "3"
else:
    model = snakemake.input["model"]
    diet = snakemake.input["diet"]
    cnst = snakemake.input["cnst_rxns"]
    out = snakemake.output["out"]
    ncores = snakemake.threads
    rxn_ix = snakemake.params["rxn"]

mod = cb.io.read_sbml_model(model)
diet = pd.read_csv(diet, index_col= 0)
cnst_rxn = pd.read_csv(cnst, index_col=0)

# adjust the diet
diet.loc[np.isnan(diet.iloc[:,0]),:] = 0
for rxn in mod.exchanges:
    if rxn.id in diet.index:
        rxn.lower_bound = -1*diet.loc[rxn.id,:][0]
    else:
        rxn.lower_bound = 0


# get the consistent model
cnst_rxns = cnst_rxn.loc[:,"rxn"].to_list()
noncon_rxns = [x.id for x in mod.reactions if x.id not in cnst_rxns]
mod2 = mod.copy()
mod2.remove_reactions(noncon_rxns)

# get the reaction to extract
rxn = cnst_rxns[int(rxn_ix)]

# create a small function to make use of parallel computation
def testReaction(model, rxn):
    fastcore = simpleFastcore(model = model, core_set = [rxn])
    try:
        fastcore.fastcore()
        return(None)
    except Exception as e:
        print(rxn)
        print(e)
        return(rxn)


# run fastcore on every reaction seperately and store those which will not work in fastcore with standard settings
out_rxn = testReaction(model = mod2, rxn = rxn)
if out_rxn == None:
    out_rxn = ""
with open(out, "w") as out_file:
    out_file.write(out_rxn)




# parallel version which causes an error for large data sets
#invalRxn = joblib.Parallel(n_jobs = ncores, verbose = 10)(joblib.delayed(testReaction)(model = mod2, rxn = r) for r in [x.id for x in mod2.reactions])
#invalDF = pd.DataFrame({"rxn" : [x for x in invalRxn if x != None]})
#invalDF.to_csv(out)
