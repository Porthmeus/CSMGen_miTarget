# Porthmeus
# 01.06.21

# calculate a FBA for each model from the extraction
dbg = False
if not dbg:
    dbg = False
    sys.stdout = sys.stderr = open(snakemake.log[0], 'w')
else:
    print("in debug mode")
    dbg = True

import cobra
import pandas as pd
import numpy as np
import warnings


if dbg:
    model = cobra.io.read_sbml_model("resources/models/colormore22.xml")
    dr = "results/data/extractedModels/MatjesAbsorption.colormore22-D5247_S53_L006.csv"
    #extract_files = os.listdir(dr)[1:20]
    #extract_files = [os.path.join(dr,x) for x in extract_files]
    extract_file = dr
    diet = pd.read_csv("resources/diets/MatjesAbsorption.csv")
    output = "temp/temp.csv"

    # set the options
    threads =  3
    frac_opt = 1 
    setMaxBounds = True 

else:
    # load the data
    model = cobra.io.read_sbml_model(snakemake.input["model"])
    extract_file = snakemake.input["extract"]
    diet = pd.read_csv(snakemake.input["diet"])
    output = snakemake.output["out"]

    # set the options
    threads = snakemake.threads
    frac_opt = snakemake.params["frac_opt"]
    setMaxBounds = snakemake.params["setMaxBounds"]

# create a data frame which holds the final fluxes
#allFBA = pd.DataFrame(0, columns= extract_files, index = [x.id for x in model.reactions])

extract = pd.read_csv(extract_file)

# reduce model
modRed = model.copy()
rxns = [x.id for x in modRed.reactions]

rm_rxns = [x for x in rxns if x not in extract.iloc[:,"rxn"].to_list()]
# security check
if (extract.shape[0] + len(rm_rxns) != len(rxns)) or len(rm_rxns) == 0 or extract.shape[0] == 0:
    raise Exception("Something went wrong with the model extraction, please check the extracted reactions file and the model, compare namespaces!")

rm_rxns = [modRed.reactions.get_by_id(x) for x in rm_rxns]
modRed.remove_reactions(rm_rxns)

diet = diet.set_index("ex_rxn",drop = False)
# adjust diet
for ex_rxn in modRed.exchanges:
    ex_rxn.lower_bound = 0
    if ex_rxn.id in list(diet.index):
        if not pd.isna(diet.loc[ex_rxn.id][1]):
            ex_rxn.lower_bound = -1*diet.loc[ex_rxn.id][1]

# if toggled, set the boundaries to max 1000/-1000
inf = float("Inf")
if setMaxBounds:
    for rxn in modRed.reactions:
        if rxn.upper_bound == inf:
            rxn.upper_bound = 1000
        if rxn.lower_bound == -inf:
            rxn.lower_bound = -1000

# run the pFBA
try:
    pFBA = cobra.flux_analysis.pfba(modRed,
            fraction_of_optimum= frac_opt)
    pFBA_fluxes = pd.DataFrame(pFBA.fluxes)
    #allFBA.loc[pFBA.fluxes.index,extrct] = pFBA.fluxes
except Exception as optErr:
    err = optErr
    pFBA_fluxes = pd.DataFrame({"fluxes": [np.nan]*len(modRed.reactions)},
            index = [x.id for x in modRed.reactions])
    #allFBA.loc[[x.id for x in modRed.reactions],extrct] = np.nan
    warnings.warn(str(optErr) + "\n Will return a table with only NaN")


# remove irrelevant data from the frame
#sel = np.apply_along_axis(sum,1,np.array(allFBA)) != 0
#allFBA = allFBA.loc[sel,]

# rename the columns to match the samplenames
#smplenames = [os.path.basename(x).split("-")[1].replace(".csv","") for x in allFBA.columns]
#allFBA.columns = smplenames
    
# write the output
pFBA_fluxes.to_csv(output)


