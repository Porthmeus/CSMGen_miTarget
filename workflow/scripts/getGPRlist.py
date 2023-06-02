# Porthmeus
# 24.05.21

import cobra as cb
import pandas as pd

try:
    dbg
except NameError:
    print("Not in debug mode")
else:
    if dbg == True :



mod = cb.io.read_sbml_model("resources/models/colormore22.xml")
df = pd.DataFrame({"rxn" : [x.id for x in mod.reactions], "GPR" : [x.gene_reaction_rule for x in mod.reactions]})
df

df.to_csv("resources/colormore22_GPR.csv",index=False)

[x.build_reaction_string() for x in mod.reactions if x.id =="r1030"]
"r1030" in list(df.rxn)
