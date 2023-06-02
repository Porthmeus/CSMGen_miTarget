# Porthmeus
# 22.01.21

# This is an auxiliary script to correct the subsystems in the colormore22 model, as CobraToolbox skrews them up, when exporting models to a SBML file.

import cobra
import pandas as pd
import re
import os

# read the data
mod = cobra.io.read_sbml_model(snakemake.input["model"])
subsys = pd.read_csv(snakemake.input["subsystems"])
subsys_indx = subsys.set_index("subSys")

# remove the present groups
while len(mod.groups) > 0:
    groups = mod.groups
    mod.remove_groups(groups)

# and add the ones in the table


for group in subsys.subSys.unique():

    # find the reactions from the rxns column of the pd.frame, find the
    # corresponding reaction object in the model and add it back to the
    # group
    if not pd.isna(group):
        rxns_names = subsys_indx.loc[group].rxns.to_list()
        rxns_obj = [mod.reactions.get_by_id(x) for x in rxns_names]
        
        # crosscheck, wheter all reactions were found
        if not (rxns_names == [x.id for x in rxns_obj]):
            raise ValueError("rxns in list and the ones found in the model are not the same")
        group = cobra.core.Group(group, name = group, members= rxns_obj)
        mod.add_groups([group])




# write the model to disk
outfile = snakemake.output["outMod"]
cobra.io.write_sbml_model(mod, outfile+".tmp")

# correct the fbc labels of the exported model, such that matlab and R can read the gene names correctly again
with open(outfile, "w") as outf:
    outf.write("")

with open(outfile+".tmp", "r") as sbml:
    for line in sbml:
        label = line.find("fbc:label")
        if label != -1:
            lineSub = line[label+10:]
            codedAscii = re.findall(r"__\d*__", lineSub)
            corAscii = [chr(int(x.strip("_"))) for x in codedAscii]
            for character in range(len(codedAscii)):
                line = re.sub('(fbc:label.*)G_(.*)' + codedAscii[character], r'\1\2' + corAscii[character], line)

        with open(outfile, "a") as outf:
            outf.write(line)


# remove the temporary file
os.remove(outfile+".tmp")


