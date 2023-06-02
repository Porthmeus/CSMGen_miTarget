# Porthmeus
# 18.08.21

# calculate different thresholds for genes which should be used to infer the tissue specific metabolic model

import numpy as np
import pandas as pd
import re
import warnings


def getCoreGenes(array,
        global_lower = 0,
        global_upper = None,
        local = None,
        subset = None):
    '''calculate the genes which are considered expressed in a sample
    @array = pandas.DataFrame - contains the TPM values for each gene (row) and sample (column) - should contain only expression values and gene names as index and sample names as columnnames.
    @global_lower = float[0-100] - defines the global lower threshold as percentile from the whole data set, genes which have an expression lower than this threshold are considered inactive.
    @global_upper = float[0-100, > global_lower] - defines the global upper threshold for genes expression as percentile of the whole data set, genes which have an expression higher than this value are considered always active. If = None, this threshold will be not employed (genes are considered active either depending on the local and/or on the global_lower threshold). If the global_upper threshold is employed, the local threshold can not = None  -> will be automatically set to 50.
    @local = float[0-100, global_lower < local < global_upper] - defines a local expression threshold as percentile which is calculated individually for each gene across the data set. Genes which have an higher expression will be considered active, if expression > global_lower, conversely genes with expression < local will be considered inactive, if exprssion < global_upper.
    @subset = list - defines either the indeces or index names of the array rows to subset the data set to only those genes in the list (to subset only metabolic active genes for example)

    Value: a pandas.DataFrame of the same dimension than the input array (or the subset of it) with containing 0 if the gene in that sample is inactive or 1 if its active. Additionally a string containing a summary of thresholds applied.
    '''
    
    # sanity check
    if global_upper != None:
        if global_upper < global_lower:
            warnings.warn("Global lower threshold is higher than global upper threshold - will use only global lower threshold and discard local and global upper thresholds")
            global_upper = None
            local = None
            
        
    
    # subset data if necessary
    if subset != None:
        if type(subset[0]) == str:
            subset = [x for x in subset if x in array.index]
            array = array.loc[subset]
        elif type(subset[0]) == int:
            subset = [x for x in subset if x < array.shape[0]]
            array = array.iloc[subset]
        elif type(subset[0]) == bool:
            subset = subset[0:array.shape[0]]
            array = array[subset]

    # test the lower threshold
    glt = np.percentile(np.array(array), global_lower)
    resDF = array > glt   
    
    # if there is a local treshold, test for it
    if local != None:
        local_string = "L"+re.sub("\.0$","",str(local))
        for gene in array.index:
            localt = np.percentile(array.loc[gene], local)
            colVals = [True if all([x,y]) else False for x,y in zip(array.loc[gene] > localt,resDF.loc[gene])]
            resDF.loc[gene,:] = colVals
    else:
        local_string = ""

    # if there is a global upper threshold, test for it
    if global_upper != None:
        gu_string = "GU" + re.sub("\.0$","",str(global_upper))
        gut = np.percentile(np.array(array), global_upper)
        upperArray = array > gut
        resDF[upperArray] = True
    else:
        gu_string = ""
    
    # create the string
    out_string = re.sub("\|+$","","|".join(["GL" +  re.sub("\.0$","",str(global_lower)) ,local_string, gu_string]))
    
    return resDF, out_string


if __name__ == "__main__":

# start the actual rule

    dbg = False
    if dbg:
        TPM = "results/data/TPM_emed_future.csv"
        conversionTable = "resources/modelGenes2ensembl/colormore22Genes2ensembl.csv"
        threshold_string = "GL25|L50"
        out = "results/data/coreGenes/coreGenes." + threshold_string + "_colormore22.csv"
    else:
        # save stuff to the log file
        sys.stdout = sys.stderr = open(snakemake.log[0], 'w')

        TPM = snakemake.input["TPM"]
        conversionTable = snakemake.input["conversionTable"]
        threshold_string = snakemake.params["thresholds"]
        out = snakemake.output["CoreGeneMatrix"]

    # extract the thresholds
    thresholds = [float(re.sub("\D","",x)) if re.sub("\D","",x) != '' else None for x in threshold_string.split("|")]
    while len(thresholds) < 3:
        thresholds.append(None)
    gl,local,gu = thresholds

    # load data
    array = pd.read_csv(TPM, index_col = 0)
    convTab = pd.read_csv(conversionTable)

    # get the conversion table right - there might be genes missing in the TPM matrix, or there might be ambiguity between the two columns in the conversion table
    subset = [i for i in range(convTab.shape[0]) if convTab.iloc[i,1] in array.index]
    convTab = convTab.iloc[subset]
    subset = list(convTab.iloc[:,1])


    # get the core genes
    coreDF, thr = getCoreGenes(array = array, global_upper = gu, global_lower = gl, local = local, subset = subset)
    coreDF.index = list(convTab.iloc[:,0])

    # save the matrix
    coreDF.to_csv(out)
