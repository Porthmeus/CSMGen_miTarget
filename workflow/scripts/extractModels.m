function [mod_tiss] = extractModels(model_in, model_out, diet, consistent_Rxns, TPM_matrix, conversionTable, thrsld, cobraPath, objective)
    % read the model, extract the genes from the tpm matrix and calculate a
    % transcriptpome specific model per column in the TPM_matrix. Write these
    % models to an SBML file

    % @ model_in - a file which can be read to a CBT model by readCbModel()
    % @ model_out - the file where the results should be stored as a csv file for the list of reactions in the model, which are tissue specific
    % @ diet - a csv file with constrains for the diet - the first column should contain the rxns names, the second column values for the fluxes of these reaction. The function will set the lower bounds to all exchange reaction not in that list. For the reactions in the list the lower bound will be set to -1*value. If such a file is not given, it is assumed, that the model is correctly constrained.
    % @ consistent_Rxns - a table of the reaction which make up the consistent model. First column should contain the reaction-IDs of the consistent model. If no file is given the model is assumed to be consistent.
    % @ TPM_matrix - a csv file with one column with ENSEMBL-ids and a second column, containing the TPM values for the tissue which should be extracted. The values will be mapped to the gene names of the model
    % @ conversionTable - a csv file for the expression to gene mapping. Use a column called Gene_stable_ID for ENSEMBL identifier and HGNC_ID for HUGO genes names in the model. The values in these columns can be adjusted to need. A table for these conversions can be easily obtained from BIOMART
    % @ thrsld - a quantile threshold, above this threshold the reaction is considered active - Defaults to 75.
    % @ cobraPath - the location of the CobraToolbox
    % @ biomass - boolean, whether or not include the objective function(s) in the reaction list of active reactions

       
    try
        cobraPath = char(cobraPath);
        % load the CobraToolbox
        run(fullfile(cobraPath, 'initCobraToolbox(false)'));
        disp('initialized CobraToolbox');


        if ~(exist('model_in','var'))
            model_in = "/home/taube/Work/miTarget/eMed_FUTURE/Pipeline_MM/resources/models/colormore22.xml";
        end
        if ~(exist('model_out','var'))
            model_out = 'extractedRxn.csv';
        end
        if ~(exist('TPM_matrix','var'))
            TPM_matrix = "/home/taube/Work/miTarget/eMed_FUTURE/Pipeline_MM/results/data/SPLIT_TPM_emed_future/SPLIT_F02343_L1_S37_L005.csv"
        end
        if ~(exist('conversionTable','var'))
            conversionTable ='~/Work/miTarget/eMed_FUTURE/Pipeline_MM/resources/ensemble2HGNC.csv';
        end
        if ~(exist('thrsld','var'))
            thrsld = 75;
        end
        if ~(exist('objective','var'))
            objective = false;
        end
        
        
        % convert the input into character in case its string, and display the
        % values for the log file
        model_in = char(model_in)
        model_out = char(model_out)
        dietF = char(dietF)
        consistent_Rxns = char(consistent_Rxns)
        TPM_matrix = char(TPM_matrix)
        conversionTable = char(conversionTable)
        thrsld = string(thrsld)
        cobraPath = char(cobraPath)
        disp(objective)

 
        % load the model
        mod = readCbModel(model_in);
        
        if exist('diet','var')
            % load the diet
            diet = readtable(dietF);
            if(~all(class(diet{1:end,2} == 'double')))
                diet.val = str2double(diet{1:end,2});
            else
                diet.val = diet{1:end,2};
            end
            diet.rxnID = findRxnIDs(mod,diet{:,1});
            diet = diet(~(diet.rxnID == 0),:);
            diet.val(isnan(diet.val)) = 0;


            % reconstrain model
            ex_rxn = findExcRxns(mod);
            EX_rxn = startsWith(mod.rxns,'EX');
            ex = ex_rxn & EX_rxn;
            % remove possible fluxes on the exchange reactions
            mod = changeRxnBounds(mod, mod.rxns(ex), 0, 'l');
            mod = changeRxnBounds(mod, mod.rxns(diet.rxnID),-1*diet.val,'l');
        end

        if (exist('consistent_Rxns','var'))
            % remove inconsistent rxns
            consRxn = readtable(consistent_Rxns);
            rxnRm = mod.rxns(~ismember(mod.rxns, consRxn.rxn));
            mod = removeRxns(mod, rxnRm);
        end
        
        % it might happen if the model was created with cobraPy, that gene
        % naming is screwed after loading the model with matlab. Thus I will
        % fix, this issue here now, but in a very specific way for the current
        % pipeline... I hope this will not end up in to much confusion later
        % on, but otherwise I would have to change the in the source of cobraPy
        % or cobraToolbox, which are not good options either... sorry for the
        % ones reading this after encountering problems.
        mod.geneNames = replace(mod.geneNames, ...
            "G_HGNC__58__", "HGNC:");

        % load the conversion table to match 
        tbl_ensemble2HGNC = readtable(conversionTable);
        opt = detectImportOptions(TPM_matrix); % I don't know why, but is required to load the table correctly
        tbl_FutureReadCounts = readtable(TPM_matrix, opt);

        % extract relevant information from the transcriptome data
        tbl_ensemble2HGNC_modOnly = tbl_ensemble2HGNC(...
            ismember(tbl_ensemble2HGNC.HGNC_ID, ...
            mod.geneNames),:);
        [dummy, ia, dummy] = unique(tbl_ensemble2HGNC_modOnly.HGNC_ID);
        tbl_ensemble2HGNC_modOnly = tbl_ensemble2HGNC_modOnly(ia,:);
        tbl_FutureReadCounts_modOnly = tbl_FutureReadCounts( ismember(tbl_FutureReadCounts{:,1}, ...
            tbl_ensemble2HGNC_modOnly.Gene_stable_ID),:);
        
        % set the model gene names to the rownames, delete the original identifier
        tbl_FutureReadCounts_MONewNames = join(tbl_FutureReadCounts_modOnly, ...
            tbl_ensemble2HGNC_modOnly(:,[1,5]), ...
            'LeftKey', 1, ...
            'RightKey', 'Gene_stable_ID');
        tbl_FutureReadCounts_MONewNames.Row = tbl_FutureReadCounts_MONewNames.HGNC_ID;
        tbl_FutureReadCounts_MONewNames(:,1) = [];
        tbl_FutureReadCounts_MONewNames.HGNC_ID = [];

     
        % transform the data to reactions of the model
        data = struct;
        data.gene = tbl_FutureReadCounts_MONewNames.Row;
        data.value = tbl_FutureReadCounts_MONewNames{:,1};
        
        disp("Number of core genes")
        disp(data)

        reactionExpression = mapExpressionToReactions(mod, data);
        reaction_IDs = find(reactionExpression >= prctile(reactionExpression, str2double(thrsld)));

        
        % another annoying patch, but for some reason fastcore does not run, if
        % the reaction 'r0364' is in the core set. I will simply remove it from
        % the core set if it is in it.
        %id2remove = find(ismember(mod.rxns(reaction_IDs), ['r0364']));
        %reaction_IDs(id2remove) = [];
        if objective
            obj_rxn = find(mod.c);
            obj_rxn = obj_rxn(~ismember(obj_rxn, reaction_IDs));
            reaction_IDs(length(reaction_IDs) + (1:length(obj_rxn))) = obj_rxn;
        end

        % extract the model with fastcore and write it an sbml-file
        opt= struct('core', reaction_IDs,'solver','fastCore');
        tic
        mod_tiss = createTissueSpecificModel(mod, opt);
        toc
        % write the reactions to a csv file
        consRxn = struct;
        consRxn.rxn = mod_tiss.rxns;
        consRxn.id = findRxnIDs(mod, consRxn.rxn);
        writetable(struct2table(consRxn), model_out);

        %writeCbModel(mod_tiss, 'format', 'sbml','fileName', model_out);
        

    catch ME
        disp(getReport(ME));
        %exit;
    end
        
